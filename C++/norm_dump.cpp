#include "norm_vectors.h"
#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace {
const char *usage =
    "Usage: straw dump-norms <input.hic> --output-dir <directory> [--norm <TYPE> ...]\n"
    "   or: straw dump-norms <input.hic> --norm <TYPE> --output <vectors.txt>\n"
    "Writes HIC_NORM_VECTORS 1 files for v10 addnorm --norm-file.\n"
    "Default selection: every stored normalization, including custom names.\n";
void check(bool ok, const std::string &message) {
    if (!ok) throw std::runtime_error(message);
}
void directory(const std::string &path) {
    check(!path.empty(), "empty output directory");
    struct stat st;
    if (stat(path.c_str(), &st) == 0) {
        check(S_ISDIR(st.st_mode), "output directory is not a directory: " + path);
        return;
    }
    auto slash = path.find_last_of('/');
    if (slash != std::string::npos && slash > 0) directory(path.substr(0, slash));
    check(mkdir(path.c_str(), 0777) == 0 || (errno == EEXIST && stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode)),
          "cannot create output directory: " + path);
}
std::string filename(const std::string &name) {
    std::string out;
    const char *hex = "0123456789ABCDEF";
    for (unsigned char b : name) {
        if ((b >= 'a' && b <= 'z') || (b >= 'A' && b <= 'Z') ||
            (b >= '0' && b <= '9') || b == '_' || b == '-') out.push_back(char(b));
        else { out += '%'; out += hex[b >> 4]; out += hex[b & 15]; }
    }
    return out + ".norm.txt";
}
void textName(const std::string &name) {
    check(!name.empty() && name.find_first_of("\r\n") == std::string::npos,
          "normalization/chromosome name cannot contain line breaks in text output");
}
struct Output {
    std::string path, temporary;
    std::ofstream stream;
    bool published = false, committed = false;
    explicit Output(const std::string &target) : path(target) {
        struct stat st;
        check(lstat(path.c_str(), &st) != 0 && errno == ENOENT,
              "output already exists or cannot be accessed: " + path);
        std::string pattern = path + ".tmp.XXXXXX";
        std::vector<char> bytes(pattern.begin(), pattern.end()); bytes.push_back(0);
        int fd = mkstemp(bytes.data());
        check(fd >= 0, "cannot create output: " + path);
        close(fd); temporary = bytes.data();
        stream.open(temporary);
        if (!stream) { unlink(temporary.c_str()); throw std::runtime_error("cannot open output: " + path); }
    }
    ~Output() {
        if (!temporary.empty()) unlink(temporary.c_str());
        if (published && !committed) unlink(path.c_str());
    }
    void finish() {
        stream.close(); check(bool(stream), "error writing output: " + path);
    }
    void publish() {
        // Link rather than rename: an output created concurrently must never
        // be overwritten. Both paths are in the same directory/filesystem.
        check(link(temporary.c_str(), path.c_str()) == 0,
              "cannot publish output (existing files are not overwritten): " + path);
        published = true; unlink(temporary.c_str()); temporary.clear();
    }
};
}
int dumpNormsMain(int argc, char **argv) {
    if (argc == 3 && (std::string(argv[2]) == "--help" || std::string(argv[2]) == "-h")) {
        std::cout << usage; return 0;
    }
    check(argc >= 5, usage);
    std::string input = argv[2], output, outputDir;
    std::set<std::string> selected;
    for (int i = 3; i < argc; ++i) {
        std::string option = argv[i];
        check(i+1 < argc, "missing value for " + option);
        std::string value = argv[++i];
        if (option == "--norm") {
            check(!value.empty() && value != "NONE", "--norm requires a stored normalization type (not NONE)");
            check(selected.insert(value).second, "duplicate --norm " + value);
        } else if (option == "--output-dir" && outputDir.empty()) outputDir = value;
        else if ((option == "--output" || option == "-o") && output.empty()) output = value;
        else throw std::runtime_error("unknown or duplicate option " + option);
    }
    check(output.empty() != outputDir.empty(), "specify exactly one of --output-dir or --output");
    check(output.empty() || selected.size() == 1, "--output requires exactly one --norm");
    StrawNormalizationVectors vectors(input);
    std::map<std::string, std::vector<size_t>> groups;
    for (size_t i = 0; i < vectors.entries().size(); ++i) {
        const auto &e = vectors.entries()[i];
        if (e.normalization == "NONE" || e.chromosome == "ALL" || e.chromosome == "All" || e.chromosome == "all") continue;
        if (selected.empty() || selected.count(e.normalization)) groups[e.normalization].push_back(i);
    }
    for (const auto &name : selected) check(groups.count(name), "no stored vectors for normalization " + name);
    check(!groups.empty(), "file contains no exportable normalization vectors");
    if (!outputDir.empty()) directory(outputDir);
    std::vector<std::unique_ptr<Output>> files;
    std::set<std::string> names;
    for (const auto &group : groups) {
        std::string base = filename(group.first), folded = base;
        std::transform(folded.begin(), folded.end(), folded.begin(), [](unsigned char c) {
            return c >= 'A' && c <= 'Z' ? char(c+32) : char(c);
        });
        // Disambiguate case-only names even on case-insensitive filesystems.
        unsigned suffix = 0;
        while (!names.insert(folded).second) {
            auto prefix = std::to_string(++suffix) + "-";
            base = prefix + filename(group.first);
            folded = prefix + filename(group.first);
            std::transform(folded.begin(), folded.end(), folded.begin(), [](unsigned char c) {
                return c >= 'A' && c <= 'Z' ? char(c+32) : char(c);
            });
        }
        auto target = output.empty() ? outputDir + "/" + base : output;
        files.emplace_back(new Output(target));
        auto &file = *files.back(); auto &out = file.stream;
        out << "HIC_NORM_VECTORS 1\n# Exact float32 divisors; normalized = raw / (N[i] * N[j]).\n";
        for (size_t index : group.second) {
            const auto &e = vectors.entries()[index];
            textName(e.normalization); textName(e.chromosome);
            out << "vector " << std::quoted(e.normalization) << ' ' << std::quoted(e.chromosome)
                << ' ' << e.unit << ' ' << e.resolution << '\n';
            if (e.count != e.binCount) out << "source-length " << e.count << '\n';
            vectors.stream(index, [&](uint32_t bits) {
                out << "bits:" << std::hex << std::setw(8) << std::setfill('0') << bits << std::dec << '\n';
            });
            out << "end\n";
            check(bool(out), "error writing normalization output: " + target);
        }
        file.finish();
    }
    // Stage every selected type before exposing any completed output.
    for (auto &file : files) file->publish();
    for (auto &file : files) file->committed = true;
    size_t i = 0;
    for (const auto &group : groups) {
        std::cerr << "Exported " << group.first << ": " << group.second.size()
                  << " chromosome/resolution vectors to " << files[i++]->path << '\n';
    }
    return 0;
}
