#include "hbs.h"
#include <cstdio>
#include <set>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

bool isHbsPath(const std::string& path) {
    return path.size() >= 7 && path.compare(path.size() - 7, 7, ".hbs.gz") == 0;
}

void HbsWriter::bytes(const void* data, unsigned size) {
    if (gzwrite(file_, data, size) != static_cast<int>(size))
        throw std::runtime_error("HBS: failed writing compressed output");
}

void HbsWriter::number(uint64_t value, unsigned size) {
    unsigned char b[8];
    for (unsigned i = 0; i < size; ++i) b[i] = static_cast<unsigned char>(value >> (8 * i));
    bytes(b, size);
}

HbsWriter::HbsWriter(const std::string& input, const std::string& output,
                     uint32_t resolution, const std::vector<chromosome>& chromosomes)
    : output_(output), resolution_(resolution) {
    struct stat in{}, out{};
    if (input == output || (stat(input.c_str(), &in) == 0 && stat(output.c_str(), &out) == 0 &&
                           in.st_dev == out.st_dev && in.st_ino == out.st_ino))
        throw std::runtime_error("HBS: input and output must differ");
    if (!resolution || resolution > INT32_MAX || chromosomes.size() > 65536)
        throw std::runtime_error("HBS: invalid resolution or too many chromosomes");
    std::set<std::string> names;
    size_t headerBytes = 20;
    for (const auto& c : chromosomes) {
        headerBytes += 10 + c.name.size();
        if (c.name.empty() || c.name.size() > 4096 || c.name.find('\0') != std::string::npos ||
            c.length <= 0 || !names.insert(c.name).second || headerBytes > 16 * 1024 * 1024 ||
            !ids_.emplace(c.index, static_cast<uint16_t>(ids_.size())).second)
            throw std::runtime_error("HBS: invalid chromosome table");
    }
    std::string pattern = output + ".tmp-XXXXXX";
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back(0);
    int fd = mkstemp(name.data());
    if (fd < 0) throw std::runtime_error("HBS: cannot create output beside " + output);
    temporary_ = name.data();
    file_ = gzdopen(fd, "wb");
    if (!file_) {
        close(fd);
        std::remove(temporary_.c_str());
        throw std::runtime_error("HBS: cannot open gzip output");
    }
    try {
        gzbuffer(file_, 128 * 1024);
        bytes("HICBS\0\r\n", 8);
        number(1, 2); // version
        number(0, 2); // reserved flags
        number(resolution, 4);
        number(chromosomes.size(), 4);
        for (const auto& c : chromosomes) {
            number(c.name.size(), 2);
            bytes(c.name.data(), static_cast<unsigned>(c.name.size()));
            number(c.length, 8);
        }
    } catch (...) {
        gzclose(file_);
        file_ = nullptr;
        std::remove(temporary_.c_str());
        throw;
    }
}

HbsWriter::~HbsWriter() {
    if (file_) gzclose(file_);
    if (!temporary_.empty()) std::remove(temporary_.c_str());
}

void HbsWriter::record(const chromosome& a, uint64_t x, const chromosome& b, uint64_t y, uint64_t count) {
    if (!count) return;
    if (x >= static_cast<uint64_t>(a.length) || y >= static_cast<uint64_t>(b.length) ||
        x % resolution_ || y % resolution_ || x / resolution_ > UINT32_MAX || y / resolution_ > UINT32_MAX)
        throw std::runtime_error("HBS: bin start outside chromosome or coordinate range");
    unsigned char data[22];
    unsigned size = 0;
    auto append = [&](uint64_t value, unsigned width) {
        for (unsigned i = 0; i < width; ++i) data[size++] = static_cast<unsigned char>(value >> (8 * i));
    };
    append(ids_.at(a.index), 2);
    append(x / resolution_, 4);
    append(ids_.at(b.index), 2);
    append(y / resolution_, 4);
    append(count < 65535 ? count : 65535, 2);
    if (count >= 65535) append(count, 8);
    bytes(data, size);
}

void HbsWriter::finish() {
    int status = gzclose(file_);
    file_ = nullptr;
    if (status != Z_OK) throw std::runtime_error("HBS: failed closing gzip output");
    if (std::rename(temporary_.c_str(), output_.c_str()) != 0)
        throw std::runtime_error("HBS: failed publishing " + output_);
    temporary_.clear();
}
