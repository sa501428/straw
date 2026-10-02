#pragma once
// Bounded local/HTTP byte ranges shared by the V9 vector and V10 readers.
#include "v10_binary.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <mutex>
#include <curl/curl.h>
namespace straw_file {
using straw_v10::Bytes;
using straw_v10::require;
using straw_v10::add;
using straw_v10::allocationLimit;
struct Source {
    std::ifstream file;
    CURL *curl = nullptr;
    uint64_t size = 0, requested = 0, responseStart = UINT64_MAX, responseEnd = 0;
    Bytes response;
    static size_t write(char *p, size_t a, size_t b, void *context) {
        auto &s = *static_cast<Source *>(context);
        if (b && a > SIZE_MAX / b)
            return 0;
        size_t n = a * b;
        if (n > s.requested - s.response.size())
            return 0;
        try {
            s.response.insert(s.response.end(), p, p + n);
        } catch (...) {
            return 0;
        }
        return n;
    }
    static size_t header(char *p, size_t a, size_t b, void *context) {
        auto &s = *static_cast<Source *>(context);
        size_t n = a * b;
        try {
            std::string line(p, n), lower = line;
            std::transform(lower.begin(), lower.end(), lower.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (lower.compare(0, 14, "content-range:") == 0) {
                unsigned long long first, last, total;
                if (std::sscanf(lower.c_str(), "content-range: bytes %llu-%llu/%llu", &first, &last,
                                &total) == 3) {
                    s.responseStart = first;
                    s.responseEnd = last;
                    s.size = total;
                }
            }
        } catch (...) {
            return 0;
        }
        return n;
    }
    explicit Source(const std::string &path) {
        if (path.compare(0, 7, "http://") == 0 || path.compare(0, 8, "https://") == 0) {
            static std::once_flag once;
            std::call_once(once, [] {
                require(curl_global_init(CURL_GLOBAL_DEFAULT) == CURLE_OK,
                        "curl initialization failed");
            });
            curl = curl_easy_init();
            require(curl != nullptr, "curl initialization failed");
            curl_easy_setopt(curl, CURLOPT_URL, path.c_str());
            curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            curl_easy_setopt(curl, CURLOPT_MAXREDIRS, 10L);
            curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
            curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 30L);
            curl_easy_setopt(curl, CURLOPT_LOW_SPEED_LIMIT, 1L);
            curl_easy_setopt(curl, CURLOPT_LOW_SPEED_TIME, 60L);
            curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write);
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, this);
            curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, header);
            curl_easy_setopt(curl, CURLOPT_HEADERDATA, this);
            curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw-v10");
        } else {
            file.open(path, std::ios::binary);
            require(bool(file), "cannot open " + path);
            file.seekg(0, std::ios::end);
            auto end = file.tellg();
            require(end >= 0, "cannot stat input");
            size = uint64_t(end);
        }
    }
    ~Source() {
        if (curl)
            curl_easy_cleanup(curl);
    }
    void interval(uint64_t pos, uint64_t len) const {
        require(len > 0 && pos <= size && len <= size - pos, "file interval out of bounds");
    }
    Bytes read(uint64_t pos, uint64_t len) {
        require(len > 0 && len <= allocationLimit, "record exceeds allocation limit");
        if (size)
            interval(pos, len);
        if (!curl) {
            Bytes out(static_cast<size_t>(len));
            file.clear();
            file.seekg(static_cast<std::streamoff>(pos));
            file.read(reinterpret_cast<char *>(out.data()), static_cast<std::streamsize>(len));
            require(bool(file), "short file read");
            return out;
        }
        requested = len;
        response.clear();
        responseStart = UINT64_MAX;
        std::string range = std::to_string(pos) + "-" + std::to_string(add(pos, len) - 1);
        curl_easy_setopt(curl, CURLOPT_RANGE, range.c_str());
        CURLcode rc = curl_easy_perform(curl);
        long status = 0;
        curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &status);
        require(rc == CURLE_OK, std::string("HTTP range read failed: ") + curl_easy_strerror(rc));
        require(status == 206 && responseStart == pos && responseEnd == pos + len - 1 &&
                    response.size() == len,
                "server must return the exact requested HTTP byte range (206)");
        interval(pos, len);
        return std::move(response);
    }
};
}
