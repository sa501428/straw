#include "../../C++/straw.cpp"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

class PythonMatrixZoomData {
  public:
    PythonMatrixZoomData(std::string file, chromosome first, chromosome second,
                         std::string type, std::string normalization,
                         std::string matrixUnit, int32_t binSize)
        : fileName(std::move(file)), chr1(std::move(first)), chr2(std::move(second)),
          matrixType(std::move(type)), norm(std::move(normalization)),
          unit(std::move(matrixUnit)), resolution(binSize) {}

    std::vector<contactRecord> getRecords(int64_t x0, int64_t x1, int64_t y0, int64_t y1) const {
        return straw(matrixType, norm, fileName, location(chr1.name, x0, x1),
                     location(chr2.name, y0, y1), unit, resolution);
    }
    py::array getRecordsAsMatrix(int64_t x0, int64_t x1, int64_t y0, int64_t y1) const {
        return py::array(py::cast(strawAsMatrix(matrixType, norm, fileName,
                                                location(chr1.name, x0, x1),
                                                location(chr2.name, y0, y1), unit, resolution)));
    }
    std::vector<double> getNormVector(int32_t chromosomeIndex) const {
        std::vector<double> values;
        const chromosome *selected = chromosomeIndex == chr1.index ? &chr1 :
                                     chromosomeIndex == chr2.index ? &chr2 : nullptr;
        if (!selected || norm == "NONE") return values;
        if (!getNormalizationVectorForFile(fileName, selected->name, resolution, norm, values))
            throw std::runtime_error("normalization vector is not available");
        return values;
    }
    std::vector<double> getExpectedValues() const {
        std::vector<double> values;
        if (!getExpectedVectorForFile(fileName, chr1.name, resolution, norm, values))
            throw std::runtime_error("expected-value vector is not available");
        return values;
    }

  private:
    static std::string location(const std::string &name, int64_t start, int64_t end) {
        return name + ":" + std::to_string(start) + ":" + std::to_string(end);
    }
    std::string fileName;
    chromosome chr1, chr2;
    std::string matrixType, norm, unit;
    int32_t resolution;
};

class PythonHiCFile {
  public:
    explicit PythonHiCFile(std::string path) : fileName(std::move(path)) {
        if (straw_v10::isV10(fileName)) {
            straw_v10::File file(fileName);
            genomeID = file.genome();
        } else {
            legacy.reset(new HiCFile(fileName));
            genomeID = legacy->getGenomeID();
        }
        chromosomes = getChromosomesForFile(fileName);
        for (const auto &c : chromosomes) chromosomeByName[c.name] = c;
    }
    std::string getGenomeID() const { return genomeID; }
    std::vector<int32_t> getResolutions() const { return getResolutionsForFile(fileName); }
    std::vector<chromosome> getChromosomes() const { return chromosomes; }
    PythonMatrixZoomData *getMatrixZoomData(const std::string &first, const std::string &second,
                                            const std::string &type, const std::string &normalization,
                                            const std::string &unit, int32_t resolution) const {
        auto a = chromosomeByName.find(first), b = chromosomeByName.find(second);
        if (a == chromosomeByName.end() || b == chromosomeByName.end())
            throw std::runtime_error("chromosome not found in the file");
        return new PythonMatrixZoomData(fileName, a->second, b->second, type, normalization,
                                        unit, resolution);
    }

  private:
    std::string fileName;
    std::string genomeID;
    std::unique_ptr<HiCFile> legacy;
    std::vector<chromosome> chromosomes;
    std::map<std::string, chromosome> chromosomeByName;
};

PYBIND11_MODULE(hicstraw, m) {
    m.doc() = "Fast tool for reading .hic files; see https://github.com/aidenlab/straw";
    m.def("strawC", &straw, "get contact records");
    m.def("straw", &straw, "get contact records");
    m.def("strawAsMatrix", [](const std::string &type, const std::string &norm,
                              const std::string &file, const std::string &chr1,
                              const std::string &chr2, const std::string &unit, int32_t resolution) {
        return py::array(py::cast(strawAsMatrix(type, norm, file, chr1, chr2, unit, resolution)));
    }, "get contact records in a numpy matrix");
    py::class_<contactRecord>(m, "contactRecord")
        .def(py::init<>()).def_readwrite("binX", &contactRecord::binX)
        .def_readwrite("binY", &contactRecord::binY).def_readwrite("counts", &contactRecord::counts);
    py::class_<chromosome>(m, "chromosome")
        .def(py::init<>()).def_readwrite("name", &chromosome::name)
        .def_readwrite("index", &chromosome::index).def_readwrite("length", &chromosome::length);
    py::class_<PythonMatrixZoomData>(m, "MatrixZoomData")
        .def("getRecords", &PythonMatrixZoomData::getRecords)
        .def("getRecordsAsMatrix", &PythonMatrixZoomData::getRecordsAsMatrix)
        .def("getNormVector", &PythonMatrixZoomData::getNormVector)
        .def("getExpectedValues", &PythonMatrixZoomData::getExpectedValues);
    py::class_<PythonHiCFile>(m, "HiCFile")
        .def(py::init<std::string>()).def("getChromosomes", &PythonHiCFile::getChromosomes)
        .def("getResolutions", &PythonHiCFile::getResolutions).def("getGenomeID", &PythonHiCFile::getGenomeID)
        .def("getMatrixZoomData", &PythonHiCFile::getMatrixZoomData, py::return_value_policy::take_ownership);
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
