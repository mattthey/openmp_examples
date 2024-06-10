#include <iostream>
#include <vector>
#include <fstream>
#include <omp.h>
#include <string>
#include <sstream>
#include <ctime>
#include <iomanip>
#include <functional>

double calculatePolygonArea(const std::vector<double>& x, const std::vector<double>& y) {
    const unsigned long n = x.size();
    double sum = 0.0;

    omp_set_num_threads(2);
    #pragma omp parallel for reduction(+:sum) default(shared)
    for (int i = 0; i < n; ++i) {
        unsigned long j = (i + 1) % n; // Next vertex index (wraps around using modulo)
        sum += (x[i] * y[j] - y[i] * x[j]);
    }

    return std::abs(sum) * 0.5;
}

double calculatePolygonArea2(const std::vector<double>& x, const std::vector<double>& y) {
    const unsigned long n = x.size();
    double sum = 0.0;

    omp_set_num_threads(2);
    #pragma omp parallel
    {
        double local_sum = 0.0;

        #pragma omp for nowait
        for (unsigned long i = 0; i < n; ++i) {
            unsigned long j = (i + 1) % n; // Next vertex index (wraps around using modulo)
            local_sum += (x[i] * y[j] - y[i] * x[j]);
        }

        #pragma omp atomic
        sum += local_sum;
    }

    return std::abs(sum) * 0.5;
}

double calculatePolygonAreaWithoutOMP(const std::vector<double>& x, const std::vector<double>& y) {
    const unsigned long n = x.size();
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        unsigned long j = (i + 1) % n; // Next vertex index (wraps around using modulo)
        sum += (x[i] * y[j] - y[i] * x[j]);
    }

    return std::abs(sum) * 0.5;
}

void executeAndPrintDuration(
        std::vector<double>& x,
        std::vector<double>& y,
        std::string methodName,
        const std::function<double(std::vector<double>, std::vector<double>)>& lambda) {
    clock_t start = clock();

    double area = lambda(x, y);

    clock_t endTime = clock();
    double elapsed_seconds = static_cast<double>(endTime - start) / CLOCKS_PER_SEC;

    std::cout << "Elapsed time for calculated area with " << methodName << ": "
              << elapsed_seconds << " sec"
              << " calculated area " << area
              << std::endl;
}

void slowReadFile(char* fileName,
                  std::vector<double>& x,
                  std::vector<double>& y
) {
    std::cout << "read data from " << fileName << std::endl;
    std::ifstream inputFile(fileName);

    clock_t start = clock();
    if (!inputFile) {
        std::cerr << "Unable to open file input.txt";
        exit(1);
    }

    // Read file into a buffer
    inputFile.seekg(0, std::ios::end);
    long fileSize = inputFile.tellg();
    inputFile.seekg(0, std::ios::beg);

    std::vector<char> buffer(fileSize);
    inputFile.read(buffer.data(), fileSize);
    inputFile.close();

    char* data = buffer.data();
    char* end = data + fileSize;

    while (data < end) {
        double xCoord, yCoord;
        xCoord = std::strtod(data, &data);
        yCoord = std::strtod(data, &data);

        x.emplace_back(xCoord);
        y.emplace_back(yCoord);

        // Skip any trailing whitespace or newlines
        while (data < end && std::isspace(*data)) {
            ++data;
        }
    }

    inputFile.close();

    clock_t endTime = clock();
    double elapsed_seconds = static_cast<double>(endTime - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time for read input data: " << std::fixed << std::setprecision(9) << elapsed_seconds << " seconds." << std::endl;
}

int main(int argc, char* argv[]) {
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(1e9);
    y.reserve(1e9);

    slowReadFile(argv[1], x, y);

    executeAndPrintDuration(x, y,
                            "old openMP",
                            [](std::vector<double> xx, std::vector<double> yy) { return calculatePolygonArea(xx, yy); });
    executeAndPrintDuration(x, y,
                            "new openMP",
                            [](std::vector<double> xx, std::vector<double> yy) { return calculatePolygonArea2(xx, yy); });
    executeAndPrintDuration(x, y,
                            "without openMP",
                            [](std::vector<double> xx, std::vector<double> yy) { return calculatePolygonAreaWithoutOMP(xx, yy); });

    return 0;
}
