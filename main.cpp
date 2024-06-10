#include <iostream>
#include <functional>
#include <vector>
#include <fstream>
#include <omp.h>
#include <string>
#include <sstream>
#include <ctime>
#include <iomanip>

double calculatePolygonArea(const std::vector<double>& x, const std::vector<double>& y, const int countThread) {
    const unsigned long n = x.size();
    double sum = 0.0;

    omp_set_num_threads(countThread);
    #pragma omp parallel for reduction(+:sum) default(shared)
    for (int i = 0; i < n; ++i) {
        unsigned long j = (i + 1) % n; // Next vertex index (wraps around using modulo)
        sum += (x[i] * y[j] - y[i] * x[j]);
    }

    return std::abs(sum) * 0.5;
}

double calculatePolygonArea2(const std::vector<double>& x, const std::vector<double>& y, const int countThread) {
    const unsigned long n = x.size();
    double sum = 0.0;

    omp_set_num_threads(countThread);
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

double calculatePolygonAreaWithoutOMP(const std::vector<double>& x, const std::vector<double>& y, const int countThread) {
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
        const int countThread,
        std::string methodName,
        const std::function<double(std::vector<double>, std::vector<double>, const int)>& lambda) {
    clock_t start = clock();

    double area = lambda(x, y, countThread);

    clock_t endTime = clock();
    double elapsed_seconds = static_cast<double>(endTime - start) / CLOCKS_PER_SEC;

    std::cout << "Elapsed time for calculated area with " << methodName << ": "
              << elapsed_seconds << " sec"
              << " count threads " << countThread
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

    std::string line;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        double xCoord, yCoord;
        if (!(iss >> xCoord >> yCoord)) {
            std::cerr << "Failed to parse line\n";
            exit(1);
        }
        x.push_back(xCoord);
        y.push_back(yCoord);
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

    for (int i = 2; i < 17; i *= 2) {
        executeAndPrintDuration(x, y, i,
                            "old openMP",
                            [](std::vector<double> xx, std::vector<double> yy, const int ii) { return calculatePolygonArea(xx, yy, ii); });
        executeAndPrintDuration(x, y, i,
                                "new openMP",
                                [](std::vector<double> xx, std::vector<double> yy, const int ii) { return calculatePolygonArea2(xx, yy, ii); });
    }
    executeAndPrintDuration(x, y, 1,
                            "without openMP",
                            [](std::vector<double> xx, std::vector<double> yy, const int ii) { return calculatePolygonAreaWithoutOMP(xx, yy, ii); });

    return 0;
}
