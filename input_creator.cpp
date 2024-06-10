#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip> // for std::setprecision

using namespace std;

struct Point {
    long double x, y;
};

void generate_polygon_coordinates(long vertices_count, const string& filename = "input4.txt") {
    if (vertices_count < 3) {
        throw invalid_argument("A polygon must have at least 3 vertices");
    }

    long double radius = 100000.0;
    long double center_x = 0.0, center_y = 0.0;
    long double angle_step = 2.0 * M_PI / vertices_count;

    ofstream file(filename);
    for (long i = 0; i < vertices_count; ++i) {
        long double angle = i * angle_step;
        long double x = center_x + radius * cos(angle);
        long double y = center_y + radius * sin(angle);

        Point point = {roundl(x * 100000) / 100000, roundl(y * 100000) / 100000};

        file << fixed << setprecision(5) << point.x << " " << point.y << "\n";
    }
    file.close();
}

int main(int argc, char* argv[]) {
    try {
        long vertices_count;
        cout << "Enter the number of vertices for the polygon: ";
        cin >> vertices_count;
        cout << "OK " << vertices_count << endl;

        generate_polygon_coordinates(vertices_count, argv[1]);

        cout << "Coordinates of the polygon with " << vertices_count << " vertices have been written to input.txt" << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
