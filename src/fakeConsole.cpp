//
// Created by Allen Li on 2019-02-17.
//
#include <iostream>
#include "../shared/eigen/Eigen/Core"
#include "../shared/eigen/Eigen/Dense"
#include "../shared/eigen/Eigen/Geometry"

using namespace std;
using namespace Eigen;

int main() {
    MatrixXd b(2,2);
    b << 2, 3,
            1, 4;
    Vector2i temp(3,4);

    cout << "Hello, World!" << endl;
    cout << temp.norm() << endl;

    cout << b(1,0) << endl;

    cout << "========================" << endl;
    cout << b.rows() << endl;
    cout << b.cols() << endl;


    return 0;
}