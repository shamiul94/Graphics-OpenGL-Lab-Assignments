//
// Created by shamiul93 on 7/13/19.
//

#include <bits/stdc++.h>

#define pi acos(-1)

using namespace std;

class Matrix;

class HomogeneousPointVector {
public:
    double x, y, z, w;

    HomogeneousPointVector() {
        x = y = z = 0;
//        w = 1;
    }

    HomogeneousPointVector(double xx, double yy, double zz) {
        x = xx;
        y = yy;
        z = zz;
//        w = 1;
    }

    double getModulus() {
        return sqrt(x * x + y * y + z * z);
    }

    void makeUnitVector() { // normalizing
        double m = getModulus();
        x /= m;
        y /= m;
        z /= m;
    }

    HomogeneousPointVector operator+(HomogeneousPointVector const &obj) {
        HomogeneousPointVector ret;
        ret.x = x + obj.x;
        ret.y = y + obj.y;
        ret.z = z + obj.z;
//        ret.w = w + obj.w;

        return ret;
    }


    HomogeneousPointVector operator-(HomogeneousPointVector const &obj) {
        HomogeneousPointVector ret;
        ret.x = x - obj.x;
        ret.y = y - obj.y;
        ret.z = z - obj.z;
//        ret.w = w - obj.w;

        return ret;
    }

    HomogeneousPointVector operator*(HomogeneousPointVector const &obj) { // cross multiplication
        HomogeneousPointVector ret;
        ret.x = y * obj.z - z * obj.y;
        ret.y = z * obj.x - x * obj.z;
        ret.z = x * obj.y - y * obj.x;
        return ret;
    }

    HomogeneousPointVector operator*(double const &obj) { // scalar multiplication
        HomogeneousPointVector ret;
        ret.x = x * obj;
        ret.y = y * obj;
        ret.z = z * obj;
        return ret;
    }

    double dotMultiplication(HomogeneousPointVector vec) {
        return (x * vec.x + y * vec.y + z * vec.z);
    }
};


class Matrix {
public:
    vector<vector<double>> matrix;

    void init() {
        for (int i = 0; i < 4; i++) {
            matrix.emplace_back(4, 0); // value == 0
        }
    }

    Matrix() {
        init();
    }

    Matrix(double a) { // useful for making identity matrix
        init();

        matrix = vector<vector<double>>{
                {a, 0, 0, 0},
                {0, a, 0, 0},
                {0, 0, a, 0},
                {0, 0, 0, 1}
        };
    }

    Matrix(int row, int colm) { // useful for making custom row column matrix
        init();
        for (int i = 0; i < row; i++) {
            matrix.emplace_back(colm, 0); // value == 0
        }
    }

    Matrix(double a, double b, double c) { // useful for making scaling matrix
        init();
        matrix = vector<vector<double>>{
                {a, 0, 0, 0},
                {0, b, 0, 0},
                {0, 0, c, 0},
                {0, 0, 0, 1}
        };
    }

    Matrix(HomogeneousPointVector vec) { // useful for making translation matrix
        init();

        matrix = vector<vector<double>>{
                {1, 0, 0, vec.x},
                {0, 1, 0, vec.y},
                {0, 0, 1, vec.z},
                {0, 0, 0, 1}
        };
    }

    Matrix(bool convertPointVecToMatrix, HomogeneousPointVector vec) { // useful to convert vector to matrix
        for (int i = 0; i < 4; i++) {
            matrix.emplace_back(1, 0);
        }

        matrix = vector<vector<double>>{
                {vec.x},
                {vec.y},
                {vec.z},
                {1}
        };
    }

    Matrix(HomogeneousPointVector l, HomogeneousPointVector r,
           HomogeneousPointVector u) { // useful for making rotation matrix
        init();
        matrix = vector<vector<double>>{
                {r.x,  r.y,  r.z,  0},
                {u.x,  u.y,  u.z,  0},
                {-l.x, -l.y, -l.z, 0},
                {0,    0,    0,    1}
        };
    }

    Matrix(double near, double far, double r, double t) { // used for projection matrix
        init();
        matrix = vector<vector<double>>{
                {near / r, 0,        0,                            0},
                {0,        near / t, 0,                            0},
                {0,        0,        -(far + near) / (far - near), -(2 * far * near) / (far - near)},
                {0,        0,        -1,                           0}
        };
    }

    Matrix operator*(Matrix const &m) {
        int row1 = static_cast<int>(matrix.size()), col1 = static_cast<int>(matrix[0].size());
        int row2 = static_cast<int>(m.matrix.size()), col2 = static_cast<int>(m.matrix[0].size());

        if (col1 != row2) {
            cout << "Can not multiply" << endl;
            exit(0);
        }

        Matrix returnMatrix(row1, col2);

        for (int i = 0; i < row1; i++) {
            for (int j = 0; j < col2; j++) {
                returnMatrix.matrix[i][j] = 0;

                for (int k = 0; k < col1; k++) {
                    returnMatrix.matrix[i][j] += matrix[i][k] * m.matrix[k][j];
                }
            }
        }
        return returnMatrix;
    }

    Matrix operator*(double const &obj) {
        auto row = static_cast<int>(matrix.size()), col = static_cast<int>(matrix[0].size());

        Matrix ret(row, col);

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                ret.matrix[i][j] = matrix[i][j] * obj;
            }
        }

        return ret;
    }


    void printMatrix() {
        for (auto row: matrix) {
            for (auto cell: row) {
                cout << cell << ' ';
            }
            cout << endl;
        }
    }
};


Matrix convertPointToMatrix(HomogeneousPointVector vec) {
    return Matrix(true, vec);
}

HomogeneousPointVector transformPoint(Matrix matrix, HomogeneousPointVector vec) {
    Matrix columnVectorMatrix = Matrix(true, vec);

    Matrix transformedFinalMatrix = matrix * columnVectorMatrix;

    HomogeneousPointVector transformedVector;

    double w = transformedFinalMatrix.matrix[3][0];

    transformedVector.x = transformedFinalMatrix.matrix[0][0] / w;
    transformedVector.y = transformedFinalMatrix.matrix[1][0] / w;
    transformedVector.z = transformedFinalMatrix.matrix[2][0] / w;

    return transformedVector;
}


HomogeneousPointVector
rodriguesFormula(HomogeneousPointVector X, HomogeneousPointVector A, double angle) {

    double cosValue = cos(angle * pi / 180);
    double sinValue = sin(angle * pi / 180);

    HomogeneousPointVector returnVector;

    returnVector = X * cosValue;

    //(1-cos theta) (a.x)
    double dotProductValue = (1 - cosValue) * A.dotMultiplication(X);

    returnVector = returnVector + (A * dotProductValue);

    HomogeneousPointVector crossProductVector = A * X;

    returnVector = returnVector + (crossProductVector * sinValue);

    return returnVector;
}


stack<Matrix> smallStack;
stack<stack<Matrix>> bigStack;

HomogeneousPointVector eye, look, up;
HomogeneousPointVector l, u, r;
double fovY, aspectRatio, near, far;


int main() {
    ifstream scin("scene.txt");
    ofstream outStage1("stage1.txt");
    ofstream outStage2("stage2.txt");
    ofstream outStage3("stage3.txt");

    scin >> eye.x >> eye.y >> eye.z;
    scin >> look.x >> look.y >> look.z;
    scin >> up.x >> up.y >> up.z;
    eye.w = 1, look.w = 1, up.w = 1;

    l = look - eye;
    l.makeUnitVector();
    r = l * up; // cross product
    r.makeUnitVector();
    u = r * l;
    u.makeUnitVector();

    Matrix translationMatrixForEye(eye * -1);
    Matrix rotationMatrixForEye(l, r, u);

    Matrix V = rotationMatrixForEye * translationMatrixForEye; // V = RT

    scin >> fovY >> aspectRatio >> near >> far;


    double fovX = fovY * aspectRatio;


    double t = near * tan((fovY / 2) * pi / 180.0);

    double r = near * tan((fovX / 2) * pi / 180.0);

    Matrix projectionMatrix(near, far, r, t);

    Matrix I(1);
    smallStack.push(I);

    string command;
    //while true
    while (true) {
        //input command
        scin >> command;


        //if command = “triangle”
        if (command == "triangle") {

            //input three trianglePoints
            HomogeneousPointVector trianglePoints[3];
            scin >> trianglePoints[0].x >> trianglePoints[0].y >> trianglePoints[0].z;
            scin >> trianglePoints[1].x >> trianglePoints[1].y >> trianglePoints[1].z;
            scin >> trianglePoints[2].x >> trianglePoints[2].y >> trianglePoints[2].z;

            //for each three point P
            for (const auto &trianglePoint : trianglePoints) {

                //P’ <- transformPoint(S.top,P)
                HomogeneousPointVector model = transformPoint(smallStack.top(), trianglePoint);
                HomogeneousPointVector view = transformPoint(V, model);
                HomogeneousPointVector projection = transformPoint(projectionMatrix, view);

                //output P’
                outStage1 << fixed << model.x << " " << model.y << " " << model.z << endl;

                outStage2 << fixed << view.x << " " << view.y << " " << view.z << endl;

                outStage3 << fixed << projection.x << " " << projection.y << " " << projection.z << endl;
            }

            outStage1 << endl;
            outStage2 << endl;
            outStage3 << endl;

        }
            //else if command = “translate”
        else if (command == "translate") {

            //input translation amounts
            double tx, ty, tz;
            scin >> tx >> ty >> tz;

            //generate the corresponding translation matrix T
            HomogeneousPointVector tem(tx, ty, tz);
            Matrix translationMatrix(tem);

            //S.push(product(S.top,T))
            smallStack.push(smallStack.top() * translationMatrix);
        }
            //else if command = “scale”
        else if (command == "scale") {

            //input scaling factors
            double sx, sy, sz;
            scin >> sx >> sy >> sz;

            //generate the corresponding scaling matrix T
            Matrix scalingMatrix(sx, sy, sz);

            //S.push(product(S.top,T))
            smallStack.push(smallStack.top() * scalingMatrix);
        }
            //else if command = “rotate”
        else if (command == "rotate") {

            //input rotation angle and axis
            double angle;
            HomogeneousPointVector axis;
            scin >> angle >> axis.x >> axis.y >> axis.z;

            //generate the corresponding rotation matrix T
            axis.makeUnitVector();


            HomogeneousPointVector i(1, 0, 0);
            HomogeneousPointVector j(0, 1, 0);
            HomogeneousPointVector k(0, 0, 1);

            HomogeneousPointVector c1 = rodriguesFormula(i, axis, angle);
            HomogeneousPointVector c2 = rodriguesFormula(j, axis, angle);
            HomogeneousPointVector c3 = rodriguesFormula(k, axis, angle);

            Matrix rotationMatrix;
            rotationMatrix.matrix = vector<vector<double>>{
                    {c1.x, c2.x, c3.x, 0},
                    {c1.y, c2.y, c3.y, 0},
                    {c1.z, c2.z, c3.z, 0},
                    {0,    0,    0,    1}
            };

            //S.push(product(S.top,T))
            smallStack.push((smallStack.top() * rotationMatrix));
        }
            //else if command = “push”
        else if (command == "push") {
            //do it yourself
            bigStack.push(smallStack);
        }
            //else if command = “pop”
        else if (command == "pop") {
            //do it yourself

            smallStack = bigStack.top();
            bigStack.pop();
        }
            //else if command = “end”
        else if (command == "end") {
            //break
            break;
        }
    }

    outStage1.close();
    outStage2.close();
    outStage3.close();

    return 0;
}