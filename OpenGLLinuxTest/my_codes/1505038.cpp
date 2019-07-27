/**
 * @author: Rumman (1505038)
 * @details: Mighty Ray Tracing
 * @todo: Ray Tracing.
 */

#include <bits/stdc++.h>
#include "bitmap_image.hpp"

#ifdef __linux__

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#else

#include <windows.h>
#include <gl.h>
#include <glu.h>
#include <glut.h>

#endif

#define pi acos(-1.0)
#define INF INT_MAX
#define NEG_INF (-100)
#define EPSILON 0.001

using namespace std;

typedef unsigned char color_t;

enum object_code {
    sphere,
    triangle,
    pyramid
};


object_code getObjectCode(string const &inString) {
    if (inString == "sphere") return sphere;
    if (inString == "pyramid") return pyramid;
    if (inString == "triangle") return triangle;
}


class Vector;

class Object;

class Ray;

class Colour;



// Vector class sometimes represents a point, sometimes a vector. Depending on the context.
// i.e. cameraPosition is a point but cameraLook, up, right are vectors.

class Vector {
public:
    double x, y, z, w;

    Vector() {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector(double X, double Y, double Z) {
        x = X;
        y = Y;
        z = Z;
    }

    Vector(const Vector &p) : x(p.x), y(p.y), z(p.z) {}

    double absoluteValue() {
        return sqrt(x * x + y * y + z * z);
    }

    void normalize() {
        double abs = absoluteValue();
        x /= abs;
        y /= abs;
        z /= abs;
    }

    Vector operator*(const Vector &vec) {
        Vector ret;
        ret.x = y * vec.z - z * vec.y;
        ret.y = -(x * vec.z - z * vec.x);
        ret.z = x * vec.y - y * vec.x;
        return ret;
    }

    double dot(const Vector &vec) {
        return (x * vec.x + y * vec.y + z * vec.z);
    }

    Vector operator*(const double &a) {
        Vector ret;
        ret.x = x * a;
        ret.y = y * a;
        ret.z = z * a;
        return ret;
    }

    Vector operator+(const Vector &v) {
        Vector ret;
        ret.x = x + v.x;
        ret.y = y + v.y;
        ret.z = z + v.z;
        return ret;
    }

    Vector operator-(const Vector &v) {
        Vector ret;
        ret.x = x - v.x;
        ret.y = y - v.y;
        ret.z = z - v.z;
        return ret;
    }

    void operator=(const Vector &v) {
        x = v.x;
        y = v.y;
        z = v.z;
    }
};

int howManyTimesCaptured = 0;
Vector cameraPos, cameraUp, cameraLook, cameraRight;
double cameraHeight;
double cameraAngle;
double cameraDistance = 10;
double aspectRatio = 1.0;
int drawaxes;
int drawgrid;
double angle;
GLdouble fovY = 90;
int reflectionLevel = 3;
vector<Vector> Lights;
vector<Object *> objectList;

int windowWidth = 500, windowHeight = 500;
int imageWidth = 768, imageHeight = 768;

void drawAxes() {
    if (drawaxes == 1) {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glColor3f(1, 0, 0);
            glVertex3f(100, 0, 0);
            glVertex3f(-100, 0, 0);

            glColor3f(0, 1, 0);
            glVertex3f(0, -100, 0);
            glVertex3f(0, 100, 0);

            glColor3f(0, 0, 1);
            glVertex3f(0, 0, 100);
            glVertex3f(0, 0, -100);
        }
        glEnd();
    }
}


void drawGrid() {
    int i;
    if (drawgrid == 1) {
        glColor3f(0.6, 0.6, 0.6);    //grey
        glBegin(GL_LINES);
        {
            for (i = -8; i <= 8; i++) {

                if (i == 0)
                    continue;    //SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i * 10, -100, 0);
                glVertex3f(i * 10, 100, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i * 10, 0);
                glVertex3f(90, i * 10, 0);
            }
        }
        glEnd();
    }
}


void drawSphere(double radius, int slices, int stacks) {
    struct point {
        double x, y, z;
    };

    struct point points[100][100];
    int i, j;
    double h, r;
    //generate points
    for (i = 0; i <= stacks; i++) {
        h = radius * sin(((double) i / (double) stacks) * (pi / 2));
        r = radius * cos(((double) i / (double) stacks) * (pi / 2));
        for (j = 0; j <= slices; j++) {
            points[i][j].x = r * cos(((double) j / (double) slices) * 2 * pi);
            points[i][j].y = r * sin(((double) j / (double) slices) * 2 * pi);
            points[i][j].z = h;
        }
    }
    //draw quads using generated points
    for (i = 0; i < stacks; i++) {
        for (j = 0; j < slices; j++) {
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(static_cast<GLfloat>(points[i][j].x), static_cast<GLfloat>(points[i][j].y),
                           static_cast<GLfloat>(points[i][j].z));
                glVertex3f(static_cast<GLfloat>(points[i][j + 1].x), static_cast<GLfloat>(points[i][j + 1].y),
                           static_cast<GLfloat>(points[i][j + 1].z));
                glVertex3f(static_cast<GLfloat>(points[i + 1][j + 1].x), static_cast<GLfloat>(points[i + 1][j + 1].y),
                           static_cast<GLfloat>(points[i + 1][j + 1].z));
                glVertex3f(static_cast<GLfloat>(points[i + 1][j].x), static_cast<GLfloat>(points[i + 1][j].y),
                           static_cast<GLfloat>(points[i + 1][j].z));

                //lower hemisphere
                glVertex3f(static_cast<GLfloat>(points[i][j].x), static_cast<GLfloat>(points[i][j].y),
                           static_cast<GLfloat>(-points[i][j].z));
                glVertex3f(static_cast<GLfloat>(points[i][j + 1].x), static_cast<GLfloat>(points[i][j + 1].y),
                           static_cast<GLfloat>(-points[i][j + 1].z));
                glVertex3f(static_cast<GLfloat>(points[i + 1][j + 1].x), static_cast<GLfloat>(points[i + 1][j + 1].y),
                           static_cast<GLfloat>(-points[i + 1][j + 1].z));
                glVertex3f(static_cast<GLfloat>(points[i + 1][j].x), static_cast<GLfloat>(points[i + 1][j].y),
                           static_cast<GLfloat>(-points[i + 1][j].z));
            }
            glEnd();
        }
    }
}


class Colour {
public:
    color_t r;
    color_t g;
    color_t b;

    Colour() {
        r = 0;
        g = 0;
        b = 0;
    }

    Colour(const color_t R, const color_t G, const color_t B) {
        r = R;
        g = G;
        b = B;
    }

    Colour operator*(const double &coefficient) {
        Colour result;
        result.r = static_cast<color_t>(r * coefficient);
        result.g = static_cast<color_t>(g * coefficient);
        result.b = static_cast<color_t>(b * coefficient);
        return result;
    }

    void operator*=(const double &coefficient) {
        r = static_cast<color_t>(r * coefficient);
        g = static_cast<color_t>(g * coefficient);
        b = static_cast<color_t>(b * coefficient);
    }


    void operator=(const Colour &colour) {
        r = colour.r;
        g = colour.g;
        b = colour.b;
    }

    Colour operator+(const Colour &colour) {
        Colour result;
        result.r = r + colour.r;
        result.g = g + colour.g;
        result.b = b + colour.b;
        return result;
    }

    void operator+=(const Colour &colour) {
        r = r + colour.r;
        g = g + colour.g;
        b = b + colour.b;
    }

    void setRGB(color_t rr, color_t gg, color_t bb) {
        r = rr;
        g = gg;
        b = bb;
    }

    void print() {
        cout << "r " << (int) r << " g " << (int) g << " b " << (int) b << endl;
    }
};


class Ray {
public:
    Vector startVector;
    Vector directionUnitVector;
    double t_value;
    Colour colour;

    Ray() {
        t_value = -1;
    }

    Ray(const Vector &start, const Vector &direction) {
        startVector = start;
        directionUnitVector = direction;
        t_value = -1;
    }
};

class Object {
public:
    Colour colour;
    int shininessValue;
    double ambientCoefficient, diffuseCoefficient, specularCoefficient, reflectionCoefficient;

    Object() {}

    virtual void drawObject() {}

    virtual double getIntersectionParameterT(Ray ray) {
        return -1;
    }

    Ray intersectAndIlluminate(Ray ray, int reflectionLevel) {
        Ray ret;
        Colour finalColor;
        Vector intersectionPoint;

        ret.t_value = getIntersectionParameterT(ray);
        intersectionPoint = ray.startVector + ray.directionUnitVector * ret.t_value;

        ret.colour = getColourAtAPoint(intersectionPoint) * ambientCoefficient;

        if (ret.t_value < 0) return ret;
        if (reflectionLevel < 1) return ret;

        finalColor = Illuminate(ray, intersectionPoint, ret.t_value, reflectionLevel);
        ret.colour = finalColor;

        return ret;
    }

    Colour Illuminate(Ray mainRay, Vector IntersectionPoint, double ParameterT, int reflectionLevel) {
        Colour resultColour = getColourAtAPoint(IntersectionPoint) * ambientCoefficient;

        Vector normalAtIntersectionPoint = getNormal(IntersectionPoint);

        double dotValue = normalAtIntersectionPoint.dot(mainRay.directionUnitVector);
        if (dotValue > 0) {
            normalAtIntersectionPoint = normalAtIntersectionPoint * (-1.0);
        }

        for (auto &Light : Lights) {
            Vector lightRayDirection = Light - IntersectionPoint; // not yet normalized
            double lightToThisObjectDistance = lightRayDirection.absoluteValue();
            lightRayDirection.normalize();

            Vector lightRayStart = IntersectionPoint + lightRayDirection * 1;
            Ray lightRayTowardsThisObject(lightRayStart, lightRayDirection);

            bool interceptedByAnotherObjectBefore = false;

            for (auto obj : objectList) {
                Ray interceptionResult = obj->intersectAndIlluminate(lightRayTowardsThisObject, 0);
                if (interceptionResult.t_value > 0 && interceptionResult.t_value < lightToThisObjectDistance) {
                    interceptedByAnotherObjectBefore = true;
                    break;
                }
            }

            if (!interceptedByAnotherObjectBefore) {
                double lightFactor = 1;
                double lambert = max(lightRayDirection.dot(normalAtIntersectionPoint), 0.0);
                Vector R = normalAtIntersectionPoint * 2.0 * lightRayDirection.dot(normalAtIntersectionPoint) -
                           lightRayDirection; // R = 2(L.N)N-L;
                double Phong = max(mainRay.directionUnitVector.dot(R), 0.0);
                resultColour += getColourAtAPoint(IntersectionPoint) * lightFactor * lambert * diffuseCoefficient;

                resultColour += Colour(255, 255, 255) * lightFactor * pow(Phong, shininessValue) *
                                specularCoefficient;
            }
        }
//        return  resultColour;
        if (reflectionLevel > 0) {
            Vector reflectionRayDirection = getReflectedVectorDirection(mainRay.directionUnitVector,
                                                                        normalAtIntersectionPoint);

            Vector reflectedRayStart = IntersectionPoint + reflectionRayDirection * 1;
            Ray reflectedRay(reflectedRayStart, reflectionRayDirection);

            int minimumObstacleIndex = NEG_INF;
            minimumObstacleIndex = getMinimumObstacleIdx(reflectedRay, minimumObstacleIndex);

            if (minimumObstacleIndex != NEG_INF) {
                Ray nextLevel = objectList[minimumObstacleIndex]->intersectAndIlluminate(reflectedRay,
                                                                                         reflectionLevel - 1);
                resultColour += nextLevel.colour * 1.0 * reflectionCoefficient;
            }
        }
        return resultColour;
    }

    int getMinimumObstacleIdx(const Ray &reflectedRay, int minimumObstacleIndex) {

        minimumObstacleIndex = NEG_INF;
        double minimumValueOfParameterT = INF;

        for (int i = 0; i < objectList.size(); i++) {
            Object *obj = objectList[i];

            // intersectAndIlluminate(reflectedRay, 0) -> 0 because just check if intersect or not
            Ray intersectOrNot = obj->intersectAndIlluminate(reflectedRay, 0);

            if (intersectOrNot.t_value > 0 && intersectOrNot.t_value < minimumValueOfParameterT) {
                minimumObstacleIndex = i;
                minimumValueOfParameterT = intersectOrNot.t_value;
            }
        }
        return minimumObstacleIndex;
    }


    virtual Vector getNormal(Vector intersectionPoint) {
        return Vector(0, 0, 1);
    }

    Vector getReflectedVectorDirection(Vector mainRay, Vector Normal) {
        Vector reflectionRay = mainRay - Normal * (2.0 * mainRay.dot(Normal));
        reflectionRay.normalize();
        return reflectionRay;
    }

    void setCoefficients(double amb, double diffuse, double spec, double reflect) {
        ambientCoefficient = amb;
        diffuseCoefficient = diffuse;
        specularCoefficient = spec;
        reflectionCoefficient = reflect;
    }

    void setShineValue(int shine) {
        shininessValue = shine;
    }

    void setColor(color_t r, color_t g, color_t b) {
        colour.setRGB(r, g, b);
    }

    virtual Colour getColourAtAPoint(Vector IntersectingPoint) {
        return colour;
    }
};


class Sphere : public Object {
public:
    Vector center;
    double radius;

    Sphere() {
        center = Vector(0, 0, 0);
        radius = 0.0;
    }

    Sphere(const Vector &Centre, double r) {
        center = Centre;
        radius = r;
    }

    void drawObject() override {
        glPushMatrix();
        {
            glColor3f((GLfloat) (colour.r / 255.0), (GLfloat) (colour.g / 255.0), (GLfloat) (colour.b / 255.0));
            glTranslatef((GLfloat) center.x, (GLfloat) center.y, (GLfloat) center.z);
            drawSphere(radius, 30, 30);
        }
        glPopMatrix();
    }

    double getIntersectionParameterT(Ray ray) override {
        double a, b, c, d, tem, t, t1, t2, r;
        Vector R0, Rd;

        r = radius;
        R0 = ray.startVector - center;
        Rd = ray.directionUnitVector;

        a = Rd.dot(Rd);
        b = 2.0 * Rd.dot(R0);
        c = R0.dot(R0) - r * r;

        tem = b * b - 4.0 * a * c;

        if (tem < 0) {
            t = -1.0;
        } else {
            d = sqrt(tem);
            t1 = (-b + d) / (2.0 * a);
            t2 = (-b - d) / (2.0 * a);

            if (t1 >= 0 && t2 >= 0) {
                t = min(t1, t2);
            } else {
                t = max(t1, t2);
            }
        }

        return t;
    }

    Vector getNormal(Vector intersectionPoint) override {
        Vector ret;
        ret = intersectionPoint - center;
        ret.normalize();
        return ret;
    }

    Colour getColourAtAPoint(Vector IntersectingPoint) override {
        return colour;
    }
};


class Triangle : public Object {
public:
    Vector A, B, C;

    Triangle() {
        A = Vector(0, 0, 0);
        B = Vector(0, 0, 0);
        C = Vector(0, 0, 0);
    }

    Triangle(const Vector &nA, const Vector &nB, const Vector &nC) {
        A = nA;
        B = nB;
        C = nC;
    }

    void drawObject() override {
        glPushMatrix();
        {
            glColor3f((GLfloat) (colour.r / 255.0), (GLfloat) (colour.g / 255.0), (GLfloat) (colour.b / 255.0));
            glBegin(GL_TRIANGLES);
            {
                glVertex3f((GLfloat) A.x, (GLfloat) A.y, (GLfloat) A.z);
                glVertex3f((GLfloat) B.x, (GLfloat) B.y, (GLfloat) B.z);
                glVertex3f((GLfloat) C.x, (GLfloat) C.y, (GLfloat) C.z);
            }
            glEnd();
        }
        glPopMatrix();
    }

    Colour getColourAtAPoint(Vector IntersectingPoint) override {
        return colour;
    }

    Vector getNormal(Vector intersectionPoint) override {
        Vector side1 = B - A;
        Vector side2 = C - A;
        Vector normal = side1 * side2;
        normal.normalize();
        return normal;
    }

    double getTriangleArea() {
        // parallelogram area = |AB * AC|, so triangle area = |AB*AC|/2.0 ;
        Vector side1 = B - A;
        Vector side2 = C - A;
        Vector crossVal = side1 * side2;
        return crossVal.absoluteValue() / 2.0;
    }

    double getDeterminant(vector<vector<double>> matrix) {
        double det, d1, d2, d3;
        d1 = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[2][1] * matrix[1][2]);
        d2 = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[2][0] * matrix[1][2]);
        d3 = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[2][0] * matrix[1][1]);
        return (d1 - d2 + d3);
    }


    double getBaryCentricT(Ray ray) {
        double Area = getTriangleArea();
        double beta, gamma, t;
        Vector o, d; // o -> origin and d -> direction

        o = ray.startVector;
        d = ray.directionUnitVector;

        beta = getDeterminant(vector<vector<double>>{
                {A.x - o.x, A.x - C.x, d.x},
                {A.y - o.y, A.y - C.y, d.y},
                {A.z - o.z, A.z - C.z, d.z}
        });

        gamma = getDeterminant(vector<vector<double>>{
                {A.x - B.x, A.x - o.x, d.x},
                {A.y - B.y, A.y - o.y, d.y},
                {A.z - B.z, A.z - o.z, d.z}
        });

        t = getDeterminant(vector<vector<double>>{
                {A.x - B.x, A.x - C.x, A.x - o.x},
                {A.y - B.y, A.y - C.y, A.y - o.y},
                {A.z - B.z, A.z - C.z, A.z - o.z}
        });

        beta /= Area;
        gamma /= Area;
        t /= Area;

        if (gamma < 0 || gamma > 1) return -1;
        if (beta < 0 || beta > 1 - gamma) return -1;
        if (t < 0) return -1;

        return t;
    }

    double getIntersectionParameterT(Ray ray) override {
        Vector side1, side2, normalOfTriangle, q, r;
        double determinant;
        double inverseDeterminant;
        double beta;
        double gamma;
        double t;

        side1 = B - A;
        side2 = C - A;

        normalOfTriangle = ray.directionUnitVector * side2;

        determinant = side1.dot(normalOfTriangle);

        if (determinant > (-EPSILON) && determinant < EPSILON) {
            return -1; // This ray is parallel to the triangle.
        }

        inverseDeterminant = 1 / determinant;

        r = ray.startVector - A;

        beta = r.dot(normalOfTriangle) * inverseDeterminant;

        if (beta < 0.0 || beta > 1.0) {
            return -1;
        }

        q = r * side1;
        gamma = ray.directionUnitVector.dot(q) * inverseDeterminant;

        if (gamma < 0.0 || (beta + gamma) > 1.0) {
            return -1;
        }
        t = side2.dot(q) * inverseDeterminant;

        if (t > EPSILON) { //ray intersection
            return t;
        }
        return -1;
    }
};


class Floor : public Object {
public:
    double totalHeight;
    double totalWidth;
    double tileHeight;
    double tileWidth;
    double initX, initY;

    Floor() {
        totalHeight = 500;
        totalWidth = 500;
        tileHeight = 3;
        tileWidth = 3;
        initX = -totalWidth / 2.0;
        initY = -totalHeight / 2.0;
    }

    Floor(double total_height, double total_width, double tile_height, double tile_width) {
        totalHeight = total_height;
        totalWidth = total_width;
        tileHeight = tile_height;
        tileWidth = tile_width;
        initX = -totalWidth / 2.0;
        initY = -totalHeight / 2.0;
        setShineValue(1);
    }

    double getIntersectionParameterT(Ray ray) override {
        return -(ray.startVector.z / ray.directionUnitVector.z);
    }

    Vector getNormal(Vector intersectionPoint) override {
        return Vector(0, 0, 1);
    }

    int isOdd(int n) {
        // 0 == even, 1 == odd
        return n % 2;
    }

    void setColourAtThisTile(int i, int j) {
        if ((i % 2 == 0 && j % 2 == 0) || (i % 2 == 1 && j % 2 == 1)) {
            glColor3f(1.0, 1.0, 1.0);
        } else {
            glColor3f(0, 0, 0);
        }
    }

    void drawFloor(double tempX, double tempY) {
        glPushMatrix();
        {
            glBegin(GL_QUADS);
            glVertex3f((GLfloat) tempX, (GLfloat) tempY, 0);
            glVertex3f((GLfloat) (tempX + tileWidth), (GLfloat) tempY, 0);
            glVertex3f((GLfloat) (tempX + tileWidth), (GLfloat) (tempY + tileHeight), 0);
            glVertex3f((GLfloat) tempX, (GLfloat) (tempY + tileHeight), 0);
            glEnd();
        }
        glPopMatrix();
    }

    void drawObject() override {
        auto tileInARow = static_cast<int>(totalWidth / tileWidth);
        auto tileInAColumn = static_cast<int>(totalHeight / tileHeight);
        tileInARow += isOdd(tileInARow); // making it an even number
        tileInAColumn += isOdd(tileInAColumn);

        for (int i = 0; i < tileInAColumn; i++) {
            for (int j = 0; j < tileInARow; j++) {
                double currX, currY;
                Colour currColour;
                currX = initX + i * tileWidth;
                currY = initY + j * tileHeight;
                setColourAtThisTile(i, j);

                drawFloor(currX, currY);
            }
        }
    }

    Colour getColourAtAPoint(Vector IntersectingPoint) override {
        auto dx = static_cast<int>((static_cast<double>(IntersectingPoint.x + totalWidth / 2.0)) / tileWidth);
        auto dy = static_cast<int>((static_cast<double>(IntersectingPoint.y + totalHeight / 2.0)) / tileHeight);

        int remI = dx % 2;
        int remJ = dy % 2;

        if (remI < 0) {
            remI *= -1;
        }
        if (remJ < 0) {
            remJ *= -1;
        }

        if ((remI == 0 && remJ == 0)) {
            return Colour(255, 255, 255);
        } else if ((remI == 1 && remJ == 1)) {
            return Colour(255, 255, 255);
        } else {
            return Colour(0, 0, 0);
        }
    }
};


void populatePixelsWithColour(const vector<vector<Colour>> &frame) {
    bitmap_image image(static_cast<const unsigned int>(imageWidth), static_cast<const unsigned int>(imageHeight));
    for (int i = 0; i < imageHeight; i++) {
        for (int j = 0; j < imageWidth; j++) {
            image.set_pixel(static_cast<const unsigned int>(j), static_cast<const unsigned int>(i),
                            frame[i][j].r, frame[i][j].g, frame[i][j].b);
        }
    }
    image.save_image("1505038.bmp");
}


Colour getPixelColour(const Ray &mainRay) {
    Colour resultColour(0, 0, 0);
    int closestObstacleIndex = NEG_INF;
    double minimumValueOfParameterT = INF;

    for (int i = 0; i < objectList.size(); i++) {
        Ray intersectOrNot = objectList[i]->intersectAndIlluminate(mainRay, 0);
        if (intersectOrNot.t_value > 0 && intersectOrNot.t_value < minimumValueOfParameterT) {
            closestObstacleIndex = i;
            minimumValueOfParameterT = intersectOrNot.t_value;
        }
    }
    if (closestObstacleIndex != NEG_INF) {
        Ray finalRayAfterAllLevelReflections = objectList[closestObstacleIndex]->intersectAndIlluminate(mainRay,
                                                                                                        reflectionLevel);
        resultColour = finalRayAfterAllLevelReflections.colour;
    }
    return resultColour;
}


void Capture() {
    vector<vector<Colour>> frameAtNearPlane; // it's correct. I cross chacked. Don't change.

    Vector currLeftCorner, eyeToPixelDirection, eyeToPixelRayStart;

    double planeDistanceFromCamera = (windowHeight / 2.0) / tan(fovY / 2.0 * (pi / 180.0));

    Vector dx, dy, dz;
    dy = cameraLook * planeDistanceFromCamera;
    dx = cameraRight * (windowWidth / 2.0);
    dz = cameraUp * (windowHeight / 2.0);

    Vector topLeftCornerOfWholeFrame = cameraPos + dy - dx + dz;

    double changeThroughRow = (windowWidth * 1.0) / imageWidth;
    double changeThroughColumn = (windowHeight * 1.0) / imageHeight;

    for (int i = 0; i < imageHeight; i++) {
        frameAtNearPlane.emplace_back(); // => frame.push_back(vector<Colour>());

        for (int j = 0; j < imageWidth; j++) {
            currLeftCorner = topLeftCornerOfWholeFrame + cameraRight * (j * changeThroughRow) -
                             cameraUp * (i * changeThroughColumn);

            eyeToPixelDirection = currLeftCorner - cameraPos;
            eyeToPixelDirection.normalize();
            eyeToPixelRayStart = cameraPos;

            Ray eyeToPixelRay(eyeToPixelRayStart, eyeToPixelDirection);

            frameAtNearPlane[i].push_back(getPixelColour(eyeToPixelRay));
        }
    }

    populatePixelsWithColour(frameAtNearPlane);
}


void setTestData() {
    windowWidth = 500;
    windowHeight = 500;
    cameraDistance = 25;
    aspectRatio = (windowWidth * 1.0) / windowHeight;
    fovY = 90;
    imageWidth = 1024;
    imageHeight = 1024;
    reflectionLevel = 10;

    Object *temp;

    Vector center(50, 35, 10);
    double radius = 10.0;
    temp = new Sphere(center, radius);
    temp->setColor(250, 250, 0);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);
    temp->setCoefficients(0, 0, 0, 1);
    temp->setShineValue(10);

    objectList.push_back(temp);
////////////////////////////////////////////////////
    center.x = 20;
    center.y = 20;
    center.z = 20;
    radius = 8.0;
    temp = new Sphere(center, radius);
    temp->setColor(0, 0, 250);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);
    temp->setCoefficients(0, 0, 0, 1);

    temp->setShineValue(15);

    objectList.push_back(temp);
/////////////////////////////////////////////////////
    center.x = -5;
    center.y = 30;
    center.z = 40;
    radius = 3.0;
    temp = new Sphere(center, radius);
    temp->setColor(0, 250, 250);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);
    temp->setCoefficients(0, 0, 0, 1);

    temp->setShineValue(20);
    objectList.push_back(temp);
/////////////////////////////////////////////////////

    center.x = -5;
    center.y = 35;
    center.z = 26;
    radius = 5;
    temp = new Sphere(center, radius);
    temp->setColor(0, 250, 0);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);
    temp->setCoefficients(0.3, 0.1, 0.2, 0.4);

    temp->setShineValue(40);
    objectList.push_back(temp);
/////////////////////////////////////////////////////

    center.x = -25;
    center.y = -35;
    center.z = 26;
    radius = 5;
    temp = new Sphere(center, radius);
    temp->setColor(0, 0, 250);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);
    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);

    temp->setShineValue(40);
    objectList.push_back(temp);
///////////////////////////////////////////////////////////////

    temp = new Triangle(Vector(-50, 0, 0), Vector(0, 0, 50), Vector(50, 0, 0));
    temp->setColor(200, 50, 30);
    temp->setCoefficients(0.3, 0.1, 0.2, 0.4);

    temp->setShineValue(1);
    objectList.push_back(temp);
///////////////////////////////////////////////////////////////

//    temp = new Triangle(Vector(-500, 500, 0), Vector(500, 500, 0), Vector(-500, -500, 0));
//    temp->setColor(255, 255, 255);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);
//
//    temp->setShineValue(1);
//
//    objectList.push_back(temp);
////////////////////////////////////////////////////////////////

//    temp = new Triangle(Vector(500, 500, 0), Vector(-500, -500, 0), Vector(500, -500, 0));
//    temp->setColor(255, 255, 255);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.2);
//
//    temp->setShineValue(1);
//    objectList.push_back(temp);

////////////////////////////////////////////////////////////////

    temp = new Floor(200, 200, 10, 10);
    temp->setCoefficients(0.3, 0.2, 0.1, 0.4);
    temp->setShineValue(1);
    objectList.push_back(temp);


    Vector light1(70, 70, 70);
//    Vector light2(100, 0, 0);
//    Vector light3(-70, 70, 70);
    Lights.push_back(light1);
//    Lights.push_back(light2);
//    Lights.push_back(light3);
}

// T1 = vector or point , T2 = Colour

template<typename T1, typename T2>
void generateTriangle(T1 A, T1 B, T1 C, T2 color, double a, double d, double s, double r, double shine) {
    Object *obj = new Triangle(A, B, C);
    obj->setColor(color.r, color.g, color.b);
    obj->setCoefficients(a, d, s, r);
    obj->setShineValue(static_cast<int>(shine));
    objectList.push_back(obj);
}

void inputDataFromFile() {
    windowWidth = 500;
    windowHeight = 500;
    cameraDistance = 25;
    aspectRatio = (windowWidth * 1.0) / windowHeight;
    fovY = 90;
    imageWidth = 1024;
    imageHeight = 1024;
    reflectionLevel = 10;

    int numOfObjects;

    Object *obj;

    ifstream input;
    input.open("description.txt");

    input >> reflectionLevel;

    input >> imageWidth;
    imageHeight = imageWidth;

    input >> numOfObjects;

    for (int i = 0; i < numOfObjects; i++) {
        string objectName;
        input >> objectName;

        switch (getObjectCode(objectName)) {
            case sphere: {
                double r, g, b;
                double x, y, z;
                double radius;
                double a, d, s, ref, shine;
                input >> x >> y >> z;
                input >> radius;
                input >> r >> g >> b;
                input >> a >> d >> s >> ref;
                input >> shine;
                obj = new Sphere(Vector(x, y, z), radius);
                obj->setColor(static_cast<color_t>(r * 255.0), static_cast<color_t>(g * 255.0),
                              static_cast<color_t>(b * 255.0));
                obj->setCoefficients(a, d, s, ref);
                obj->setShineValue(static_cast<int>(shine));
                objectList.push_back(obj);

                break;
            }
            case pyramid: {
                double ox, oy, oz;
                double d, h;
                double r, g, b;
                double a, diff, s, ref;
                double shine;

                input >> ox >> oy >> oz;
                input >> d >> h;
                input >> r >> g >> b;
                input >> a >> diff >> s >> ref;
                input >> shine;

                // 1
                generateTriangle(
                        Vector(ox - d / 2.0, oy + d / 2.0, oz + 0),
                        Vector(ox + d / 2.0, oy + d / 2.0, oz + 0),
                        Vector(ox - d / 2.0, oy - d / 2.0, oz + 0),
                        Colour(static_cast<const color_t>(r * 255.0),
                               static_cast<const color_t>(g * 255.0),
                               static_cast<const color_t>(b * 255.0)),
                        a, diff, s, ref, shine
                );

                // 2
                generateTriangle(
                        Vector(ox - d / 2.0, oy - d / 2.0, oz + 0),
                        Vector(ox + d / 2.0, oy - d / 2.0, oz + 0),
                        Vector(ox + d / 2.0, oy + d / 2.0, oz + 0),
                        Colour(static_cast<const color_t>(r * 255.0),
                               static_cast<const color_t>(g * 255.0),
                               static_cast<const color_t>(b * 255.0)),
                        a, diff, s, ref, shine
                );

                // 3
                generateTriangle(
                        Vector(ox - d / 2.0, oy - d / 2.0, oz + 0),
                        Vector(ox - d / 2.0, oy + d / 2.0, oz + 0),
                        Vector(ox + 0, oy + 0, oz + h),
                        Colour(static_cast<const color_t>(r * 255.0),
                               static_cast<const color_t>(g * 255.0),
                               static_cast<const color_t>(b * 255.0)),
                        a, diff, s, ref, shine
                );

                // 4
                generateTriangle(
                        Vector(ox + d / 2.0, oy + d / 2.0, oz + 0),
                        Vector(ox - d / 2.0, oy + d / 2.0, oz + 0),
                        Vector(ox + 0, oy + 0, oz + h),
                        Colour(static_cast<const color_t>(r * 255.0),
                               static_cast<const color_t>(g * 255.0),
                               static_cast<const color_t>(b * 255.0)),
                        a, diff, s, ref, shine
                );

                // 5
                generateTriangle(
                        Vector(ox + d / 2.0, oy + d / 2.0, oz + 0),
                        Vector(ox + d / 2.0, oy - d / 2.0, oz + 0),
                        Vector(ox + 0, oy + 0, oz + h),
                        Colour(static_cast<const color_t>(r * 255.0),
                               static_cast<const color_t>(g * 255.0),
                               static_cast<const color_t>(b * 255.0)),
                        a, diff, s, ref, shine
                );

                // 6
                generateTriangle(
                        Vector(ox - d / 2.0, oy - d / 2.0, oz + 0),
                        Vector(ox + d / 2.0, oy - d / 2.0, oz + 0),
                        Vector(ox + 0, oy + 0, oz + h),
                        Colour(static_cast<const color_t>(r * 255.0),
                               static_cast<const color_t>(g * 255.0),
                               static_cast<const color_t>(b * 255.0)),
                        a, diff, s, ref, shine
                );

                break;
            }
            case triangle:
                break;
            default:
                break;
        }
    }

    obj = new Floor(1000, 1000, 10, 10);
    obj->setCoefficients(0.4, 0.1, 0.2, 0.3);
    obj->setShineValue(1);
    objectList.push_back(obj);

    int numOfLights;
    input >> numOfLights;
    for (int i = 0; i < numOfLights; i++) {
        double lx, ly, lz;
        input >> lx >> ly >> lz;
        Vector light(lx, ly, lz);
        Lights.push_back(light);
    }

    input.close();
    cout << "Data Load Done. Press 0 to capture." << endl;

}

void keyboardListener(unsigned char key, int x, int y) {
    Vector l = cameraLook, r = cameraRight, u = cameraUp;
    double ang = .05;

    switch (key) {

        case '1':   //rotate right. rotate r and l WRT u;
            ang *= -1;
            cameraRight.x = r.x * cos(ang) + l.x * sin(ang);
            cameraRight.y = r.y * cos(ang) + l.y * sin(ang);
            cameraRight.z = r.z * cos(ang) + l.z * sin(ang);
            cameraRight.normalize();

            cameraLook = cameraUp * cameraRight;
            cameraLook.normalize();

//                cameraUp.print();
            break;
        case '2':
            cameraRight.x = r.x * cos(ang) + l.x * sin(ang);
            cameraRight.y = r.y * cos(ang) + l.y * sin(ang);
            cameraRight.z = r.z * cos(ang) + l.z * sin(ang);
            cameraRight.normalize();

            cameraLook = cameraUp * cameraRight;
            cameraLook.normalize();
            break;
        case '3':       //rotate up. rotate l and u vectors WRT r.
            cameraLook.x = l.x * cos(ang) + u.x * sin(ang);
            cameraLook.y = l.y * cos(ang) + u.y * sin(ang);
            cameraLook.z = l.z * cos(ang) + u.z * sin(ang);
            cameraLook.normalize();

            cameraUp = cameraRight * cameraLook;
            cameraUp.normalize();
            break;
        case '4':       //rotate down. rotate l and u vectors WRT r.
            ang *= -1;
            cameraLook.x = l.x * cos(ang) + u.x * sin(ang);
            cameraLook.y = l.y * cos(ang) + u.y * sin(ang);
            cameraLook.z = l.z * cos(ang) + u.z * sin(ang);
            cameraLook.normalize();

            cameraUp = cameraRight * cameraLook;
            cameraUp.normalize();
            break;

        case '5':       //tilt camera clockwise.rotate r and u WRT l.
            cameraUp.x = u.x * cos(ang) + r.x * sin(ang);
            cameraUp.y = u.y * cos(ang) + r.y * sin(ang);
            cameraUp.z = u.z * cos(ang) + r.z * sin(ang);
            cameraUp.normalize();

            cameraRight = cameraLook * cameraUp;
            cameraRight.normalize();
            break;
        case '6':
            ang *= -1;
            cameraUp.x = u.x * cos(ang) + r.x * sin(ang);
            cameraUp.y = u.y * cos(ang) + r.y * sin(ang);
            cameraUp.z = u.z * cos(ang) + r.z * sin(ang);
            cameraUp.normalize();

            cameraRight = cameraLook * cameraUp;
            cameraRight.normalize();
            break;

        case '0':
            howManyTimesCaptured++;
            Capture();
            cout << "Image Captured. Capture No. " << howManyTimesCaptured << endl;
            break;

        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    double mov = 3;
    switch (key) {
        case GLUT_KEY_PAGE_DOWN:        //down arrow key
            cameraPos.x -= mov * cameraUp.x;
            cameraPos.y -= mov * cameraUp.y;
            cameraPos.z -= mov * cameraUp.z;
            break;
        case GLUT_KEY_PAGE_UP:        // up arrow key
            cameraPos.x += mov * cameraUp.x;
            cameraPos.y += mov * cameraUp.y;
            cameraPos.z += mov * cameraUp.z;
            break;

        case GLUT_KEY_RIGHT:
            cameraPos.x += mov * cameraRight.x;
            cameraPos.y += mov * cameraRight.y;
            cameraPos.z += mov * cameraRight.z;
            break;
        case GLUT_KEY_LEFT:
            cameraPos.x -= mov * cameraRight.x;
            cameraPos.y -= mov * cameraRight.y;
            cameraPos.z -= mov * cameraRight.z;
            break;

        case GLUT_KEY_UP:
            cameraPos.x += mov * cameraLook.x; // shamne piche.
            cameraPos.y += mov * cameraLook.y;
            cameraPos.z += mov * cameraLook.z;
            break;
        case GLUT_KEY_DOWN:
            cameraPos.x -= mov * cameraLook.x;
            cameraPos.y -= mov * cameraLook.y;
            cameraPos.z -= mov * cameraLook.z;
            break;

        default:
            break;
    }
}


void mouseListener(int button, int state, int x, int y)    //x, y is the x-y of the screen (2D)
{
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if (state == GLUT_DOWN)        // 2 times?? in ONE click? -- solution is checking DOWN or UP
            {
                drawaxes = 1 - drawaxes;
            }
            break;

        case GLUT_RIGHT_BUTTON:
            //........
            break;

        case GLUT_MIDDLE_BUTTON:
            //........
            break;

        default:
            break;
    }
}


void display() {

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);    //color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /*********************
    # set-up camera here #
    *********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    gluLookAt(cameraPos.x, cameraPos.y, cameraPos.z,
              cameraPos.x + cameraDistance * cameraLook.x, cameraPos.y + cameraDistance * cameraLook.y,
              cameraPos.z + cameraDistance * cameraLook.z,
              cameraUp.x, cameraUp.y, cameraUp.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/

    glPushMatrix();
    {
        //add objects

        drawAxes();
        drawGrid();
        for (auto &i : objectList) {
            i->drawObject();
        }
        glColor3f(1, 1, 1);
        for (auto &Light : Lights) {
            glPushMatrix();
            {
                glTranslatef((GLfloat) Light.x, (GLfloat) Light.y, (GLfloat) Light.z);
                glutSolidSphere(0.5, 7, 7);
            }
            glPopMatrix();
        }
    }
    glPopMatrix();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate() {
    angle += 0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init() {
    //codes for initialization
    drawgrid = 0;
    drawaxes = 1;
    cameraHeight = 150.0;
    cameraAngle = 1.0;
    angle = 0;

    /** Camera initialization **/
    cameraPos = Vector(-150, -150, 0);
    cameraLook = Vector(1 / sqrt(2.0), 1 / sqrt(2.0), 0);
    cameraUp = Vector(0, 0, 1);
    cameraRight = Vector(1 / sqrt(2.0), -1 / sqrt(2.0), 0);
    drawgrid = 0;
    drawaxes = 1;
    cameraHeight = 80;
    cameraAngle = pi / 4;;

    //clear the screen
    glClearColor(0, .5, .5, 0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(fovY, 1, 1, 1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}


int main(int argc, char **argv) {
    glutInit(&argc, argv);

//    setTestData();
    inputDataFromFile();
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(300, 100);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);    //Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);    //enable Depth Testing

    glutDisplayFunc(display);    //display callback function
    glutIdleFunc(animate);        //what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();        //The main loop of OpenGL

    return 0;
}