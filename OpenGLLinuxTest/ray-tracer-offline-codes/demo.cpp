#include <iostream>
#include <fstream>

#ifdef __linux__

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#else

#include "windows.h"
#include <gl.h>
#include <glu.h>
#include <glut.h>

#endif

#include <cmath>
#include <vector>
#include "bitmap_image.hpp"

#define PI acos(-1.0)

using namespace std;

class Vector;

class Object;

class Camera;

int drawaxes = 1;
int drawgrid = 1;

int windowWidth = 500;
int windowHeight = 500;
double cameraDistance = 25;
GLdouble fovY = 90;
GLdouble aspectRatio = 1.0;
GLdouble zNear = 1;
GLdouble zFar = 1000.0;
unsigned int imageWidth = 768;
unsigned int imageHeight = 768;
int rayTraceLevel = 10;
int refraction = 1;

vector<Vector> lights;
vector<Object *> objects;

class Vector {
public:
    double x;
    double y;
    double z;

    Vector() {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector(double Nx, double Ny, double Nz) {
        x = Nx;
        y = Ny;
        z = Nz;
    }

    void operator=(const Vector &a) {
        x = a.x;
        y = a.y;
        z = a.z;
    }

    Vector operator+(const Vector &a) {
        Vector result;
        result.x = x + a.x;
        result.y = y + a.y;
        result.z = z + a.z;
        return result;
    }

    friend Vector operator+(double scalar, const Vector &a) {
        Vector result;
        result.x = a.x + scalar;
        result.y = a.y + scalar;
        result.z = a.z + scalar;
        return result;
    }

    friend Vector operator+(const Vector &a, double scalar) {
        Vector result;
        result.x = a.x + scalar;
        result.y = a.y + scalar;
        result.z = a.z + scalar;
        return result;
    }

    Vector operator-(const Vector &a) {
        Vector result;
        result.x = x - a.x;
        result.y = y - a.y;
        result.z = z - a.z;
        return result;
    }

    friend Vector operator-(double scalar, const Vector &a) {
        Vector result;
        result.x = scalar - a.x;
        result.y = scalar - a.y;
        result.z = scalar - a.z;
        return result;
    }

    friend Vector operator-(const Vector &a, double scalar) {
        Vector result;
        result.x = a.x - scalar;
        result.y = a.y - scalar;
        result.z = a.z - scalar;
        return result;
    }

    double dotMultiply(const Vector &a) {
        double result;
        result = x * a.x + y * a.y + z * a.z;
        return result;
    }

    Vector operator*(const Vector &a) {
        Vector result;
        result.x = y * a.z - z * a.y;
        result.y = z * a.x - x * a.z;
        result.z = x * a.y - y * a.x;
        return result;
    }

    friend Vector operator*(const Vector &a, double scalar) {
        Vector result;
        result.x = a.x * scalar;
        result.y = a.y * scalar;
        result.z = a.z * scalar;
        return result;
    }

    friend Vector operator*(double scalar, const Vector &a) {
        Vector result;
        result.x = a.x * scalar;
        result.y = a.y * scalar;
        result.z = a.z * scalar;
        return result;
    }

    void operator+=(const Vector &a) {
        x = x + a.x;
        y = y + a.y;
        z = z + a.z;
    }

    void operator-=(const Vector &a) {
        x = x - a.x;
        y = y - a.y;
        z = z - a.z;
    }

    void operator*=(const Vector &a) {
        Vector result;

        result.x = y * a.z - z * a.y;
        result.y = z * a.x - x * a.z;
        result.z = x * a.y - y * a.x;

        x = result.x;
        y = result.y;
        z = result.z;
    }

    void operator*=(double scalar) {
        x = x * scalar;
        y = y * scalar;
        z = z * scalar;
    }

    void normalize() {
        double absoluteValue = sqrt(x * x + y * y + z * z);
        x = x / absoluteValue;
        y = y / absoluteValue;
        z = z / absoluteValue;
    }

    void print() {
        cout << x << " " << y << " " << z << endl;
    }
};

class Colour {
public:
    unsigned char r;
    unsigned char g;
    unsigned char b;

    Colour() {
        r = 0;
        g = 0;
        b = 0;
    }

    Colour(unsigned char inputR, unsigned char inputG, unsigned char inputB) {
        r = inputR;
        g = inputG;
        b = inputB;
    }

    friend Colour operator*(const Colour &colour, double coefficient) {
        Colour result;
        result.r = (unsigned char) (colour.r * coefficient);
        result.g = (unsigned char) (colour.g * coefficient);
        result.b = (unsigned char) (colour.b * coefficient);
        return result;
    }

    friend Colour operator*(double coefficient, const Colour &colour) {
        Colour result;
        result.r = (unsigned char) (colour.r * coefficient);
        result.g = (unsigned char) (colour.g * coefficient);
        result.b = (unsigned char) (colour.b * coefficient);
        return result;
    }

    void operator=(const Colour &colour) {
        r = colour.r;
        g = colour.g;
        b = colour.b;
    }

    void operator+=(const Colour &colour) {
        r = r + colour.r;
        g = g + colour.g;
        b = b + colour.b;
    }

    void print() {
        cout << "r " << (int) r << " g " << (int) g << " b " << (int) b << endl;
    }

};

class Ray {
public:
    Vector start;
    Vector direction;
    double t;
    Colour colour;

    Ray() {
        t = -1;
    }

    Ray(Vector startInput, Vector directionInput) {
        start = startInput;
        direction = directionInput;
        t = -1;
    }
};

class Object {
public:
    Vector referencePoint;
    int shine;
    Colour colour;
    double coefficients[4];
    double length, width, height;
    double refractiveCoefficient, refractiveIndex;

    Object() {
        length = 0;
        width = 0;
        height = 0;
        shine = 0;
    }

    virtual void draw() {}

    virtual Colour getColourAt(Vector point) {
        return colour;
    }

    virtual double getIntersectionT(Ray ray) {
        return -1;
    }

    virtual Vector getNormal(Vector intersectionPoint) {
        Vector normal(1, 0, 0);
        return normal;
    }

    Vector getReflection(Vector incidentRay, Vector normal) {
        Vector reflectedRay = incidentRay - 2 * incidentRay.dotMultiply(normal) * normal;
        reflectedRay.normalize();
        return reflectedRay;
    }

    Ray getRefraction(Vector incidentRay, Vector normal, int cycle) {
        Ray refractedRayFinal;
        Vector refractedRay;
        double nDotI, k, actualRefractiveIndex;
        if (cycle) {
            actualRefractiveIndex = 1.0 / refractiveIndex;
        } else {
            actualRefractiveIndex = refractiveIndex;
        }

        nDotI = normal.dotMultiply(incidentRay);
        k = 1.0 - actualRefractiveIndex * actualRefractiveIndex * (1.0 - nDotI * nDotI);
        if (k < 0.0) {
            refractedRayFinal.t = -1;
            return refractedRayFinal;
        } else {
            refractedRay = actualRefractiveIndex * incidentRay - (actualRefractiveIndex * nDotI + sqrt(k)) * normal;
            refractedRay.normalize();
            refractedRayFinal.direction = refractedRay;
            refractedRayFinal.t = 1;
            return refractedRayFinal;
        }
    }

    Ray intersect(Ray ray, int level) {
        int i, j;
        Ray result;
        Vector intersectionPoint;

        result.t = getIntersectionT(ray);
        intersectionPoint = ray.start + ray.direction * result.t;
        result.colour = getColourAt(intersectionPoint) * coefficients[0];

        if (result.t < 0 || level < 1) {
            return result;
        }


        Vector normal = getNormal(intersectionPoint);
        if (normal.dotMultiply(ray.direction) > 0) {
            normal = Vector(0, 0, 0) - normal;
        }
        for (i = 0; i < lights.size(); ++i) {
            bool flag = true;
            Vector direction = lights[i] - intersectionPoint;
            double maximumT = sqrt(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z);
            direction.normalize();

            Vector start = intersectionPoint + direction * 1;
            Ray shineRay(start, direction);

            for (j = 0; j < objects.size(); ++j) {
                Ray obscureResult = objects[j]->intersect(shineRay, 0); // intersect can be bool or double
                if (obscureResult.t > 0) {
                    if (obscureResult.t < maximumT) {
                        flag = false;
                        break;
                    }
                }
            }

            if (flag) {
                double intensity = maximumT * maximumT;
                double lightFactor = intensity / (maximumT * maximumT);
                double lambert = max(direction.dotMultiply(normal), 0.0);
                Vector r = 2 * direction.dotMultiply(normal) * normal - direction;
                double phong = max(ray.direction.dotMultiply(r), 0.0);
                result.colour += lightFactor * lambert * coefficients[1] * getColourAt(intersectionPoint);
                result.colour += lightFactor * pow(phong, shine) * coefficients[2] * getColourAt(intersectionPoint);
            }
        }
        if (level > 0) {
            Vector reflection = getReflection(ray.direction, normal);
            Vector start = intersectionPoint + reflection * 1;
            Ray reflectedRay(start, reflection);

            int nearest = -1;

            double minimumT = INFINITY;

            for (i = 0; i < objects.size(); ++i) {
                Ray checkResult = objects[i]->intersect(reflectedRay, 0);
                if (checkResult.t <= 0) {
                    continue;
                }
                if (checkResult.t < minimumT) {
                    nearest = i;
                    minimumT = result.t;
//                    minimumT = checkResult.t; // rumman
                }
            }
            if (nearest > -1) {
                Ray resultNext = objects[nearest]->intersect(reflectedRay, level - 1);
                result.colour += resultNext.colour * coefficients[3];
            }
        }

        return result;
    }

    void setReferencePoint(const Vector &a) {
        referencePoint = a;
    }

    void setShine(int shineInput) {
        shine = shineInput;
    }

    void setColour(unsigned char r, unsigned char g, unsigned char b) {
        colour.r = r;
        colour.g = g;
        colour.b = b;
    }

    void setCoefficients(double ambient, double diffuse, double specular, double reflection) {
        coefficients[0] = ambient;
        coefficients[1] = diffuse;
        coefficients[2] = specular;
        coefficients[3] = reflection;
    }

    void setDimensions(double lengthInput, double widthInput, double heightInput) {
        length = lengthInput;
        width = widthInput;
        height = heightInput;
    }

    void setRefraction(double coefficient, double index) {
        refractiveCoefficient = coefficient;
        refractiveIndex = index;
    }
};

class Sphere : public Object {
    double radius;
public:
    Sphere(Vector center, double radiusInput) {
        setReferencePoint(center);
        radius = radiusInput;
    }

    void draw() override {
        glPushMatrix();
        {
            glColor3f((GLfloat) (colour.r / 255.0), (GLfloat) (colour.g / 255.0), (GLfloat) (colour.b / 255.0));
            glTranslatef((GLfloat) referencePoint.x, (GLfloat) referencePoint.y, (GLfloat) referencePoint.z);
            glutSolidSphere(radius, 100, 100);
        }
        glPopMatrix();
    }

    Colour getColourAt(Vector point) override {
        return colour;
    }

    double getIntersectionT(Ray ray) override {
        double a, b, c, d, t, t1, t2;
        a = ray.direction.dotMultiply(ray.direction);
        b = 2 * (ray.direction.dotMultiply(ray.start - referencePoint));
        c = (ray.start - referencePoint).dotMultiply(ray.start - referencePoint) - radius * radius;
        d = b * b - 4 * a * c;
        if (d >= 0) {
            d = sqrt(d);
            t1 = (-b + d) / (2 * a);
            t2 = (-b - d) / (2 * a);
            if (t1 < 0 || t2 < 0) {
                t = max(t1, t2);
            } else {
                t = min(t1, t2);
            }
        } else {
            t = -1;
        }
        return t;
    }

    Vector getNormal(Vector intersectionPoint) override {
        Vector result;
        result = intersectionPoint - referencePoint;
        result.normalize();
        return result;
    }
};

class Floor : public Object {
    double height, width, tileHeight, tileWidth;
public:
    Floor(double floorHeight, double floorWidth, double tileHeightInput, double tileWidthInput) {
        height = floorHeight;
        width = floorWidth;
        setShine(1);
        tileHeight = tileHeightInput;
        tileWidth = tileWidthInput;
    }

    void draw() override {
        int i, j, numberOfI, numberOfJ;
        double x = -width / 2, y = -height / 2, tempX, tempY;
        numberOfI = (int) (height / tileHeight);
        numberOfJ = (int) (width / tileWidth);

        if (numberOfI % 2) {
            numberOfI++;
        }
        if (numberOfJ % 2) {
            numberOfJ++;
        }

        for (i = 0; i < numberOfI; i++) {
            for (j = 0; j < numberOfJ; j++) {
                tempX = x + i * tileWidth;
                tempY = y + j * tileHeight;
                if ((i % 2 == 0 && j % 2 == 0) || (i % 2 == 1 && j % 2 == 1)) {
                    glColor3f(1, 1, 1);
                } else {
                    glColor3f(0, 0, 0);
                }
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
        }
    }

    Colour getColourAt(Vector point) override {
        double xDistance = point.x + width / 2.0;
        double yDistance = point.y + height / 2.0;
        int i = (int) (xDistance / tileWidth);
        int j = (int) (yDistance / tileHeight);
        int remainderI = i % 2;
        int remainderJ = j % 2;
        if (remainderI < 0) {
            remainderI *= -1;
        }
        if (remainderJ < 0) {
            remainderJ *= -1;
        }
        if ((remainderI == 0 && remainderJ == 0) || (remainderI == 1 && remainderJ == 1)) {
            return Colour(255, 255, 255);
        } else {
            return Colour(0, 0, 0);
        }
    }

    double getIntersectionT(Ray ray) override {
        double t = -(ray.start.z / ray.direction.z);
        return t;
    }

    Vector getNormal(Vector intersectionPoint) override {
        Vector result(0, 0, 1);
        return result;
    }
};

class Triangle : public Object {
    Vector a;
    Vector b;
    Vector c;
public:
    Triangle(Vector p, Vector q, Vector r) {
        a = p;
        b = q;
        c = r;
    }

    void draw() override {
        glPushMatrix();
        {
            glColor3f((GLfloat) (colour.r / 255.0), (GLfloat) (colour.g / 255.0), (GLfloat) (colour.b / 255.0));
            glBegin(GL_TRIANGLES);
            {
                glVertex3f((GLfloat) a.x, (GLfloat) a.y, (GLfloat) a.z);
                glVertex3f((GLfloat) b.x, (GLfloat) b.y, (GLfloat) b.z);
                glVertex3f((GLfloat) c.x, (GLfloat) c.y, (GLfloat) c.z);
            }
            glEnd();
        }
        glPopMatrix();
    }

    Colour getColourAt(Vector point) override {
        return colour;
    }

    double getIntersectionT(Ray ray) override {
        double epsilon = 0.001, determinant, inverseDeterminant, u, v, t;
        Vector edge1, edge2, p, q, r;

        //Find vectors for two edges sharing V1
        edge1 = b - a;
        edge2 = c - a;
        //Begin calculating determinant - also used to calculate u parameter
        p = ray.direction * edge2;
        //if determinant is near zero, ray lies in plane of triangle or ray is parallel to plane of triangle
        determinant = edge1.dotMultiply(p);
        //NOT CULLING
        if (determinant > (-epsilon) && determinant < epsilon) {
            return -1;
        }
        inverseDeterminant = 1 / determinant;
        //calculate distance from V1 to ray origin
        r = ray.start - a;
        //Calculate u parameter and test bound
        u = r.dotMultiply(p) * inverseDeterminant;
        //The intersection lies outside of the triangle
        if (u < 0.0 || u > 1.0) {
            return -1;
        }
        //Prepare to test v parameter
        q = r * edge1;
        //Calculate V parameter and test bound
        v = ray.direction.dotMultiply(q) * inverseDeterminant;
        //The intersection lies outside of the triangle
        if (v < 0.0 || (u + v) > 1.0) {
            return -1;
        }
        t = edge2.dotMultiply(q) * inverseDeterminant;
        if (t > epsilon) { //ray intersection
            return t;
        }
        // No hit, no win
        return -1;
    }

    Vector getNormal(Vector intersectionPoint) override {
        Vector edge1 = b - a;
        Vector edge2 = c - a;
        Vector normal = edge1 * edge2;
        normal.normalize();
        return normal;
    }
};

class GeneralQuadratic : public Object {
    double A, B, C, D, E, F, G, H, I, J;
public:
    GeneralQuadratic(double a, double b, double c, double d, double e, double f, double g, double h, double i,
                     double j) {
        A = a;
        B = b;
        C = c;
        D = d;
        E = e;
        F = f;
        G = g;
        H = h;
        I = i;
        J = j;
    }

    void draw() override {}

    Colour getColourAt(Vector point) override {
        return colour;
    }

    double getIntersectionT(Ray ray) override {
        double a, b, c, d, t, t1, t2;

        a = A * ray.direction.x * ray.direction.x + B * ray.direction.y * ray.direction.y +
            C * ray.direction.z * ray.direction.z +
            D * ray.direction.x * ray.direction.y + E * ray.direction.y * ray.direction.z +
            F * ray.direction.z * ray.direction.x;

        b = A * 2 * ray.start.x * ray.direction.x + B * 2 * ray.start.y * ray.direction.y +
            C * 2 * ray.start.z * ray.direction.z +
            D * (ray.start.x * ray.direction.y + ray.start.y * ray.direction.x) +
            E * (ray.start.y * ray.direction.z + ray.start.z * ray.direction.y) +
            F * (ray.start.z * ray.direction.x + ray.start.x * ray.direction.z) +
            G * ray.direction.x + H * ray.direction.y + I * ray.direction.z;

        c = A * ray.start.x * ray.start.x + B * ray.start.y * ray.start.y + C * ray.start.z * ray.start.z +
            D * ray.start.x * ray.start.y + E * ray.start.y * ray.start.z + F * ray.start.z * ray.start.x +
            G * ray.start.x + H * ray.start.y + I * ray.start.z + J;
        d = b * b - 4 * a * c;
        if (d >= 0) {
            d = sqrt(d);
            t1 = (-b + d) / (2 * a);
            t2 = (-b - d) / (2 * a);

            Vector intersectionPoint;
            double xDistance, yDistance, zDistance;

            intersectionPoint = ray.start + t1 * ray.direction;
            xDistance = intersectionPoint.x - referencePoint.x;
            yDistance = intersectionPoint.y - referencePoint.y;
            zDistance = intersectionPoint.z - referencePoint.z;
            if (length) {
                if (intersectionPoint.x < referencePoint.x || xDistance > length) {
                    t1 = -1;
                }
            }
            if (width) {
                if (intersectionPoint.y < referencePoint.y || yDistance > width) {
                    t1 = -1;
                }
            }
            if (height) {
                if (intersectionPoint.z < referencePoint.z || zDistance > height) {
                    t1 = -1;
                }
            }

            intersectionPoint = ray.start + t2 * ray.direction;
            xDistance = intersectionPoint.x - referencePoint.x;
            yDistance = intersectionPoint.y - referencePoint.y;
            zDistance = intersectionPoint.z - referencePoint.z;
            if (length) {
                if (intersectionPoint.x < referencePoint.x || xDistance > length) {
                    t2 = -1;
                }
            }
            if (width) {
                if (intersectionPoint.y < referencePoint.y || yDistance > width) {
                    t2 = -1;
                }
            }
            if (height) {
                if (intersectionPoint.z < referencePoint.z || zDistance > height) {
                    t2 = -1;
                }
            }

            if (t1 < 0 || t2 < 0) {
                t = max(t1, t2);
            } else {
                t = min(t1, t2);
            }
        } else {
            t = -1;
        }
        return t;
    }

    Vector getNormal(Vector intersectionPoint) override {
        Vector normal(2 * A * intersectionPoint.x + D * intersectionPoint.y + F * intersectionPoint.z + G,
                      2 * B * intersectionPoint.y + D * intersectionPoint.x + E * intersectionPoint.z + H,
                      2 * C * intersectionPoint.z + E * intersectionPoint.y + F * intersectionPoint.x + I);
        normal.normalize();
        return normal;
    }
};

class Camera {
    double unitDistance;
    double unitAngle;
public:
    Vector cameraPosition, cameraLook, cameraRight, cameraUp;

    Camera() {
        unitDistance = 2;
        unitAngle = 0.05;
    }

    void button1() {
        cameraRight = cameraRight * cos(-unitAngle) + cameraLook * sin(-unitAngle);
        cameraRight.normalize();
        cameraLook = cameraUp * cameraRight;
        cameraLook.normalize();
    }

    void button2() {
        cameraRight = cameraRight * cos(unitAngle) + cameraLook * sin(unitAngle);
        cameraRight.normalize();
        cameraLook = cameraUp * cameraRight;
        cameraLook.normalize();
    }

    void button3() {
        cameraLook = cameraLook * cos(unitAngle) + cameraUp * sin(unitAngle);
        cameraLook.normalize();
        cameraUp = cameraRight * cameraLook;
        cameraUp.normalize();
    }

    void button4() {
        cameraLook = cameraLook * cos(-unitAngle) + cameraUp * sin(-unitAngle);
        cameraLook.normalize();
        cameraUp = cameraRight * cameraLook;
        cameraUp.normalize();
    }

    void button5() {
        cameraUp = cameraUp * cos(unitAngle) + cameraRight * sin(unitAngle);
        cameraUp.normalize();
        cameraRight = cameraLook * cameraUp;
        cameraRight.normalize();
    }

    void button6() {
        cameraUp = cameraUp * cos(-unitAngle) + cameraRight * sin(-unitAngle);
        cameraUp.normalize();
        cameraRight = cameraLook * cameraUp;
        cameraRight.normalize();
    }

    void buttonUpArrow() {
        cameraPosition += cameraLook * unitDistance;
    }

    void buttonDownArrow() {
        cameraPosition -= cameraLook * unitDistance;
    }

    void buttonLeftArrow() {
        cameraPosition -= cameraRight * unitDistance;
    }

    void buttonRightArrow() {
        cameraPosition += cameraRight * unitDistance;
    }

    void buttonPageDown() {
        cameraPosition -= cameraUp * unitDistance;
    }

    void buttonPageUp() {
        cameraPosition += cameraUp * unitDistance;
    }
};

Camera camera;

void capture() {
    unsigned int i, j, k;
    int nearest;
    double planeDistance;
    Vector topLeft, corner, direction;
    double du, dv, minimumT;

    Colour **frameBuffer;
    frameBuffer = new Colour *[imageHeight];
    for (i = 0; i < imageHeight; i++) {
        frameBuffer[i] = new Colour[imageWidth];
    }

    planeDistance = (windowHeight / 2) / tan(fovY / 2 * (PI / 180));
    topLeft = camera.cameraPosition
              + camera.cameraLook * planeDistance
              - camera.cameraRight * (windowWidth / 2)
              + camera.cameraUp * (windowHeight / 2);
    du = (windowWidth * 1.0) / imageWidth;
    dv = (windowHeight * 1.0) / imageHeight;

    for (i = 0; i < imageHeight; i++) {
        for (j = 0; j < imageWidth; j++) {
            corner = topLeft + camera.cameraRight * (j * du) - camera.cameraUp * (i * dv);
            direction = corner - camera.cameraPosition;
            direction.normalize();

            Ray ray(camera.cameraPosition, direction);

            nearest = -1;
            minimumT = INFINITY;
            for (k = 0; k < objects.size(); ++k) {
                Ray checkResult = objects[k]->intersect(ray, 0);
                if (checkResult.t <= 0) {
                    continue;
                }
                if (checkResult.t < minimumT) {
                    nearest = k;
                    minimumT = checkResult.t;
                }
            }
            if (nearest > -1) {
                Ray result = objects[nearest]->intersect(ray, rayTraceLevel);
                frameBuffer[i][j] = result.colour;
            }
        }
    }

    bitmap_image image(imageWidth, imageHeight);
    for (i = 0; i < imageHeight; i++) {
        for (j = 0; j < imageWidth; j++) {
            image.set_pixel(j, i, frameBuffer[i][j].r, frameBuffer[i][j].g, frameBuffer[i][j].b);
        }
    }

    image.save_image("output.bmp");
//    image.save_image("output.bmp");
}

/***Keyboard and mouse listeners***/
void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {
        case '1':
            camera.button1();
            break;
        case '2':
            camera.button2();
            break;
        case '3':
            camera.button3();
            break;
        case '4':
            camera.button4();
            break;

        case '5':
            camera.button5();
            break;
        case '6':
            camera.button6();
            break;
        case '0':
            capture();
            cout << "Image Captured " << endl;
            break;
        case '9':
            refraction = 1 - refraction;
            if (refraction) {
                cout << "Refraction turned on." << endl;
            } else {
                cout << "Refraction turned off." << endl;
            }
        default:
            break;
    }

}

void specialKeyListener(int key, int x, int y) {

    switch (key) {
        case GLUT_KEY_HOME:
            break;
        case GLUT_KEY_END:
            break;
        case GLUT_KEY_PAGE_DOWN:
            camera.buttonPageDown();
            break;
        case GLUT_KEY_PAGE_UP:
            camera.buttonPageUp();
            break;
        case GLUT_KEY_RIGHT:
            camera.buttonRightArrow();
            break;
        case GLUT_KEY_LEFT:
            camera.buttonLeftArrow();
            break;
        case GLUT_KEY_UP:
            camera.buttonUpArrow();
            break;
        case GLUT_KEY_DOWN:
            camera.buttonDownArrow();
            break;
        default:
            break;
    }
}

void mouseListener(int button, int state, int x, int y) {
    //x, y is the x-y of the screen (2D)
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

        case 3:             //scroll up
            break;
        case 4:             //scroll down
            break;
        default:
            break;
    }
}

void drawAxes() {
    glPushMatrix();
    {
        if (drawaxes == 1) {
            glBegin(GL_LINES);
            {
                glColor3f(1.0, 0, 0);
                glVertex3f(100, 0, 0);
                glVertex3f(-100, 0, 0);

                glColor3f(0, 1.0, 0);
                glVertex3f(0, -100, 0);
                glVertex3f(0, 100, 0);

                glColor3f(0, 0, 1.0);
                glVertex3f(0, 0, 100);
                glVertex3f(0, 0, -100);
            }
            glEnd();
        }
    }
    glPopMatrix();
}

void drawGrid() {
    int i;
    if (drawgrid == 1) {
        glColor3f(0.6, 0.6, 0.6);    //grey
        double a = 30;
        glBegin(GL_LINES);
        {
            for (i = int(-a); i <= a; i++) {
                //lines parallel to Y-axis
                glVertex3f(i * 10, GLfloat(-10 * a), 0);
                glVertex3f(i * 10, GLfloat(10 * a), 0);

                //lines parallel to X-axis
                glVertex3f(GLfloat(-10 * a), i * 10, 0);
                glVertex3f(GLfloat(10 * a), i * 10, 0);
            }
        }
        glEnd();
    }
}


/***main display function***/
void display() {

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);    //color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);

    gluLookAt(camera.cameraPosition.x, camera.cameraPosition.y, camera.cameraPosition.z,
              camera.cameraPosition.x + cameraDistance * camera.cameraLook.x,
              camera.cameraPosition.y + cameraDistance * camera.cameraLook.y,
              camera.cameraPosition.z + cameraDistance * camera.cameraLook.z,
              camera.cameraUp.x, camera.cameraUp.y, camera.cameraUp.z);

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
        for (int i = 0; i < objects.size(); ++i) {
            objects[i]->draw();
        }
        glColor3f(1, 1, 1);
        for (int j = 0; j < lights.size(); ++j) {
            glPushMatrix();
            {
                glTranslatef((GLfloat) lights[j].x, (GLfloat) lights[j].y, (GLfloat) lights[j].z);
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
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init() {
    //codes for initialization

    /** Camera initialization **/
    Vector origin;
    camera.cameraPosition = Vector(100, 100, 0);
    camera.cameraLook = origin - camera.cameraPosition;
    camera.cameraLook.normalize();
    camera.cameraUp = Vector(0, 0, 1);
    camera.cameraRight = camera.cameraLook * camera.cameraUp;
    drawgrid = 1;
    drawaxes = 1;

    //clear the screen
    glClearColor(0, 0, 0, 0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(fovY, aspectRatio, zNear, zFar);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

void loadTestData() {
    windowWidth = 500;
    windowHeight = 500;
    cameraDistance = 25;
    aspectRatio = (windowWidth * 1.0) / windowHeight;
    fovY = 90;
    zNear = 1;
    zFar = 1000.0;
    imageWidth = 1024;
    imageHeight = 1024;
    rayTraceLevel = 10;
    refraction = 1;

    Object *temp;
    Vector center(50, 50, 10);
    double radius = 10.0;
    temp = new Sphere(center, radius);
    temp->setColour(250, 0, 0);
    temp->setCoefficients(0.4, 0.2, 0.2, 0.6);
    temp->setRefraction(0.5, 1.2);
    temp->setShine(1);

    objects.push_back(temp);

//    center.x = 20;
//    center.y = 20;
//    center.z = 20;
//    radius = 8.0;
//    temp = new Sphere(center, radius);
//    temp->setColour(0, 0, 250);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.5);
//    temp->setRefraction(0.5, 1.2);
//    temp->setShine(1);
//
//    objects.push_back(temp);

    temp = new Floor(200, 200, 10, 10);
    temp->setCoefficients(0.4, 0.2, 0.2, 0.5);
    temp->setRefraction(0.1, 1.2);
    temp->setShine(1);
    objects.push_back(temp);

//    center.x = -5;
//    center.y = -30;
//    center.z = 40;
//    radius = 3.0;
//    temp = new Sphere(center, radius);
//    temp->setColour(0, 250, 250);
//    temp->setCoefficients(0.4, 0.2, 0.2, 0.5);
//    temp->setRefraction(0.2, 1.2);
//    temp->setShine(1);
//    objects.push_back(temp);

//    center.x = -5;
//    center.y = 5;
//    center.z = 26;
//    radius = 5;
//    temp = new Sphere(center, radius);
//    temp->setColour(0, 250, 0);
//    temp->setCoefficients(0.4, 0.3, 0.3, 0.4);
//    temp->setRefraction(0.3, 1.2);
//    temp->setShine(1);
//    objects.push_back(temp);



    temp = new Triangle(Vector(-50, 0, 0), Vector(0, 0, 50), Vector(50, 0, 0));
    temp->setColour(200, 50, 30);
    temp->setCoefficients(0.3, 0.1, 0.2, 0.4);

    temp->setShine(1);
    objects.push_back(temp);

//    temp = new GeneralQuadratic(144, 64, -36, 0, 0, 0, 0, 0, 0, 0);
//    temp->setReferencePoint(Vector(0, 0, 0));
//    temp->setDimensions(0, 0, 10);
//    temp->setColour(250, 0, 250);
//    temp->setCoefficients(0.3, 0.9, 0.9, 0.5);
//    temp->setRefraction(0.5, 1.2);
//    temp->setShine(1);
//    objects.push_back(temp);

    Vector light1(15, 15, 60);
    lights.push_back(light1);
}

void loadActualData() {
    int numberOfObjects;
    int numberOfLights;
    string objectType;
    windowWidth = 500;
    windowHeight = 500;
    cameraDistance = 25;
    aspectRatio = (windowWidth * 1.0) / windowHeight;
    fovY = 90;
    zNear = 1;
    zFar = 1000.0;
    imageWidth = 1024;
    imageHeight = 1024;
    rayTraceLevel = 10;
    refraction = 0;
    Object *temp;

    ifstream inputFile;
    inputFile.open("scene.txt");
    inputFile >> rayTraceLevel;
    cout << "raytracelevel: " << rayTraceLevel << endl << endl;
    inputFile >> imageWidth;
    imageHeight = imageWidth;
    inputFile >> numberOfObjects;
    for (int i = 0; i < numberOfObjects; ++i) {
        inputFile >> objectType;
        if (objectType == "sphere") {
            double x, y, z, radius, r, g, b, ambient, diffuse, specular, reflection;
            int shine;
            inputFile >> x >> y >> z;
            inputFile >> radius;
            inputFile >> r >> g >> b;
            inputFile >> ambient >> diffuse >> specular >> reflection;
            inputFile >> shine;
            temp = new Sphere(Vector(x, y, z), radius);
            temp->setColour((unsigned char) (r * 255.0), (unsigned char) (g * 255.0), (unsigned char) (b * 255.0));
            temp->setCoefficients(ambient, diffuse, specular, reflection);
            temp->setShine(shine);
            temp->setRefraction(0.2, 1.3);
            objects.push_back(temp);
        }
        if (objectType == "triangle") {
            double x1, y1, z1, x2, y2, z2, x3, y3, z3, r, g, b, ambient, diffuse, specular, reflection;
            int shine;
            inputFile >> x1 >> y1 >> z1;
            inputFile >> x2 >> y2 >> z2;
            inputFile >> x3 >> y3 >> z3;
            inputFile >> r >> g >> b;
            inputFile >> ambient >> diffuse >> specular >> reflection;
            inputFile >> shine;
            Vector v1(x1, y1, z1), v2(x2, y2, z2), v3(x3, y3, z3);
            temp = new Triangle(v1, v2, v3);
            temp->setColour((unsigned char) (r * 255.0), (unsigned char) (g * 255.0), (unsigned char) (b * 255.0));
            temp->setCoefficients(ambient, diffuse, specular, reflection);
            temp->setShine(shine);
            temp->setRefraction(0.4, 1.5);
            objects.push_back(temp);
        }
        if (objectType == "general") {
            double A, B, C, D, E, F, G, H, I, J, x, y, z, length, width, height, r, g, b, ambient, diffuse, specular, reflection;
            int shine;
            inputFile >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J;
            inputFile >> x >> y >> z >> length >> width >> height;
            inputFile >> r >> g >> b;
            inputFile >> ambient >> diffuse >> specular >> reflection;
            inputFile >> shine;
            temp = new GeneralQuadratic(A, B, C, D, E, F, G, H, I, J);
            temp->setReferencePoint(Vector(x, y, z));
            temp->setDimensions(length, width, height);
            temp->setColour((unsigned char) (r * 255.0), (unsigned char) (g * 255.0), (unsigned char) (b * 255.0));
            temp->setCoefficients(ambient, diffuse, specular, reflection);
            temp->setShine(shine);
            temp->setRefraction(0.4, 1.5);
            objects.push_back(temp);
        }
    }


    temp = new Floor(1000, 1000, 10, 10);
    temp->setCoefficients(0.4, 0.2, 0.2, 0.4);
    temp->setShine(1);
    temp->setRefraction(0.0, 1.5);
    objects.push_back(temp);

    inputFile >> numberOfLights;
    for (int j = 0; j < numberOfLights; ++j) {
        double x, y, z;
        inputFile >> x >> y >> z;
        Vector light(x, y, z);
        lights.push_back(light);
    }
    cout << "Data load complete." << endl << "Press 9 to toggle refraction." << endl;
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
//    loadTestData();
    loadActualData();
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);    //Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracer");

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
