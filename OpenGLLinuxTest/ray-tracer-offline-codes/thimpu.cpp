#include <bits/stdc++.h>

//Ray Tracing

#include <vector>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <vector>
#include <math.h>
#include <iostream>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "bitmap_image.hpp"

using namespace std;

#define AMBIENT 0
#define DIFFUSE 1
#define SPECULAR 2
#define REFLECTION 3

#define EPSILON 0.000001

extern int recursionLevel;


struct Point {
    double x, y, z;

    Point() {
        x = 0;
        y = 0;
        z = 0;
    }

    Point(double xN, double yN, double zN) {
        x = xN;
        y = yN;
        z = zN;
    }

    Point operator+(Point pt) {
        Point temp(this->x + pt.x, this->y + pt.y, this->z + pt.z);
        return temp;
    }

    Point operator-(Point pt) {
        Point temp(this->x - pt.x, this->y - pt.y, this->z - pt.z);
        return temp;
    }

    Point operator*(Point pt) {
        Point temp(this->x * pt.x, this->y * pt.y, this->z * pt.z);
        return temp;
    }

    Point operator*(double val) {
        Point temp(this->x * val, this->y * val, this->z * val);
        return temp;
    }

    void normalize() {
        double dividend = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        this->x = this->x / dividend;
        this->y = this->y / dividend;
        this->z = this->z / dividend;
    }

    void printPoint() {
        cout << x << " " << y << " " << z << endl;
    }
};

double dotProduct(Point p1, Point p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

Point crossProduct(Point p1, Point p2) {
    Point crossProd(p1.y * p2.z - p2.y * p1.z, -(p1.x * p2.z - p2.x * p1.z), p1.x * p2.y - p2.x * p1.y);
    return crossProd;
}

class Ray {
public:
    Point start;
    Point dir;

    Ray(Point start, Point dir) {
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
};


class Object {
protected:
    Point reference_Point;
    double height, width, length;
    int Shine;
    double color[3];
    double co_efficients[4];
    double source_factor = 1.0;
    double refractionIndex = 1.5;
public:
    Object() {};

    virtual void draw() {};

    virtual double getIntersectingT(Ray *r) {
        return -1;
    }

    virtual double intersect(Ray *r, double *current_color, int level) {
        return -1;
    }

    void setColor(double color1, double color2, double color3) {
        color[0] = color1;
        color[1] = color2;
        color[2] = color3;
    }

    void setShine(int shineValue) {
        Shine = shineValue;
    }

    void setCoEfficients(double val1, double val2, double val3, double val4) {
        co_efficients[0] = val1;
        co_efficients[1] = val2;
        co_efficients[2] = val3;
        co_efficients[3] = val4;
    }

    Point getReflection(Ray *ray, Point normal) {
        Point temp = ray->dir - normal * dotProduct(ray->dir, normal) * 2.0;
        temp.normalize();
        return temp;
    }

    Point getRefraction(Ray *ray, Point normal) {

        Point refraction(0, 0, 0);

        double dot = dotProduct(normal, ray->dir);
        double k = 1.0 - refractionIndex * refractionIndex * (1.0 - dot * dot);

        if (k >= 0) {
            refraction = ray->dir * refractionIndex - normal * (refractionIndex * dot + sqrt(k));
            refraction.normalize();
        }

        return refraction;
    }

    void setCurrentColor(double *current_color) {
        for (int i = 0; i < 3; i++) {
            current_color[i] = this->color[i] * this->co_efficients[AMBIENT];
        }
    }

    void ClipCurrentColor(double *current_color) {
        for (int i = 0; i < 3; i++) {
            if (current_color[i] < 0) current_color[i] = 0;
            if (current_color[i] > 1) current_color[i] = 1;
        }
    }


};


extern vector<Object *> objects;
extern vector<Point> lights;


class Sphere : public Object {
public:
    Sphere(Point center, double radius) {
        reference_Point = center;
        length = radius;
    }

    void draw() {
        glPushMatrix();
        glColor3f(color[0], color[1], color[2]);
        glTranslatef(reference_Point.x, reference_Point.y, reference_Point.z);
        glutSolidSphere(length, 100, 100);
        glPopMatrix();
    }

    double getIntersectingT(Ray *r) {
        Point origin = r->start - reference_Point;

        double a = dotProduct(r->dir, r->dir);
        double b = 2 * dotProduct(r->dir, origin);
        double c = dotProduct(origin, origin) - length * length;

        double d = b * b - 4 * a * c;

        if (d < 0) return -1;

        double t1 = (double) (-b + sqrt(d)) / (double) (2 * a);
        double t2 = (double) (-b - sqrt(d)) / (double) (2 * a);

        if (t1 < t2)return t1;
        else return t2;
    }

    Point getNormal(Point intersectionPoint) {
        Point temp = intersectionPoint - reference_Point;
        temp.normalize();
        return temp;
    }

    double intersect(Ray *r, double *current_color, int level) {

        double t = getIntersectingT(r);

        if (t <= 0)return -1;
        if (level == 0)return t;


        Point intersectionPoint = r->start + r->dir * t;

        setCurrentColor(current_color);

        Point normal = getNormal(intersectionPoint);
        Point reflection = getReflection(r, normal);
        Point refraction = getRefraction(r, normal);

        for (int i = 0; i < lights.size(); i++) {
            Point direction = lights[i] - intersectionPoint;
            double rayLength = sqrt(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z);
            direction.normalize();
            Point start = intersectionPoint + direction * 1.0;

            Ray *tempRay = new Ray(start, direction);

            bool clearPath = true;

            for (int j = 0; j < objects.size(); j++) {
                double tempT = objects[j]->getIntersectingT(tempRay);
                if (tempT >= 0 && tempT <= rayLength) {
                    clearPath = false;
                    break;

                }
            }

            if (clearPath) {
                double lambertValue = dotProduct(normal, tempRay->dir);
                double phongValue = pow(dotProduct(reflection, r->dir), Shine);

                if (lambertValue < 0)lambertValue = 0;
                if (phongValue < 0)phongValue = 0;

                for (int k = 0; k < 3; k++) {
                    current_color[k] += source_factor * lambertValue * co_efficients[DIFFUSE] * this->color[k];
                    current_color[k] += source_factor * phongValue * co_efficients[SPECULAR] * this->color[k];
                }

            }

            if (level < recursionLevel) {

                struct Point start = intersectionPoint + reflection * 1.0;
                Ray *reflectionRay = new Ray(start, reflection);

                int nearest = -1;
                double reflectedcolor[3];
                double minimum_t = 999999999;

                for (int k = 0; k < objects.size(); k++) {
                    double t = objects[k]->intersect(reflectionRay, reflectedcolor, 0);
                    if (t <= 0)continue;

                    if (t < minimum_t) {
                        minimum_t = t;
                        nearest = k;
                    }
                }
                if (nearest != -1) {
                    double t = objects[nearest]->intersect(reflectionRay, reflectedcolor, level + 1);
                    for (int p = 0; p < 3; p++) {
                        current_color[p] += reflectedcolor[p] * co_efficients[REFLECTION];
                        if (current_color[p] < 0)current_color[p] = 0;
                        if (current_color[p] > 1)current_color[p] = 1;
                    }
                }

                //Refraction Part starts here

                start = intersectionPoint + refraction * 1.0;

                Ray *refractionRay = new Ray(start, refraction);

                nearest = -1;
                minimum_t = 9999999;
                double refracted_color[3];

                for (int k = 0; k < objects.size(); k++) {

                    double tk = objects[k]->getIntersectingT(refractionRay);

                    if (tk <= 0) {
                        continue;
                    } else if (tk < minimum_t) {
                        minimum_t = tk;
                        nearest = k;
                    }

                    //cout<<tk<<endl;
                }

                if (nearest != -1) {

                    objects[nearest]->intersect(refractionRay, refracted_color, level + 1);

                    for (int p = 0; p < 3; p++) {
                        current_color[p] += refracted_color[p] * refractionIndex;
                    }
                }

            }
        }

        ClipCurrentColor(current_color);

        return t;

    }

};

class Floor : public Object {
    bitmap_image texture;
    double textureHeight, textureWidth;
    Point reference_Point;
    double length;
public:
    Floor(double FloorWidth, double TileWidth) {
        reference_Point = Point(-FloorWidth / 2, -FloorWidth / 2, 0);
        length = TileWidth;
        texture = bitmap_image("texture.bmp");
        textureHeight = texture.height() / FloorWidth;
        textureWidth = texture.width() / FloorWidth;
    }

    void draw() {

        int numOfTiles = abs(reference_Point.x * 2 / length);
        int seq = 0;

        for (int i = 0; i < numOfTiles; i++) {
            for (int j = 0; j < numOfTiles; j++) {

                if ((i + j) % 2) {
                    glColor3f(0, 0, 0);
                } else {
                    glColor3f(1, 1, 1);
                }
                glRecti(reference_Point.x + i * length, reference_Point.y + j * length,
                        reference_Point.x + (i + 1) * length, reference_Point.y + (j + 1) * length);
            }
        }
    }

    Point getNormal(Point intersection) {
        Point temp(0, 0, 1);
        temp.normalize();
        return temp;
    }

    double getIntersectingT(Ray *ray) {

        Point normal = getNormal(reference_Point);

        double t = dotProduct(normal, ray->start) * (-1) / dotProduct(normal, ray->dir);

        return t;
    }

    double intersect(Ray *r, double *current_color, int level) {
        double t = getIntersectingT(r);

        if (t <= 0)return -1;


        Point intersectionPoint = r->start + r->dir * t;
        double xMin = reference_Point.x;
        double xMax = xMin * (-1);

        double yMin = reference_Point.y;
        double yMax = yMin * (-1);


        if (xMin > intersectionPoint.x || intersectionPoint.x > xMax ||
            yMin > intersectionPoint.y || intersectionPoint.y > yMax) {
            return -1;
        }

        if (level == 0)return t;

        int xCord = intersectionPoint.x / length;
        int yCord = intersectionPoint.y / length;

        if ((xCord + yCord) % 2) {
            this->color[0] = this->color[1] = this->color[2] = 0;
        } else {
            this->color[0] = this->color[1] = this->color[2] = 1;
        }

        // setCurrentColor(current_color);
        unsigned char red, green, blue;
        int x = (intersectionPoint.x + abs(reference_Point.x)) * textureWidth;
        int y = (intersectionPoint.y + abs(reference_Point.y)) * textureHeight;

        //cout<<x<<','<<y<<endl;

        texture.get_pixel(x, y, red, green, blue);

        double rgb[] = {red, green, blue};
        //cout<<rgb[0];

        for (int i = 0; i < 3; i++) {
            current_color[i] = color[i] * co_efficients[AMBIENT] * rgb[i] / 255.0;
        }

        Point normal = getNormal(intersectionPoint);
        Point reflection = getReflection(r, normal);

        for (int i = 0; i < lights.size(); i++) {
            Point direction = lights[i] - intersectionPoint;
            double rayLength = sqrt(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z);
            direction.normalize();
            Point start = intersectionPoint + direction * 1.0;

            Ray *tempRay = new Ray(start, direction);

            bool clearPath = true;

            for (int j = 0; j < objects.size(); j++) {
                double tempT = objects[j]->getIntersectingT(tempRay);
                if (tempT >= 0 && tempT <= rayLength) {
                    clearPath = false;
                    break;
                }
            }

            if (clearPath) {
                double lambertValue = dotProduct(normal, tempRay->dir);
                double phongValue = pow(dotProduct(reflection, r->dir), Shine);

                if (lambertValue < 0)lambertValue = 0;
                if (phongValue < 0)phongValue = 0;

                for (int k = 0; k < 3; k++) {
                    current_color[k] += source_factor * lambertValue * co_efficients[DIFFUSE] * this->color[k];
                    current_color[k] += source_factor * phongValue * co_efficients[SPECULAR] * this->color[k];
                }

            }


            if (level < recursionLevel) {

                struct Point start = intersectionPoint + reflection * 1.0;
                Ray *reflectionRay = new Ray(start, reflection);

                int nearest = -1;
                double reflectedcolor[3];
                double minimum_t = 999999999;

                for (int k = 0; k < objects.size(); k++) {
                    double t = objects[k]->intersect(reflectionRay, reflectedcolor, 0);
                    if (t <= 0)continue;

                    if (t < minimum_t) {
                        minimum_t = t;
                        nearest = k;
                    }
                }
                if (nearest != -1) {
                    double t = objects[nearest]->intersect(reflectionRay, reflectedcolor, level + 1);
                    for (int p = 0; p < 3; p++) {
                        current_color[p] += reflectedcolor[p] * co_efficients[REFLECTION];
                    }
                }


            }

        }
        ClipCurrentColor(current_color);

        return t;
    }


};


class Triangle : public Object {

public:

    Point vertex1, vertex2, vertex3;

    Triangle(Point a, Point b, Point c) {
        this->vertex1 = a;
        this->vertex2 = b;
        this->vertex3 = c;
    }

    void draw() {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(vertex1.x, vertex1.y, vertex1.z);
            glVertex3f(vertex2.x, vertex2.y, vertex2.z);
            glVertex3f(vertex3.x, vertex3.y, vertex3.z);
        }
        glEnd();
    }

    Point getNormal(Point intersection) {

        Point u = vertex2 - vertex1;
        Point v = vertex3 - vertex1;

        Point temp = crossProduct(u, v);
        temp.normalize();
        return temp;
    }


    double getIntersectingT(Ray *ray) {
        Point edge1 = vertex2 - vertex1;
        Point edge2 = vertex3 - vertex1;

        Point h = crossProduct(ray->dir, edge2);
        double a = dotProduct(edge1, h);

        if (a > -EPSILON && a < EPSILON) {
            return -1;
        }

        double f = 1.0 / a;

        Point s = ray->start - vertex1;

        double u = f * dotProduct(s, h);

        if ((u < 0.0) || (u > 1.0)) {
            return -1;
        }

        Point q = crossProduct(s, edge1);

        double v = f * dotProduct(ray->dir, q);

        if ((v < 0.0) || ((u + v) > 1.0)) {
            return -1;
        }

        double t = f * dotProduct(edge2, q);

        if (t > EPSILON) { //ray intersection
            return t;
        }

        return -1;
    }


    double intersect(Ray *r, double *current_color, int level) {

        double t = getIntersectingT(r);

        if (t <= 0)return -1;
        if (level == 0)return t;


        Point intersectionPoint = r->start + r->dir * t;


        setCurrentColor(current_color);

        Point normal = getNormal(intersectionPoint);
        Point reflection = getReflection(r, normal);

        for (int i = 0; i < lights.size(); i++) {
            Point direction = lights[i] - intersectionPoint;
            double rayLength = sqrt(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z);
            direction.normalize();
            Point start = intersectionPoint + direction * 1.0;

            Ray *tempRay = new Ray(start, direction);

            bool clearPath = true;

            for (int j = 0; j < objects.size(); j++) {
                double tempT = objects[j]->getIntersectingT(tempRay);
                if (tempT >= 0 && tempT <= rayLength) {
                    clearPath = false;
                    break;
                }
            }

            if (clearPath) {
                double lambertValue = dotProduct(normal, tempRay->dir);
                double phongValue = pow(dotProduct(reflection, r->dir), Shine);

                if (lambertValue < 0)lambertValue = 0;
                if (phongValue < 0)phongValue = 0;

                for (int k = 0; k < 3; k++) {
                    current_color[k] += source_factor * lambertValue * co_efficients[DIFFUSE] * this->color[k];
                    current_color[k] += source_factor * phongValue * co_efficients[SPECULAR] * this->color[k];
                }

            }

            if (level < recursionLevel) {

                struct Point start = intersectionPoint + reflection * 1.0;
                Ray *reflectionRay = new Ray(start, reflection);

                int nearest = -1;
                double reflectedcolor[3];
                double minimum_t = 999999999;

                for (int k = 0; k < objects.size(); k++) {
                    double t = objects[k]->intersect(reflectionRay, reflectedcolor, 0);
                    if (t <= 0)continue;

                    if (t < minimum_t) {
                        minimum_t = t;
                        nearest = k;
                    }
                }
                if (nearest != -1) {
                    double t = objects[nearest]->intersect(reflectionRay, reflectedcolor, level + 1);
                    for (int p = 0; p < 3; p++) {
                        current_color[p] += reflectedcolor[p] * co_efficients[REFLECTION];
                    }
                }


            }
        }
        ClipCurrentColor(current_color);

        return t;

    }
};

class GeneralQuadratic : public Object {

public:

    double A, B, C, D, E, F, G, H, I, J;

    GeneralQuadratic(double coeff[10], Point reff, double length, double width, double height) {
        this->A = coeff[0];
        this->B = coeff[1];
        this->C = coeff[2];
        this->D = coeff[3];
        this->E = coeff[4];
        this->F = coeff[5];
        this->G = coeff[6];
        this->H = coeff[7];
        this->I = coeff[8];
        this->J = coeff[9];
        this->reference_Point = reff;
        this->height = height;
        this->width = width;
        this->length = length;
    }

    void draw() {}

    Point getNormal(Point intersection) {

        double u = 2 * A * intersection.x + D * intersection.y + F * intersection.z + G;
        double v = 2 * B * intersection.y + D * intersection.x + E * intersection.z + H;
        double z = 2 * C * intersection.z + E * intersection.y + F * intersection.x + I;

        Point normal(u, v, z);
        normal.normalize();

        return normal;
    }

    double getIntersectingT(Ray *ray) {

        double a = A * ray->dir.x * ray->dir.x + B * ray->dir.y * ray->dir.y + C * ray->dir.z * ray->dir.z;
        double b = 2 * (A * ray->start.x * ray->dir.x + B * ray->start.y * ray->dir.y + C * ray->start.z * ray->dir.z);
        double c = A * ray->start.x * ray->start.x + B * ray->start.y * ray->start.y + C * ray->start.z * ray->start.z;

        a += D * ray->dir.x * ray->dir.y + E * ray->dir.y * ray->dir.z + F * ray->dir.z * ray->dir.x;
        b += D * (ray->start.x * ray->dir.y + ray->dir.x * ray->start.y)
             + E * (ray->start.y * ray->dir.z + ray->dir.y * ray->start.z)
             + F * (ray->start.z * ray->dir.x + ray->dir.z * ray->start.x);
        c += D * ray->start.x * ray->start.y + E * ray->start.y * ray->start.z + F * ray->start.z * ray->start.x;

        b += G * ray->dir.x + H * ray->dir.y + I * ray->dir.z;
        c += G * ray->start.x + H * ray->start.y + I * ray->start.z + J;


        double d = b * b - 4 * a * c;


        if (d < 0) {
            return -1;
        }

        double t1 = (-b + sqrt(d)) / (2.0 * a);
        double t2 = (-b - sqrt(d)) / (2.0 * a);


        Point intersectionPoint1 = ray->start + ray->dir * t1;
        Point intersectionPoint2 = ray->start + ray->dir * t2;

        double xMin = reference_Point.x;
        double xMax = xMin + length;

        double yMin = reference_Point.y;
        double yMax = yMin + width;

        double zMin = reference_Point.z;
        double zMax = zMin + height;


        bool flag1 = (length > 0 && (xMin > intersectionPoint1.x || intersectionPoint1.x > xMax) ||
                      width > 0 && (yMin > intersectionPoint1.y || intersectionPoint1.y > yMax) ||
                      height > 0 && (zMin > intersectionPoint1.z || intersectionPoint1.z > zMax));

        bool flag2 = (length > 0 && (xMin > intersectionPoint2.x || intersectionPoint2.x > xMax) ||
                      width > 0 && (yMin > intersectionPoint2.y || intersectionPoint2.y > yMax) ||
                      height > 0 && (zMin > intersectionPoint2.z || intersectionPoint2.z > zMax));


        if (flag1 && flag2) {
            return -1;
        } else if (flag1) {
            return t2;
        } else if (flag2) {
            return t1;
        } else {
            return min(t1, t2);
        }
    }

    double intersect(Ray *r, double *current_color, int level) {

        double t = getIntersectingT(r);

        if (t <= 0)return -1;
        if (level == 0)return t;


        Point intersectionPoint = r->start + r->dir * t;

        setCurrentColor(current_color);

        Point normal = getNormal(intersectionPoint);
        Point reflection = getReflection(r, normal);

        for (int i = 0; i < lights.size(); i++) {
            Point direction = lights[i] - intersectionPoint;
            double rayLength = sqrt(direction.x * direction.x + direction.y * direction.y + direction.z * direction.z);
            direction.normalize();
            Point start = intersectionPoint + direction * 1.0;

            Ray *tempRay = new Ray(start, direction);

            bool clearPath = true;

            for (int j = 0; j < objects.size(); j++) {
                double tempT = objects[j]->getIntersectingT(tempRay);
                if (tempT >= 0 && tempT <= rayLength) {
                    clearPath = false;
                    break;
                }
            }

            if (clearPath) {
                double lambertValue = dotProduct(normal, tempRay->dir);
                double phongValue = pow(dotProduct(reflection, r->dir), Shine);

                if (lambertValue < 0)lambertValue = 0;
                if (phongValue < 0)phongValue = 0;

                for (int k = 0; k < 3; k++) {
                    current_color[k] += source_factor * lambertValue * co_efficients[DIFFUSE] * this->color[k];
                    current_color[k] += source_factor * phongValue * co_efficients[SPECULAR] * this->color[k];
                }

            }

            if (level < recursionLevel) {

                struct Point start = intersectionPoint + reflection * 1.0;
                Ray *reflectionRay = new Ray(start, reflection);

                int nearest = -1;
                double reflectedcolor[3];
                double minimum_t = 999999999;

                for (int k = 0; k < objects.size(); k++) {
                    double t = objects[k]->intersect(reflectionRay, reflectedcolor, 0);
                    if (t <= 0)continue;

                    if (t < minimum_t) {
                        minimum_t = t;
                        nearest = k;
                    }
                }
                if (nearest != -1) {
                    double t = objects[nearest]->intersect(reflectionRay, reflectedcolor, level + 1);
                    for (int p = 0; p < 3; p++) {
                        current_color[p] += reflectedcolor[p] * co_efficients[REFLECTION];
                    }
                }


            }
        }
        ClipCurrentColor(current_color);

        return t;

    }


};

#include "bitmap_image.hpp"

using namespace std;

#define pi (2 * acos(0.0))
#define VIEW_ANGLE 80

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

//Ray Tracing
struct Point eye;
struct Point u;
struct Point r;
struct Point l;

double windowWidth = 500;
double windowHeight = 500;

int imageHeight, imageWidth;

vector<Object *> objects;
vector<Point> lights;

int recursionLevel;

//Ray Tracing

void Capture() {
    bitmap_image image(imageHeight, imageWidth);
    double plane_distance =
            (double) ((double) windowHeight / (double) 2.0) / (double) tan((double) (VIEW_ANGLE * pi) / (double) 360.0);

    Point topleft = eye + l * plane_distance - r * (windowWidth / 2.0) + u * (windowHeight / 2.0);

    double du = (double) windowWidth / (double) imageWidth;
    double dv = (double) windowHeight / (double) imageWidth;

    for (int i = 0; i < imageWidth; i++) {
        for (int j = 0; j < imageWidth; j++) {
            Point corner = topleft + r * j * du - u * i * dv;

            Ray *ray = new Ray(eye, corner - eye);

            int nearest = -1;
            double color[3];
            double minimum_t = 999999999;

            for (int k = 0; k < objects.size(); k++) {
                double t = objects[k]->intersect(ray, color, 0);
                if (t <= 0)
                    continue;

                if (t < minimum_t) {
                    minimum_t = t;
                    nearest = k;
                }
            }
            if (nearest != -1) {
                double t = objects[nearest]->intersect(ray, color, 1);
                image.set_pixel(j, i, color[0] * 255, color[1] * 255, color[2] * 255);
            }
        }
    }

    image.save_image("output.bmp");
    cout << "Image Generated" << endl;
}

void drawAxes() {
    if (drawaxes == 1) {
        glColor3f(1.0, 0.6, 0.8);
        glBegin(GL_LINES);
        {
            glVertex3f(100, 0, 0);
            glVertex3f(-100, 0, 0);

            glVertex3f(0, -100, 0);
            glVertex3f(0, 100, 0);

            glVertex3f(0, 0, 100);
            glVertex3f(0, 0, -100);
        }
        glEnd();
    }
}

void drawSquare(double a) {
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 0);
        glVertex3f(a, -a, 0);
        glVertex3f(-a, -a, 0);
        glVertex3f(-a, a, 0);
    }
    glEnd();
}

void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {
        case '1': {
            double angle = -(double) 3 / (double) 360 * 2 * pi;

            struct Point lCrossu = crossProduct(l, u);
            l = l * cos(angle) + lCrossu * sin(angle);

            struct Point rCrossu = crossProduct(r, u);
            r = r * cos(angle) + rCrossu * sin(angle);

            break;
        }
        case '2': {
            double angle = (double) 3 / (double) 360 * 2 * pi;

            struct Point lCrossu = crossProduct(l, u);
            l = l * cos(angle) + lCrossu * sin(angle);

            struct Point rCrossu = crossProduct(r, u);
            r = r * cos(angle) + rCrossu * sin(angle);

            break;
        }
        case '3': {
            double angle = -(double) 3 / (double) 360 * 2 * pi;

            struct Point lCrossr = crossProduct(l, r);
            l = l * cos(angle) + lCrossr * sin(angle);

            struct Point uCrossr = crossProduct(u, r);
            u = u * cos(angle) + uCrossr * sin(angle);

            break;
        }

        case '4': {
            double angle = (double) 3 / (double) 360 * 2 * pi;

            struct Point lCrossr = crossProduct(l, r);
            l = l * cos(angle) + lCrossr * sin(angle);

            struct Point uCrossr = crossProduct(u, r);
            u = u * cos(angle) + uCrossr * sin(angle);

            break;
        }
        case '5': {
            double angle = -(double) 3 / (double) 360 * 2 * pi;

            struct Point rCrossl = crossProduct(r, l);
            r = r * cos(angle) + rCrossl * sin(angle);

            struct Point uCrossl = crossProduct(u, l);
            u = u * cos(angle) + uCrossl * sin(angle);

            break;
        }

        case '6': {
            double angle = (double) 3 / (double) 360 * 2 * pi;

            struct Point rCrossl = crossProduct(r, l);
            r = r * cos(angle) + rCrossl * sin(angle);

            struct Point uCrossl = crossProduct(u, l);
            u = u * cos(angle) + uCrossl * sin(angle);

            break;
        }
        case '0': {
            Capture();
        }

        default:
            break;
    }
}

void specialKeyListener(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_DOWN: //down arrow key
            eye.x -= l.x;
            eye.y -= l.y;
            eye.y -= l.z;
            //cameraHeight -= 3.0;
            break;
        case GLUT_KEY_UP: // up arrow key
            eye.x += l.x;
            eye.y += l.y;
            eye.y += l.z;
            //cameraHeight += 3.0;
            break;

        case GLUT_KEY_RIGHT:
            //cameraAngle += 0.03;
            eye.x += r.x;
            eye.y += r.y;
            eye.z += r.z;
            break;
        case GLUT_KEY_LEFT:
            //cameraAngle -= 0.03;
            eye.x -= r.x;
            eye.y -= r.y;
            eye.z -= r.z;
            break;

        case GLUT_KEY_PAGE_UP:
            eye.x += u.x;
            eye.y += u.y;
            eye.z += u.z;
            break;
        case GLUT_KEY_PAGE_DOWN:
            eye.x -= u.x;
            eye.y -= u.y;
            eye.z -= u.z;
            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            break;
        case GLUT_KEY_END:
            break;

        default:
            break;
    }
}

void mouseListener(int button, int state, int x, int y) { //x, y is the x-y of the screen (2D)
    switch (button) {
        case GLUT_LEFT_BUTTON:
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

void drawAll() {

    for (int i = 0; i < (int) objects.size(); i++) {
        objects[i]->draw();
    }

    glColor3f(1, 1, 1);
    glBegin(GL_POINTS);
    for (int i = 0; i < lights.size(); i++) {
        glVertex3f(lights[i].x, lights[i].y, lights[i].z);
    }
    glEnd();
}

void display() {

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0); //color black
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

    //gluLookAt(eye.x, eye.y, eye.z, eye.x + l.x, eye.y + l.y, eye.z + l.z, u.x, u.y, u.z);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    //gluLookAt(0,0,200,	0,0,0,	0,1,0);

    //Ray Tracing
    gluLookAt(eye.x, eye.y, eye.z, eye.x + l.x, eye.y + l.y, eye.z + l.z, u.x, u.y, u.z);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);

    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    /****************************
	/ Add your objects from here
	****************************/
    //add objects

    drawAxes();
    drawAll();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate() {
    glutPostRedisplay();
}

void init() {
    //codes for initialization
    drawgrid = 0;
    drawaxes = 1;
    cameraHeight = 150.0;
    cameraAngle = 1.0;
    angle = 0;

    eye.x = 0;
    eye.y = -100;
    eye.z = 10;

    u.x = 0;
    u.y = 0;
    u.z = 1;

    r.x = 1;
    r.y = 0;
    r.z = 0;

    l.x = 0;
    l.y = 1;
    l.z = 0;

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
    gluPerspective(VIEW_ANGLE, 1, 1, 1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

//Ray Tracing
void loadTestData() {
    recursionLevel = 3;

    Object *temp;
    Point Center(0.0, 0.0, 10.0);
    double Radius = 10;
    temp = new Sphere(Center, Radius); // Center(0,0,10), Radius 10
    temp->setColor(1, 0, 0);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);
    objects.push_back(temp);

    Point Center2(20.0, 10.0, 40.0);
    double Radius2 = 5;
    temp = new Sphere(Center2, Radius2); // Center(0,0,10), Radius 10
    temp->setColor(0, 1, 0);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.7);
    temp->setShine(1);
    objects.push_back(temp);

    Point Center3(20.0, -50.0, 0.0);
    double Radius3 = 3;
    temp = new Sphere(Center3, Radius3); // Center(0,0,10), Radius 10
    temp->setColor(0, 0, 1);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.5);
    temp->setShine(1);
    objects.push_back(temp);

    Point Center4(30.0, -20.0, 5.0);
    double Radius4 = 7;
    temp = new Sphere(Center4, Radius4); // Center(0,0,10), Radius 10
    temp->setColor(0, 0, 1);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.9);
    temp->setShine(1);
    objects.push_back(temp);

    Point Center5(20.0, -50.0, 0.0);
    double Radius5 = 8;
    temp = new Sphere(Center5, Radius5); // Center(0,0,10), Radius 10
    temp->setColor(1, 0, 1);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);
    objects.push_back(temp);

    Point Center6(4.0, -16.0, 4.0);
    double Radius6 = 4;
    temp = new Sphere(Center6, Radius6); // Center(0,0,10), Radius 10
    temp->setColor(1, 1, 0);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.6);
    temp->setShine(1);
    objects.push_back(temp);

    Point Center7(50.0, 40.0, 50.0);
    double Radius7 = 15;
    temp = new Sphere(Center7, Radius7); // Center(0,0,10), Radius 10
    temp->setColor(0, 1, 1);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.8);
    temp->setShine(1);
    objects.push_back(temp);

    Point a(50, 30, 0);
    Point b(70, 60, 0);
    Point c(50, 45, 50);

    temp = new Triangle(a, b, c);
    temp->setColor(1.0, 0.0, 0.0);
    temp->setCoEfficients(0.4, 0.2, 0.1, 0.3);
    temp->setShine(5);
    objects.push_back(temp);

    Point light1(-50, 50, 50);
    lights.push_back(light1);

    Object *floorTemp = new Floor(1000, 20);
    floorTemp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    floorTemp->setShine(1);
    objects.push_back(floorTemp);

    imageWidth = 768;
}

void loadActualData() {

    freopen("scene.txt", "r", stdin);

    cin >> recursionLevel;
    cin >> imageWidth;
    imageHeight = imageWidth;

    int numOfObjects;
    cin >> numOfObjects;

    string objType;
    double a, b, c, d;

    Object *temp;

    for (int i = 0; i < numOfObjects; i++) {
        cin >> objType;

        if (objType == "sphere") {

            cin >> a >> b >> c;
            Point center(a, b, c);

            cin >> d;
            temp = new Sphere(center, d);

            cin >> a >> b >> c;
            temp->setColor(a, b, c);

            cin >> a >> b >> c >> d;
            temp->setCoEfficients(a, b, c, d);

            cin >> a;
            temp->setShine(a);

            objects.push_back(temp);
        } else if (objType == "triangle") {

            cin >> a >> b >> c;
            Point A(a, b, c);

            cin >> a >> b >> c;
            Point B(a, b, c);

            cin >> a >> b >> c;
            Point C(a, b, c);

            temp = new Triangle(A, B, C);

            cin >> a >> b >> c;
            temp->setColor(a, b, c);

            cin >> a >> b >> c >> d;
            temp->setCoEfficients(a, b, c, d);

            cin >> a;
            temp->setShine(a);

            objects.push_back(temp);
        } else if (objType == "general") {

            double coeff[10];
            for (int i = 0; i < 10; i++) {
                cin >> coeff[i];
            }

            cin >> a >> b >> c;
            Point reff(a, b, c);

            cin >> a >> b >> c;
            temp = new GeneralQuadratic(coeff, reff, a, b, c);

            cin >> a >> b >> c;
            temp->setColor(a, b, c);

            cin >> a >> b >> c >> d;
            temp->setCoEfficients(a, b, c, d);

            cin >> a;
            temp->setShine(a);

            objects.push_back(temp);
        }
    }

    cin >> numOfObjects;
    for (int i = 0; i < numOfObjects; i++) {
        cin >> a >> b >> c;

        Point light(a, b, c);
        lights.push_back(light);
    }

    temp = new Floor(1000, 20);
    temp->setCoEfficients(0.4, 0.2, 0.2, 0.2);
    temp->setShine(1);
    objects.push_back(temp);
}

int main(int argc, char **argv) {
    //Ray Tracing
    //loadTestData();
    loadActualData();

    glutInit(&argc, argv);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); //Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST); //enable Depth Testing

    glutDisplayFunc(display); //display callback function
    glutIdleFunc(animate);    //what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop(); //The main loop of OpenGL

    return 0;
}
