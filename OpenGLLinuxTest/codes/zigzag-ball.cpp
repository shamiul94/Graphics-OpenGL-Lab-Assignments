#include<bits/stdc++.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

using namespace std;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double cylinderHeight = 20, cylinderRadius = 30, sphereRadius = 20;
double movementOnZaxes = 40, sphereZaxisMovement = 40, movementOnYaxes = 40, movementOnXaxes = 40;
double distanceCoveredPerRotation, distanceCoveredPerRotationX, distanceCoveredPerRotationY;
double distanceCoveredPerRotationZ, angleY, angleZ, angleChange, period = 60;

struct Point {
    double x, y, z;

    Point() {}

    Point(double X, double Y, double Z) {
        x = X;
        y = Y;
        z = Z;
    }

    Point(const Point &p) : x(p.x), y(p.y), z(p.z) {}
};

class Position {
public:
    double dx, dy, dz, dtheta;

    Position() {}

    Position(double x, double y, double z, double theta) : dx(x), dy(y), dz(z), dtheta(theta) {}
};

Position movement, currPos;

Point crossProduct(const Point &a, const Point &b) {
    Point ret;
    ret.x = a.y * b.z - a.z * b.y;
    ret.y = -(a.x * b.z - a.z * b.x);
    ret.z = a.x * b.y - a.y * b.x;
    return ret;
}

Point cameraPos, cameraUp, cameraLook, cameraRight;

double degreesToRadians(double angle_in_degrees) {
    return angle_in_degrees * (pi / 180.0);
}


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
    glColor3f(0.6, 0.6, 0.6);    //grey
    glBegin(GL_LINES);
    {
        for (i = -28; i <= 28; i++) {

//            if (i == 0)
//                continue;    //SKIP the MAIN axes

            //lines parallel to Y-axis
            glVertex3f(i * 10, -300, 0);
            glVertex3f(i * 10, 300, 0);

            //lines parallel to X-axis
            glVertex3f(-300, i * 10, 0);
            glVertex3f(300, i * 10, 0);
        }
    }
    glEnd();
}


void drawSphere(double radius, int slices, int stacks) {
    struct Point points[100][100];
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
    for (i = 0; i <= stacks; i += 2) {
//        glColor3f((double) i / (double) stacks, (double) i / (double) stacks, (double) i / (double) stacks);
        for (j = 0; j <= slices; j += 2) {

            glPushMatrix();
            {
                glColor3f(1, 1, 1);

                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                    glVertex3f(points[i + 2][j].x, points[i + 2][j].y, points[i + 2][j].z);
                    glVertex3f(points[i + 2][j + 2].x, points[i + 2][j + 2].y, points[i + 2][j + 2].z);
                    glVertex3f(points[i][j + 2].x, points[i][j + 2].y, points[i][j + 2].z);

                    //lower hemisphere

                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                    glVertex3f(points[i + 2][j].x, points[i + 2][j].y, -points[i + 2][j].z);
                    glVertex3f(points[i + 2][j + 2].x, points[i + 2][j + 2].y, -points[i + 2][j + 2].z);
                    glVertex3f(points[i][j + 2].x, points[i][j + 2].y, -points[i][j + 2].z);
                }
                glEnd();
            }
            glPopMatrix();

            glPushMatrix();
            {
                glColor3f(1, 0, 0);
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j + 2].x, points[i + 1][j + 2].y, points[i + 1][j + 2].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);


                    glVertex3f(points[i + 1][j + 2].x, points[i + 1][j + 2].y, points[i + 1][j + 2].z);
                    glVertex3f(points[i + 2][j + 1].x, points[i + 2][j + 1].y, points[i + 2][j + 1].z);
                    glVertex3f(points[i + 2][j].x, points[i + 2][j].y, points[i + 2][j].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);

                    //lower hemisphere

                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j + 2].x, points[i + 1][j + 2].y, -points[i + 1][j + 2].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);


                    glVertex3f(points[i + 1][j + 2].x, points[i + 1][j + 2].y, -points[i + 1][j + 2].z);
                    glVertex3f(points[i + 2][j + 1].x, points[i + 2][j + 1].y, -points[i + 2][j + 1].z);
                    glVertex3f(points[i + 2][j].x, points[i + 2][j].y, -points[i + 2][j].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                }
                glEnd();
            }
            glPopMatrix();
        }
    }
}


void keyboardListener(unsigned char key, int x, int y) {

    switch (key) {

        case 'a':
            angleZ += (2 * 180 / period);
            break;

        case 'd':
            angleZ -= (2 * 180 / period);
            break;

        case 'w':
            movement.dx -= distanceCoveredPerRotation * cos(degreesToRadians(angleZ));
            movement.dy -= distanceCoveredPerRotation * sin(degreesToRadians(angleZ));
            angleY -= (2 * 180 / period);
            break;

        case 's':
            movement.dx += distanceCoveredPerRotation * cos(degreesToRadians(angleZ));
            movement.dy += distanceCoveredPerRotation * sin(degreesToRadians(angleZ));
            angleY += (2 * 180 / period);
            break;

        default:
            break;
    }
}


void specialKeyListener(int key, int x, int y) {
    double mov = 3;
    switch (key) {

        case GLUT_KEY_RIGHT:
            cameraAngle += .05;
            break;
        case GLUT_KEY_LEFT:
            cameraAngle -= .05;
            break;

        case GLUT_KEY_UP:
            cameraHeight += 5;
            break;
        case GLUT_KEY_DOWN:
            cameraHeight -= 5;
            break;

        case GLUT_KEY_HOME:

            break;
        case GLUT_KEY_END:
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


    gluLookAt(200 * cos(cameraAngle), 200 * sin(cameraAngle), cameraHeight, 0, 0, 0, 0, 0, 1);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    glColor3f(1, 0, 0);
    drawAxes();
//    drawGrid();

//    drawWheel();
    glRotatef(90, 0, 1, 0);
    drawSphere(100, 50, 50);
    glutSwapBuffers();
}


void animate() {
    angle += 0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init() {
    //codes for initialization
    drawgrid = 1;
    drawaxes = 1;
    cameraHeight = 150.0;
    cameraAngle = 1.0;
    angle = 0;

    /*** position movement initialization ***/

    movement = Position(0, 0, 0, 0);
    currPos = Position(0, 0, 0, 0);
    distanceCoveredPerRotation = 2 * pi * cylinderRadius / period;
    angleZ = 0;
    angleY = 0;
    angleChange = 0;

    /** Camera initialization **/

    drawgrid = 0;
    drawaxes = 1;
    cameraHeight = 80;
    cameraAngle = pi / 4;

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
    gluPerspective(80, 1, 1, 1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
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