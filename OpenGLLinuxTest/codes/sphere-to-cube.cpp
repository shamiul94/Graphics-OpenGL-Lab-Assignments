#include<bits/stdc++.h>

#include <GL/glut.h>

#define pi (2*acos(0.0))

using namespace std;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

double cylinderHeight = 80, cylinderRadius = 20, sphereRadius = 20;
double movementOnZaxes = 40, sphereZaxisMovement = 40, movementOnYaxes = 40, movementOnXaxes = 40;

struct Vector {
    double x, y, z;

    Vector() {}

    Vector(double X, double Y, double Z) {
        x = X;
        y = Y;
        z = Z;
    }

    Vector(const Vector &p) : x(p.x), y(p.y), z(p.z) {}
};

Vector crossProduct(const Vector &a, const Vector &b) {
    Vector ret;
    ret.x = a.y * b.z - a.z * b.y;
    ret.y = -(a.x * b.z - a.z * b.x);
    ret.z = a.x * b.y - a.y * b.x;
    return ret;
}

Vector cameraPos, cameraUp, cameraLook, cameraRight;

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


void drawOctaSphere(double radius, int slices, int stacks) {
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
        glColor3f(1, 0, 0);
        for (j = 0; j < slices / 4; j++) {
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
            }
            glEnd();
        }
    }
}


void drawCylinder(double height, double radius, int slices) {
    struct point {
        double x, y, z;
    };

    int i;
    struct point pointsUp[1000];
    struct point pointsDown[1000];
    glColor3f(0.7, 0.7, 0.7);


    //generate points



    for (i = 0; i <= slices; i++) {
        pointsUp[i].x = pointsDown[i].x = radius * cos(((double) i / (double) slices) * 2 * pi);
        pointsUp[i].y = pointsDown[i].y = radius * sin(((double) i / (double) slices) * 2 * pi);
        pointsUp[i].z = height / 2;
        pointsDown[i].z = -1.0 * pointsUp[i].z;
    }
    //draw segments using generated points
    for (i = 0; i < slices / 4; i++) {
        glBegin(GL_LINES);
        {
            glVertex3f(static_cast<GLfloat>(pointsUp[i].x), static_cast<GLfloat>(pointsUp[i].y),
                       static_cast<GLfloat>(pointsUp[i].z));
            glVertex3f(static_cast<GLfloat>(pointsUp[i + 1].x), static_cast<GLfloat>(pointsUp[i + 1].y),
                       static_cast<GLfloat>(pointsUp[i + 1].z));

            glVertex3f(static_cast<GLfloat>(pointsDown[i].x), static_cast<GLfloat>(pointsDown[i].y),
                       static_cast<GLfloat>(pointsDown[i].z));
            glVertex3f(static_cast<GLfloat>(pointsDown[i + 1].x), static_cast<GLfloat>(pointsDown[i + 1].y),
                       static_cast<GLfloat>(pointsDown[i + 1].z));
        }
        glEnd();


        glColor3f(0, 1, 0);
        glBegin(GL_QUADS);
        {
            glVertex3f(static_cast<GLfloat>(pointsUp[i].x), static_cast<GLfloat>(pointsUp[i].y),
                       static_cast<GLfloat>(pointsUp[i].z));
            glVertex3f(static_cast<GLfloat>(pointsDown[i].x), static_cast<GLfloat>(pointsDown[i].y),
                       static_cast<GLfloat>(pointsDown[i].z));
            glVertex3f(static_cast<GLfloat>(pointsDown[i + 1].x), static_cast<GLfloat>(pointsDown[i + 1].y),
                       static_cast<GLfloat>(pointsDown[i + 1].z));
            glVertex3f(static_cast<GLfloat>(pointsUp[i + 1].x), static_cast<GLfloat>(pointsUp[i + 1].y),
                       static_cast<GLfloat>(pointsUp[i + 1].z));
        }
        glEnd();
    }
}


void drawObjectSquare() {
    glBegin(GL_QUADS);
    {
        glVertex3f(static_cast<GLfloat>(cylinderHeight / 2), 0, static_cast<GLfloat>(cylinderHeight / 2));
        glVertex3f(static_cast<GLfloat>(cylinderHeight / 2), 0, static_cast<GLfloat>(-cylinderHeight / 2));
        glVertex3f(static_cast<GLfloat>(-cylinderHeight / 2), 0, static_cast<GLfloat>(-cylinderHeight / 2));
        glVertex3f(static_cast<GLfloat>(-cylinderHeight / 2), 0, static_cast<GLfloat>(cylinderHeight / 2));
    }
    glEnd();
}


void drawShell() {
    int sphereSlices = 80, sphereStacks = 80, cylinderSlices = 20;

    glPushMatrix();
    {
        glTranslatef(static_cast<GLfloat>(movementOnXaxes), static_cast<GLfloat>(movementOnYaxes),
                     static_cast<GLfloat>((sphereZaxisMovement / 2) + cylinderRadius));
        drawOctaSphere(sphereRadius, sphereSlices, sphereStacks);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(static_cast<GLfloat>(-movementOnXaxes), static_cast<GLfloat>(movementOnYaxes),
                     static_cast<GLfloat>((sphereZaxisMovement / 2) + cylinderRadius));
        glRotatef(90, 0, 0, 1);
        drawOctaSphere(sphereRadius, sphereSlices, sphereStacks);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(static_cast<GLfloat>(movementOnXaxes), static_cast<GLfloat>(movementOnYaxes),
                     static_cast<GLfloat>(-((sphereZaxisMovement / 2) + cylinderRadius)));
        glRotatef(-90, 1, 0, 0);
        drawOctaSphere(sphereRadius, sphereSlices, sphereStacks);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(static_cast<GLfloat>(-movementOnXaxes), static_cast<GLfloat>(movementOnYaxes),
                     static_cast<GLfloat>(-((sphereZaxisMovement / 2) + cylinderRadius)));
        glRotatef(-90, 1, 0, 0);
        glRotatef(90, 0, 0, 1);
        drawOctaSphere(sphereRadius, sphereSlices, sphereStacks);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(static_cast<GLfloat>(movementOnXaxes), static_cast<GLfloat>(movementOnYaxes), 0);
        drawCylinder(cylinderHeight, cylinderRadius, cylinderSlices);
    }
    glPopMatrix();


    glPushMatrix();
    {
        glTranslatef(static_cast<GLfloat>(-1.0 * movementOnXaxes), static_cast<GLfloat>(movementOnYaxes), 0);
        glRotatef(90, 0, 0, 1);
        drawCylinder(cylinderHeight, cylinderRadius, cylinderSlices);
    }
    glPopMatrix();


    glPushMatrix();
    {
        glTranslatef(0, 0, static_cast<GLfloat>(movementOnZaxes));
        glRotatef(90, 0, 1, 0);
        glTranslatef(0, static_cast<GLfloat>(movementOnYaxes), 0);
        glRotatef(90, 0, 0, 1);
        drawCylinder(cylinderHeight, cylinderRadius, cylinderSlices);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(0, 0, static_cast<GLfloat>(-movementOnZaxes));
        glRotatef(90, 0, 1, 0);
        glTranslatef(0, static_cast<GLfloat>(movementOnYaxes), 0);
        drawCylinder(cylinderHeight, cylinderRadius, cylinderSlices);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(0, static_cast<GLfloat>(cylinderHeight / 2 + cylinderRadius), 0);
        glColor3f(1, 1, 1);
        drawObjectSquare();
    }
    glPopMatrix();

    glPushMatrix();
    {
        glRotatef(90, 1, 0, 0);
        glTranslatef(0, static_cast<GLfloat>(cylinderHeight / 2 + cylinderRadius), 0);
        glColor3f(1, 1, 1);
        drawObjectSquare();
    }
    glPopMatrix();

}

void drawFinalObject() {
    glPushMatrix();
    {
        drawShell();
    }
    glPopMatrix();

    glPushMatrix();
    {
        glRotatef(90, 0, 0, 1);
        drawShell();
    }
    glPopMatrix();

    glPushMatrix();
    {
        glRotatef(-90, 0, 0, 1);
        drawShell();
    }
    glPopMatrix();

    glPushMatrix();
    {
        glRotatef(180, 1, 0, 0);
        drawShell();
    }
    glPopMatrix();
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

            cameraLook = crossProduct(cameraUp, cameraRight);
//                cameraUp.print();
            break;
        case '2':
            cameraRight.x = r.x * cos(ang) + l.x * sin(ang);
            cameraRight.y = r.y * cos(ang) + l.y * sin(ang);
            cameraRight.z = r.z * cos(ang) + l.z * sin(ang);

            cameraLook = crossProduct(cameraUp, cameraRight);
            break;
        case '3':       //rotate up. rotate l and u vectors WRT r.
            cameraLook.x = l.x * cos(ang) + u.x * sin(ang);
            cameraLook.y = l.y * cos(ang) + u.y * sin(ang);
            cameraLook.z = l.z * cos(ang) + u.z * sin(ang);

            cameraUp = crossProduct(cameraRight, cameraLook);
            break;
        case '4':       //rotate down. rotate l and u vectors WRT r.
            ang *= -1;
            cameraLook.x = l.x * cos(ang) + u.x * sin(ang);
            cameraLook.y = l.y * cos(ang) + u.y * sin(ang);
            cameraLook.z = l.z * cos(ang) + u.z * sin(ang);

            cameraUp = crossProduct(cameraRight, cameraLook);
            break;

        case '5':       //tilt camera clockwise.rotate r and u WRT l.
            cameraUp.x = u.x * cos(ang) + r.x * sin(ang);
            cameraUp.y = u.y * cos(ang) + r.y * sin(ang);
            cameraUp.z = u.z * cos(ang) + r.z * sin(ang);

            cameraRight = crossProduct(cameraLook, cameraUp);

            break;
        case '6':
            ang *= -1;
            cameraUp.x = u.x * cos(ang) + r.x * sin(ang);
            cameraUp.y = u.y * cos(ang) + r.y * sin(ang);
            cameraUp.z = u.z * cos(ang) + r.z * sin(ang);

            cameraRight = crossProduct(cameraLook, cameraUp);
            break;


        case 'r':
            cylinderHeight = 80, cylinderRadius = 20, sphereRadius = 20;
            movementOnXaxes = 40, movementOnYaxes = 40, movementOnZaxes = 40, sphereZaxisMovement = 40;
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

        case GLUT_KEY_HOME:
            sphereRadius += 2.0;
            cylinderRadius = sphereRadius;

            cylinderHeight -= 4.0;


            if (sphereRadius > 60) {
                cylinderHeight = 0;
                cylinderRadius = sphereRadius = 60;
            } else {
                sphereZaxisMovement -= 8.0;
            }

            movementOnXaxes = cylinderHeight / 2;
            movementOnYaxes = cylinderHeight / 2;
            movementOnZaxes = cylinderHeight / 2;

            break;
        case GLUT_KEY_END:
            sphereRadius -= 2.0;
            cylinderRadius -= 2.0;

            cylinderHeight += 4.0;

            if (sphereRadius < 0) {
                cylinderHeight = 120;
                cylinderRadius = sphereRadius = 0;
            } else {
                sphereZaxisMovement += 8.0;
            }


            movementOnXaxes = cylinderHeight / 2;
            movementOnYaxes = cylinderHeight / 2;
            movementOnZaxes = cylinderHeight / 2;
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

//    gluLookAt(200 * cos(cameraAngle), 200 * sin(cameraAngle), cameraHeight, 0, 0, 0, 0, 0, 1);

    gluLookAt(cameraPos.x, cameraPos.y, cameraPos.z,
              cameraPos.x + 10 * cameraLook.x, cameraPos.y + 10 * cameraLook.y, cameraPos.z + 10 * cameraLook.z,
              cameraUp.x, cameraUp.y, cameraUp.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();

    glColor3f(1, 0, 0);

    drawFinalObject();

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
//    cameraHeight=0;
    cameraAngle = pi / 4;;
//    cameraAngle=0;

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