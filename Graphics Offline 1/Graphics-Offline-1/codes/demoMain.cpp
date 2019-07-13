#include<bits/stdc++.h>

#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

//double cylinderHeight = 80, cylinderRadius = 20, sphereRadius = 20;
double cylinderHeight = 80, cylinderRadius = 20, sphereRadius = 20;
double movementOnZaxes = 40, sphereZaxisMovement = 40;
GLfloat movementOnXaxes = 40;
GLfloat movementOnYaxes = 40;

struct point {
    double x, y, z;
};


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

void drawSquare(double a) {
    //glColor3f(1.0,0.0,0.0);
    glBegin(GL_QUADS);
    {
        glVertex3f(a, a, 2);
        glVertex3f(a, -a, 2);
        glVertex3f(-a, -a, 2);
        glVertex3f(-a, a, 2);
    }
    glEnd();
}


void drawCircle(double radius, int segments) {
    int i;
    struct point points[100];
    glColor3f(0.7, 0.7, 0.7);
    //generate points
    for (i = 0; i <= segments; i++) {
        points[i].x = radius * cos(((double) i / (double) segments) * 2 * pi);
        points[i].y = radius * sin(((double) i / (double) segments) * 2 * pi);
    }
    //draw segments using generated points
    for (i = 0; i < segments; i++) {
        glBegin(GL_LINES);
        {
            glVertex3f(points[i].x, points[i].y, 0);
            glVertex3f(points[i + 1].x, points[i + 1].y, 0);
        }
        glEnd();
    }
}

void drawCone(double radius, double height, int segments) {
    int i;
    double shade;
    struct point points[100];
    //generate points
    for (i = 0; i <= segments; i++) {
        points[i].x = radius * cos(((double) i / (double) segments) * 2 * pi);
        points[i].y = radius * sin(((double) i / (double) segments) * 2 * pi);
    }
    //draw triangles using generated points
    for (i = 0; i < segments; i++) {
        //create shading effect
        if (i < segments / 2)
            shade = 2 * (double) i / (double) segments;
        else
            shade = 2 * (1.0 - (double) i / (double) segments);
        glColor3f(shade, shade, shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0, 0, height);
            glVertex3f(points[i].x, points[i].y, 0);
            glVertex3f(points[i + 1].x, points[i + 1].y, 0);
        }
        glEnd();
    }
}


void drawOctaSphere(double radius, int slices, int stacks) {
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
//        glColor3f((double) i / (double) stacks, (double) i / (double) stacks, (double) i / (double) stacks);
        glColor3f(1, 0, 0);
        for (j = 0; j < slices / 4; j++) {
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                //lower hemisphere
//                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
//                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
//                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
//                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
            }
            glEnd();
        }
    }
}


void drawCylinder(double height, double radius, int slices) {
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
            glVertex3f(pointsUp[i].x, pointsUp[i].y, pointsUp[i].z);
            glVertex3f(pointsUp[i + 1].x, pointsUp[i + 1].y, pointsUp[i + 1].z);

            glVertex3f(pointsDown[i].x, pointsDown[i].y, pointsDown[i].z);
            glVertex3f(pointsDown[i + 1].x, pointsDown[i + 1].y, pointsDown[i + 1].z);
        }
        glEnd();

//        glColor3f((double) i / (double) slices, (double) i / (double) slices, (double) i / (double) slices);
//        glColor3f(1, 0, 0);
//        glBegin(GL_TRIANGLES);
//        {
//            glVertex3f(pointsUp[i].x, pointsUp[i].y, pointsUp[i].z);
//            glVertex3f(0, 0, pointsUp[i].z);
//            glVertex3f(pointsUp[i + 1].x, pointsUp[i + 1].y, pointsUp[i + 1].z);
//
//            glVertex3f(pointsDown[i].x, pointsDown[i].y, pointsDown[i].z);
//            glVertex3f(0, 0, pointsDown[i].z);
//            glVertex3f(pointsDown[i + 1].x, pointsDown[i + 1].y, pointsDown[i + 1].z);
//        }
//        glEnd();

        glColor3f(0, 1, 0);
        glBegin(GL_QUADS);
        {
            glVertex3f(pointsUp[i].x, pointsUp[i].y, pointsUp[i].z);
            glVertex3f(pointsDown[i].x, pointsDown[i].y, pointsDown[i].z);
            glVertex3f(pointsDown[i + 1].x, pointsDown[i + 1].y, pointsDown[i + 1].z);
            glVertex3f(pointsUp[i + 1].x, pointsUp[i + 1].y, pointsUp[i + 1].z);
        }
        glEnd();
    }
}


void drawObjectSquare() {
    glBegin(GL_QUADS);
    {
        glVertex3f(cylinderHeight / 2, 0, cylinderHeight / 2);
        glVertex3f(cylinderHeight / 2, 0, -cylinderHeight / 2);
        glVertex3f(-cylinderHeight / 2, 0, -cylinderHeight / 2);
        glVertex3f(-cylinderHeight / 2, 0, cylinderHeight / 2);
    }
    glEnd();
}


void drawShell() {

    glPushMatrix();
    {
        glTranslatef(movementOnXaxes, movementOnYaxes, (sphereZaxisMovement / 2) + cylinderRadius);
        drawOctaSphere(sphereRadius, 24, 20);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(-movementOnXaxes, movementOnYaxes, (sphereZaxisMovement / 2) + cylinderRadius);
        glRotatef(90, 0, 0, 1);
        drawOctaSphere(sphereRadius, 24, 20);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(movementOnXaxes, movementOnYaxes, -((sphereZaxisMovement / 2) + cylinderRadius));
        glRotatef(-90, 1, 0, 0);
        drawOctaSphere(sphereRadius, 24, 20);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(-movementOnXaxes, movementOnYaxes, -((sphereZaxisMovement / 2) + cylinderRadius));
        glRotatef(-90, 1, 0, 0);
        glRotatef(90, 0, 0, 1);
        drawOctaSphere(sphereRadius, 24, 20);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(movementOnXaxes, movementOnYaxes, 0);
        drawCylinder(cylinderHeight, cylinderRadius, 80);
    }
    glPopMatrix();


    glPushMatrix();
    {
        glTranslatef(-1.0 * movementOnXaxes, movementOnYaxes, 0);
        glRotatef(90, 0, 0, 1);
        drawCylinder(cylinderHeight, cylinderRadius, 80);
    }
    glPopMatrix();


    glPushMatrix();
    {
        glTranslatef(0, 0, movementOnZaxes);
        glRotatef(90, 0, 1, 0);
        glTranslatef(0, movementOnYaxes, 0);
        glRotatef(90, 0, 0, 1);
        drawCylinder(cylinderHeight, cylinderRadius, 80);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(0, 0, -movementOnZaxes);
        glRotatef(90, 0, 1, 0);
        glTranslatef(0, movementOnYaxes, 0);
        drawCylinder(cylinderHeight, cylinderRadius, 80);
    }
    glPopMatrix();

    glPushMatrix();
    {
        glTranslatef(0, cylinderHeight / 2 + cylinderRadius, 0);
        glColor3f(1, 1, 1);
        drawObjectSquare();
    }
    glPopMatrix();

    glPushMatrix();
    {
        glRotatef(90, 1, 0, 0);
        glTranslatef(0, cylinderHeight / 2 + cylinderRadius, 0);
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


//SS = Solar System.

void drawSimpleSS() {
    glColor3f(1, 0, 0);
    drawSquare(20);

/*
 * Basically we need to have clear concept about glpushmatrix() and glpopmatrix()
 *
 * what happens is -
 *
 * 1. glpushmatrix() makes a copy of current matrix (as it is right now) and push it in the stack. So, the top element
 * and second top element in that stack are same. OpenGL keeps a stack of matrices to quickly
 * apply and remove transformations. glPushMatrix copies the top matrix
 * and pushes it onto the stack, while glPopMatrix pops the top matrix off the stack.
 * All transformation functions (glScaled, etc.) function on the top matrix, and
 * the top matrix is what all rendering commands use to transform their data.
 * By pushing and popping matrices, you can control what transformations apply to which objects,
 * as well as apply transformations to groups of objects, and easily
 * reverse the transformations so that they don't affect other objects.
 *
 * 2. Whenever we do a transform on any matrix, every thing (even the origin) get transformed. So, if you wanna
 * render next shape around new axis, you don't need any glpush or pop.
 *
 *
 */


    glPushMatrix();
    {
        glRotatef(angle, 0, 0, 1); // main axis er charpashe.
        glTranslatef(110, 0, 0);
        glRotatef(2 * angle, 0, 0, 1); //nijer axis er char pashe ghure.
        glColor3f(0, 1, 0);
//        drawSquare(15);
        drawOctaSphere(10, 10, 10);
    }
    glPopMatrix();

//    drawAxes();

    glPushMatrix();
    {
        glRotatef(angle, 0, 0, 1);
        glTranslatef(60, 0, 0);
        glRotatef(2 * angle, 0, 0, 1);
        glColor3f(0, 0, 1);
//        drawSquare(10);
        drawCylinder(50, 30, 80);
    }
    glPopMatrix();


    glPushMatrix();
    {
        glRotatef(3 * angle, 0, 0, 1);
        glTranslatef(40, 0, 0);
        glRotatef(4 * angle, 0, 0, 1);
        glColor3f(1, 1, 0);
        drawSquare(5);
//    drawAxes();
    }
    glPopMatrix();
}


void drawSS() {
    glColor3f(1, 0, 0);
    drawSquare(20);

//    glRotatef(angle, 0, 1, 0); // main axis er charpashe.
//    glRotatef(angle, 1, 0, 0); // main axis er charpashe.


    glRotatef(angle, 0, 0, 1); // main axis er charpashe.
    glTranslatef(110, 0, 0);
    glRotatef(2 * angle, 0, 0, 1); //nijer axis er char pashe ghure.
    glColor3f(0, 1, 0);
    drawSquare(15);

    drawAxes();

    glPushMatrix();
    {
        glRotatef(angle, 0, 0, 1);
        glTranslatef(60, 0, 0);
        glRotatef(2 * angle, 0, 0, 1);
        glColor3f(0, 0, 1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(angle, 0, 0, 1);
    glTranslatef(60, 0, 0);
    glRotatef(2 * angle, 0, 0, 1);
    glColor3f(0, 0, 1);
    drawSquare(10);
    drawAxes();


    glRotatef(3 * angle, 0, 0, 1);
    glTranslatef(40, 0, 0);
    glRotatef(4 * angle, 0, 0, 1);
    glColor3f(1, 1, 0);
    drawSquare(5);
    drawAxes();

}

void keyboardListener(unsigned char key, int x, int y) {
    switch (key) {

        case '1':
            drawgrid = 1 - drawgrid;
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
    switch (key) {
        case GLUT_KEY_DOWN:        //down arrow key
            cameraHeight -= 3.0;
            break;
        case GLUT_KEY_UP:        // up arrow key
            cameraHeight += 3.0;
            break;

        case GLUT_KEY_RIGHT:
            cameraAngle += 0.03;
            break;
        case GLUT_KEY_LEFT:
            cameraAngle -= 0.03;
            break;

        case GLUT_KEY_PAGE_UP:

            break;
        case GLUT_KEY_PAGE_DOWN:

            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            sphereRadius += 2.0;
            cylinderRadius = sphereRadius;

            cylinderHeight -= 4.0;


            if (sphereRadius > 60) {
                cylinderHeight = 0;
                cylinderRadius = sphereRadius = 60;
            }
            else {
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
            }
            else {
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

//    gluLookAt(100,100,100,	0,0,0,	0,0,1);
    gluLookAt(200 * cos(cameraAngle), 200 * sin(cameraAngle), cameraHeight, 0, 0, 0, 0, 0, 1);
//    gluLookAt(0, 0, 200, 0, 0, 0, 0, 1, 0);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
//    drawGrid();

    glColor3f(1, 0, 0);
//    drawSquare(10);

//    drawSS();
//    drawSimpleSS();

//    drawCircle(30, 24);

//    drawCone(20, 50, 24);

//    drawOctaSphere(30, 24, 20);

//    drawCylinder(50, 20, 20, 60);

    drawFinalObject();
//    drawObject();

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
    drawgrid = 1;
    drawaxes = 1;
    cameraHeight = 150.0;
    cameraAngle = 1.0;
    angle = 0;

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
