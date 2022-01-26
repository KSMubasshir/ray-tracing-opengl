#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<GL/glut.h>
#include<bits/stdc++.h>
#include "bitmap_image.hpp"
#include "classes.h"
#define pi (2*acos(0.0))
#define ita 0.5
#define Window_width 500
#define Window_height 500
using namespace std;
inline double degreeToRadian(int degree)
{
    return degree * pi/180;
}
Point pos,l,r,u;
//camera movement functions
void moveCameraForward(int unit)
{
    pos.x+= (unit*l.x);
    pos.y+= (unit*l.y);
    pos.z+= (unit*l.z);
}
void moveCameraBackward(int unit)
{
    pos.x-= (unit*l.x);
    pos.y-= (unit*l.y);
    pos.z-= (unit*l.z);
}
void moveCameraRight(int unit)
{
    pos.x+= (unit*r.x);
    pos.y+= (unit*r.y);
    pos.z+= (unit*r.z);
}
void moveCameraLeft(int unit)
{
    pos.x-= (unit*r.x);
    pos.y-= (unit*r.y);
    pos.z-= (unit*r.z);
}
void moveCameraUp(int unit)
{
    pos.x+= (unit*u.x);
    pos.y+= (unit*u.y);
    pos.z+= (unit*u.z);
}
void moveCameraDown(int unit)
{
    pos.x-= (unit*u.x);
    pos.y-= (unit*u.y);
    pos.z-= (unit*u.z);
}
//camera rotation functions
void rotateCameraLeft(int unit)
{
    double cosA = cos(degreeToRadian(unit));
    double sinA = sin(degreeToRadian(unit));
    r.x = (r.x * cosA) + (l.x * sinA) ;
    r.y = (r.y * cosA) + (l.y * sinA) ;
    r.z = (r.z * cosA) + (l.z * sinA) ;
    l.x = u.y * r.z - r.y * u.z;
    l.y = r.x * u.z - u.x * r.z;
    l.z = u.x * r.y - r.x * u.y;
}
void rotateCameraRight(int unit)
{
    double cosA = cos(degreeToRadian(unit));
    double sinA = sin(degreeToRadian(unit));
    l.x = (l.x * cosA) + (r.x * sinA) ;
    l.y = (l.y * cosA) + (r.y * sinA) ;
    l.z = (l.z * cosA) + (r.z * sinA) ;
    r.x = l.y * u.z - u.y * l.z;
    r.y = u.x * l.z - l.x * u.z;
    r.z = l.x * u.y - u.x * l.y;
}
void rotateCameraUp(int unit)
{
    double cosA = cos(degreeToRadian(unit));
    double sinA = sin(degreeToRadian(unit));
    l.x = (l.x * cosA) + (u.x * sinA) ;
    l.y = (l.y * cosA) + (u.y * sinA) ;
    l.z = (l.z * cosA) + (u.z * sinA) ;
    u.x = r.y * l.z - l.y * r.z;
    u.y = l.x * r.z - r.x * l.z;
    u.z = r.x * l.y - l.x * r.y;
}
void rotateCameraDown(int unit)
{
    double cosA = cos(degreeToRadian(unit));
    double sinA = sin(degreeToRadian(unit));
    u.x = (u.x * cosA) + (l.x * sinA) ;
    u.y = (u.y * cosA) + (l.y * sinA) ;
    u.z = (u.z * cosA) + (l.z * sinA) ;
    l.x = u.y * r.z - r.y * u.z;
    l.y = r.x * u.z - u.x * r.z;
    l.z = u.x * r.y - r.x * u.y;
}
void tiltCameraClockwise(int unit)
{
    double cosA = cos(degreeToRadian(unit));
    double sinA = sin(degreeToRadian(unit));
    u.x = (u.x * cosA) + (r.x * sinA) ;
    u.y = (u.y * cosA) + (r.y * sinA) ;
    u.z = (u.z * cosA) + (r.z * sinA) ;
    r.x = l.y * u.z - u.y * l.z;
    r.y = u.x * l.z - l.x * u.z;
    r.z = l.x * u.y - u.x * l.y;
}
void tiltCameraAnticlockwise(int unit)
{
    double cosA = cos(degreeToRadian(unit));
    double sinA = sin(degreeToRadian(unit));
    r.x = (r.x * cosA) + (u.x * sinA) ;
    r.y = (r.y * cosA) + (u.y * sinA) ;
    r.z = (r.z * cosA) + (u.z * sinA) ;
    u.x = r.y * l.z - l.y * r.z;
    u.y = l.x * r.z - r.x * l.z;
    u.z = r.x * l.y - l.x * r.y;
}

void drawLightSource(Light pt)
{
    glColor3f(pt.color[0], pt.color[1], pt.color[2]);
    glBegin(GL_QUADS);
    {
        glVertex3f(pt.x + th, pt.y, pt.z + th);
        glVertex3f(pt.x + th, pt.y, pt.z - th);
        glVertex3f(pt.x - th, pt.y, pt.z - th);
        glVertex3f(pt.x - th, pt.y, pt.z + th);
    }
    glEnd();
}

void capture()
{
    double planeDistance, du, dv, mod;
    Point topLeft;
    Point** frameBuf;
    frameBuf = new point * [imageWidth];
    for (int i = 0; i < imageWidth; i++)
    {
        frameBuf[i] = new point[imageHeight];
        for (int j = 0; j < imageHeight; j++)
        {
            Point color;
            color.x = 0 ;
            color.y = 0 ;
            color.z = 0 ;
            frameBuf[i][j] = color;
        }
    }

    planeDistance = (Window_height / 2) / tan(fov * pi / 360);

    topLeft.x = pos.x + (l.x * planeDistance - r.x * (Window_width / 2) + u.x * (Window_height / 2));
    topLeft.y = pos.y + (l.y * planeDistance - r.y * (Window_width / 2) + u.y * (Window_height / 2));
    topLeft.z = pos.z + (l.z * planeDistance - r.z * (Window_width / 2) + u.z * (Window_height / 2));
    du = Window_width * 1.0 / imageWidth, dv = Window_height * 1.0 / imageHeight;
    int aa = 0;

    for (int i = 0; i < imageWidth; i++)
    {
        for (int j = 0; j < imageHeight; j++)
        {

            Point corner;
            corner.x = topLeft.x + r.x * j * du - u.x * i * dv;
            corner.y = topLeft.y + r.y * j * du - u.y * i * dv;
            corner.z = topLeft.z + r.z * j * du - u.z * i * dv;

            Point temp;
            temp.x = corner.x - pos.x;
            temp.y = corner.y - pos.y;
            temp.z = corner.z - pos.z;

            Ray ray;
            ray.start.x = pos.x;
            ray.start.y = pos.y;
            ray.start.z = pos.z;
            ray.dir.x = temp.x ;
            ray.dir.y = temp.y ;
            ray.dir.z = temp.z ;

            mod = sqrt(ray.dir.x*ray.dir.x + ray.dir.y*ray.dir.y + ray.dir.z*ray.dir.z);
            ray.dir.x /= mod ;
            ray.dir.y /= mod ;
            ray.dir.z /= mod ;


            int nearest = -1;
            double mint = 99999999;
            double ncolor[3];

            for (int k = 0; k < objects.size(); k++)
            {

                double t = objects[k]->intersect(&ray, ncolor, 0);

                if (t <= 0)
                    continue;

                else if (min(mint, t) == t)
                {
                    mint = t;
                    nearest = k;
                }
            }

            if (nearest != -1)
            {

                double t = objects[nearest]->intersect(&ray, ncolor, 1);

                frameBuf[i][j].x = ncolor[0];
                frameBuf[i][j].y = ncolor[1];
                frameBuf[i][j].z = ncolor[2];
            }
        }

    }


    bitmap_image img(imageWidth, imageHeight);
    for (int i = 0; i < imageWidth; i++)
    {
        for (int j = 0; j < imageHeight; j++)
        {
            img.set_pixel(j, i, frameBuf[i][j].x * 255, frameBuf[i][j].y * 255, frameBuf[i][j].z * 255);
        }
    }
    img.save_image("out.bmp");
    cout << "Image Captured" << endl;
}

void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {
    case '1':
        rotateCameraLeft(rotateUnit);
        break;
    case '2':
        rotateCameraRight(rotateUnit);
        break;
    case '3':
        rotateCameraUp(rotateUnit);
        break;
    case '4':
        rotateCameraDown(rotateUnit);
        break;
    case '5':
        tiltCameraClockwise(rotateUnit);
        break;
    case '6':
        tiltCameraAnticlockwise(rotateUnit);
        break;

    case '0':
        capture();
        break;
    default:
        break;
    }
}

void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_DOWN:
        moveCameraBackward(movementUnit);
        break;
    case GLUT_KEY_UP:
        moveCameraForward(movementUnit);
        break;

    case GLUT_KEY_RIGHT:
        moveCameraRight(movementUnit);
        break;
    case GLUT_KEY_LEFT:
        moveCameraLeft(movementUnit);
        break;

    case GLUT_KEY_PAGE_UP:
        moveCameraUp(movementUnit);
        break;
    case GLUT_KEY_PAGE_DOWN:
        moveCameraDown(movementUnit);
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

void display()
{
    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0);	//color black
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
    gluLookAt(pos.x,pos.y,pos.z,	pos.x+l.x,pos.y+l.y,pos.z+l.z,	u.x,u.y,u.z);
    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);
    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    for (int i = 0; i < objects.size(); i++)
    {
        objects[i]->draw();
    }
    for (int i = 0; i < lights.size(); i++)
    {
        drawLightSource(lights[i]);
    }
    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate()
{
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    //codes for initialization
    th = 0.5;
    ang = pi / 60;
    fov = 80;
    slices = 24;
    stacks = 20;
    movementUnit = 2*5;
    rotateUnit=3*5;
    pos.x = 100;
    pos.y = 100;
    pos.z = 20;
    u.x = 0;
    u.y = 0;
    u.z=1;
    r.x = -1/sqrt(2);
    r.y = -r.x;
    r.z=0;
    l.x = r.x;
    l.y = l.x;
    l.z =0;

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

void loadData()
{
    freopen("scene.txt", "r", stdin);
    int n;
    cin >> recursion_level >> imageWidth >> n;
    imageHeight = imageWidth;
    string obj;
    double a, b, c, d, radius, e, f, g, h, i, j;
    double x,y,z,baseLength,height;
    Object* temp;
    for (int i = 0; i < n; i++)
    {
        cin >> obj;
        if (obj == "sphere")
        {
            cin >> a >> b >> c >> radius >> d >> e >> f;
            temp = new Sphere;
            temp->center.x = a ;
            temp->center.y = b ;
            temp->center.z = c ;
            temp->radius = radius ;
            temp->color[0] = d ;
            temp->color[1] = e ;
            temp->color[2] = f ;

            cin >> a >> b >> c >> d >> g;
            temp->co_efficients[0] = a;
            temp->co_efficients[1] = b;
            temp->co_efficients[2] = c;
            temp->co_efficients[3] = d;

            temp->shine = g;
            //temp->generatePoints();
            objects.push_back(temp);
        }
        else if (obj == "triangle")
        {
            temp = new Triangle ;
            cin >> temp->A.x >> temp->A.y >> temp->A.z ;
            cin >> temp->B.x >> temp->B.y >> temp->B.z ;
            cin >> temp->C.x >> temp->C.y >> temp->C.z ;
            cin >> temp->color[0] >> temp->color[1]  >> temp->color[2]  ;
            cin >> temp->co_efficients[0] >> temp->co_efficients[1] >> temp->co_efficients[2] >> temp->co_efficients[3];
            cin >> temp->shine ;
            objects.push_back(temp);
        }
        else if (obj == "general")
        {
            cin >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j;
            cin >> a >> b >> c >> d >> e >> f;
            cin >> a >> b >> c;
            cin >> a >> b >> c >> d;
            cin >> h;
            // todo
        }
    }

    cin >> n;
    for (int i = 0; i < n; i++)
    {
        cin >> a >> b >> c;
        cin >> d >> e >> f; //todo
        Light light;
        light.x = a;
        light.y = b;
        light.z = c;
        light.color[0] = d;
        light.color[1] = e;
        light.color[2] = f;
        lights.push_back(light);
    }
    temp = new Plane(1000, 20);
    temp->co_efficients[0] = 0.4;
    temp->co_efficients[1] = 0.4;
    temp->co_efficients[2] = 0.2;
    temp->co_efficients[3] = 0.0;
    temp->shine = 5;
    objects.push_back(temp);
}

int main(int argc, char** argv)
{
    loadData();
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("My OpenGL Program");
    init();
    glEnable(GL_DEPTH_TEST);	//enable Depth Testing
    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();		//The main loop of OpenGL
    return 0;
}
