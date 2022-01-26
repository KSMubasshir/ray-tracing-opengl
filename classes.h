#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<vector>
#include<GL/glut.h>
#include<bits/stdc++.h>
#include "bitmap_image.hpp"
#define pi (2*acos(0.0))
#define ita 0.5
#define Window_width 500
#define Window_height 500
using namespace std;
class Object;
int slices;
int stacks;
int recursion_level;
int imageWidth;
int imageHeight;
int fov;
int movementUnit,rotateUnit;
double th;
double ang;
double left, up, right;

struct point
{
    double x;
    double y;
    double z;
};
struct light
{
    double x;
    double y;
    double z;
    double color[3];
};
typedef struct point Point;
typedef struct light Light;
vector <Light> lights;
vector <Object*> objects;
struct ray
{
    Point start;
    Point dir;
};
typedef ray Ray;

class Object
{
public:
    Point A,B,C;
    Point center;
    double source_factor = 1.0;
    double radius;
    double color[3];
    double co_efficients[4];
    int shine;
    virtual void draw()=0;
    //virtual void generatePoints()=0;
    virtual double intersect(Ray* ray,double ncolor[3],int level)=0;
    virtual double findClosest(Ray *ray)=0;
};

class Sphere:public Object
{
public:
    double height, width ;
    Point points[100][100];
    double h, r;

    void generatePoints()
    {
        for (int i = 0; i <= stacks; i++)
        {
            h = radius * sin(((double)i / (double)stacks) * (pi / 2));
            r = radius * cos(((double)i / (double)stacks) * (pi / 2));
            for (int j = 0; j <= slices; j++)
            {
                points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
                points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
                points[i][j].z = h;
            }
        }
    }

    void draw()
    {
        glColor3f(color[0], color[1], color[2]);
        for (int i = 0; i < stacks; i++)
        {
            h = radius * sin(((double)i / (double)stacks) * (pi / 2));
            r = radius * cos(((double)i / (double)stacks) * (pi / 2));
            for (int j = 0; j < slices; j++)
            {
                points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
                points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
                points[i][j].z = h;
                glBegin(GL_QUADS);
                {
                    glVertex3f(points[i][j].x + center.x, points[i][j].y + center.y, points[i][j].z + center.z);
                    glVertex3f(points[i][j + 1].x + center.x, points[i][j + 1].y + center.y, points[i][j + 1].z + center.z);
                    glVertex3f(points[i + 1][j + 1].x + center.x, points[i + 1][j + 1].y + center.y, points[i + 1][j + 1].z + center.z);
                    glVertex3f(points[i + 1][j].x + center.x, points[i + 1][j].y + center.y, points[i + 1][j].z + center.z);

                    glVertex3f(points[i][j].x + center.x, points[i][j].y + center.y, -points[i][j].z + center.z);
                    glVertex3f(points[i][j + 1].x + center.x, points[i][j + 1].y + center.y, -points[i][j + 1].z + center.z);
                    glVertex3f(points[i + 1][j + 1].x + center.x, points[i + 1][j + 1].y + center.y, -points[i + 1][j + 1].z + center.z);
                    glVertex3f(points[i + 1][j].x + center.x, points[i + 1][j].y + center.y, -points[i + 1][j].z + center.z);
                }
                glEnd();
            }
        }
    }

    double findClosest(Ray* r)
    {
        double t,mod;
        Point s;
        s.x = r->start.x - center.x;
        s.y = r->start.y - center.y;
        s.z = r->start.z - center.z;
        double a, b, c, d;
        a = r->dir.x * r->dir.x +  r->dir.y * r->dir.y + r->dir.z * r->dir.z ;
        b = 2 * (r->dir.x * s.x + r->dir.y * s.y + r->dir.z * s.z);
        c = (s.x * s.x) + (s.y * s.y) + (s.z * s.z)  - radius * radius;
        d = (b * b - 4 * a * c);
        if (d < 0)
            t = -1;
        else
        {
            d = sqrt(d);
            double root1 = (-b + d) / (2.0 * a);
            double root2 = (-b - d) / (2.0 * a);
            t = min(root1, root2);
        }
        return t;
    }
    double intersect(Ray* r, double ncolor[3], int level)
    {
        //finds closest t
        double t,mod;
        Point s;
        s.x = r->start.x - center.x;
        s.y = r->start.y - center.y;
        s.z = r->start.z - center.z;
        double a, b, c, d;

        a = r->dir.x * r->dir.x +  r->dir.y * r->dir.y + r->dir.z * r->dir.z ;
        b = 2 * (r->dir.x * s.x + r->dir.y * s.y + r->dir.z * s.z);
        c = (s.x * s.x) + (s.y * s.y) + (s.z * s.z)  - radius * radius;

        d = (b * b - 4 * a * c);
        if (d < 0)
            t = -1;
        else
            d = sqrt(d);

        double root1 = (-b + d) / (2.0 * a);
        double root2 = (-b - d) / (2.0 * a);

        t = min(root1, root2);

        //finds intersection
        if (t <= 0)
            return -1;
        if (level == 0)
            return t;

        for (int i = 0; i < 3; i++)
        {
            ncolor[i] = color[i] * co_efficients[0];
        }

        Point intersectionPoint;

        intersectionPoint.x = r->start.x + r->dir.x * t;
        intersectionPoint.y = r->start.y + r->dir.y * t;
        intersectionPoint.z = r->start.z + r->dir.z * t;

        Point normal ;
        normal.x = intersectionPoint.x - center.x;
        normal.y = intersectionPoint.y - center.y;
        normal.z = intersectionPoint.z - center.z;
        mod = sqrt( normal.x*normal.x + normal.y*normal.y + normal.z*normal.z );

        normal.x /= mod ;
        normal.y /= mod ;
        normal.z /= mod ;

        //finds reflection
        Point reflection;
        double dd = r->dir.x * normal.x + r->dir.y * normal.y + r->dir.z * normal.z;
        reflection.x = r->dir.x - 2.0 * normal.x * dd;
        reflection.y = r->dir.y - 2.0 * normal.y * dd;
        reflection.z = r->dir.z - 2.0 * normal.z * dd;

        mod = sqrt( reflection.x*reflection.x + reflection.y*reflection.y + reflection.z*reflection.z );

        reflection.x /= mod ;
        reflection.y /= mod ;
        reflection.z /= mod ;

        for (int i = 0; i < lights.size(); i++)
        {
            Point direction;
            direction.x = lights[i].x - intersectionPoint.x;
            direction.y = lights[i].y - intersectionPoint.y;
            direction.z = lights[i].z - intersectionPoint.z;

            double l = sqrt(pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2));
            direction.x /= l ;
            direction.y /= l ;
            direction.z /= l ;
            Point start;

            start.x = intersectionPoint.x + direction.x * 1;
            start.y = intersectionPoint.y + direction.y * 1;
            start.z = intersectionPoint.z + direction.z * 1;

            Ray L;
            L.start.x = start.x ;
            L.start.y = start.y ;
            L.start.z = start.z ;
            L.dir.x = direction.x ;
            L.dir.y = direction.y ;
            L.dir.z = direction.z ;

            int f = 1;

            for (int j = 0; j < objects.size(); j++)
            {
                double obs = objects[j]->findClosest(&L);
                if (obs > 0)
                {
                    f = 0;
                    break;
                }
                if ((obs) > l)
                {
                    f = 0;
                    break;
                }
            }

            //for floor

            if (f == 1)
            {
                double lambert, phong;
                lambert = (L.dir.x * normal.x) + (L.dir.y * normal.y) + (L.dir.z * normal.z);
                phong = (reflection.x*r->dir.x) + (reflection.y*r->dir.y) + (reflection.z*r->dir.z);
                phong = pow(phong, shine);

                if(lambert < 0)
                    lambert = 0;
                if(phong < 0)
                    phong = 0;

                for (int m = 0; m < 3; m++)
                {
                    ncolor[m] += source_factor * lambert * co_efficients[1] * color[m]  * lights[i].color[m];
                    ncolor[m] += source_factor * phong * co_efficients[2] * color[m]  * lights[i].color[m];
                }
            }

            if (level < recursion_level)
            {
                start.x = intersectionPoint.x + reflection.x * 1;
                start.y = intersectionPoint.y + reflection.y * 1;
                start.z = intersectionPoint.z + reflection.z * 1;

                Ray reflectionRay;

                reflectionRay.start.x = start.x ;
                reflectionRay.start.y = start.y ;
                reflectionRay.start.z = start.z ;
                reflectionRay.dir.x = reflection.x ;
                reflectionRay.dir.y = reflection.y ;
                reflectionRay.dir.z = reflection.z ;

                double mint = 99999;
                double refl_color[3];
                int nearest = -1;

                for (int k = 0; k < objects.size(); k++)
                {
                    double tt = objects[k]->findClosest(&reflectionRay);
                    if (tt <= 0)
                        continue;
                    else if(min(mint, tt) == tt)
                    {
                        mint = tt;
                        nearest = k;
                    }
                }

                if (nearest != -1)
                {
                    objects[nearest]->intersect(&reflectionRay, refl_color, ++level);
                    for (int k = 0; k < 3; k++)
                    {
                        ncolor[k] += refl_color[k] * co_efficients[3];
                    }
                }
            }

            for (int k = 0; k < 3; k++)
            {
                if (ncolor[k] < 0)
                    ncolor[k] = 0;
                else if (ncolor[k] > 1)
                    ncolor[k] = 1;
            }
        }
        return t;
    }
};

class Plane:public Object
{
public:
    double length ;
    Point ref_point;
    Point normal;

    Plane(double fw, double tw)
    {
        length = tw;
        ref_point.x =  -fw / 2;
        ref_point.y =  -fw / 2;
        ref_point.z =   0;
        normal.x = 0;
        normal.y = 0;
        normal.z = 1;
    }

    void draw()
    {
        int s = 0, n = abs(2 * ref_point.x / length);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if ((i + j) % 2 == 1)
                {
                    glColor3f(0, 0, 0);
                }
                else
                {
                    glColor3f(1, 1, 1);
                }
                glBegin(GL_QUADS);
                {
                    glVertex3f(ref_point.x + length * i, ref_point.y + length * j, ref_point.z);
                    glVertex3f(ref_point.x + length * (i + 1), ref_point.y + length * j, ref_point.z);
                    glVertex3f(ref_point.x + length * (i + 1), ref_point.y + length * (j + 1), ref_point.z);
                    glVertex3f(ref_point.x + length * i, ref_point.y + length * (j + 1), ref_point.z);
                }
                glEnd();
            }
        }
    }

    void generatePoints()
    {

    }

    double findClosest(Ray* ray)
    {
        double t;
        t = -((normal.x * ray->start.x) + (normal.y * ray->start.y) + (normal.z * ray->start.z))
            / ((normal.x * ray->dir.x) + (normal.y * ray->dir.y) + (normal.z * ray->dir.z));
        return t;
    }

    double intersect(Ray* r, double ncolor[3], int level)
    {
        double t;
        t = - ((normal.x * r->start.x) + (normal.y * r->start.y) + (normal.z * r->start.z)) / ((normal.x * r->dir.x) + (normal.y * r->dir.y) + (normal.z * r->dir.z));
        Point intersectionPoint;
        intersectionPoint.x = r->start.x + r->dir.x * t ;
        intersectionPoint.y = r->start.y + r->dir.y * t ;
        intersectionPoint.z = r->start.z + r->dir.z * t ;

        if (ref_point.x > intersectionPoint.x || -ref_point.x < intersectionPoint.x)
            return -1;
        if (ref_point.y > intersectionPoint.y || -ref_point.y < intersectionPoint.y)
            return -1;

        int xWhich = intersectionPoint.x / length;
        int yWhich = intersectionPoint.y / length;

        if ((xWhich + yWhich) % 2)
        {
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
        }
        else
        {
            color[0] = 1;
            color[1] = 1;
            color[2] = 1;
        }

        for (int i = 0; i < 3; i++)
        {
            ncolor[i] = color[i] * co_efficients[0];
            ncolor[i] /= 255.0;
        }

        Point reflection;
        double dd = r->dir.x * normal.x + r->dir.y * normal.y + r->dir.z * normal.z;
        reflection.x = r->dir.x - 2.0 * normal.x * dd;
        reflection.y = r->dir.y - 2.0 * normal.y * dd;
        reflection.z = r->dir.z - 2.0 * normal.z * dd;

        double mod = sqrt( reflection.x*reflection.x + reflection.y*reflection.y + reflection.z*reflection.z );

        reflection.x /= mod ;
        reflection.y /= mod ;
        reflection.z /= mod ;

        for (int i = 0; i < lights.size(); i++)
        {
            Point direction;
            direction.x = lights[i].x - intersectionPoint.x;
            direction.y = lights[i].y - intersectionPoint.y;
            direction.z = lights[i].z - intersectionPoint.z;
            double l = sqrt(pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2));
            direction.x /= l;
            direction.y /= l;
            direction.z /= l;

            Point start;
            start.x = intersectionPoint.x + direction.x * 1;
            start.y = intersectionPoint.y + direction.y * 1;
            start.z = intersectionPoint.z + direction.z * 1;

            Ray L;
            L.start.x = start.x ;
            L.start.y = start.y ;
            L.start.z = start.z ;
            L.dir.x = direction.x ;
            L.dir.y = direction.y ;
            L.dir.z = direction.z ;
            int f = 1;
            for (int j = 0; j < objects.size(); j++)
            {
                double obs = objects[j]->findClosest(&L);
                if (obs > 0)
                {
                    f = 0;
                    break;
                }
                if ((obs) > l)
                {
                    f = 0;
                    break;
                }
            }


            if (f == 1)
            {
                double lambert, phong;
                lambert = (L.dir.x * normal.x) + (L.dir.y * normal.y) + (L.dir.z * normal.z);
                phong = (reflection.x*r->dir.x) + (reflection.y*r->dir.y) + (reflection.z*r->dir.z);
                phong = pow(phong, shine);

                if (lambert < 0)
                    lambert = 0;
                if (phong < 0)
                    phong = 0;
                for (int m = 0; m < 3; m++)
                {
                    ncolor[m] += source_factor * lambert * co_efficients[1] * color[m] * lights[i].color[m];

                    ncolor[m] += source_factor * phong * co_efficients[2] * color[m] * lights[i].color[m];
                }
            }

            if (level < recursion_level)
            {

                start.x = intersectionPoint.x + reflection.x * 1;
                start.y = intersectionPoint.y + reflection.y * 1;
                start.z = intersectionPoint.z + reflection.z * 1;

                Ray reflectionRay;

                reflectionRay.start.x = start.x ;
                reflectionRay.start.y = start.y ;
                reflectionRay.start.z = start.z ;
                reflectionRay.dir.x = reflection.x ;
                reflectionRay.dir.y = reflection.y ;
                reflectionRay.dir.z = reflection.z ;

                double mint = 99999;
                double refl_color[3];
                int nearest = -1;
                for (int k = 0; k < objects.size(); k++)
                {
                    double tt = objects[k]->findClosest(&reflectionRay);
                    if (tt <= 0)
                        continue;
                    else if (min(mint, tt) == tt)
                    {
                        mint = tt;
                        nearest = k;
                    }
                }

                if (nearest != -1)
                {
                    objects[nearest]->intersect(&reflectionRay, refl_color, ++level);
                    for (int k = 0; k < 3; k++)
                    {
                        ncolor[k] += refl_color[k] * co_efficients[3];
                    }
                }
            }
            for (int k = 0; k < 3; k++)
            {
                if (ncolor[k] < 0)
                    ncolor[k] = 0;
                else if (ncolor[k] > 1)
                    ncolor[k] = 1;
            }
        }
        return t;
    }
};

class Triangle:public Object
{
public:
    Point normal ;
    Triangle(){
        double mod;
        Point p1;
        p1.x = B.x-A.x;
        p1.y = B.y-A.y;
        p1.z = B.z-A.z;
        Point p2;
        p2.x = C.x-A.x;
        p2.y = C.y-A.y;
        p2.z = C.z-A.z;

        normal.x = p1.y * p2.z - p2.y * p1.z;
        normal.y = p2.x * p1.z - p1.x * p2.z;
        normal.z = p1.x * p2.y - p2.x * p1.y;
        mod = sqrt((normal.x*normal.x)+(normal.y*normal.y)+(normal.z*normal.z));
        normal.x /= mod ;
        normal.y /= mod ;
        normal.z /= mod ;
    }

    void draw()
    {
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(A.x,A.y,A.z);
            glVertex3f(B.x,B.y,B.z);
            glVertex3f(C.x,C.y,C.z);
        }
        glEnd();
    }



    double findClosest(Ray *ray)
    {
        float EPSILON=0.0000001;
        float aa,f,u,v;
        Point edge1;
        edge1.x = B.x-A.x ;
        edge1.y = B.y-A.y ;
        edge1.z = B.z-A.z ;

        Point edge2;
        edge2.x = C.x-A.x;
        edge2.y = C.y-A.y;
        edge2.z = C.z-A.z;

        Point h;
        h.x = ray->dir.y * edge2.z - edge2.y * ray->dir.z;
        h.y = edge2.x * ray->dir.z - ray->dir.x * edge2.z;
        h.z = ray->dir.x * edge2.y - edge2.x * ray->dir.y;

        aa = (edge1.x*h.x) + (edge1.y*h.y) + (edge1.z*h.z) ;
        if (aa > -EPSILON && aa < EPSILON)
            return -1;
        f = 1.0/aa;
        Point s;
        s.x = ray->start.x - A.x;
        s.y = ray->start.y - A.y;
        s.z = ray->start.z - A.z;

        u = f * ((s.x*h.x)+(s.y*h.y)+(s.z*h.z)) ;
        if (u < 0.0 || u > 1.0)
            return -1;
        Point q ;
        q.x = s.y * edge1.z - edge1.y * s.z ;
        q.y = edge1.x * s.z - s.x * edge1.z ;
        q.z = s.x * edge1.y - edge1.x * s.y ;

        v = f * ((ray->dir.x * q.x) + (ray->dir.y * q.y) + (ray->dir.z * q.z));
        if (v < 0.0 || u + v > 1.0)
            return -1;
        float t ;
        t = f * ((edge2.x*q.x)+(edge2.y*q.y)+(edge2.z*q.z));
        if (t<= EPSILON) return -1;
        else return t;
    }


    double intersect(Ray* r, double ncolor[3], int level)
    {
        double t=findClosest(r);
        if(t<=0) return -1;
        if(level==0) return t;

        for(int i=0; i<3; i++)
            ncolor[i]=color[i]*co_efficients[0];

        Point intersectionPoint;
        intersectionPoint.x = r->start.x+ r->dir.x*t ;
        intersectionPoint.y = r->start.y+ r->dir.y*t;
        intersectionPoint.z = r->start.z+ r->dir.z*t;

        Point normal;

        Point reflection;
        double dd = r->dir.x * normal.x + r->dir.y * normal.y + r->dir.z * normal.z;
        reflection.x = r->dir.x - 2.0 * normal.x * dd;
        reflection.y = r->dir.y - 2.0 * normal.y * dd;
        reflection.z = r->dir.z - 2.0 * normal.z * dd;

        double mod = sqrt( reflection.x*reflection.x + reflection.y*reflection.y + reflection.z*reflection.z );

        reflection.x /= mod ;
        reflection.y /= mod ;
        reflection.z /= mod ;

        for(int i=0; i<lights.size(); i++)
        {
            Point direction;
            direction.x = lights[i].x - intersectionPoint.x;
            direction.y = lights[i].y - intersectionPoint.y;
            direction.z = lights[i].z - intersectionPoint.z;
            double l = sqrt(pow(direction.x, 2) + pow(direction.y, 2) + pow(direction.z, 2));
            direction.x /= l;
            direction.y /= l;
            direction.z /= l;
            Point start;
            start.x = intersectionPoint.x + direction.x*1;
            start.y = intersectionPoint.y + direction.y*1;
            start.z = intersectionPoint.z + direction.z*1;
            Ray L;
            L.start.x = start.x ;
            L.start.y = start.y ;
            L.start.z = start.z ;
            L.dir.x = direction.x ;
            L.dir.y = direction.y ;
            L.dir.z = direction.z ;
            int f=1;
            for (int j=0; j < objects.size(); j++)
            {
                double obs=objects[j]->findClosest(&L);
                if(obs>0)
                {
                    f=0;
                    break;
                }
                if((obs)>l)
                {
                    f=0;
                    break;
                }
            }

            if (f==1)
            {
                double lambert,phong;
                lambert = (L.dir.x * normal.x) + (L.dir.y * normal.y) + (L.dir.z * normal.z);
                phong = (reflection.x*r->dir.x) + (reflection.y*r->dir.y) + (reflection.z*r->dir.z);
                phong=pow(phong,shine);
                if(lambert<0) lambert=0;
                if(phong<0) phong=0;
                for (int m=0; m<3; m++)
                {
                    ncolor[m]+=source_factor* lambert*co_efficients[1]*color[m]*lights[i].color[m] ;
                    ncolor[m]+= source_factor*phong*co_efficients[2]*color[m]*lights[i].color[m];
                }
            }
            if (level<recursion_level)
            {
                start.x = intersectionPoint.x + reflection.x*1;
                start.y = intersectionPoint.y + reflection.y*1;
                start.z = intersectionPoint.z + reflection.z*1;

                Ray reflectionRay;
                reflectionRay.start.x = start.x ;
                reflectionRay.start.y = start.y ;
                reflectionRay.start.z = start.z ;
                reflectionRay.dir.x = reflection.x ;
                reflectionRay.dir.y = reflection.y ;
                reflectionRay.dir.z = reflection.z ;

                double mint = 99999;
                double refl_color[3];
                int nearest=-1;
                for (int k=0; k<objects.size(); k++)
                {
                    double tt= objects[k]->findClosest(&reflectionRay);
                    if(tt<=0)
                        continue;
                    else if(min(mint,tt)==tt)
                    {
                        mint=tt;
                        nearest=k;
                    }
                }

                if(nearest!=-1)
                {
                    objects[nearest]->intersect(&reflectionRay, refl_color, ++level);
                    for (int k=0; k<3; k++)
                        ncolor[k]+=refl_color[k]*co_efficients[3];
                }
            }
            for(int k=0; k<3; k++)
            {
                if(ncolor[k]<0)
                    ncolor[k]=0;
                else if(ncolor[k]>1)
                    ncolor[k]=1;
            }
        }
        return t;
    }


    float findClosestbarycentric(
            Point RayOrigin
            , Point RayDelta
            , Point Vertex1
            , Point Vertex2
            , Point Vertex3
            , float MinT)
    {
        const float NoIntersection = FLT_MAX;

        Point e1 ;
        Point e2 ;

        Point n ;

        float dot ;

        if (!(dot < 0.0f))
            return NoIntersection;

        float d ;

        float t = d;

        if (!(t <= 0.0f))
            return NoIntersection;

        if (!(t >= dot * MinT))
            return NoIntersection;

        t /= dot;
        if(t >= 0.0f && t <= MinT){

        }
        Point p ;

        float u0 = 0.f, u1 = 0.f, u2 = 0.f;
        float v0 = 0.f, v1 = 0.f, v2 = 0.f;

        if (std::fabs(n.x) > std::fabs(n.y))
        {
            if (std::fabs(n.x) > std::fabs(n.z))
            {
                u0 = p.y - Vertex1.y;
                u1 = Vertex2.y - Vertex1.y;
                u2 = Vertex3.y - Vertex1.y;

                v0 = p.z - Vertex1.z;
                v1 = Vertex2.z - Vertex1.z;
                v2 = Vertex3.z - Vertex1.z;
            }
            else
            {
                u0 = p.x - Vertex1.x;
                u1 = Vertex2.x - Vertex1.x;
                u2 = Vertex3.x - Vertex1.x;

                v0 = p.y - Vertex1.y;
                v1 = Vertex2.y - Vertex1.y;
                v2 = Vertex3.y - Vertex1.y;
            }
        }
        else
        {
            if (std::fabs(n.y) > std::fabs(n.z))
            {
                u0 = p.x - Vertex1.x;
                u1 = Vertex2.x - Vertex1.x;
                u2 = Vertex3.x - Vertex1.x;

                v0 = p.z - Vertex1.z;
                v1 = Vertex2.z - Vertex1.z;
                v2 = Vertex3.z - Vertex1.z;
            }
            else
            {
                u0 = p.x - Vertex1.x;
                u1 = Vertex2.x - Vertex1.x;
                u2 = Vertex3.x - Vertex1.x;

                v0 = p.y - Vertex1.y;
                v1 = Vertex2.y - Vertex1.y;
                v2 = Vertex3.y - Vertex1.y;
            }
        }

        float Denominator = u1*v2 - v1*u2;
        if (!(Denominator != 0.0f))
            return NoIntersection;
        Denominator = 1.0f / Denominator;

        float alpha = (u0*v2 - v0*u2) * Denominator;
        if (!(alpha >= 0.0f))
            return NoIntersection;

        float beta = (u1*v0 - v1*u0) * Denominator;
        if (!(beta >= 0.0f))
            return NoIntersection;

        float gamma = 1.0f - alpha - beta;
        if (!(gamma >= 0.0f))
            return NoIntersection;

        return t;
    }


};
