#define _USE_MATH_DEFINES
#include <object.h>
#include <iostream>
using namespace std;

class Plane{
    public:
        virtual double intersect(const Ray &r, Record &rec) {}

};

class flip_normals : public Plane {
    public:
        Plane *plane;
        double d;
        flip_normals(Plane *p) : plane(p) {}
        double intersect(const Ray &r, Record &rec){
            if (d = plane->intersect(r, rec)) {
                rec.normal = rec.normal*-1;
                return d;
            }
            else 
                return 0;
        }
};

class xy_plane : public Plane
{
    public:
        Refl_t refl;
        float x0, x1, y0, y1, k;
        Vec c, e, normal = Vec(0, 0, 1);
        xy_plane(float _x0, float _x1, float _y0, float _y1, float _k, Vec _e, Vec _c, Refl_t material) : 
            x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), e(_e), c(_c), refl(material) {};
        virtual double intersect(const Ray &r, Record &rec);
};

class xz_plane : public Plane
{
    public:
        Refl_t refl;
        float x0, x1, z0, z1, k;
        Vec c, e, normal = Vec(0, 1, 0);
        virtual double intersect(const Ray &r, Record &rec);        
        xz_plane(float _x0, float _x1, float _z0, float _z1, float _k, Vec _e, Vec _c, Refl_t material) : 
            x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), e(_e), c(_c), refl(material) {};
};

class yz_plane : public Plane
{
    public:
        Refl_t refl;
        float y0, y1, z0, z1, k;
        Vec c, e, normal = Vec(1, 0, 0);
        virtual double intersect(const Ray &r, Record &rec);   
        yz_plane(float _y0, float _y1, float _z0, float _z1, float _k, Vec _e, Vec _c, Refl_t material) : 
            y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), e(_e), c(_c), refl(material) {};
};

double xy_plane::intersect(const Ray &r, Record &rec)
{
/* if intersect, return distance, if not, return 0 */
/* negtive t is not allowd */
    float t = (k-r.o.z) / r.d.z;
    if (t < 0.1 || t > 100)
        return 0;
    float x = r.o.x + t*r.d.x;
    float y = r.o.y + t*r.d.y;
    if (x < x0 || x > x1 || y < y0 || y > y1)
        return 0;
    rec.e = e;
    rec.c = c;
    rec.normal = normal;
    rec.refl = refl;
    return t;
};

double xz_plane::intersect(const Ray &r, Record &rec)
{
/* if intersect, return distance, if not, return 0 */
/* negtive t is not allowd */
    float t = (k-r.o.y) / r.d.y;
    if (t < 0.1 || t > 100)
        return 0;
    float x = r.o.x + t*r.d.x;
    float z = r.o.z + t*r.d.z;
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return 0;
    rec.e = e;
    rec.c = c;
    rec.normal = normal;
    rec.refl = refl;
    return t;
};

double yz_plane::intersect(const Ray &r, Record &rec)
{
/* if intersect, return distance, if not, return 0 */
/* negtive t is not allowd */
    float t = (k-r.o.x) / r.d.x;
    if (t < 0.1 || t > 100)
        return 0;
    float y = r.o.y + t*r.d.y;
    float z = r.o.z + t*r.d.z;
    if (y < y0 || y > y1 || z < z0 || z > z1) 
        return 0;
    rec.e = e;
    rec.c = c;
    rec.normal = normal;
    rec.refl = refl;
    return t;
};