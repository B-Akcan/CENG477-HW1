#include "parser.h"
using namespace parser;

typedef struct {
    Vec3f e, d; //eye and direction (ray = e + dt)
} Ray;

typedef struct
{
    float t; // intersection time
	Vec3f point;
	Vec3f normal;
	int mat_id;
} Intersection;

Vec3f multiplyScalar(Vec3f v, float s);
Vec3f normalize(Vec3f v);
Vec3f add(Vec3f a, Vec3f b);
Vec3f subtract(Vec3f a, Vec3f b);
Vec3f cross(Vec3f a, Vec3f b);
Vec3f negateVector(Vec3f v);
float dot(Vec3f a, Vec3f b);
float determinant(float a1, float a2, float a3, float b1, float b2, float b3, float c1, float c2, float c3);
float norm(Vec3f v);
Vec3f multiplyVector(Vec3f a, Vec3f b);