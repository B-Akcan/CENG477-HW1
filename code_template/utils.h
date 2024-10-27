#include "parser.h"
using namespace parser;

typedef struct {
    Vec3f e, d; //eye and direction (ray = e + dt)
} Ray;

Vec3f multiplyScalar(Vec3f v, float s);
void normalize(Vec3f &v);
Vec3f add(Vec3f a, Vec3f b);
Vec3f cross(Vec3f a, Vec3f b);
Vec3f negateVector(Vec3f v);
float dot(Vec3f a, Vec3f b);