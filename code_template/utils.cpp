#include "parser.h"
#include <cmath>
#include "utils.h"
using namespace parser;

Vec3f multiplyScalar(Vec3f v, float s) {
    return {v.x * s, v.y * s, v.z * s};
}

void normalize(Vec3f &v) {
    Vec3f result;
    float length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    result.x = v.x / length;
    result.y = v.y / length;
    result.z = v.z / length;

    v.x = result.x;
    v.y = result.y;
    v.z = result.z;
}

Vec3f add(Vec3f a, Vec3f b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3f cross(Vec3f a, Vec3f b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

Vec3f negate(Vec3f v) {
    return multiplyScalar(v, -1);
}

