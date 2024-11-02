#include "parser.h"
#include <cmath>
#include "utils.h"
#include <iostream>
using namespace parser;

Vec3f multiplyScalar(Vec3f v, float s) {
    return {v.x * s, v.y * s, v.z * s};
}

Vec3f normalize(Vec3f v) {
    Vec3f result;
    float normOfVector = norm(v);

    result.x = v.x / normOfVector;
    result.y = v.y / normOfVector;
    result.z = v.z / normOfVector;

    return result;
}

Vec3f add(Vec3f a, Vec3f b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3f subtract(Vec3f a, Vec3f b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3f cross(Vec3f a, Vec3f b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

Vec3f negateVector(Vec3f v) {
    return multiplyScalar(v, -1);
}

float dot(Vec3f a, Vec3f b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void printVec3f(Vec3f vec) {
    std::cout << "x: " << vec.x << " y: " << vec.y << " z: " << vec.z << std::endl;
}

float determinant(float a1, float a2, float a3, float b1, float b2, float b3, float c1, float c2, float c3)
{
	return (a1*(c3*b2-b3*c2)-a2*(c3*b1-c1*b3)+a3*(c2*b1-c1*b2));
}

float norm(Vec3f v)
{
    return sqrt((v.x)*(v.x)+(v.y)*(v.y)+(v.z)*(v.z));
}

Vec3f multiplyVector(Vec3f a, Vec3f b) {
    return { a.x * b.x, a.y * b.y, a.z * b.z };
}

