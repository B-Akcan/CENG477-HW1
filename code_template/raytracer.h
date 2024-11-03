#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "utils.h"
#include <cmath>
#include <limits>

using namespace parser;
using namespace std;

Ray calculateRay(Camera cam, int i, int j);
Intersection findClosestIntersection( Ray ray, Scene scene, bool is_shadow);
Intersection intersectSphere(Ray r, Sphere s, vector<Vec3f> vertex_data, Scene scene);
Intersection intersectTriangle(Ray ray, Triangle triangle, vector<Vec3f> vertex_data, Scene scene);
Vec3f computeColor(Ray r, Scene scene);
Vec3f applyShading(Ray ray, Scene scene, Intersection intersection);
Vec3f ambientShading(Scene scene, Intersection intersection);
Vec3f diffuseShading(Scene scene, Intersection intersection, PointLight pl);
Vec3f specularShading(Scene scene, Intersection intersection, Ray ray, PointLight pl);
void colorPixel(unsigned char* &image, int pixelPosition, Vec3f color);