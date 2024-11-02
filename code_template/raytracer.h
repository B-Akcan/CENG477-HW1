#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "utils.h"
#include <cmath>
#include <limits>

using namespace parser;
using namespace std;

Ray calculateRay(Camera cam, int i, int j);
Intersection findClosestIntersection(float &min_t, Ray ray, Scene scene);
Intersection intersectSphere(Ray r, Sphere s, vector<Vec3f> vertex_data);
Intersection intersectTriangle(Ray ray, Triangle triangle, vector<Vec3f> vertex_data);
Vec3f computeColor(Ray r, Scene scene);
void colorPixel(unsigned char* &image, int pixelPosition, Vec3f color);
void calculateCameraUVector (Camera &cam);
Vec3f ambientShading(Scene scene, Intersection intersection);
Vec3f diffuseShading(Scene scene, Intersection intersection, PointLight pl);
Vec3f specularShading(Scene scene, Intersection intersection, Ray ray, PointLight pl);