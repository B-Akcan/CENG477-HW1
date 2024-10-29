#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "utils.h"
#include <cmath>

using namespace parser;
using namespace std;

Ray calculateRay(Camera cam, int i, int j);
float intersectSphere(Ray r, Sphere s, vector<Vec3f>& vertex_data);
Vec3f computeColor(Ray r, vector<Sphere> & sphereVector, vector<Vec3f>& vertex_data, vector<PointLight>& pointVector);
void colorPixel(unsigned char* &image, int pixelPosition, Vec3f color);
void calculateCameraUVector (Camera &cam);