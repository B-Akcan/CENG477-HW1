#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "utils.h"
#include <cmath>

using namespace parser;
using namespace std;

typedef unsigned char RGB[3];
Ray calculateRay(Camera cam, int i, int j);
float intersectSphere(Ray r, Sphere s, vector<Vec3f>& vertex_data);
Vec3f computeColor(Ray r, vector<Sphere> & sphereVector, vector<Vec3f>& vertex_data, vector<PointLight>& pointVector);

int main(int argc, char* argv[])
{
    parser::Scene scene;
    scene.loadFromXml(argv[1]);

    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int columns, rows;
    for (Camera cam : scene.cameras) { //for each camera in the scene, we will generate another output file
        int width = cam.image_width, height = cam.image_height;
        unsigned char* image = new unsigned char [width * height * 3]; //that will bring huge memory overhead
        //since it will be generated for each camera
        //cannot give variable-size length to c-array since size should be known at compile time
        columns = cam.image_width;
        rows = cam.image_height;
        normalize(cam.gaze);
        cam.u = cross(cam.up, negateVector(cam.gaze)); // u = v x w
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                Ray r = calculateRay(cam, i, j);
                Vec3f color = computeColor(r, scene.spheres, scene.vertex_data, scene.point_lights); //will return which object is intersected (closest one) 
                //cout << "Color: " << color.x << " " << color.y << " " << color.z << endl;
                image[ (j*height  + i)*3] = color.x * 255;
                image[ (j*height  + i)*3 + 1] = color.y * 255;
                image[ (j*height  + i)*3 + 2] = color.z * 255;
                /*
                image[i*width*3+j].x = color.x * 255;
                image[i][j].y = color.y * 255;
                image[i][j].z = color.z * 255;
                */
            }
        }

        std::string output_path = "../outputs_dev/" + cam.image_name;
        write_ppm(output_path.c_str(), image, cam.image_width, cam.image_height);
    }
    /*
    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }
    */



    
}

Ray calculateRay(Camera cam, int i, int j) {
    Ray r;
    Vec3f s; //pixel position at (i,j)
    float su = (i+0.5) * (cam.near_plane.y - cam.near_plane.x) / cam.image_width;
    float sv = (j+0.5) * (cam.near_plane.w - cam.near_plane.z) / cam.image_height;

    Vec3f m = add(cam.position, multiplyScalar(cam.gaze, cam.near_distance));
    Vec3f q = add(m, add(multiplyScalar(cam.u, cam.near_plane.x), multiplyScalar(cam.up, cam.near_plane.w)));
    s = add(q, add(multiplyScalar(cam.u, su), multiplyScalar(cam.up, -sv)));
    r.e = cam.position;
    r.d = add(s, multiplyScalar(cam.position, -1)); //direction vector */
    return r;
}

float intersectSphere(Ray r, Sphere s, vector<Vec3f>& vertex_data){
    float A,B,C; //coefficients of the quadratic equation
    float delta;
    Vec3f c;
    c = vertex_data[s.center_vertex_id];
 
    //all scalars there
    //A -> d.d  t^2
    //B -> 2d(e-c) . t (o or e , same thing)
    //C -> (e-c)^2 - r^2
    C = (r.e.x - c.x) * (r.e.x - c.x) + (r.e.y - c.y) * (r.e.y - c.y) + (r.e.z - c.z) * (r.e.z - c.z) - s.radius * s.radius;
    B = 2 * (r.e.x - c.x) * r.d.x + 2 * (r.e.y - c.y) * r.d.y + 2 * (r.e.z - c.z) * r.d.z;
    A = r.d.x * r.d.x + r.d.y * r.d.y + r.d.z * r.d.z;

    delta = B*B - 4*A*C;
    if(delta < 0) return -1; //no intersection
    else if(delta == 0) return -B / (2*A); //one intersection point (tangent)
    else {
        delta = sqrt(delta);
        A = 2*A;
        float t1 = (-B - delta) / (A);
        float t2 = (-B + delta) / (A);
        return (t1 < t2) ? t1 : t2; //return the closest intersection point
    }

}

Vec3f computeColor(Ray r, vector<Sphere> & sphereVector, vector<Vec3f>& vertex_data, vector<PointLight>& pointVector) { //we find closest intersecting object and return its color
    int i, min_i = -1; //min_i -> id of the closest object
    Vec3f c ={0,0,0}; //color
    Vec3f L, N, P; //light vector and normal vector and intersection point
    
    float min_t = 1000000; //initially infinity
    float t; //time of intersection
    int size = sphereVector.size();
    for(i = 0; i < size; i++){
        t = intersectSphere(r, sphereVector[i], vertex_data); //get intersection point with that spesific sphere
        if(t < min_t && t>= 1) { //closer and in front of the image window
            min_t = t; 
            min_i = i;
            //c = sphereVector[i].color; //can also return the id of the object
            c = {255,0,0}; //TODO: red color for trial
        }
    }
    if(min_i == -1) 
        return c; //no intersection, returning background color

    P = add(r.e, multiplyScalar(r.d, min_t)); //intersection point
    L = add(pointVector[0].position, multiplyScalar(P , -1));   //light point - intersection point
    N = add(P, multiplyScalar(vertex_data[sphereVector[min_i].center_vertex_id], -1)); //normal vector (P - center of the sphere)
    normalize(L);
    normalize(N);
    if(dot(N,L) < 0) //check whether light ray really sees the surface point
        c.x = c.y = c.z = 0; //if the light is coming from the back side of the sphere, return the color
    else
        c = multiplyScalar(c, dot(N,L)); //diffuse component
    //diffuse component -> (color) x (N . L) 

    //normally, if dot product is greater than 0 ,we should return the color
    return c;
}
