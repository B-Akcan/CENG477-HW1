#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "utils.h"
#include <cmath>

using namespace parser;

typedef unsigned char RGB[3];
Ray calculateRay(Camera cam, int i, int j);

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.

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
    //each camera has its own image plane
    
    //int columnWidth = width / 8;

    
    

    int columns, rows;
    for (Camera cam : scene.cameras) { //for each camera in the scene, we will generate another output file
        int width = scene.cameras[0].image_width, height = scene.cameras[0].image_height;
        unsigned char* image = new unsigned char [width * height * 3]; //that will bring huge memory overhead
        //since it will be generated for each camera
        //cannot give variable-size length to c-array since size should be known at compile time
        columns = cam.image_width;
        rows = cam.image_height;
        normalize(cam.gaze);
        cam.u = cross(cam.up, negate(cam.gaze)); // u = v x w
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                Ray r = calculateRay(cam, i, j);

            }
        }
        write_ppm(cam.image_name.c_str(), image, cam.image_width, cam.image_height);
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

    /* Vec3f m = add(cam.position, multiplyScalar(cam.gaze, cam.near_distance));
    Vec3f q = add(m, add(multiplyScalar(u, left), multiplyScalar(v, top)));
    s = add(q, add(multiplyScalar(u, su), multiplyScalar(v, -sv)));
    r.e = eye;
    r.d = add(s, multiplyScalar(eye, -1)); //direction vector */
    return r;
}

