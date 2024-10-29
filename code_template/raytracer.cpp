#include "raytracer.h"

int main(int argc, char* argv[])
{
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    for (Camera cam : scene.cameras) {
        int width = cam.image_width;
        int height = cam.image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        normalize(cam.gaze);
        calculateCameraUVector(cam);

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                Ray r = calculateRay(cam, i, j);
                Vec3f color = computeColor(r, scene.spheres, scene.vertex_data, scene.point_lights);
                int pixelPosition = (i * height + j) * 3;
                colorPixel(image, pixelPosition, color);
            }
        }

        std::string output_path = "../outputs_dev/" + cam.image_name;
        write_ppm(output_path.c_str(), image, cam.image_width, cam.image_height);
    }
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

void colorPixel(unsigned char* &image, int pixelPosition, Vec3f color) {
    image[ pixelPosition] = color.x * 255; // R
    image[ pixelPosition + 1] = color.y * 255; // G
    image[ pixelPosition + 2] = color.z * 255; // B
}

void calculateCameraUVector (Camera &cam) {
    cam.u = cross(cam.up, negateVector(cam.gaze));
}
