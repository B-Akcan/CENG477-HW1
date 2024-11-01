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

        int index = 0;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                Ray r = calculateRay(cam, i, j);
                Vec3f color = computeColor(r, scene);
                int pixelPosition = index * 3;
                colorPixel(image, pixelPosition, color);
                index++;
            }
        }

        std::string output_path = "../outputs_dev/" + cam.image_name;
        write_ppm(output_path.c_str(), image, cam.image_width, cam.image_height);
    }
}

Ray calculateRay(Camera cam, int i, int j) {
    Ray r;
    Vec3f s; //pixel position at (i-th row, j-th column)
    float su = (j+0.5) * (cam.near_plane.y - cam.near_plane.x) / cam.image_width; //su is horizontal distance
    float sv = (i+0.5) * (cam.near_plane.w - cam.near_plane.z) / cam.image_height; //sv is vertical distance

    Vec3f m = add(cam.position, multiplyScalar(cam.gaze, cam.near_distance));
    Vec3f q = add(m, add(multiplyScalar(cam.u, cam.near_plane.x), multiplyScalar(cam.up, cam.near_plane.w)));
    s = add(q, add(multiplyScalar(cam.u, su), multiplyScalar(cam.up, -sv)));
    r.e = cam.position;
    r.d = add(s, multiplyScalar(cam.position, -1)); //direction vector */
    return r;
}

Vec3f computeColor(Ray ray, Scene scene) { //we find closest intersecting object and return its color
    Vec3f color = { scene.background_color.x / 255., scene.background_color.y / 255., scene.background_color.z / 255. };
    Intersection intersection;
    int min_i = -1; //min_i -> id of the closest object
    float min_t = __FLT_MAX__;
    
    for(int i = 0; i < scene.spheres.size(); i++){
        intersection = intersectSphere(ray, scene.spheres[i], scene.vertex_data); //get intersection point with that spesific sphere
        if(intersection.t < min_t && intersection.t >= 1) { //closer and in front of the image window
            min_t = intersection.t; 
            min_i = i;
            color = {1.0, 0.0, 0.0}; //TODO: red color for trial
        }
    }
    for (int i = 0; i < scene.triangles.size(); i++) {
        intersection = intersectTriangle(ray, scene.triangles[i], scene.vertex_data);
        if (intersection.t < min_t && intersection.t >= 1) {
            min_t = intersection.t; 
            min_i = i;
            color = {0.0, 1.0, 0.0}; //TODO: red color for trial
        }
    }
    for (int i = 0; i < scene.meshes.size(); i++) {
        for (int j = 0; j < scene.meshes[i].faces.size(); j++) {
            Triangle triangle = { scene.meshes[i].material_id, scene.meshes[i].faces[j] };
            intersection = intersectTriangle(ray, triangle, scene.vertex_data);
            if (intersection.t < min_t && intersection.t >= 1) {
                min_t = intersection.t; 
                min_i = i;
                color = {1.0, 1.0, 1.0}; //TODO: red color for trial
            }
        }
    }

    if (min_i == -1)
        return color; //no intersection, returning background color

    color = diffuseShading(color, scene, intersection);

    return color;
}

Intersection intersectSphere(Ray r, Sphere s, vector<Vec3f> vertex_data){
    Intersection intersection;
    float A,B,C; //coefficients of the quadratic equation
    float delta;
    Vec3f center;
    center = vertex_data[s.center_vertex_id - 1];
 
    //all scalars there
    //A -> d.d  t^2
    //B -> 2d(e-c) . t (o or e , same thing)
    //C -> (e-c)^2 - r^2
    C = (r.e.x - center.x) * (r.e.x - center.x) + (r.e.y - center.y) * (r.e.y - center.y) + (r.e.z - center.z) * (r.e.z - center.z) - s.radius * s.radius;
    B = 2 * (r.e.x - center.x) * r.d.x + 2 * (r.e.y - center.y) * r.d.y + 2 * (r.e.z - center.z) * r.d.z;
    A = r.d.x * r.d.x + r.d.y * r.d.y + r.d.z * r.d.z;

    delta = B*B - 4*A*C;
    if (delta < 0)
        intersection.t = -1; //no intersection
    else if (delta == 0)
        intersection.t = -B / (2*A); //one intersection point (tangent)
    else {
        delta = sqrt(delta);
        A = 2*A;
        float t1 = (-B - delta) / (A);
        float t2 = (-B + delta) / (A);
        intersection.t = (t1 < t2) ? t1 : t2; //return the closest intersection point
    }

    intersection.point = add(r.e, multiplyScalar(r.d, intersection.t));
    intersection.normal = subtract(intersection.point, center);
    normalize(intersection.normal);
    intersection.mat_id = s.material_id;

    return intersection;
}

Intersection intersectTriangle(Ray ray, Triangle triangle, vector<Vec3f> vertex_data) {
    Intersection intersection;
    float t, alpha, beta, gama, detA;
    Vec3f v1 = vertex_data[triangle.indices.v0_id - 1];
    Vec3f v2 = vertex_data[triangle.indices.v1_id - 1];
    Vec3f v3 = vertex_data[triangle.indices.v2_id - 1];
    
    detA = determinant( v1.x-v2.x, v1.x-v3.x, ray.d.x,
                        v1.y-v2.y, v1.y-v3.y, ray.d.y,
                        v1.z-v2.z, v1.z-v3.z, ray.d.z );
    beta = determinant( v1.x-ray.e.x, v1.x-v3.x, ray.d.x,
                        v1.y-ray.e.y, v1.y-v3.y, ray.d.y,
                        v1.z-ray.e.z, v1.z-v3.z, ray.d.z ) / detA;
    gama = determinant( v1.x-v2.x, v1.x-ray.e.x, ray.d.x,
                        v1.y-v2.y, v1.y-ray.e.y, ray.d.y,
                        v1.z-v2.z, v1.z-ray.e.z, ray.d.z ) / detA;
    t = determinant(v1.x-v2.x, v1.x-v3.x, v1.x-ray.e.x,
                    v1.y-v2.y, v1.y-v3.y, v1.y-ray.e.y,
                    v1.z-v2.z, v1.z-v3.z, v1.z-ray.e.z ) / detA;

    if (beta + gama <= 1 && beta >= 0 && gama >= 0 && t > 0) {
        intersection.t = t;
        intersection.point = add(ray.e, multiplyScalar(ray.d, t));
        intersection.normal = cross(subtract(v3, v2), subtract(v1, v2));
        normalize(intersection.normal);
        intersection.mat_id = triangle.material_id;
    }
    else {
        intersection.t = -1;
    }

    return intersection;
}

void colorPixel(unsigned char* &image, int pixelPosition, Vec3f color) {
    image[ pixelPosition ] = color.x * 255; // R
    image[ pixelPosition + 1] = color.y * 255; // G
    image[ pixelPosition + 2] = color.z * 255; // B
}

void calculateCameraUVector (Camera &cam) {
    cam.u = cross(cam.up, negateVector(cam.gaze));
}

Vec3f diffuseShading(Vec3f color, Scene scene, Intersection intersection) {
    Vec3f L = add(scene.point_lights[0].position, multiplyScalar(intersection.point, -1));   //light point - intersection point
    normalize(L);

    if (dot(intersection.normal, L) < 0) //check whether light ray really sees the surface point
        color = { 0, 0, 0 }; //if the light is coming from the back side of the sphere, return the color
    else
        color = multiplyScalar(color, dot(intersection.normal, L)); //diffuse component

    return color;
}
