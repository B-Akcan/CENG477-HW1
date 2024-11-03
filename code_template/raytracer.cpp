#include "raytracer.h"

int main(int argc, char* argv[])
{
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    for (Camera cam : scene.cameras) {
        int width = cam.image_width;
        int height = cam.image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        cam.gaze = normalize(cam.gaze);
        calculateCameraUVector(cam);

        int index = 0;
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                Ray r = calculateRay(cam, i, j);
                Vec3f color = color = computeColor(r, scene);
                if((i %100 == 0)  && j == 0){
                    cout << color.x << " " << color.y << " " << color.z << "at line: " << i << endl;
                    color = computeColor(r, scene);
                }
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

Vec3f computeColor(Ray ray, Scene scene) {
    float min_t = __FLT_MAX__;
    bool is_shadow = false;
    Intersection intersection = findClosestIntersection( ray, scene, is_shadow);

    if (intersection.t == -1)
        return { (float) scene.background_color.x, (float) scene.background_color.y, (float) scene.background_color.z }; //no intersection, returning background color

    Vec3f color = { 0.0, 0.0, 0.0 };;
    color = add(color, ambientShading(scene, intersection));

    Vec3f pointPlusEpsilon = add(intersection.point, multiplyScalar(intersection.normal, scene.shadow_ray_epsilon));
    for (PointLight pl : scene.point_lights) {
        Vec3f lightDirection = subtract(pl.position, intersection.point);
        Vec3f normalizedLightDirection = normalize(lightDirection);
        Ray shadowRay = { pointPlusEpsilon, normalizedLightDirection };
        float temp = __FLT_MAX__;
        is_shadow = true;
        Intersection shadowIntersection = findClosestIntersection(shadowRay, scene, is_shadow);

        if (shadowIntersection.t >= (norm(lightDirection)-scene.shadow_ray_epsilon) || shadowIntersection.t == -1) {
            color = add(color, diffuseShading(scene, intersection, pl));
            color = add(color, specularShading(scene, intersection, ray, pl));
        }
    }

    return color;
}


Intersection findClosestIntersection(Ray ray, Scene scene, bool is_shadow) {
    float min_t = __FLT_MAX__; //initially infinite
    Intersection temp_intersection;
    Intersection intersection;
    intersection.t = -1;

    for(int i = 0; i < scene.spheres.size(); i++) {
        temp_intersection = intersectSphere(ray, scene.spheres[i], scene.vertex_data, scene); //get intersection point with that spesific sphere
        if(is_shadow){
            if (temp_intersection.t < min_t && temp_intersection.t >= -scene.shadow_ray_epsilon){
                min_t = temp_intersection.t;
                intersection = temp_intersection;
            } 
        }
        else{
            if (temp_intersection.t < min_t && temp_intersection.t >= 1-scene.shadow_ray_epsilon) { //closer and in front of the image window
                min_t = temp_intersection.t;
                intersection = temp_intersection;
            }
        }

    }
    for (int i = 0; i < scene.triangles.size(); i++) {
        temp_intersection = intersectTriangle(ray, scene.triangles[i], scene.vertex_data, scene);
        if(is_shadow){
            if (temp_intersection.t < min_t && temp_intersection.t >= -scene.shadow_ray_epsilon){
                min_t = temp_intersection.t;
                intersection = temp_intersection;
            } 
        }
        else{
            if (temp_intersection.t < min_t && temp_intersection.t >= 1-scene.shadow_ray_epsilon) { //closer and in front of the image window
                min_t = temp_intersection.t;
                intersection = temp_intersection;
            }
        }
    }
    for (int i = 0; i < scene.meshes.size(); i++) {
        for (int j = 0; j < scene.meshes[i].faces.size(); j++) {
            Triangle triangle = { scene.meshes[i].material_id, scene.meshes[i].faces[j] };
            temp_intersection = intersectTriangle(ray, triangle, scene.vertex_data, scene);
            if(is_shadow){
                if (temp_intersection.t < min_t && temp_intersection.t >= -scene.shadow_ray_epsilon){
                    min_t = temp_intersection.t;
                    intersection = temp_intersection;
                } 
            }
            else{
                if (temp_intersection.t < min_t && temp_intersection.t >= 1-scene.shadow_ray_epsilon) { //closer and in front of the image window
                    min_t = temp_intersection.t;
                    intersection = temp_intersection;
                }
            }
        }
    }

    return intersection;
}

Intersection intersectSphere(Ray r, Sphere s, vector<Vec3f> vertex_data, Scene scene) { 
    //if t becomes negative, even though delta is positive, it is handled at caller function as expected
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
    if (delta < scene.shadow_ray_epsilon)
        intersection.t = -1; //no intersection

    else{
        float t1 = (-B - sqrt(delta)) / (2*A);
        float t2 = (-B + sqrt(delta)) / (2*A);
        intersection.t = (t1 < t2) ? t1 : t2; //finding closes point
        if(intersection.t < scene.shadow_ray_epsilon) {
            intersection.t = (t1 > t2) ? t1 : t2; //if the closest is negative, checking the other
            if(intersection.t < scene.shadow_ray_epsilon)
                intersection.t = -1; //if both are negative, no desired intersection
        }
    }



    intersection.point = add(r.e, multiplyScalar(r.d, intersection.t));
    intersection.normal = normalize(subtract(intersection.point, center));
    intersection.mat_id = s.material_id;

    return intersection;
}



Intersection intersectTriangle(Ray ray, Triangle triangle, vector<Vec3f> vertex_data, Scene scene) {
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

    if (beta + gama <= (1+scene.shadow_ray_epsilon) && beta >= -scene.shadow_ray_epsilon && gama >= -scene.shadow_ray_epsilon && t > 0) {
        intersection.t = t;
        intersection.point = add(ray.e, multiplyScalar(ray.d, t));
        intersection.normal = normalize(cross(subtract(v3, v2), subtract(v1, v2)));
        intersection.mat_id = triangle.material_id;
    }
    else {
        intersection.t = -1;
    }

    return intersection;
}

void colorPixel(unsigned char* &image, int pixelPosition, Vec3f color) {
    image[pixelPosition] = (int) round(min(max(0.f, color.x), 255.f));
    image[pixelPosition + 1] = (int) round(min(max(0.f, color.y), 255.f));
    image[pixelPosition + 2] =  (int) round(min(max(0.f, color.z), 255.f));
}

void calculateCameraUVector (Camera &cam) {
    cam.u = cross(cam.up, negateVector(cam.gaze));
}

Vec3f ambientShading(Scene scene, Intersection intersection) {
    return multiplyVector(scene.ambient_light, scene.materials[intersection.mat_id - 1].ambient);
}

Vec3f diffuseShading(Scene scene, Intersection intersection, PointLight pl) {
    Vec3f L = subtract(pl.position, intersection.point); // light vector
    float cos_theta = max(0.f, dot(normalize(L), intersection.normal));
    float dist_squared = dot(L, L);
    Vec3f colorForPl = multiplyVector(pl.intensity, multiplyScalar(scene.materials[intersection.mat_id - 1].diffuse, cos_theta/dist_squared));

    return colorForPl;
}

Vec3f specularShading(Scene scene, Intersection intersection, Ray ray, PointLight pl) {
    Vec3f L = subtract(pl.position, intersection.point); // light vector

    float cos_theta = max(0.f, dot(normalize(L), intersection.normal));
    float theta = acos(cos_theta);
    if (theta >= 90)
        return { 0.f, 0.f, 0.f };

    Vec3f halfVector = multiplyScalar(add(normalize(L), normalize(subtract(ray.e, intersection.point))), 0.5);
    float cos_alpha_with_phong = pow(max(0.f, dot(normalize(halfVector), intersection.normal)), scene.materials[intersection.mat_id - 1].phong_exponent);
    float dist_squared = dot(L, L);

    Vec3f colorForPl = multiplyVector(scene.materials[intersection.mat_id - 1].specular, multiplyScalar(pl.intensity, cos_alpha_with_phong / dist_squared));

    return colorForPl;
}
