#include "raytracer.h"

vector<ExtendedTriangle> extendedTriangleVector; //triangle with normal
vector<ExtendedMesh> extendedMeshVector; //mesh with normal

int main(int argc, char* argv[])
{
    parser::Scene scene;
    scene.loadFromXml(argv[1]);

    for(Triangle &t: scene.triangles){
        Vec3f v1 = scene.vertex_data[t.indices.v0_id - 1];
        Vec3f v2 = scene.vertex_data[t.indices.v1_id - 1];
        Vec3f v3 = scene.vertex_data[t.indices.v2_id - 1];
        Vec3f normal = normalize(cross(subtract(v3, v2), subtract(v1, v2)));
        extendedTriangleVector.push_back({t.material_id, t.indices, normal});
    }

    for(Mesh &m: scene.meshes){
        ExtendedMesh em;
        em.material_id = m.material_id;
        em.faces = m.faces;
        for(Face &f: m.faces){
            Vec3f v1 = scene.vertex_data[f.v0_id - 1];
            Vec3f v2 = scene.vertex_data[f.v1_id - 1];
            Vec3f v3 = scene.vertex_data[f.v2_id - 1];
            Vec3f normal = normalize(cross(subtract(v3, v2), subtract(v1, v2)));
            em.normals.push_back(normal); 
        }
        extendedMeshVector.push_back(em);
    }

    for (Camera cam : scene.cameras) {
        int width = cam.image_width;
        int height = cam.image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        cam.gaze = normalize(cam.gaze);
        
        thread t1(rayTracing, ref(scene), ref(cam), ref(image), 0, height/2, 0, width/2);
        thread t2(rayTracing, ref(scene), ref(cam), ref(image), 0, height/2, width/2, width);
        thread t3(rayTracing, ref(scene), ref(cam), ref(image), height/2, height, 0, width/2);
        thread t4(rayTracing, ref(scene), ref(cam), ref(image), height/2, height, width/2, width);
        t1.join();
        t2.join();
        t3.join();
        t4.join();

        std::string output_path = "../outputs_dev/" + cam.image_name;
        write_ppm(output_path.c_str(), image, cam.image_width, cam.image_height);
        delete[] image;
    }
}

void rayTracing(Scene &scene, Camera &cam, unsigned char* &image, int heightStart, int heightEnd, int widthStart, int widthEnd) {
    for (int i = heightStart; i < heightEnd; i++) {
        for (int j = widthStart; j < widthEnd; j++) {
            int index = (i * cam.image_width + j) * 3;
            Ray ray = calculateRay(cam, i, j);
            Vec3f color = computeColor(ray, scene, false);
            colorPixel(image, index, color);
        }
    }
}

Ray calculateRay(Camera &cam, int i, int j) {
    Ray r;
    Vec3f cameraUVector = cross(cam.up, negateVector(cam.gaze));
    float su = (j+0.5) * (cam.near_plane.y - cam.near_plane.x) / cam.image_width; //su is horizontal distance
    float sv = (i+0.5) * (cam.near_plane.w - cam.near_plane.z) / cam.image_height; //sv is vertical distance
    Vec3f m = add(cam.position, multiplyScalar(cam.gaze, cam.near_distance));
    Vec3f q = add(m, add(multiplyScalar(cameraUVector, cam.near_plane.x), multiplyScalar(cam.up, cam.near_plane.w)));
    Vec3f s = add(q, add(multiplyScalar(cameraUVector, su), multiplyScalar(cam.up, -sv))); //pixel position at (i-th row, j-th column)

    r.e = cam.position;
    r.d = add(s, multiplyScalar(cam.position, -1)); //direction vector */
    r.depth = 0;

    return r;
}

Intersection findClosestIntersection(Ray ray, Scene &scene, bool is_shadow) {
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
    for (int i = 0; i < extendedTriangleVector.size(); i++) {
        temp_intersection = intersectTriangle(ray, extendedTriangleVector[i], scene.vertex_data, scene);
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
    for (int i = 0; i < extendedMeshVector.size(); i++) {
        for (int j = 0; j < extendedMeshVector[i].faces.size(); j++) {
            ExtendedTriangle et = { extendedMeshVector[i].material_id, extendedMeshVector[i].faces[j], extendedMeshVector[i].normals[j] };
            temp_intersection = intersectTriangle(ray, et, scene.vertex_data, scene);
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

Intersection intersectSphere(Ray r, Sphere s, vector<Vec3f> &vertex_data, Scene &scene) { 
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
    if (delta < -scene.shadow_ray_epsilon)
        intersection.t = -1; //no intersection
    else if (delta > -scene.shadow_ray_epsilon && delta < scene.shadow_ray_epsilon)
        intersection.t = -B / (2*A);
    else {
        delta = sqrt(delta);
        A = 2*A;
        float t1 = (-B - delta) / (A);
        float t2 = (-B + delta) / (A);
        intersection.t = (t1 < t2) ? t1 : t2; //return the closest intersection point
    }

    intersection.point = add(r.e, multiplyScalar(r.d, intersection.t));
    intersection.normal = normalize(subtract(intersection.point, center));
    intersection.mat_id = s.material_id;

    return intersection;
}

Intersection intersectTriangle(Ray ray, ExtendedTriangle triangle, vector<Vec3f> &vertex_data, Scene &scene) {
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
        intersection.normal = triangle.normal; //NORMAL
        intersection.mat_id = triangle.material_id;
    }
    else {
        intersection.t = -1;
    }

    return intersection;
}

Vec3f computeColor(Ray ray, Scene &scene, bool is_shadow_or_reflection) {
    if (ray.depth > scene.max_recursion_depth)
        return { 0.f, 0.f, 0.f };

    Intersection intersection = findClosestIntersection(ray, scene, is_shadow_or_reflection);

    if (intersection.t != -1) {
        return applyShading(ray, scene, intersection);
    }
    else if (ray.depth == 0) {
        return { (float) scene.background_color.x, (float) scene.background_color.y, (float) scene.background_color.z };
    }
    else {
        return { 0.f, 0.f, 0.f };
    }
}

Vec3f applyShading(Ray ray, Scene &scene, Intersection intersection) {
    Vec3f color = ambientShading(scene, intersection);
    Vec3f pointPlusEpsilon = add(intersection.point, multiplyScalar(intersection.normal, scene.shadow_ray_epsilon));
    Material material = scene.materials[intersection.mat_id - 1];

    for (PointLight pl : scene.point_lights) {
        Vec3f lightDirection = subtract(pl.position, intersection.point);
        Vec3f normalizedLightDirection = normalize(lightDirection);
        Ray shadowRay = { pointPlusEpsilon, normalizedLightDirection, 0 };
        Intersection shadowIntersection = findClosestIntersection(shadowRay, scene, true);
        bool notInShadow = shadowIntersection.t >= (norm(lightDirection)-scene.shadow_ray_epsilon) || shadowIntersection.t == -1;

        if (notInShadow) {
            color = add(color, diffuseShading(scene, intersection, pl));
            color = add(color, specularShading(scene, intersection, ray, pl));
        }
    }

    if (material.is_mirror) {
        Vec3f reflectionDirection = subtract(multiplyScalar(intersection.normal, 2 * dot(normalize(subtract(ray.e, intersection.point)), intersection.normal)), normalize(subtract(ray.e, intersection.point)));
        Ray reflectionRay = {pointPlusEpsilon, reflectionDirection, ray.depth + 1};
        color = add(color, multiplyVector(computeColor(reflectionRay, scene, true), material.mirror));
    }

    return color;
}

Vec3f ambientShading(Scene &scene, Intersection intersection) {
    return multiplyVector(scene.ambient_light, scene.materials[intersection.mat_id - 1].ambient);
}

Vec3f diffuseShading(Scene &scene, Intersection intersection, PointLight pl) {
    Vec3f L = subtract(pl.position, intersection.point); // light vector
    float cos_theta = max(0.f, dot(normalize(L), intersection.normal));
    float dist_squared = dot(L, L);
    Vec3f colorForPl = multiplyVector(pl.intensity, multiplyScalar(scene.materials[intersection.mat_id - 1].diffuse, cos_theta/dist_squared));

    return colorForPl;
}

Vec3f specularShading(Scene &scene, Intersection intersection, Ray ray, PointLight pl) {
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

void colorPixel(unsigned char* &image, int pixelPosition, Vec3f color) {
    image[pixelPosition] = (int) round(min(max(0.f, color.x), 255.f));
    image[pixelPosition + 1] = (int) round(min(max(0.f, color.y), 255.f));
    image[pixelPosition + 2] =  (int) round(min(max(0.f, color.z), 255.f));
}