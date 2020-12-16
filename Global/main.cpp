#include <iostream>
#include "../glm/glm/glm.hpp"
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <string>
#include <ctime>

using namespace std;
using namespace glm;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCDFAInspection"


const double indirectLight = 0.1;
const double directLight = 1;
const int cutoff_value = 200;
const int nr_shadowRays = 10;


struct ColorDbl{

    double r, g, b;
    ColorDbl(): r(0.0), g(0.0), b(0.0) {} //default constructor
    ColorDbl(double intensity): r(intensity), g(intensity), b(intensity){}
    ColorDbl(double red, double green, double blue) : r(red), g(green), b(blue){}
    [[maybe_unused]] ColorDbl(vec3 color)
    {
        r = color.x;
        g = color.y;
        b = color.z;
    }

    ColorDbl& operator+=(ColorDbl const& other){
        r += other.r;
        g += other.g;
        b += other.b;
        return *this;
    }

    ColorDbl& operator*(ColorDbl const& other){
        r *= other.r;
        g *= other.g;
        b *= other.b;
        return *this;
    }
    ColorDbl& operator*=(ColorDbl const& other){
        r *= other.r;
        g *= other.g;
        b *= other.b;
        return *this;
    }

    ColorDbl& operator/=(ColorDbl const& other){
        r /= other.r;
        g /= other.g;
        b /= other.b;
        return *this;
    }

    ColorDbl& operator/=(int val){
        r /= (double) val;
        g /= (double) val;
        b /= (double) val;
        return *this;
    }

    ColorDbl& operator*=(double scalar){
        r *= scalar;
        g *= scalar;
        b *= scalar;
        return *this;
    }

    bool operator==(const ColorDbl& other) const {
        return (r == other.r && g == other.g && b == other.b);
    }

};
#pragma clang diagnostic pop

[[maybe_unused]] string to_string(ColorDbl color){
    ostringstream res{};
    res << "R: " << color.r << ", G: " << color.g << ", B: " << color.b;
    return res.str();
}

string to_string(vec4 vector){
    ostringstream res{};
    res << "x: " << vector.x << ", y: " << vector.y << ", z: " << vector.z << ", w: " << vector.w;
    return res.str();
}

[[maybe_unused]] string to_string(vec3 vector){
    ostringstream res{};
    res << "x: " << vector.x << ", y: " << vector.y << ", z: " << vector.z;
    return res.str();
}

ColorDbl sqrt(ColorDbl input){
    ColorDbl res;
    res.r = sqrt(input.r);
    res.g = sqrt(input.g);
    res.b = sqrt(input.b);
    return res;

}

/*ColorDbl operator*(ColorDbl const& color, double scalar){
    ColorDbl ret{};
    ret.r = color.r*scalar;
    ret.g = color.g*scalar;
    ret.b = color.b*scalar;
    return ret;
}
*/
vec3 operator* (double const& scalar,vec3 vector){
    vec3 tmp;
    tmp.x = scalar * vector.x;
    tmp.y = scalar * vector.y;
    tmp.z = scalar * vector.z;
    return tmp;
}


//Image generation from https://stackoverflow.com/a/62946358
struct image{
    image(int width, int height): w(width), h(height), data(w*h*3)
    {}
    uint8_t & r(int x, int y) {return data[(x + y*w)*3 + 2];}
    uint8_t & g(int x, int y) {return data[(x + y*w)*3 + 1];}
    uint8_t & b(int x, int y) {return data[(x + y*w)*3 + 0];}

    int w,h;
    vector<uint8_t> data;
};

template<class Stream>
Stream & operator<<(Stream & out, image const& img){
    uint32_t w, h = img.h;
    w = img.w;
    uint32_t pad = w * -3 & 3; // NOLINT(hicpp-signed-bitwise)
    uint32_t total = 54 + 3*w*h + pad*h;
    uint32_t head[13] = {total, 0, 54, 40, w, h, (24<<16)|1}; // NOLINT(hicpp-signed-bitwise)
    char const* data = (char const*)img.data.data();

    out.write("BM", 2);
    out.write((char*)head, 52);

    for(int i=0; i<h ;i++){
        out.write(data + (3*w * i), 3 * w);
        out.write((char*)&pad, pad);
    }/*
    for(int i=h-1; i>=0;i--){
        out.write(data + (3*w * i), 3 * w);
        out.write((char*)&pad, pad);
    }*/
    return out;
}

struct Vertex {
    vec4 position;
    Vertex():position(vec3(0.0,0.0,0.0), DBL_MAX){}   //Default constructor
    Vertex(double x, double y, double z, double w) : position(vec4(x, y, z, w)){}
    Vertex(double x, double y, double z) : position(vec4(x, y, z, DBL_MAX)){}
    Vertex(vec3 in): position(in, DBL_MAX){}

    Vertex& operator+=(Vertex const& other){
        position.x += other.position.x;
        position.y += other.position.y;
        position.z += other.position.z;
        position.w += other.position.w;
        return *this;
    }

    glm::vec3 operator- (Vertex const& obj) const {
        vec3 res;
        res.x = position.x - obj.position.x;
        res.y = position.y - obj.position.y;
        res.z = position.z - obj.position.z;
        return res;
    }

    Vertex& operator*=(Vertex const& other){
        position.x *= other.position.x;
        position.y *= other.position.y;
        position.z *= other.position.z;
        position.w *= other.position.w;
        return *this;
    }

};

struct Direction {
    vec3 direction;
    Direction(): direction(vec3{}) {}
    Direction(Vertex in) : direction(in.position){}

};

struct Triangle;

struct Ray {
    Ray() : start{}, end{}, endTriangle{}, sphereIntersection(false), direction(){}
    Ray(Vertex vertex1, Vertex vertex2) : direction(vertex2 - vertex1){
        start = vertex1;
        end = vertex2;
        endTriangle = nullptr;
        sphereIntersection = false;
    }

    Ray(Vertex startP, Direction dir ) : direction(dir)
    {
        start = startP;
        endTriangle = nullptr;
        sphereIntersection = false;
    }

    Vertex start, end;
    ColorDbl color;
    Triangle* endTriangle;
    Vertex intersectionPoint{vec4(DBL_MAX,DBL_MAX,DBL_MAX,DBL_MAX)};
    Direction direction;
    int bounces = 20;
    bool sphereIntersection;

    void calculateDirection() {
        direction.direction = end-start;
    }
};

enum Type {
    LAMBERTIAN, SPECULAR, LIGHTSOURCE
};

struct Material{
    Type typ_;
    ColorDbl color_;
    Material() : color_(ColorDbl{}), typ_(Type{}){}

    Material(ColorDbl color, Type typ){
        color_ = color;
        typ_ = typ;
    }
};


struct Triangle {
    Vertex vec0, vec1, vec2;
    vec3 normal{};
    Material mat;

   // Triangle() : Triangle(Vertex(vec3{}),Vertex(vec3{}),Vertex(vec3{})){}
    Triangle() = default;
    Triangle(Vertex v1, Vertex v2, Vertex v3, Material mat_){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = normalize(cross((vec3)(v2.position-v1.position),(vec3)(v3.position-v1.position)));
        mat = mat_;
    }

    Triangle(Vertex v1, Vertex v2, Vertex v3, Vertex v4, Material mat_){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = normalize(cross((vec3)(v2.position-v1.position),(vec3)(v3.position-v1.position)));
        mat = mat_;
        if(dot(normal, (vec3)v4.position) > 0){
            normal = normalize(-normal);
        }
    }

    bool rayIntersection(Ray& intersectingRay, bool print = false) {

        //Möller-Trumbore
        vec3 T, E_1, E_2, D, P, Q;
        T = (vec3) (intersectingRay.start.position - vec0.position);
        E_1 = (vec3) (vec1.position - vec0.position);
        E_2 = (vec3) (vec2.position - vec0.position);
        D = intersectingRay.direction.direction;
        P = cross(D, E_2);
        Q = cross(T, E_1);
        vec3 tuv = vec3(dot(Q, E_2), dot(P, T), dot(Q, D)) / dot(P, E_1);
        if(print) {
            cout << "t: " << tuv.x << ", u: " << tuv.y << ", v: " << tuv.z << ", u+v = " << tuv.y + tuv.z << endl;
            cout << "vec0 = x: " << vec0.position.x << " y: " << vec0.position.y << " z: " << vec0.position.z << endl;
            cout << "Intersecting ray intersection point:" << endl << "x = "
                 << intersectingRay.intersectionPoint.position.x <<
                 ", y = " << intersectingRay.intersectionPoint.position.y <<
                 ", z = " << intersectingRay.intersectionPoint.position.z <<
                 ", w = " << intersectingRay.intersectionPoint.position.w << endl;

        }
        if (tuv.x > 0 && tuv.y >= 0 && tuv.z >= 0 && tuv.y + tuv.z <= (1.000005)) {
            if (intersectingRay.intersectionPoint.position.w >= tuv.x) {
                //Update intersectingRay here
                intersectingRay.endTriangle = this;
                intersectingRay.intersectionPoint = Vertex(tuv);
                intersectingRay.intersectionPoint.position.w = tuv.x;
                intersectingRay.color = this->mat.color_;
                intersectingRay.sphereIntersection = false;
            }
            return true;
        } else {
            return false;
        }
    }

    Triangle& operator= (Triangle const &Tri) = default;

};


struct Tetrahedron {
    vector<Triangle> triangles;
    Material mat = Material(ColorDbl(0,0,0), SPECULAR);

    Tetrahedron() = default;
    Tetrahedron(Vertex v1, Vertex v2, Vertex v3, Vertex v4){
        triangles = {Triangle(v1, v2, v3, v4, mat), Triangle(v1, v2, v4, v3,mat), Triangle(v1, v3, v4, v2,mat), Triangle(v3, v2, v4, v1, mat)};
    }

    bool rayIntersection(Ray& intersectingRay, bool print = false) {
        bool collision = false;
        for(auto & triangle : triangles) {
            //Triangle* tri = new Triangle(triangles[i]);
            if(triangle.rayIntersection(intersectingRay, print)){
                //intersectingRay.color = mat.color_;
                collision = true;
            }
            //delete(tri);
        }
        return collision;
    }
};


struct Sphere{
    vec3 centerOfSphere{};
    double rad{};
    Material mat;

    Sphere()= default;

    Sphere(double radius, vec3 centerSphere){
        centerOfSphere = centerSphere;
        rad = radius;
        mat = Material(ColorDbl(1.0,0,0),LAMBERTIAN);
    }

    //inspired by
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
    bool rayIntersection(Ray& ray, bool print = false) const {
        vec3 direction = vec3(ray.direction.direction);
        vec3 normDirection = normalize(direction);
        vec3 L = centerOfSphere - (vec3) ray.start.position;

        double tca = dot(L, normDirection);

        if (tca < DBL_EPSILON) return false;

        double radiusSquared = (rad * rad);
        double dsquared = dot(L, L) - (tca * tca);

        if (radiusSquared < dsquared) return false;     //check if ray hits sphere

        double thc = sqrt(radiusSquared - dsquared);

        double t0 = tca - thc;
        double t1 = tca + thc;

        if (t0 > t1) swap(t0, t1);

        if (t0 < DBL_EPSILON) {
            t0 = t1;
            if (t0 < DBL_EPSILON) {
                return false;
            }
        }
        if (ray.intersectionPoint.position.x > t0) {
            //Hit!
            ray.intersectionPoint = vec3(ray.start.position.x, ray.start.position.y, ray.start.position.z) + t0 * normDirection;
            ray.intersectionPoint.position.w = t0;
            ray.sphereIntersection = true;
            return true;
        }
        return false;

    }

    bool rayShadowIntersection(Ray& ray, bool print = false) const {
        vec3 direction = ray.direction.direction;
        vec3 normDirection = normalize(direction);
        vec3 L = centerOfSphere - (vec3) ray.start.position;

        double tca = dot(L, normDirection);

        if (tca < DBL_EPSILON) return false;

        double radiusSquared = (rad * rad);
        double dsquared = dot(L, L) - (tca * tca);

        if (radiusSquared < dsquared) return false;     //check if ray hits sphere

        double thc = sqrt(radiusSquared - dsquared);

        double t0 = tca - thc;
        double t1 = tca + thc;

        if (t0 > t1) swap(t0, t1);

        if (t0 < DBL_EPSILON) {
            t0 = t1;
            if (t0 < DBL_EPSILON) {
                return false;
            }
        }
        t0 = t0/length(direction);
        if (ray.intersectionPoint.position.x > t0) {
            //Hit!

            ray.intersectionPoint = vec3(ray.start.position.x, ray.start.position.y, ray.start.position.z) + t0 * normDirection;
            ray.intersectionPoint.position.w = t0;
            ray.color = mat.color_;
            ray.sphereIntersection = true;
            return true;
        }
        return false;

    }

    [[nodiscard]] vec3 getSphereNormal (vec4 middle) const{
        return normalize((vec3)middle-centerOfSphere);
    }

};

struct LightSource{
    vec3 position;
    vec3 normal{};
    vector<Triangle> triangles;
    LightSource(): position(vec3{}){
        Vertex corner[4];
        position = vec3(5,0,5);
        corner[0] = Vertex(position.x-0.5,position.y-0.5,position.z);
        corner[1] = Vertex(position.x-0.5,position.y+0.5,position.z);
        corner[2] = Vertex(position.x+0.5,position.y-0.5,position.z);
        corner[3] = Vertex(position.x+0.5,position.y+0.5,position.z);
        triangles.emplace_back(corner[0],corner[1],corner[2],Material(ColorDbl(1.0),LIGHTSOURCE));
        triangles.emplace_back(corner[3],corner[2],corner[1],Material(ColorDbl(1.0),LIGHTSOURCE));
        normal = triangles[0].normal;
    }
    bool rayIntersection(Ray& intersectingRay, bool print = false) {
        bool collision = false;
        for(auto & triangle : triangles) {
            if(triangle.rayIntersection(intersectingRay, print)){
                collision = true;
            }
        }
        return collision;
    }

};


struct Scene {
    Triangle walls[24]{};
    Tetrahedron tetras{};
    Sphere spheres{};
    LightSource areaLight{};
    random_device rd;
    mt19937 gen;
    bool randomRoof = false;
    Scene() {
        gen = (mt19937(rd()));
    }

    static Vertex getShadowrayStart(Ray& intersectingRay){
        double u = intersectingRay.intersectionPoint.position.y;
        double v = intersectingRay.intersectionPoint.position.z;
        Vertex start = Vertex{(1 - u - v) * intersectingRay.endTriangle->vec0.position +
                              u * intersectingRay.endTriangle->vec1.position +
                              v * intersectingRay.endTriangle->vec2.position};
        start.position = vec4((vec3) start.position -
                              0.0001 * intersectingRay.direction.direction, DBL_MAX);
        return start;
    }

    double getShadeValue(Ray& intersectingRay){
        uniform_real_distribution<> shadowDistributionX(areaLight.triangles[0].vec0.position.x,areaLight.triangles[1].vec0.position.x);
        uniform_real_distribution<> shadowDistributionY(areaLight.triangles[0].vec0.position.y,areaLight.triangles[1].vec0.position.y);
        Vertex start = getShadowrayStart(intersectingRay);
        int litShadowrays = 0;
        for (int i = 0; i < nr_shadowRays; i++) {
            Vertex end = Vertex(shadowDistributionX(gen), shadowDistributionY(gen), areaLight.position.z);
            Ray shadowRay{start, end};
            rayShadowIntersection(shadowRay);
            if (!(shadowRay.intersectionPoint.position.w <= 1.0 - DBL_EPSILON &&
                  shadowRay.intersectionPoint.position.w >= 0)) {
                litShadowrays++;
            }
        }
        return (double)((double) litShadowrays / (double)nr_shadowRays);
    }

    void rayIntersection(Ray& intersectingRay, bool print = false) { // NOLINT(misc-no-recursion)
        uniform_int_distribution<> distrib(0,1000);
        int russian = distrib(gen);
        intersectingRay.bounces--;

        //Cut ray if we end it on russian roulette      TODO: NOPE!
        if(russian < cutoff_value || intersectingRay.bounces == 0){
            if(intersectingRay.bounces < 10) cout << "Many bounces" << endl;
            return;
        }
        //Check if the ray intersects the area light source. Terminates the ray if it's the closest hit.
        if (areaLight.rayIntersection(intersectingRay)) {
            intersectingRay.color = ColorDbl(1.0);       // TODO: Better calculation here initial value L_0 doesn't matter? 1000? >1000?
        //If hit on the area light source, don't bother checking walls
        //If missing ALS, check walls
        } else {
            for (auto & wall : walls) {
                //If wall is hit, calculate end point through multiple shadowrays.
                if (wall.rayIntersection(intersectingRay)) {

                    Vertex start = getShadowrayStart(intersectingRay);
                    double fractionLight = getShadeValue(intersectingRay);
                    intersectingRay.color = intersectingRay.color * fractionLight;
                    double lightFactor = dot(intersectingRay.endTriangle->normal,normalize(vec3(areaLight.position-(vec3)start.position)));
                    //Clamp lightFactor
                    if(lightFactor < 0){
                        lightFactor = 0;
                    }
                    intersectingRay.color *= lightFactor;

                    break;
                }
            }
            //Walls have been checked, check for closer hit on Tetrahedron
            if (tetras.rayIntersection(intersectingRay)) {
                Ray *reflectedRay;
                reflectedRay = reflectRay(intersectingRay, intersectingRay.endTriangle->normal);
                double percentShadow = getShadeValue(intersectingRay);
                intersectingRay.color = intersectingRay.color*percentShadow;

                Vertex start = getShadowrayStart(intersectingRay);
                double lightFactor = dot(intersectingRay.endTriangle->normal,
                                         normalize(vec3(areaLight.position - (vec3) start.position)));
                intersectingRay.color = reflectedRay->color * lightFactor;
                delete (reflectedRay);

            }
            //Same as for Tetrahedron, walls are hit, find closer sphere hit
            if (spheres.rayIntersection(intersectingRay, print)) {
                vec3 normal = spheres.getSphereNormal(intersectingRay.intersectionPoint.position);
                Ray *reflectedRay;
                reflectedRay = reflectRay(intersectingRay,
                                          spheres.getSphereNormal(intersectingRay.intersectionPoint.position));
                double lightFactor = -dot(spheres.getSphereNormal(intersectingRay.intersectionPoint.position),
                                          normalize(intersectingRay.direction.direction));
                intersectingRay.color = reflectedRay->color * lightFactor;
                delete (reflectedRay);
            }
        }
    }


    Ray* reflectRay(Ray incomingRay, vec3 normal){ // NOLINT(misc-no-recursion)
        vec3 incomingDirection = incomingRay.end.position-incomingRay.start.position;
        vec4 start = (vec4(incomingRay.direction.direction, 1.0)*incomingRay.intersectionPoint.position.w)+incomingRay.start.position;
        start = vec4((vec3) start - 0.00001 * incomingRay.direction.direction, 1.0);
        vec3 reflectedDirection = reflect(incomingDirection, normal);
        Ray* reflectedRay = new Ray(Vertex(start),Direction(reflectedDirection));
        rayIntersection(*reflectedRay);
        return reflectedRay;
    }


    void rayShadowIntersection(Ray& intersectingRay, bool print = false) {
        if(areaLight.rayIntersection(intersectingRay)){
            double a = 0;
        } else {
            for (auto &wall : walls) {
                if (wall.rayIntersection(intersectingRay)) { break; }
            }
        }
        tetras.rayIntersection(intersectingRay);
        spheres.rayShadowIntersection(intersectingRay, print);
    }

    void setRandomRoof(bool b) {
        randomRoof = b;
    }
};

struct Pixel{
    Pixel() : color(ColorDbl{}){}
    Pixel(ColorDbl rgb) : color(rgb){}
    ColorDbl color;
    [[maybe_unused]] Ray intersectingRay;        //Prepare for array/list/pointer
};

struct Camera{
    Camera(Vertex left, Vertex right, int width, int height): leftEye(left), rightEye(right), image(width*height*3), perspective(false){}
    Vertex leftEye, rightEye;
    vector<Pixel> image;
    bool perspective;

    void setPerspective(bool p){
        perspective = p;
    }

    [[maybe_unused]] void togglePerspective(){
        perspective = !perspective;
    }

    [[nodiscard]] Vertex getEye() const {
        if(perspective){return leftEye;}
        else{return rightEye;}
    }
};


void createScene(Scene *world){

    //Vertex points r = roof, f = floor
    Vertex vrtx0r = Vertex(5.0, 0.0, 5.0, 1.0); //vrtx0r mitten toppen
    Vertex vrtx1r = Vertex(-3.0, 0.0, 5.0, 1.0); //vrtx1r
    Vertex vrtx2r = Vertex(0.0, 6.0, 5.0, 1.0); //vrtx2r
    Vertex vrtx3r = Vertex(10.0, 6.0, 5.0, 1.0); //vrtx3r
    Vertex vrtx4r = Vertex(13.0, 0.0, 5.0, 1.0); //vrtx4r
    Vertex vrtx5r = Vertex(10.0, -6.0, 5.0, 1.0); //vrtx5r
    Vertex vrtx6r = Vertex(0.0, -6.0, 5.0, 1.0); //vrtx6r

    Vertex vrtx0f = Vertex(5.0, 0.0, -5.0, 1.0); //vrtx0f mitten botten
    Vertex vrtx1f = Vertex(-3.0, 0.0, -5.0, 1.0); //vrtx1f
    Vertex vrtx2f = Vertex(0.0, 6.0, -5.0, 1.0); //vrtx2f
    Vertex vrtx3f = Vertex(10.0, 6.0, -5.0, 1.0); //vrtx3f
    Vertex vrtx4f = Vertex(13.0, 0.0, -5.0, 1.0); //vrtx4f
    Vertex vrtx5f = Vertex(10.0, -6.0, -5.0, 1.0); //vrtx5f
    Vertex vrtx6f = Vertex(0.0, -6.0, -5.0, 1.0); //vrtx6f


    Material white_lam = Material(ColorDbl(1.0,1.0,1.0), LAMBERTIAN);
    Material red_lam = Material(ColorDbl(1.0,0,0), LAMBERTIAN);
    Material blue_lam = Material(ColorDbl(0,0,1.0), LAMBERTIAN);
    Material green_lam = Material(ColorDbl(0,1.0,0), LAMBERTIAN);
    Material purple_lam = Material(ColorDbl(1.0,0,1.0), LAMBERTIAN);
    Material grey_lam = Material(ColorDbl(1.0,1.0,1.0), LAMBERTIAN);
    Material yellow_lam = Material(ColorDbl(1.0,1.0,0), LAMBERTIAN);
    Material teal_lam = Material(ColorDbl(0,1.0,1.0), LAMBERTIAN);

    uniform_int_distribution<> distrib2(50,150);
    double val = 1.0;
    if(world->randomRoof){
        val = (double)distrib2 (world->gen)/255.99;
    }
    cout << "Grey value: " << val << endl;
    grey_lam.color_ = ColorDbl(val);

    //Roof = Grey/white
    world->walls[0] = Triangle(vrtx0r, vrtx1r, vrtx2r, grey_lam);
    world->walls[1] = Triangle(vrtx0r, vrtx2r, vrtx3r, grey_lam);
    world->walls[2] = Triangle(vrtx0r, vrtx3r, vrtx4r, grey_lam);
    world->walls[3] = Triangle(vrtx0r, vrtx4r, vrtx5r, grey_lam);
    world->walls[4] = Triangle(vrtx0r, vrtx5r, vrtx6r, grey_lam);
    world->walls[5] = Triangle(vrtx0r, vrtx6r, vrtx1r, grey_lam);

    //Floor = White
    world->walls[6] = Triangle(vrtx0f, vrtx2f, vrtx1f, white_lam);
    world->walls[7] = Triangle(vrtx0f, vrtx3f, vrtx2f, white_lam);
    world->walls[8] = Triangle(vrtx0f, vrtx4f, vrtx3f, white_lam);
    world->walls[9] =  Triangle(vrtx0f, vrtx5f, vrtx4f, white_lam);
    world->walls[10] = Triangle(vrtx0f, vrtx6f, vrtx5f, white_lam);
    world->walls[11] = Triangle(vrtx0f, vrtx1f, vrtx6f, white_lam);

    //Walls
    //Wall 1 = Red
    world->walls[12] = Triangle(vrtx2r, vrtx1r, vrtx2f, red_lam);
    world->walls[13] = Triangle(vrtx1f, vrtx2f, vrtx1r, red_lam);

    //Wall 2 = Yellow
    world->walls[14] = Triangle(vrtx3r, vrtx2r, vrtx3f, yellow_lam);
    world->walls[15] = Triangle(vrtx2f, vrtx3f, vrtx2r, yellow_lam);

    //Wall 3 = Green
    world->walls[16] = Triangle(vrtx4r, vrtx3r, vrtx4f, green_lam);
    world->walls[17] = Triangle(vrtx3f, vrtx4f, vrtx3r, green_lam);

    //Wall 4 = Teal
    world->walls[18] = Triangle(vrtx5r, vrtx4r, vrtx4f, teal_lam);
    world->walls[19] = Triangle(vrtx5r, vrtx4f, vrtx5f, teal_lam);

    //Wall 5 = Blue
    world->walls[20] = Triangle(vrtx6r, vrtx5r, vrtx6f, blue_lam);
    world->walls[21] = Triangle(vrtx5f, vrtx6f, vrtx5r, blue_lam);

    //Wall 6 = Purple
    world->walls[22] = Triangle(vrtx1r, vrtx6r, vrtx6f, purple_lam);
    world->walls[23] = Triangle(vrtx1r, vrtx6f, vrtx1f, purple_lam);
/*
    for (int i = 0; i < 24; i++)
    {
        world->walls[i].mat.color_ = blue_lam.color_;
    }*/

    //create Tetrahedron
    vec3 tetraCenter = vec3(8,-2,-1);
    double size = 0.8;
    world->tetras = Tetrahedron(vec3(tetraCenter.x+0,tetraCenter.y+0,tetraCenter.z+2*size),
                                vec3(tetraCenter.x-size,tetraCenter.y-size,tetraCenter.z-3*size),
                                vec3(tetraCenter.x+2*size,tetraCenter.y+3*size,tetraCenter.z-3*size),
                                vec3(tetraCenter.x+2*size,tetraCenter.y-3*size,tetraCenter.z-3*size));
    //Create Sphere
    world->spheres = Sphere(1.8,vec3(4,2,0));

}

/*  TODO:
 *  1. Reflection rays      check?
 *  2. Perfect reflection objects   check?
 *  3. Physically correct light calculation
*/

int main() {
    cout << "This is the start of our RayTracer!" << endl;
    const int xWidth = 1920;
    const int yWidth = 1080;


    Scene world;
    const int subPixelsPerAxis = 1;        //Anti-aliasing level, subPixelsPerAxis = k, k = [1,2,3,...];
    const int raysPerPixel = subPixelsPerAxis*subPixelsPerAxis;
    const int width = 800;
    const int height = 800;
    const bool randomRays = true;            //Select if ray directions are randomized/dithered
    const bool brightSpots = true;           //Are there bright spots in scene?

    const double pixelSize = 2.0 / width;
    Camera cam{Vertex(-2, 0, 0), Vertex(-1, 0, 0), width, height};
    cam.setPerspective(false);
    world.setRandomRoof(false);      //Use this to select if roof color randoms between runs.

    image img{width, height};

    cout << "Creating Scene" << endl;
    createScene(&world);

    double maxIntensity = 0;
    random_device randomDevice;
    mt19937 gen(randomDevice());
    uniform_real_distribution<> distrib(0, 1.0);
    time_t timer;


    double seconds = time(&timer);
    for (int i = 0; i < width; i++) {
        double random = distrib(gen);
        if(random > 0.95){
            double progress = (double)(i*height)/(width*height);
            int progressPercent = (int) (progress*100);
            cout << "Rendering ..." << progressPercent << "%." << endl;
        }
        for (int j = 0; j < height; j++) {
            cam.image[i * width + j] = ColorDbl(0.4, 0.4, 0.4);
            ColorDbl pixelAvg{};
            for (int r = 0; r < raysPerPixel; r++) {
                double dy = (randomRays) ? distrib(gen) : 0;
                double dz = (randomRays) ? distrib(gen) : 0;

                Vertex end = Vertex(0,
                                     (i - floor(width/2) + (double) (r % subPixelsPerAxis) / (double) subPixelsPerAxis +
                                      dy / subPixelsPerAxis) * pixelSize,
                                     (j - floor(height/2) + floor(r / subPixelsPerAxis) / subPixelsPerAxis +
                                      dz / subPixelsPerAxis) * pixelSize);
                Ray* current = new Ray(cam.getEye(),end);
                world.rayIntersection(*current);

                if(current != NULL) {
                    pixelAvg += (current->color);
                    if (current->color.r > maxIntensity) { maxIntensity = current->color.r; }
                    if (current->color.g > maxIntensity) { maxIntensity = current->color.g; }
                    if (current->color.b > maxIntensity) { maxIntensity = current->color.b; }
                    delete (current);
                }
            }
            pixelAvg /= raysPerPixel;
            if(brightSpots) pixelAvg = sqrt(pixelAvg);
            cam.image[i * width + j] = Pixel(pixelAvg);
        }
    }
    seconds = time(&timer) - seconds;
    cout << "Rendering ...100%" << ", rendering took " << seconds << " seconds." << endl;
    cout << "MaxIntensity: " << maxIntensity << endl;

    if(brightSpots) maxIntensity = sqrt(maxIntensity);
    //write cam.image to output img
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            img.r(i, j) = 255.99 * cam.image[(height - 1 - i) * width + j].color.r/maxIntensity;
            img.g(i, j) = 255.99 * cam.image[(height - 1 - i) * width + j].color.g/maxIntensity;
            img.b(i, j) = 255.99 * cam.image[(height - 1 - i) * width + j].color.b/maxIntensity;
        }
    }

    cout << "Generated an image!" << endl;
    string filename = "Scene_aa_" + to_string((raysPerPixel)) + ".bmp";
    ofstream(filename, ios_base::out | ios_base::binary) << img;
    cout << "Wrote file " << filename << endl;

    return 0;
}
