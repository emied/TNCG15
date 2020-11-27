#include <iostream>
#include "../glm/glm/glm.hpp"
#include <vector>
#include <fstream>
#include <random>
#include <list>
#include <string>
#include <ctime>

using namespace std;
using namespace glm;
long double EPSILON = 0.01;

struct ColorDbl{

    double r, g, b;
    ColorDbl(): r(0.0), g(0.0), b(0.0) {} //default constructor
    ColorDbl(double red, double green, double blue) : r(red), g(green), b(blue){}
    ColorDbl(vec3 color)
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
};

ColorDbl operator*(ColorDbl const& color, double scalar){
    ColorDbl ret{};
    ret.r = color.r*scalar;
    ret.g = color.g*scalar;
    ret.b = color.b*scalar;
    return ret;
}

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
    uint32_t w = img.w, h = img.h;
    uint32_t pad = w * -3 & 3;
    uint32_t total = 54 + 3*w*h + pad*h;
    uint32_t head[13] = {total, 0, 54, 40, w, h, (24<<16)|1};
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
    Vertex():position(vec3(0.0,0.0,0.0), 1.0){}   //Default constructor
    Vertex(double x, double y, double z, double w) : position(vec4(x, y, z, w)){}
    Vertex(double x, double y, double z) : position(vec4(x, y, z, 1.0)){}
    Vertex(vec3 in): position(in, 1.0){}

    Vertex& operator+=(Vertex const& other){
        position.x += other.position.x;
        position.y += other.position.y;
        position.z += other.position.z;
        position.w += other.position.w;
        return *this;
    }

    glm::vec3 operator- (Vertex const& obj) {
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
    glm::vec3 direction;
    Direction(): direction(glm::vec3{}) {}
    Direction(Vertex in) : direction(in.position){}

};

struct Triangle;

struct Ray {
    Ray() : start{}, end{}, endTriangle{} {}
    Ray(Vertex vertex1, Vertex vertex2) {
        start = vertex1;
        end = vertex2;
        endTriangle = nullptr;
        direction.direction = normalize(vertex2 - vertex1);
    }

    Ray(Vertex startP, Direction dir )
    {
        start = startP;
        endTriangle = nullptr;
        direction = dir;
    }

    Vertex start, end;
    ColorDbl color;
    Triangle* endTriangle;
    Vertex intersectionPoint{vec3(DBL_MAX,DBL_MAX,DBL_MAX)};
    Direction direction;
};

enum Type {
    LAMBERTIAN, SPECULAR, LIGHTSOURCE, ORENNAYAR, GLASS
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
    vec3 normal;
    Material mat;

   // Triangle() : Triangle(Vertex(vec3{}),Vertex(vec3{}),Vertex(vec3{})){}
    Triangle() {}
    Triangle(Vertex v1, Vertex v2, Vertex v3, Material mat_){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = normalize(cross((vec3)(v2.position-v1.position),(vec3)(v3.position-v1.position)));
        double d = dot(normal, (vec3)vec0.position);
        mat = mat_;
    }

    Triangle(Vertex v1, Vertex v2, Vertex v3, Vertex v4, Material mat_){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = cross((vec3)(v2.position-v1.position),(vec3)(v3.position-v1.position));
        double d = dot(normal, (vec3)vec0.position);
        mat = mat_;
        if(dot(normal, (vec3)v4.position) > 0){
            normal = -normal;
        }
    }

    bool rayIntersection(Ray& intersectingRay, bool print = false) {

        //Möller-Trumbore
        vec3 T, E_1, E_2, D, P, Q;
        T = (vec3) (intersectingRay.start.position - vec0.position);
        E_1 = (vec3) (vec1.position - vec0.position);
        E_2 = (vec3) (vec2.position - vec0.position);
        D = intersectingRay.end.position - intersectingRay.start.position;
        P = cross(D, E_2);
        Q = cross(T, E_1);
        vec3 tuv = vec3(dot(Q, E_2), dot(P, T), dot(Q, D)) / dot(P, E_1);
        if(print){
            cout << "t: " << tuv.x << ", u: " << tuv.y << ", v: " << tuv.z << ", u+v = " << tuv.y + tuv.z << endl;
            cout << "vec0 = x: " << vec0.position.x << " y: " << vec0.position.y << " z: " << vec0.position.z << endl;
            cout << "Intersecting ray intersection point:" << endl << "x = " << intersectingRay.intersectionPoint.position.x <<
                                                                ", y = " << intersectingRay.intersectionPoint.position.y <<
                                                                ", z = " << intersectingRay.intersectionPoint.position.z << endl;

        }
        if (tuv.x > 0 && tuv.y >= 0 && tuv.z >= 0 && tuv.y + tuv.z <= (1.000005)) {
            if (intersectingRay.intersectionPoint.position.x >= tuv.x) {
                //Update intersectingRay here
                intersectingRay.endTriangle = this;
                intersectingRay.intersectionPoint = Vertex(tuv);
                intersectingRay.color = this->mat.color_;
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
    Material mat = Material(ColorDbl(255,80,255), LAMBERTIAN);

    Tetrahedron() {}
    Tetrahedron(Vertex v1, Vertex v2, Vertex v3, Vertex v4){
        triangles = {Triangle(v1, v2, v3, v4, mat), Triangle(v1, v2, v4, v3,mat), Triangle(v1, v3, v4, v2,mat), Triangle(v3, v2, v4, v1, mat)};
    }

    bool rayIntersection(Ray& intersectingRay, bool print = false) {
        bool collision = false;
        for(int i=0; i < triangles.size(); i++) {
            Triangle* tri = new Triangle(triangles[i]);
            if(tri->rayIntersection(intersectingRay, print)){
                intersectingRay.color = mat.color_;
                collision = true;
            }
        }
        return collision;
    }
};

struct Sphere{
    vec3 centerOfSphere;
    double rad;
    Material mat;

    Sphere(){}

    Sphere(double radius, vec3 centerSphere){
        centerOfSphere = centerSphere;
        rad = radius;
        mat = Material(ColorDbl(255,0,0),LAMBERTIAN);
    }

    //inspired by
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
    bool sphereRayIntersection(Ray& ray){
        double t = DBL_MAX;
        double r = rad;
        vec3 direction = normalize(vec3(ray.end - ray.start));
        vec3 L = centerOfSphere - direction;

        double tca = dot((L), direction);

        if(tca < EPSILON) return false;

        double tcaSquared = (tca*tca);
        double radiusSquared = (r*r);
        float dsquared = dot(L, L) - tcaSquared;

        if(radiusSquared < dsquared) return false;     //check if ray hits sphere

        double thc = sqrt(radiusSquared - dsquared);

        double t0 = tca-thc;
        double t1 = tca+thc;

        if(t0 > t1) swap(t0,t1);

        if(t0 < EPSILON && t1 < EPSILON){
            return false;
        }

        ray.intersectionPoint = (ray.start.position.x,ray.start.position.y,ray.start.position.z) + direction;
        //hit
        t = t0;
        float f = (float) t;
        ray.color = mat.color_;
        return true;

    }

    vec3 getSphereNormal (Vertex middle){
        return normalize(middle-centerOfSphere);
    }

};

struct Scene;

struct Scene {
    Triangle walls[24]{};
    Tetrahedron tetras{};
    Sphere spheres{};
    bool randomRoof = false;
    Scene() = default;


    void rayIntersection(Ray& intersectingRay, bool print = false){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> distrib(0, 1.0);
        double random;
        for(int i = 0; i < 24 ;i++){

            if ((i == 18 || i == 19 || i == 9) && print) cout << "Checking wall i = " << i << endl;
            if(walls[i].rayIntersection(intersectingRay,((i == 18 || i == 19 || i == 9) && print))){break;}
        }
        //tetras.rayIntersection(intersectingRay, random > 0.9999);
        tetras.rayIntersection(intersectingRay,print);
        //spheres.sphereRayIntersection(intersectingRay);
    }

    void setRandomRoof(bool b) {
        randomRoof = b;
    }
};

struct Pixel{
    Pixel() : color(ColorDbl{}){}
    Pixel(ColorDbl rgb) : color(rgb){}
    ColorDbl color;
    Ray intersectingRay;        //Prepare for array/list/pointer
};

struct Camera{
    Camera(Vertex left, Vertex right, int width, int height): leftEye(left), rightEye(right), image(width*height*3), perspective(false){}
    Vertex leftEye, rightEye;
    vector<Pixel> image;
    bool perspective;

    void setPerspective(bool p){
        perspective = p;
    }

    void togglePerspective(){
        perspective = !perspective;
    }

    Vertex getEye() {
        if(perspective){return leftEye;}
        else{return rightEye;}
    }
};

struct LightSource{
    vec3 position;
    ColorDbl color;
    LightSource(): position(vec3{}), color(vec3{}){}
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


    Material white_lam = Material(ColorDbl(255,255,255), LAMBERTIAN);
    Material red_lam = Material(ColorDbl(255,0,0), LAMBERTIAN);
    Material blue_lam = Material(ColorDbl(0,0,255), LAMBERTIAN);
    Material green_lam = Material(ColorDbl(0,255,0), LAMBERTIAN);
    Material purple_lam = Material(ColorDbl(255,0,255), LAMBERTIAN);
    Material grey_lam = Material(ColorDbl(254,254,254), LAMBERTIAN);
    Material yellow_lam = Material(ColorDbl(255,255,0), LAMBERTIAN);
    Material teal_lam = Material(ColorDbl(0,255,255), LAMBERTIAN);
/*
    random_device rd;
    mt19937  gen(rd());
    uniform_int_distribution<> distrib(50,200);
    int val = 50;
    if(world->randomRoof){
        val = distrib ( gen);
    }
    cout << "Grey value: " << val << endl;
     colors[0] = ColorDbl(val,val,val); */

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
    world->walls[13] = Triangle(vrtx1f, vrtx1r, vrtx2f, red_lam);

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
    world->tetras = Tetrahedron(vec3(6,0,-3),
                                vec3(8,2,-3),
                                vec3(8,-2,-3),
                                vec3(7,0,0));
    //Create Sphere
    //world->spheres = Sphere(1,vec3(7,3,0));

}


int main() {
    cout << "This is the start of our RayTracer!" << endl;
    const int xWidth = 1920;
    const int yWidth = 1080;


    Scene world;
    const int subPixelsPerAxis = 1;        //Anti-aliasing level, subPixelsPerAxis = k, k = [1,2,3,...];
    const int raysPerPixel = subPixelsPerAxis*subPixelsPerAxis;
    const int width = 800;
    const int height = 800;
    const double pixelSize = 2.0 / width;
    Camera cam{Vertex(-2, 0, 0), Vertex(-1, 0, 0), width, height};
    cam.setPerspective(false);
    world.setRandomRoof(true);      //Use this to select if roof color randoms between runs.

    image img{width, height};

    cout << "Creating Scene" << endl;
    createScene(&world);
    double maxIntensity = 0;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> distrib(0, 1.0);
    time_t timer;

    //Add point light
    LightSource light;
    light.color = vec3{1.0,1.0,1.0};
    light.position = vec3{5,0,1};            //behind light
    //light.position = vec3{3,-1,1};              //front light
    //light.position = cam.getEye().position;

    double seconds = time(&timer);
    for (int i = 0; i < width; i++) {
        double random = distrib(gen);
        if(random > 0.95){
            double progress = (double)(i*width)/(width*height);
            int progressPercent = progress*100;
            cout << "Rendering ..." << progressPercent << "%." << endl;
        }
        for (int j = 0; j < height; j++) {
            cam.image[i * width + j] = ColorDbl(100, 100, 100);
            ColorDbl pixelAvg{};
            for (int r = 0; r < raysPerPixel; r++) {
                Ray current{};
                double dy = distrib(gen);
                double dz = distrib(gen);
                //dy = dz = 0;      //disable random rays
                current.start = cam.getEye();

                current.end = Vertex(0,
                                     (i - 400 + (double)(r%subPixelsPerAxis)/(double)subPixelsPerAxis + dy/subPixelsPerAxis) * pixelSize,
                                     (j - 400 + floor(r/subPixelsPerAxis)/subPixelsPerAxis + dz/subPixelsPerAxis) * pixelSize);
                world.rayIntersection(current, false);

                Ray shadow{};
                double u , v;
                ColorDbl shadowOrNot = {1.0,1.0,1.0};

                //bool sph = world.spheres.sphereRayIntersection(current);
                bool sph = false;
                //Shadow for triangle objects
                double distanceToLight;

                if (!sph) {
                    u = current.intersectionPoint.position.y;
                    v = current.intersectionPoint.position.z;
                    shadow.start = Vertex{(1 - u - v) * current.endTriangle->vec0.position +
                                          u * current.endTriangle->vec1.position +
                                          v * current.endTriangle->vec2.position};
                    //Move start out of object
                    shadow.start.position = vec4((vec3) shadow.start.position - 0.00001 * normalize(current.end.position - current.start.position), 1.0);

                    /*random = distrib(gen);
                    if (random > 0.999) {
                        cout << "Ray index: (" << i << "," << j << ")." << endl;
                        cout << "Current ray end triangle vec0:" << endl <<
                        "x: " << current.endTriangle->vec0.position.x <<
                        " y: " << current.endTriangle->vec0.position.y <<
                        " z: " << current.endTriangle->vec0.position.z <<
                        " w: " << current.endTriangle->vec0.position.w << endl;
                        cout << "Shadow start pos:" << endl << "x: " << shadow.start.position.x << " y: " << shadow.start.position.y
                             << " z: " << shadow.start.position.z << " w: " << shadow.start.position.w << endl << endl;
                    }*/
                }

                    //Shadow for implicit objects
                else if(sph) {
                    random = distrib(gen);
                    shadow.start = Vertex{(vec3) current.intersectionPoint.position + world.spheres.centerOfSphere};
                    /*if (random > 0.99) {
                        cout << "Ray index: (" << i << "," << j << ")." << endl;
                        cout << "Shadow start pos, x: " << shadow.start.position.x << " y: " << shadow.start.position.y
                             << " z: " << shadow.start.position.z << endl;
                    }*/
                }

                shadow.end = light.position;
                world.rayIntersection(shadow);
                bool shadedRay;
                if(!sph){shadedRay = (shadow.intersectionPoint.position.x <= 1.0 + DBL_EPSILON && shadow.intersectionPoint.position.x >= 0) ;}
                if(sph){shadedRay = true;}
                vec3 dist = shadow.end - shadow.start;
                distanceToLight = sqrt(pow(dist.x, 2) + pow(dist.y, 2) + pow(dist.z, 2));

                shadowOrNot = (shadedRay) ? ColorDbl{.1,.1,.1} : ColorDbl{1.0, 1.0, 1.0};

                shadowOrNot *= (1.0/ distanceToLight);
                pixelAvg += (current.color) * shadowOrNot;

                random = distrib(gen);
                if(random > 1.9999){
                    cout << "Ray index: (" << i << "," << j << ")." << endl;
                    cout << "Shadow start pos:" << endl << "x: " << shadow.start.position.x <<
                            " y: " << shadow.start.position.y <<
                            " z: " << shadow.start.position.z <<
                            " w: " << shadow.start.position.w << endl;
                    cout << "Shaded ray t-value: " << shadow.intersectionPoint.position.x << ", shaded? = "<< ((shadedRay)?"Yes":"No") << endl;
                    cout << "shadowOrNot  " << shadowOrNot.r << " " << shadowOrNot.g << " " << shadowOrNot.b << endl;
                    cout << "Current triangle color: " << endl << "r: " << current.endTriangle->mat.color_.r <<
                    ", g: " << current.endTriangle->mat.color_.g << ", b: " << current.endTriangle->mat.color_.b << endl;
                    cout << "Collision triangle  normal: " << endl << "x: " << current.endTriangle->normal.x <<
                    ", y: " << current.endTriangle->normal.y << ", z: " << current.endTriangle->normal.z << endl << endl;
                }
                if(!sph) {
                    if (current.color.r > maxIntensity) { maxIntensity = current.color.r; }
                    if (current.color.g > maxIntensity) { maxIntensity = current.color.g; }
                    if (current.color.b > maxIntensity) { maxIntensity = current.color.b; }
                }
            }
            pixelAvg /= raysPerPixel;
            cam.image[i * width + j] = Pixel(pixelAvg);
        }
    }
    seconds = time(&timer) - seconds;
    cout << "Rendering ...100%" << ", rendering took " << seconds << " seconds." << endl;

    //write cam.image to output img
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
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
