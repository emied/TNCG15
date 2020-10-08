#include <iostream>
#include "../glm/glm/glm.hpp"
#include <vector>
#include <fstream>
#include <random>
#include <list>
#include <string>

using namespace std;
using namespace glm;

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
};

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

void generateGradientImage(image *image);

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
    }
    return out;
}

struct Vertex {
    vec4 position;
    Vertex():position(vec3(0.0,0.0,0.0), 1.0){}   //Default constructor
    Vertex(double x, double y, double z, double w) : position(vec4(x, y, z, w)){}
    Vertex(double x, double y, double z) : position(vec4(x, y, z, 1.0)){}
    Vertex(vec3 in): position(in, 1.0){}
};

struct Direction {
    vec3 direction;
    Direction(): direction(vec3{}){}
    Direction(Vertex in) : direction(in.position){}

};

struct Triangle;

struct Ray {
    //KAOS
    Ray() : start{}, end{} {}
    Ray(Vertex vertex1, Vertex vertex2) {
        start = vertex1;
        end = vertex2;
    }

    Vertex start, end;
    ColorDbl color;
    Triangle* endTriangle;
    Vertex intersectionPoint{vec3(DBL_MAX,DBL_MAX,DBL_MAX)};


};

//KAOS
struct Intersections{
    Vertex pos;
    //och sluttriangel
};

struct Triangle {

    ColorDbl color;
    Vertex vec0, vec1, vec2;
    vec3 normal{};
    double d{};


    Triangle() : Triangle(Vertex(vec3{}),Vertex(vec3{}),Vertex(vec3{})){}

    Triangle(Vertex v1, Vertex v2, Vertex v3){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = cross((vec3)(v2.position-v1.position),(vec3)(v3.position-v1.position));
        d = dot(normal, (vec3)vec0.position);
        color = ColorDbl();
    }

    Triangle(Vertex v1, Vertex v2, Vertex v3, Vertex v4){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = cross((vec3)(v2.position-v1.position),(vec3)(v3.position-v1.position));
        d = dot(normal, (vec3)vec0.position);
        color = ColorDbl();
        if(dot(normal, (vec3)v4.position) > 0){
            normal = -normal;
        }
    }


    bool rayIntersection(Ray& intersectingRay) {
        //Möller-Trumbore
        vec3 T, E_1, E_2, D, P, Q;
        T = (vec3) (intersectingRay.start.position - vec0.position);
        E_1 = (vec3) (vec1.position - vec0.position);
        E_2 = (vec3) (vec2.position - vec0.position);
        D = intersectingRay.end.position - intersectingRay.start.position;
        P = cross(D, E_2);
        Q = cross(T, E_1);
        vec3 tuv = vec3(dot(Q, E_2), dot(P, T), dot(Q, D)) / dot(P, E_1);

        //KAOS
        std::list<Intersections> IntersectionsTmp = {};
        Intersections inter;

        if(tuv.x > 0 && tuv.y > 0 && tuv.z > 0 && tuv.y+tuv.z <=1.0){
            /*if(rand() % 100 > 98){
                cout << "Vector print; t: " << tuv.x << ", u: " << tuv.y << ", v: " << tuv.z << endl;
            }*/
            if(intersectingRay.intersectionPoint.position.x > tuv.x) {
                //Update intersectingRay here
                intersectingRay.endTriangle = this;
                intersectingRay.intersectionPoint = Vertex(tuv);
                intersectingRay.color = this->color;

                //KAOS
                inter.pos = Vertex(tuv);
                IntersectionsTmp.push_back(inter);

            }
            return true;
        } else {
            return false;
        }
    }

};


struct Tetrahedron {
    ColorDbl color;
    vector<Triangle> triangles;
    Tetrahedron() {
        color = ColorDbl{};
        //triangles = vector<Triangle>{};

    }
    Tetrahedron(Vertex v1, Vertex v2, Vertex v3, Vertex v4,ColorDbl c){
        color = c;
        ColorDbl c1,c2,c3,c4;
        triangles = {Triangle(v1, v2, v3, v4), Triangle(v1, v2, v4, v3), Triangle(v1, v3, v4, v2), Triangle(v2, v3, v4, v1)};
    }

    Tetrahedron(Vertex v1, Vertex v2, Vertex v3, Vertex v4) : Tetrahedron(v1,v2,v3,v4, ColorDbl{}){}

    bool rayIntersection(Ray& intersectingRay) {
        bool collision = false;
        for(Triangle tri : triangles) {
            if(tri.rayIntersection(intersectingRay)){
                intersectingRay.color = color;
                collision = true;

            }
        }
        return collision;
    }
};

struct Scene;

void createScene(Scene *world);

struct Scene {
    Triangle triangles[24]{};
    Tetrahedron tetras{};
    //Vertex vertices[14];
    ColorDbl colors[8];
    bool randomRoof;
    Scene() = default;//{
        //createScene(this);
    //}


    void rayIntersection(Ray& intersectingRay){
        for(int i = 0; i < 24 ;i++){
            if(triangles[i].rayIntersection(intersectingRay)){break;}
        }
        tetras.rayIntersection(intersectingRay);

    }

    void setRandomRoof(bool b) {
        randomRoof = b;
    }
};

//KAOS
//    vec3 CastShadowRay(Scene scen, vec3 hitSurface, vec3 lightSource){
//
//        Vertex startingPoint = Vertex(hitSurface);
//        Vertex lightPoint = lightSource;
//        Ray ShadowRay = Ray(startingPoint, lightPoint);
//
//        //skicka shadow ray
//        scen.rayIntersection(ShadowRay);
//        //lägg in värden i en lista
//        //list<Ray> intersections =
//
//
//
//
//        double distanceLight = distance(hitSurface, lightSource);
//        double distanceIntersection;
//    }


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

    //Roof = Random Grey
    random_device rd;
    mt19937  gen(rd());
    uniform_int_distribution<> distrib(50,200);
    int val = 50;
    if(world->randomRoof){
        val = distrib ( gen);
    }
    cout << "Grey value: " << val << endl;
    world->colors[0] = ColorDbl(val,val,val);
    world->triangles[0] = Triangle(vrtx0r, vrtx1r, vrtx2r);
    world->triangles[1] = Triangle(vrtx0r, vrtx2r, vrtx3r);
    world->triangles[2] = Triangle(vrtx0r, vrtx3r, vrtx4r);
    world->triangles[3] = Triangle(vrtx0r, vrtx4r, vrtx5r);
    world->triangles[4] = Triangle(vrtx0r, vrtx5r, vrtx6r);
    world->triangles[5] = Triangle(vrtx0r, vrtx6r, vrtx1r);

    //Floor = White
    world->colors[1] = ColorDbl(255,255,255);
    world->triangles[6] = Triangle(vrtx0f, vrtx2f, vrtx1f);
    world->triangles[7] = Triangle(vrtx0f, vrtx3f, vrtx2f);
    world->triangles[8] = Triangle(vrtx0f, vrtx4f, vrtx3f);
    world->triangles[9] =  Triangle(vrtx0f, vrtx5f, vrtx4f);
    world->triangles[10] = Triangle(vrtx0f, vrtx6f, vrtx5f);
    world->triangles[11] = Triangle(vrtx0f, vrtx1f, vrtx6f);


    //Walls
    //Wall 1 = Red
    world->colors[2] = ColorDbl(255,0,0);
    world->triangles[12] = Triangle(vrtx2r, vrtx1r, vrtx2f);// ska vara 2,1,8, blir dock knas då
    world->triangles[13] = Triangle(vrtx1f, vrtx1r, vrtx2f);

    //Wall 2 = Yellow
    world->colors[3] = ColorDbl(255,255,0);
    world->triangles[14] = Triangle(vrtx3r, vrtx2r, vrtx3f);
    world->triangles[15] = Triangle(vrtx2f, vrtx2r, vrtx3f);

    //Wall 3 = Green
    world->colors[4] = ColorDbl(0,255,0);
    world->triangles[16] = Triangle(vrtx4r, vrtx3r, vrtx4f);
    world->triangles[17] = Triangle(vrtx3f, vrtx4f, vrtx3r);

    //Wall 4 = Teal
    world->colors[5] = ColorDbl(0,255,255);
    world->triangles[18] = Triangle(vrtx5r, vrtx4r, vrtx4f);
    world->triangles[19] = Triangle(vrtx5r, vrtx4f, vrtx5f);

    //Wall 5 = Blue
    world->colors[6] = ColorDbl(0,0,255);
    world->triangles[20] = Triangle(vrtx6r, vrtx5r, vrtx6f);
    world->triangles[21] = Triangle(vrtx5f, vrtx6f, vrtx5r);

    //Wall 6 = Purple
    world->colors[7] = ColorDbl(255,0,255);
    world->triangles[22] = Triangle(vrtx1r, vrtx6r, vrtx6f);
    world->triangles[23] = Triangle(vrtx1r, vrtx6f, vrtx1f);

    //Actually assign colours
    for(int j = 0; j < 24; j++){
        if(j < 12) {
            world->triangles[j].color = world->colors[j/6];
        } else {
            world->triangles[j].color = world->colors[(int)floor(j/2) - 4];
        }
    }
    //create Tetrahedron
    world->tetras = Tetrahedron(vec3(8,0,2),
                                vec3(7,0,-3),
                                vec3(9,3,-1),
                                vec3(9,-2,-1),
                                ColorDbl(0,150,0));


    //Add point light
    LightSource().color = vec3{1.0,1.0,1.0};
    LightSource().position = vec3{5.0,5.0,5.0};

}


int main() {
    cout << "This is the start of our RayTracer!" << endl;
    const int xWidth = 1920;
    const int yWidth = 1080;


    Scene world;
    const int subPixelsPerAxis = 3;        //Anti-aliasing level, subPixelsPerAxis = k, k = [1,2,3,...];
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
                current.start = cam.getEye();
                current.end = Vertex(0,
                                     (i - 401 + r%subPixelsPerAxis + dy/subPixelsPerAxis) * pixelSize,
                                     (j - 401 + floor(r/subPixelsPerAxis)/subPixelsPerAxis + dz/subPixelsPerAxis) * pixelSize);
                world.rayIntersection(current);
                pixelAvg += current.color;
                if (current.color.r > maxIntensity) { maxIntensity = current.color.r; }
                if (current.color.g > maxIntensity) { maxIntensity = current.color.g; }
                if (current.color.b > maxIntensity) { maxIntensity = current.color.b; }
            }
            pixelAvg /= raysPerPixel;
            cam.image[i * width + j] = Pixel(pixelAvg);
        }
    }
    cout << "Rendering ...100%" << endl;

    //write cam.image to output img
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            img.r(i, j) = cam.image[(height - 1 - i) * width + j].color.r;
            img.g(i, j) = cam.image[(height - 1 - i) * width + j].color.g;
            img.b(i, j) = cam.image[(height - 1 - i) * width + j].color.b;
        }
    }

    cout << "Generated an image!" << endl;
    string filename = "Scene_aa_" + to_string((raysPerPixel)) + ".bmp";
    ofstream(filename, ios_base::out | ios_base::binary) << img;
    cout << "Wrote file " << filename << endl;

    return 0;
}




//to be removed
void generateGradientImage(image *image) {
    int width = image->w, height = image->h;


    for(int i = 0; i < height; i++){
        double normI = ((double)i)/((double)height);
        for(int j = 0; j < width; j++){
            double normJ = ((double)j)/((double)width);
            image->r(j,i) = (255*normI);
            if(normJ >= normI) {
                image->g(j,i) = (255 * normI / normJ);
            } else {
                image->g(j,i) = (255 * normJ / normI);
            }
            image->b(j,i) = (255*normJ);
        }
    }
}
