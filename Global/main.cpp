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
long double EPSILON = 0.00000001;

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

    ColorDbl& operator*=(float scalar){
        r *= scalar;
        g *= scalar;
        b *= scalar;
        return *this;
    }
};
ColorDbl operator*(ColorDbl const& color, float scalar){
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
    Direction(): direction(glm::vec3{}) {
        direction = glm::vec3(0.0,0.0,0.0);
    }
    Direction(Vertex in) : direction(in.position){}

};

struct Triangle;

struct Ray {
    Ray() : start{}, end{} {}
    Ray(Vertex vertex1, Vertex vertex2) {
        start = vertex1;
        end = vertex2;
        direction.direction = normalize(vertex2 - vertex1);
    }

    Ray(Vertex startP, Direction dir )
    {
        start = startP;
        direction = dir;
    }

    Vertex start, end;
    ColorDbl color;
    Triangle* endTriangle;
    Vertex intersectionPoint{vec3(DBL_MAX,DBL_MAX,DBL_MAX)};
    Direction direction;


};


struct Triangle {

    ColorDbl color;
    Vertex vec0, vec1, vec2;
    vec3 normal{};
    double d{};
    Vertex PointOfIntersection;
    double Tout;


    Triangle() : Triangle(Vertex(vec3{}),Vertex(vec3{}),Vertex(vec3{})){}

    Triangle(Vertex v1, Vertex v2, Vertex v3){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = normalize(cross((vec3)(v2.position-v1.position),(vec3)(v3.position-v1.position)));
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

    bool rayIntersection(Ray& intersectingRay, bool print = false, glm::vec3 point = vec3{}) {

         //Möller-Trumbore
         vec3 T, E_1, E_2, D, P, Q;
         T = (vec3) (intersectingRay.start.position - vec0.position);
         E_1 = (vec3) (vec1.position - vec0.position);
         E_2 = (vec3) (vec2.position - vec0.position);
         D = intersectingRay.end.position - intersectingRay.start.position;
         P = cross(D, E_2);
         Q = cross(T, E_1);
         vec3 tuv = vec3(dot(Q, E_2), dot(P, T), dot(Q, D)) / dot(P, E_1);

        if(tuv.x > 0 && tuv.y > 0 && tuv.z > 0 && tuv.y+tuv.z <=1.0){
            if(intersectingRay.intersectionPoint.position.x > tuv.x) {
                //Update intersectingRay here
                if(print){
                    cout << "vec0 = x: " << vec0.position.x << " y: " << vec0.position.y << " z: " << vec0.position.z << endl;
                }
                intersectingRay.endTriangle = this;
                intersectingRay.intersectionPoint = Vertex(tuv);
                intersectingRay.color = this->color;
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
        triangles = {Triangle(v1, v2, v3, v4), Triangle(v1, v2, v4, v3), Triangle(v1, v3, v4, v2), Triangle(v2, v3, v4, v1)};
    }

    //Tetrahedron(Vertex v1, Vertex v2, Vertex v3, Vertex v4) : Tetrahedron(v1,v2,v3,v4, ColorDbl{}){}

    bool rayIntersection(Ray& intersectingRay, bool print = false) {
        bool collision = false;
        glm::vec3 tmp{}; // tmp är ny här och i ifen 2 rader under
        for(Triangle tri : triangles) {
            if(tri.rayIntersection(intersectingRay, print)){
                intersectingRay.color = color;
                collision = true;

            }
        }
        return collision;
    }
};
/*
struct IntersectionPoint{
    Triangle tri;
    glm::vec3 PointTriangleIntersection;
    float w;
};
*/

struct Sphere{
    ColorDbl color;
    vec3 centerOfSphere;
    double rad;


    Sphere(){}

    Sphere(double radius, vec3 centerSphere, ColorDbl c){
    centerOfSphere = centerSphere;
    rad = radius;
    color = ColorDbl(255,0,0);
    }

    //inspired by
    //https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection
    bool sphereRayIntersection(Ray& ray){
        double t = 100000;
        double r = rad;
        vec3 direction = normalize(vec3(ray.end - ray.start));
        vec3 sphereCenter = centerOfSphere;
        vec3 L = sphereCenter - direction;

        double tca = dot((L), direction);

        if(tca < EPSILON) return false;

        double tcaSquared = (tca*tca);
        double radiusSquared = (r*r);
        float d = dot(L,L) - tcaSquared;

        if(radiusSquared < d) return false;

        double thc = sqrt(radiusSquared-d);

        double t0 = tca-thc;
        double t1 = tca+thc;

        if(t0 > t1) swap(t0,t1);

        if(t0 < EPSILON){
            t0 = t1;
            if (t0 < EPSILON){
                return false;
            }
        }

      //  ray.intersectionPoint = (ray.start.position.x,ray.start.position.y,ray.start.position.z) + direction*f;
        //hit
        t = t0;
        float f = (float) t;
        ray.color = color;
        cout << "color: " << ray.color.r << " t" << t << endl;
        return true;

    }

    vec3 getSphereNormal (Vertex middle){
        return normalize(middle-centerOfSphere);
    }

};




struct Scene;

struct Scene {
    Triangle triangles[24]{};
    Tetrahedron tetras{};
    Sphere spheres{};
    //Vertex vertices[14];
    bool randomRoof;
    Scene() = default;


    void rayIntersection(Ray& intersectingRay, bool print = false){
        random_device rd;
        mt19937 gen(rd());
        uniform_real_distribution<> distrib(0, 1.0);
        double random = distrib(gen);
        for(int i = 0; i < 24 ;i++){
            random = distrib(gen);
            if(triangles[i].rayIntersection(intersectingRay)){break;}
        }
        random = distrib(gen);
        tetras.rayIntersection(intersectingRay, random > 0.9999);
        spheres.sphereRayIntersection(intersectingRay);
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

    ColorDbl colors[8];
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
    colors[0] = ColorDbl(val,val,val);
    world->triangles[0] = Triangle(vrtx0r, vrtx1r, vrtx2r);
    world->triangles[1] = Triangle(vrtx0r, vrtx2r, vrtx3r);
    world->triangles[2] = Triangle(vrtx0r, vrtx3r, vrtx4r);
    world->triangles[3] = Triangle(vrtx0r, vrtx4r, vrtx5r);
    world->triangles[4] = Triangle(vrtx0r, vrtx5r, vrtx6r);
    world->triangles[5] = Triangle(vrtx0r, vrtx6r, vrtx1r);

    //Floor = White
    colors[1] = ColorDbl(255,255,255);
    world->triangles[6] = Triangle(vrtx0f, vrtx2f, vrtx1f);
    world->triangles[7] = Triangle(vrtx0f, vrtx3f, vrtx2f);
    world->triangles[8] = Triangle(vrtx0f, vrtx4f, vrtx3f);
    world->triangles[9] =  Triangle(vrtx0f, vrtx5f, vrtx4f);
    world->triangles[10] = Triangle(vrtx0f, vrtx6f, vrtx5f);
    world->triangles[11] = Triangle(vrtx0f, vrtx1f, vrtx6f);


    //Walls
    //Wall 1 = Red
    colors[2] = ColorDbl(255,0,0);
    world->triangles[12] = Triangle(vrtx2r, vrtx1r, vrtx2f);
    world->triangles[13] = Triangle(vrtx1f, vrtx1r, vrtx2f);

    //Wall 2 = Yellow
    colors[3] = ColorDbl(255,255,0);
    world->triangles[14] = Triangle(vrtx3r, vrtx2r, vrtx3f);
    world->triangles[15] = Triangle(vrtx2f, vrtx3f, vrtx2r);

    //Wall 3 = Green
    colors[4] = ColorDbl(0,255,0);
    world->triangles[16] = Triangle(vrtx4r, vrtx3r, vrtx4f);
    world->triangles[17] = Triangle(vrtx3f, vrtx4f, vrtx3r);

    //Wall 4 = Teal
    colors[5] = ColorDbl(0,255,255);
    world->triangles[18] = Triangle(vrtx5r, vrtx4r, vrtx4f);
    world->triangles[19] = Triangle(vrtx5r, vrtx4f, vrtx5f);

    //Wall 5 = Blue
    colors[6] = ColorDbl(0,0,255);
    world->triangles[20] = Triangle(vrtx6r, vrtx5r, vrtx6f);
    world->triangles[21] = Triangle(vrtx5f, vrtx6f, vrtx5r);

    //Wall 6 = Purple
    colors[7] = ColorDbl(255,0,255);
    world->triangles[22] = Triangle(vrtx1r, vrtx6r, vrtx6f);
    world->triangles[23] = Triangle(vrtx1r, vrtx6f, vrtx1f);

    //Actually assign colours
    for(int j = 0; j < 24; j++){
        if(j < 12) {
            world->triangles[j].color = colors[j/6];
        } else {
            world->triangles[j].color = colors[(int)floor(j/2) - 4];
        }
    }
    //create Tetrahedron
    /*world->tetras = Tetrahedron(vec3(8,0,2),
                                vec3(7,0,-3),
                                vec3(9,3,-1),
                                vec3(9,-2,-1),
                                ColorDbl(0,150,0));*/
    world->tetras = Tetrahedron(vec3(6,0,-3),
                                vec3(8,2,-3),
                                vec3(7.5,-2,-3),
                                vec3(7,0,0),
                                ColorDbl(0,150,0));
    //Create Sphere
    world->spheres = Sphere(1,vec3(6,0,0),ColorDbl(255, 0, 0));

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
    light.position = vec3{3.0,-1.0,1.0};

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
                current.start = cam.getEye();

                current.end = Vertex(0,
                                     (i - 400 + (double)(r%subPixelsPerAxis)/(double)subPixelsPerAxis + dy/subPixelsPerAxis) * pixelSize,
                                     (j - 400 + floor(r/subPixelsPerAxis)/subPixelsPerAxis + dz/subPixelsPerAxis) * pixelSize);
                world.rayIntersection(current);

                Ray shadow{};
                double u , v;
                u = current.intersectionPoint.position.y;
                v = current.intersectionPoint.position.z;

                shadow.start = Vertex{(1-u-v)*(vec3)current.endTriangle->vec0.position + u*current.endTriangle->vec1.position + v*current.endTriangle->vec2.position};
                //Move start out of object
                shadow.start = (vec3)shadow.start.position + 0.1*(current.start.position - current.end.position);
                shadow.end = light.position;


                random = distrib(gen);
                world.rayIntersection(shadow);
                bool shadedRay = false;
                if(shadow.intersectionPoint.position.x <= 1 && shadow.intersectionPoint.position.x >= 0){
                    shadedRay = true;
                }


                ColorDbl shadowOrNot = (shadedRay)?ColorDbl{}:ColorDbl{1.0,1.0,1.0};

                pixelAvg += (current.color) * shadowOrNot;
                //KAOS, bör inte vara 0 hela tiden
                if(random > 0.9999){
                    //cout << "Shadow start pos, x: " << shadow.start.position.x << " y: " << shadow.start.position.y << " z: " << shadow.start.position.z << endl;
                    //cout << "Shaded ray t-value: " << shadow.intersectionPoint.position.x << ", shaded? = "<< shadedRay <<endl;
                    //cout << "shadowOrNot  " << shadowOrNot.r << " " << shadowOrNot.g << " " << shadowOrNot.b << endl;
                    /*cout << "Current triangle color: " << endl << "r: " << current.endTriangle->color.r <<
                    ", g: " << current.endTriangle->color.g << ", b: " << current.endTriangle->color.b << endl;*/
                    /*cout << "Collision triangle  normal: " << endl << "x: " << current.endTriangle->normal.x <<
                    ", y: " << current.endTriangle->normal.y << ", z: " << current.endTriangle->normal.z << endl;*/
                }
                if (current.color.r > maxIntensity) { maxIntensity = current.color.r; }
                if (current.color.g > maxIntensity) { maxIntensity = current.color.g; }
                if (current.color.b > maxIntensity) { maxIntensity = current.color.b; }
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