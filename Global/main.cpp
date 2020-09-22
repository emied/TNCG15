#include <iostream>
#include "../glm/glm/glm.hpp"
#include <vector>
#include <fstream>

using namespace std;
using namespace glm;

struct ColorDbl{

    double r, g, b;
    glm::vec3 rgb;
    ColorDbl(): r(0.0), g(0.0), b(0.0) {} //default constructor

    ColorDbl(glm::vec3 color)
    {
      r = color.x;
      g = color.y;
      b = color.z;
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
    for(int i=0; i<h;i++){
        out.write(data + (3*w * i), 3 * w);
        out.write((char*)&pad, pad);
    }
    return out;
}

struct Vertex {
    glm::vec4 position;
    //double x, y, z, w;
    Vertex():position(glm::vec3(0.0,0.0,0.0), 1.0){}   //Default constructor
    Vertex(double x, double y, double z, double w) : position(glm::vec4(x, y, z, w)){}
    Vertex(double x, double y, double z) : position(glm::vec4(x, y, z, 1.0)){}
    Vertex(glm::vec3 in): position(in, 1.0){}
};

struct Direction {
    glm::vec3 direction;
    Direction(): direction(glm::vec3{}){}
    Direction(Vertex in) : direction(in.position){}

};

struct Triangle;

struct Ray {
    Vertex start, end;
    ColorDbl color;
    Triangle* endTriangle;
    Vertex intersectionPoint;

};

struct Triangle {

    ColorDbl color;
    Vertex vec0, vec1, vec2;
    glm::vec3 normal{};
    double d{};


    Triangle() : Triangle(Vertex(glm::vec3{}),Vertex(glm::vec3{}),Vertex(glm::vec3{})){}

    Triangle(Vertex v1, Vertex v2, Vertex v3){
        vec0 = v1;
        vec1 = v2;
        vec2 = v3;
        normal = glm::cross((glm::vec3)(v2.position-v1.position),(glm::vec3)(v3.position-v1.position));
        d = glm::dot(normal, (glm::vec3)vec0.position);
        color = ColorDbl();
    }


    bool rayIntersection(Ray& intersectingRay) {
        //Möller-Trumbore
        glm::vec3 T, E_1, E_2, D, P, Q;
        T = (glm::vec3) (intersectingRay.start.position - vec0.position);
        E_1 = (glm::vec3) (vec1.position - vec0.position);
        E_2 = (glm::vec3) (vec2.position - vec0.position);
        D = intersectingRay.end.position - intersectingRay.start.position;
        P = glm::cross(D, E_2);
        Q = glm::cross(T, E_1);
        glm::vec3 tuv = glm::vec3(glm::dot(Q, E_2), glm::dot(P, T), glm::dot(Q, D)) / glm::dot(P, E_1);
        if(tuv.y > 0 && tuv.z > 0 && tuv.y+tuv.z <=1.0){
            //Update intersectingRay here
            intersectingRay.endTriangle = this;
            intersectingRay.intersectionPoint = Vertex(tuv);
            intersectingRay.color = this->color;
            return true;
        } else {
            return false;
        }
    }

};

struct Scene;

void createScene(Scene *world);

struct Scene {
    Triangle triangles[24]{};
    //Vertex vertices[14];
    ColorDbl colors[8];
    Scene(){
        createScene(this);
    }


    void rayIntersection(Ray& intersectingRay){
        for(int i = 0; i < 24 ;i++){
            if(triangles[i].rayIntersection(intersectingRay)){break;}
        }
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

void createScene(Scene *world){

    //Vertex points r = roof, f = floor
    Vertex vrtx0 = Vertex(5.0, 0.0, 5.0, 1.0); //vrtx0r mitten toppen
    Vertex vrtx1 = Vertex(-3.0, 0.0, 5.0, 1.0); //vrtx1r
    Vertex vrtx2 = Vertex(0.0, 6.0, 5.0, 1.0); //vrtx2r
    Vertex vrtx3 = Vertex(10.0, 6.0, 5.0, 1.0); //vrtx3r
    Vertex vrtx4 = Vertex(13.0, 0.0, 5.0, 1.0); //vrtx4r
    Vertex vrtx5 = Vertex(10.0, -6.0, 5.0, 1.0); //vrtx5r
    Vertex vrtx6 = Vertex(0.0, -6.0, 5.0, 1.0); //vrtx6r

    Vertex vrtx7 = Vertex(5.0, 0.0, -5.0, 1.0); //vrtx0f mitten botten
    Vertex vrtx8 = Vertex(-3.0, 0.0, -5.0, 1.0); //vrtx1f
    Vertex vrtx9 = Vertex(0.0, 6.0, -5.0, 1.0); //vrtx2f
    Vertex vrtx10 = Vertex(10.0, 6.0, -5.0, 1.0); //vrtx3f
    Vertex vrtx11 = Vertex(13.0, 0.0, -5.0, 1.0); //vrtx4f
    Vertex vrtx12 = Vertex(10.0, -6.0, -5.0, 1.0); //vrtx5f
    Vertex vrtx13 = Vertex(0.0, -6.0, -5.0, 1.0); //vrtx6f

    //Roof = Red
    Triangle tri1= Triangle(vrtx0, vrtx1, vrtx2);
    Triangle tri2= Triangle(vrtx0, vrtx2, vrtx3);
    Triangle tri3= Triangle(vrtx0, vrtx3, vrtx4);
    Triangle tri4= Triangle(vrtx0, vrtx4, vrtx5);
    Triangle tri5= Triangle(vrtx0, vrtx5, vrtx6);
    Triangle tri6= Triangle(vrtx0, vrtx6, vrtx1);

    world->colors[0] = ColorDbl(glm::vec3(255,0,0));
    world->triangles[0] = tri1;
    world->triangles[1] = tri2;
    world->triangles[2] = tri3;
    world->triangles[3] = tri4;
    world->triangles[4] = tri5;
    world->triangles[5] = tri6;
    int j = 0;
    for (int i = j*6; i < (j+1)*6; i++){
        world->triangles[i].color = world->colors[j];
    }
    j++;

    //Floor = Green
    Triangle tri7= Triangle(vrtx7, vrtx9, vrtx8);
    Triangle tri8= Triangle(vrtx7, vrtx10, vrtx9);
    Triangle tri9= Triangle(vrtx7, vrtx11, vrtx10);
    Triangle tri10= Triangle(vrtx7, vrtx12, vrtx11);
    Triangle tri11= Triangle(vrtx7, vrtx13, vrtx12);
    Triangle tri12= Triangle(vrtx7, vrtx8, vrtx13);

    world->colors[1] = ColorDbl(vec3(0,255,0));
    world->triangles[6] = tri7;
    world->triangles[7] = tri8;
    world->triangles[8] = tri9;
    world->triangles[9] = tri10;
    world->triangles[10] = tri11;
    world->triangles[11] = tri12;
    for (int i = j*6; i < (j+1)*6; i++){
        world->triangles[i].color = world->colors[j];
    }
    j++;


    //Walls
    //Wall 1 = Blue
    Triangle tri13= Triangle(vrtx2, vrtx1, vrtx8);
    Triangle tri14= Triangle(vrtx2, vrtx8, vrtx9);
    world->colors[2] = ColorDbl(vec3(0,0,255));
    world->triangles[12] = tri13;
    world->triangles[13] = tri14;
    for (int i = (j+4)*2; i < (j+4+1)*2; i++){
        world->triangles[i].color = world->colors[j];
    }
    j++;
    //Wall 2 = Teal
    Triangle tri15= Triangle(vrtx3, vrtx2, vrtx9);
    Triangle tri16= Triangle(vrtx3, vrtx9, vrtx10);
    world->colors[3] = ColorDbl(vec3(0,255,255));
    world->triangles[14] = tri15;
    world->triangles[15] = tri16;
    for (int i = (j+4)*2; i < (j+4+1)*2; i++){
        world->triangles[i].color = world->colors[j];
    }
    j++;
    //Wall 3 = Yellow
    Triangle tri17= Triangle(vrtx4, vrtx3, vrtx10);
    Triangle tri18= Triangle(vrtx4, vrtx10, vrtx11);
    world->colors[4] = ColorDbl(vec3(255,255,0));
    world->triangles[16] = tri17;
    world->triangles[17] = tri18;
    for (int i = (j+4)*2; i < (j+4+1)*2; i++){
        world->triangles[i].color = world->colors[j];
    }
    j++;
    //Wall 4 = Purple
    Triangle tri19= Triangle(vrtx5, vrtx4, vrtx11);
    Triangle tri20= Triangle(vrtx5, vrtx11, vrtx12);
    world->colors[5] = ColorDbl(vec3(255,0,255));
    world->triangles[18] = tri19;
    world->triangles[19] = tri20;
    for (int i = (j+4)*2; i < (j+4+1)*2; i++){
        world->triangles[i].color = world->colors[j];
    }
    j++;
    //Wall 5 = Black
    Triangle tri21= Triangle(vrtx6, vrtx5, vrtx12);
    Triangle tri22= Triangle(vrtx6, vrtx12, vrtx13);
    world->colors[6] = ColorDbl(vec3(0,0,0));
    world->triangles[20] = tri21;
    world->triangles[21] = tri22;
    for (int i = (j+4)*2; i < (j+4+1)*2; i++){
        world->triangles[i].color = world->colors[j];
    }
    j++;
    //Wall 6 = White
    Triangle tri23= Triangle(vrtx1, vrtx6, vrtx13);
    Triangle tri24= Triangle(vrtx1, vrtx13, vrtx8);
    world->colors[7] = ColorDbl(vec3(255,255,255));
    world->triangles[22] = tri23;
    world->triangles[23] = tri24;
    for (int i = (j+4)*2; i < (j+4+1)*2; i++){
        world->triangles[i].color = world->colors[j];
    }
}


int main() {
    cout << "This is the start of our RayTracer!" << endl;
    const int xWidth = 1920;
    const int yWidth = 1080;
    image img{xWidth,yWidth};
    generateGradientImage(&img);

    cout << "Generated an image!" << endl;
    ofstream("generated.bmp", ios_base::out | ios_base::binary) << img;
    cout << "Wrote file generated.bmp" << endl;

    Scene world;
    const int width = 800;
    const int height = 800;
    const double pixelSize = 2.0/width;
    Camera cam{Vertex(-2,0,0), Vertex(-1,0,0), width, height};
    cam.setPerspective(false);


    createScene(&world);
    double maxIntensity = 0;
    double dy = 0.5;            //Random values, chosen by creator. Guaranteed to be random.
    double dz = 0.5;            //Random values, chosen by creator. Guaranteed to be random.
    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            Ray current{};
            current.start = cam.getEye();
            current.end = Vertex(0,(i-401 + dy)*pixelSize,(j-401 + dz)*pixelSize);
            world.rayIntersection(current);
            cam.image[i*width+j] = Pixel(current.color);
            if(current.color.r > maxIntensity){maxIntensity = current.color.r;}
            if(current.color.g > maxIntensity){maxIntensity = current.color.g;}
            if(current.color.b > maxIntensity){maxIntensity = current.color.b;}
        }
    }

    return 0;
}

void generateGradientImage(image *image) {
    int width = image->w, height = image->h;


    for(int i = 0; i < height; i++){
        double normI = ((double)i)/((double)height);
        for(int j = 0; j < width; j++){
            double normJ = ((double)j)/((double)width);
            image->r(j,i) = (255*normI);
            if(normJ >= normI) {
                image->g(j, i) = (255 * normI / normJ);
            } else {
                image->g(j,i) = (255 * normJ / normI);
            }
            image->b(j,i) = (255*normJ);
        }
    }
}
