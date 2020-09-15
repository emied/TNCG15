#include <iostream>
#include "../glm/glm/glm.hpp"
#include <vector>
#include <fstream>

using namespace std;

struct ColorDbl{
    double r, g, b;
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
/*
struct Vertex {
    double x, y, z, w;
};

struct Direction {
    double x, y, z;
};

struct Triangle {
    Vertex a, b, c; //names?
    ColorDbl color;
    Direction normal;
};

struct Scene {
    Triangle triangles[24];
    ColorDbl colors[6];
};
struct Ray {
    Vertex start, end;
    ColorDbl color;
};

struct Pixel{
    ColorDbl r, g, b;
    Ray intersectingRay;        //Prepare for array/list/pointer
};

struct Camera{
    Vertex leftEye, rightEye;
    Pixel image[800][800];

};
*/


int main() {
    cout << "This is the start of our RayTracer!" << endl;
    const int xWidth = 1920;
    const int yWidth = 1080;
    image img{xWidth,yWidth};
    generateGradientImage(&img);

    cout << "Generated an image!" << endl;
    ofstream(".\\generated.bmp", ios_base::out | ios_base::binary) << img;
    cout << "Wrote file generated.bmp" << endl;
    return 0;
}

void generateGradientImage(image *image) {
    int width = image->w, height = image->h;


    for(int i = 0; i < height; i++){
        double normI = ((double)i)/((double)height);
        for(int j = 0; j < width; j++){
            double normJ = ((double)j)/((double)width);
            image->r(j,i) = (255*normI);
            image->g(j,i) = (255*normJ);
            image->b(j,i) = (255*((double)(i+j)/(double)(width+height)));
        }
    }
}



