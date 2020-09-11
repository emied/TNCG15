#include <iostream>

using namespace std;

struct ColorDbl{
    double r, g, b;
};

struct Vertex {
    double x, y, z, w;
};

struct Direction {
    double x, y, z;
};

struct Triangle {           //Convert to class
    Vertex a, b, c; //names?
    ColorDbl color;
    Direction normal;
};

struct Scene {          //Convert to class
    Triangle triangles[24];
    ColorDbl colors[6];
};
struct Ray {
    Vertex start, end;
    ColorDbl color;
};

struct Pixel{       //Possibly class
    ColorDbl r, g, b;
    Ray intersectingRay;        //Prepare for array/list/pointer
};

struct Camera{      //Convert to class
    Vertex leftEye, rightEye;
    Pixel image[800][800];

};


int main() {
    cout << "This is the start of our RayTracer!" << endl;
    const int xWidth = 1920;
    const int yWidth = 1080;
    string filename = "Test.bmp";
    ColorDbl* image = new ColorDbl[yWidth*xWidth];
    for(int i = 0; i < yWidth; i++){
        for(int j = 0; j < xWidth; j++){
            image[i*xWidth+j] = ColorDbl{(double)(255*i/yWidth), (double)(255*j/xWidth), (double)(255*(i+j)/(xWidth+yWidth))};
        }
    }
    cout << "Generated an image!" << endl;
    return 0;
}