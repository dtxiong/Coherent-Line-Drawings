#include <iostream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include "a9.cuh"
#include "a9.h"

using namespace std;

int main()
{
    // Test your intermediate functions
    cout << "Starting main(): " << endl;

    
    // Image flower = Image("Input/flower.png");
    // ETF(flower, 0).debug_write();
    // ETF(flower, 3, 3).debug_write();
    // Image ETFflower = ETF(flower, 3);
    // ETFflower.write("Output/ETFflower.png");
    // lineConstructionCUDA(flower, 5.0f, 3.0f, 0.99f, 1.0f, 5.0f, 3, 1.0f).write("Output/flowerLine.png");
    // Image flowerLine = lineConstructionCUDA(flower, 3.0f, 2.0f, 0.99f, 1.0f, 5.0f, 3, 1.0f);
    // flowerLine.write("Output/flower.png");


    // Image circle = Image("Input/circle.png");
    //ETF(circle, 0).debug_write();
    // ETF(circle, 5, 3).debug_write();
    // ETF(circle, 10, 7).debug_write();
    // ETF(circle, 20, 3).debug_write();
    // ETF(circle, 10, 3, 2).debug_write();
    // ETF(circle, 2, 3).debug_write();
    // Image ETFcircle = ETF(circle, 5);
    // ETFcircle.write("Output/ETFcircle.png");
    // lineConstruction(circle, 5.0f, 3.0f, 0.99f, 1.0f, 3.0f, 3, 1.0f).write("Output/circleLine.png");
    // lineConstruction(circle, 5.0f, 3.0f, 0.99f, 1.0f, 10.0f, 3, 1.0f).write("Output/circleLine1.png");

    
    // Image baboon = Image("Input/baboon.png");
    // ETF(baboon, 0).debug_write();
    // ETF(baboon, 3, 3).debug_write();
    // Image ETFbaboon = ETF(baboon, 3);
    // ETFbaboon.write("Output/ETFbaboon.png");
    // lineConstructionCUDA(baboon, 3.0f, 2.0f, 0.99f, 1.0f, 3.0f, 3, 1.0f).write("Output/baboonLine.png");
    
    /*
    Image circleETF = Image(100, 100, 3);
    for (size_t i = 0; i < 100; i++)
    {
        for (size_t j = 0; j < 100; j++)
        {
            vector<float> xy = {i-49.5f, j-49.5f};
            vector<float> tan = {-xy[1], xy[0]};
            float mag = sqrt(xy[0] * xy[0] + xy[1] * xy[1]); 
            tan = {tan[0]/mag, tan[1]/mag};
            circleETF(i, j, 0) = (tan[0] + 1)/2;
            circleETF(i, j, 1) = (tan[1] + 1)/2;
        }
        
    }
    circleETF.write("Input/circleETF.png");
    cout << "test" << circleETF(25, 75, 0) << circleETF(25, 75, 1) << endl;
    Image circleETF1 = Image(100, 100, 3);
    cout << circleETF1(25, 75, 0) << circleETF1(25, 75, 1) << endl;
    */
    
    // Image dtxiong = Image("Input/dtxiong.png");
    // Image ETFdtxiong = ETF(dtxiong, 3);
    // ETFdtxiong.write("Output/ETFdtxiong.png");
    // lineConstructionCUDA(dtxiong).write("Output/dtxiongLine.png");
    
    /*
    auto start = chrono::high_resolution_clock::now();
    Image C0261 = Image("Input/C0261/0030.png");
    Image C0261Line = lineConstructionCUDA(C0261, 3.0f, 2.0f, 0.99f, 1.0f, 3.0f, 3, 1.0f);
    C0261Line.write("Output/C0261.png");
    //lineConstruction(montauk, 3.0f, 2.0f, 0.99f, 1.0f, 3.0f, 3, 1.0f).debug_write();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    cout << "Time taken by frame: "
        << duration.count() << " microseconds" << endl;
    */
    //video sequence
    
    for (int i = 1; i <= 369; i = i+1)
    {
        auto start = chrono::high_resolution_clock::now();

        cout << "frame " << i << endl;

        string directory = "Input/C0261/";
        stringstream ss;
        ss << setw(4) << setfill('0') << i;
        string frame = ss.str();
        string filename = directory + frame + ".png";
        Image frameImage = Image(filename);
        lineConstructionCUDA(frameImage, 3.0f, 2.0f, 0.99f, 1.0f, 3.0f, 3, 1.0f).write("Output/C0261/" + frame + ".png");
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

        cout << "Time taken by frame " << i << ": "
            << duration.count() << " microseconds" << endl;
    }
    
    
    return EXIT_SUCCESS;
}
