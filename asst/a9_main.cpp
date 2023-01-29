#include <iostream>

#include "a9.h"

using namespace std;

int main()
{
    // Test your intermediate functions
    cout << "Starting testing: " << endl;
    Image flower = Image("Input/flower.png");
    ETF(flower, 0).debug_write();
    ETF(flower, 3, 0).debug_write();
    Image ETFflower = ETF(flower, 3);
    ETFflower.write("Output/ETFflower.png");
    lineConstruction(flower).write("Output/flowerLine.png");

    Image dtxiong = Image("Input/dtxiong.png");
    Image ETFdtxiong = ETF(dtxiong, 3);
    ETFdtxiong.write("Output/ETFdtxiong.png");
    lineConstruction(dtxiong).write("Output/dtxiongLine.png");

    testGauss(5, 15).debug_write();

    return EXIT_SUCCESS;
}
