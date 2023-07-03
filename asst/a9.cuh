#ifndef A9_H_CUDA
#define A9_H_CUDA

#include "Image.h"
#include "basicImageManipulation.h"
#include "filtering.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <math.h>

using namespace std;


Image computeETFCUDA(const Image &im, const Image &tcur, float radius = 3.0f, float eta = 1.0f);

Image ETFCUDA(const Image &im, float radius = 3.0f, int n = 1, float eta = 1.0f);

__global__ void lineKernel(float* lumiData, float* outputData, float* etfData, int width, int height,
    float sigmam, float sigmac, float tau, 
    float* gaussC, float* gaussS, float* gaussM, float* f);

__device__ float interpolateLinCUDA(float* im, float x, float y, int z, int width, int height, bool clamp = false);

__device__ float smartAccessorCUDA(float* im, int x, int y, int z, int width, int height, bool clamp = false);

Image lineConstructionCUDA(const Image &im, float sigmam = 3.0f, float sigmac = 2.0f,float rho = 0.99f, float tau = 1.0f, float radius = 3.0f, int n = 1, float eta = 1.0f);

vector<vector<float>> getStreamlineCUDA(Image tangent, float x, float y, int n, float stepsize);

vector<float> calculateGaussValuesCUDA(float sigma, float n);

Image medianFilter(const Image &im, int radius = 1);

#endif /* end of include guard: A10_H_PHUDVTKB */

