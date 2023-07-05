#ifndef A9_H
#define A9_H

#include "Image.h"
#include "basicImageManipulation.h"
#include "filtering.h"
#include <iostream>
#include <math.h>

using namespace std;

// Write your declarations here, or extend the Makefile if you add source
// files
Image computeETF(const Image &im, const Image &tcur, float radius = 3.0f, float eta = 1.0f);

Image ETF(const Image &im, float radius = 3.0f, int n = 1, float eta = 1.0f);

Image lineConstruction(const Image &im, float sigmam = 3.0f, float sigmac = 2.0f,float rho = 0.99f, float tau = 1.0f, float radius = 3.0f, int n = 1, float eta = 1.0f);

vector<vector<float>> getStreamline(Image tangent, float x, float y, int n, float stepsize, bool clamp);

vector<float> calculateGaussValues(float sigma, float n);

Image medianFilter(const Image& im, int radius = 1);

#endif /* end of include guard: A10_H_PHUDVTKB */

