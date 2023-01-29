#include <iostream>

#include "a9.h"

using namespace std;

// Write your implementations here, or extend the Makefile if you add source
// files

// Implements http://www.umsl.edu/cmpsci/faculty-sites/kang/publications/2007/npar07/kang_npar07_hi.pdf


Image computeETF(const Image &im, const Image &tcur, float radius, float eta) {
    // Computes tnew from tcur
    // Uses the red channel for x direction, green channel for y direction. 
    Image lumi = lumiChromi(im)[0];
    Image gradMag = gradientMagnitude(lumi);
    gradMag = gradMag / gradMag.max(); //normalized
    Image tnew = Image(im.width(), im.height(), 3);
    for (int x = 0; x < im.width(); x++)
    {
        for (int y = 0; y < im.height(); y++)
        {
            float tnewsumx = 0.0f;
            float tnewsumy = 0.0f;
            for (int i = -ceil(radius); i <= ceil(radius); i++)
            {
                for (int j = -ceil(radius); j <= ceil(radius); j++)
                {
                    // Equations 1-5
                    float ws = (float) i*i + j*j < radius*radius;
                    float wm = 0.5f * (1.0f + tanh(eta * (gradMag.smartAccessor(x + i, y + j, 0) - gradMag(x, y, 0)))); //normalized grad?
                    float tdot = tcur(x, y, 0) * tcur.smartAccessor(x+i, y+j, 0) + tcur(x, y, 1) * tcur.smartAccessor(x+i, y+j, 1);
                    float wd = abs(tdot);
                    int sign = tdot > 0 ? 1 : -1;
                    tnewsumx += tcur(x, y, 0) * ws * wm * wd * sign;
                    tnewsumy += tcur(x, y, 1) * ws * wm * wd * sign;
                }
                
            }
            float k = sqrt(tnewsumx*tnewsumx + tnewsumy*tnewsumy);
            if (k != 0) {        
                tnew(x, y, 0) = tnewsumx/k;
                tnew(x, y, 1) = tnewsumy/k;
            }
        }
        
    }
    
    return tnew;
}

Image ETF(const Image &im, float radius, int n, float eta) {
    //First compute t0
    Image lumi = lumiChromi(im)[0];
    Image gradX = gradientX(lumi);
    Image gradY = gradientY(lumi);
    Image t0 = Image(im.width(), im.height(), 3);

    //Calculate t0
    for (int i = 0; i < im.stride(2); i++)
    {
        float gradx = gradX(i);
        float grady = gradY(i);
        float mag = sqrt(gradx * gradx + grady * grady); 
        //perp to grad, counterclockwise
        if (mag != 0) {
            t0(i) = -grady/mag; 
            t0(i + im.stride(2)) = gradx/mag;
        }
    }
    Image tnew = t0;
    for (size_t i = 0; i < n; i++)
    {
        tnew = computeETF(im, t0, radius);
        t0 = tnew;
    }
    //tnew = computeETF(im, t0, radius);
    return tnew;
}

Image lineConstruction(const Image &im, float sigmam, float sigmac,float rho, float tau, float radius, int n, float eta) {
    //Draws the lines given t from ETF
    Image lumi = lumiChromi(im)[0];
    Image etf = ETF(im, radius, n, eta);
    Image output = Image(im.width(), im.height(), 1);

    float deltam = 1.0f;
    float deltan = 1.0f;
    int p = 15;
    int q = 15;
    int pm = (p-1)/2;
    int qm = (q-1)/2;
    float sigmas = 1.6 * sigmac;

    // Precompute gauss, f values
    vector<float> gaussC = calculateGaussValues(sigmac, q);
    vector<float> gaussS = calculateGaussValues(sigmas, q);
    vector<float> gaussM = calculateGaussValues(sigmam, p);
    vector<float> f(q, 0.0f);
    for (int i = 0; i < q; i++)
    {
        f[i] = (gaussC[i] - rho * gaussS[i]); //DoG, Equation 7
        //cout << f[i] << endl;
    }
    
    for (int x = 0; x < im.width(); x++)
    {
        for (int y = 0; y < im.height(); y++)
        {
            //Find region around x, y according to ETF
            float tx = etf(x, y, 0); //tangent direction
            float ty = etf(x, y, 1);
            vector<vector<float>> cx(p);
            cx[pm] = {(float)x, (float)y}; 
            //cout << x << " " << y << " " << tx << " " << ty<< endl;
            // Calculate cx values
            for (int i = 1; i <= pm; i++)
            {
                vector<float> z = cx[pm + i - 1];
                vector<float> tz = {interpolateLin(etf, z[0], z[1], 0), interpolateLin(etf, z[0], z[1], 1)};
                cx[pm + i] = {z[0] + deltam * tz[0], z[1] + deltam * tz[1]};
                z = cx[pm - i + 1];
                tz = {interpolateLin(etf, z[0], z[1], 0), interpolateLin(etf, z[0], z[1], 1)};
                cx[pm - i] = {z[0] - deltam * tz[0], z[1] - deltam * tz[1]};
            }
            // For cx,
            float H = 0.0f;
            for (int i = 0; i < p; i++)
            {
                // Calculate q values
                //cout << x << " " << y << " " << i << " " << cx[i][0] << " " << cx[i][1] << endl;
                vector<vector<float>> ls(q);
                ls[qm] = cx[i]; 
                for (int j = 1; j <= qm; j++)
                {
                    vector<float> z = ls[qm + j - 1];
                    vector<float> gz = {interpolateLin(etf, z[0], z[1], 1), interpolateLin(etf, z[0], z[1], 0)}; //grad, perp to t
                    ls[qm + j] = {z[0] + deltan * gz[0], z[1] + deltan * gz[1]};
                    z = ls[qm - j + 1];
                    gz = {interpolateLin(etf, z[0], z[1], 1), interpolateLin(etf, z[0], z[1], 0)};
                    ls[qm - j] = {z[0] - deltan * gz[0], z[1] - deltan * gz[1]};
                }
                float Fs = 0.0f;
                for (int j = 0; j < q; j++)
                {  
                    //Integral, Equation 6
                    //cout << j << " " << ls[j][0] << " " << ls[j][1] << " " << f[j] << endl;
                    Fs += interpolateLin(lumi, ls[j][0], ls[j][1], 0) * f[j] * deltan;
                }
                //cout << Fs << endl;
                //Integral, Equation 9
                H += gaussM[i] * Fs * deltam;
            }
            // cout << x << " " << y << " " << H << endl;
            output(x, y) = (float) !((H < 0) && (1 + tanh(H) < tau));
        }
    }
    return output;
}

vector<float> calculateGaussValues(float sigma, float n) {
        //assume n is odd
        int offset     = (n-1)/2;
        int filterSize = n;
        vector <float> fData (filterSize, 0.0f);
        for( int i = 0; i < filterSize; i++) {
            fData[i] = (1.0f/ sigma / sqrt(2 * M_PI)) * exp( -(i - offset)*(i - offset) / (2.0f *sigma*sigma) );
        }
        return fData;
    }


Image testGauss(float sigma, float n) {
    vector<float> values = calculateGaussValues(sigma, n);
    Image output = Image(n, 10, 1);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < 10; j++)
        {
            output(i, j) = values[i];
        }
        
    }
    output = output / output.max();
    return output;
    
}