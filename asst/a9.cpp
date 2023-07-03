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
                    float ws = (float) (i*i + j*j ) <= radius*radius;
                    float wm = 0.5f * (1.0f + tanh(eta * (gradMag.smartAccessor(x + i, y + j, 0, true) - gradMag(x, y, 0)))); //normalized grad?
                    float tdot = tcur(x, y, 0) * tcur.smartAccessor(x+i, y+j, 0, true) + tcur(x, y, 1) * tcur.smartAccessor(x+i, y+j, 1, true);
                    float wd = abs(tdot);
                    int sign = tdot > 0 ? 1 : -1;
                    tnewsumx += tcur.smartAccessor(x+i, y+j, 0, true) * ws * wm * wd * sign;
                    tnewsumy += tcur.smartAccessor(x+i, y+j, 1, true) * ws * wm * wd * sign;
                }
                
            }
            float k = sqrt(tnewsumx*tnewsumx + tnewsumy*tnewsumy);
            if (k > 0.00001) {        
                tnew(x, y, 0) = tnewsumx/k;
                tnew(x, y, 1) = tnewsumy/k;
            } else {
                tnew(x, y, 0) = 0;
                tnew(x, y, 1) = 0;
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
        //cout << gradx << " " << grady << endl;
        float mag = sqrt(gradx * gradx + grady * grady); 
        //perp to grad, counterclockwise
        if (mag >= 0.00001) {
            t0(i) = -grady/mag; 
            t0(i + im.stride(2)) = gradx/mag;
        }
        else {
            t0(i) = 0; 
            t0(i + im.stride(2)) = 0;
        }
    }
    Image tnew = t0;
    //tnew.debug_write();
    for (int i = 0; i < n; i++)
    {
        tnew = computeETF(im, t0, radius, eta);
        t0 = tnew;
    }
    //tnew.debug_write();
    Image output = (tnew + 1)/2;
    return output;
}

Image lineConstruction(const Image &im, float sigmam, float sigmac,float rho, float tau, float radius, int n, float eta) {
    //Draws the lines given t from ETF
    Image lumi = lumiChromi(im)[0];
    cout << "Starting ETF" << endl;
    Image etf = ETF(im, radius, n, eta);
    //etf = Image("Input/circleETF.png");
    etf = etf * 2 - 1;
    Image output = Image(im.width(), im.height(), 1);
    cout << "Starting Line Construction" << endl;
    float deltam = 1.0f;
    float deltan = 1.0f;
    float sigmas = 1.6 * sigmac;
    int pm = ceil(sigmam * 3);
    int qm = ceil(sigmas * 3);
    int p = 2 * pm + 1;
    int q = 2 * qm + 1;

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
            //get streamline
            //Find region around x, y according to ETF
            vector<vector<float>> cx(p);
            cx = getStreamline(etf, x, y, p, deltam, true);
            // For cx,
            float H = 0.0f;
            for (int i = 0; i < p; i++)
            {
                // Calculate q values
                vector<vector<float>> ls(q);
                ls[qm] = cx[i]; 
                float gx = etf(x, y, 1); //grad direction
                float gy = -etf(x, y, 0);
                for (int j = 1; j <= qm; j++)
                {
                    ls[qm + j] = {cx[i][0] + j * deltan * gx, cx[i][1] + j * deltan * gy};
                    ls[qm - j] = {cx[i][0] - j * deltan * gx, cx[i][1] - j * deltan * gy};
                }
                float Fs = 0.0f;
                for (int j = 0; j < q; j++)
                {  
                    //Integral, Equation 6
                    Fs += interpolateLin(lumi, ls[j][0], ls[j][1], 0, true) * f[j] * deltan;
                }
                //cout << Fs << endl;
                //Integral, Equation 9
                H += gaussM[i] * Fs * deltam;
            }
            //cout << x << y << endl;
            output(x, y) = (float) !((H < 0) && (1 + tanh(H) < tau));
        }
    }
    return output;
}

vector<vector<float>> getStreamline(Image tangent, float x, float y, int n, float stepsize, bool clamp) {
    vector<vector<float>> streamline(n);
    int nm = (n-1)/2;
    streamline[nm] = {x, y}; 
    // Calculate cx values
    for (int i = 1; i <= nm; i++)
    {
        vector<float> z;
        vector<float> tz; //tangent at z
        //forward
        z = streamline[nm + i - 1]; //previous coord
        tz = {interpolateLin(tangent, z[0], z[1], 0, clamp), interpolateLin(tangent, z[0], z[1], 1, clamp)}; 
        streamline[nm + i] = {z[0] + stepsize * tz[0], z[1] + stepsize * tz[1]};
        //backward
        z = streamline[nm - i + 1];
        tz = {interpolateLin(tangent, z[0], z[1], 0, clamp), interpolateLin(tangent, z[0], z[1], 1, clamp)};
        streamline[nm - i] = {z[0] - stepsize * tz[0], z[1] - stepsize * tz[1]};
    }
    return streamline;
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

