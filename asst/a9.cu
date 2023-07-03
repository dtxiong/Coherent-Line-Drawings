#include <iostream>
#include "a9.cuh"


using namespace std;

// Write your implementations here, or extend the Makefile if you add source
// files

// Implements http://www.umsl.edu/cmpsci/faculty-sites/kang/publications/2007/npar07/kang_npar07_hi.pdf


Image computeETFCUDA(const Image &im, const Image &tcur, float radius, float eta) {
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

Image ETFCUDA(const Image &im, float radius, int n, float eta) {
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
        tnew = computeETFCUDA(im, t0, radius, eta);
        t0 = tnew;
    }
    //tnew.debug_write();
    Image output = (tnew + 1)/2;
    return output;
}

__global__ void lineKernel(float* lumiData, float* outputData, float* etfData, int width, int height, 
                           float sigmam, float sigmac, float tau, 
                           float* gaussC, float* gaussS, float* gaussM, float* f) {
    
    int x = (blockIdx.x) * blockDim.x + threadIdx.x;
    int y = (blockIdx.y) * blockDim.y + threadIdx.y;
    int index = x + y * width;
    
    float deltam = 1.0f;
    float deltan = 1.0f;
    float sigmas = 1.6 * sigmac;
    const int pm = ceil(sigmam * 3);
    const int qm = ceil(sigmas * 3);
    const int p = 2 * pm + 1; 
    const int q = 2 * qm + 1;

    float* streamline = (float*)malloc(p * 2 * sizeof(float));
    float* ls = (float*) malloc(q * 2 * sizeof(float));
    //BANDAID FIX. ls and streamline sometimes change to nullptr, and have address 0x0
    while (ls == nullptr) {
        free(ls);
        ls = (float*)malloc(q * 2 * sizeof(float));
    }
    while (streamline == nullptr) {
        free(streamline);
        streamline = (float*)malloc(p * 2 * sizeof(float));
    }
    if (x >= 0 && x < width && y >= 0 && y < height) {
        //Calculate Streamline
        streamline[2 * pm] = (float)x;
        streamline[2 * pm + 1] = (float)y;
        // Calculate cx values
        for (int i = 1; i <= pm; i++)
        {
            float tz[2]; //tangent at z
            float z[2];

            //forward
            z[0] = streamline[2 * (pm + i - 1) + 0]; //previous coord
            z[1] = streamline[2 * (pm + i - 1) + 1];
            tz[0] = interpolateLinCUDA(etfData, z[0], z[1], 0, width, height, true);
            tz[1] = interpolateLinCUDA(etfData, z[0], z[1], 1, width, height, true);
            streamline[2 * (pm + i) + 0] = z[0] + deltam * tz[0];
            streamline[2 * (pm + i) + 1] = z[1] + deltam * tz[1];

            //backward
            z[0] = streamline[2 * (pm - i + 1) + 0];
            z[1] = streamline[2 * (pm - i + 1) + 1];
            tz[0] = interpolateLinCUDA(etfData, z[0], z[1], 0, width, height, true);
            tz[1] = interpolateLinCUDA(etfData, z[0], z[1], 1, width, height, true);
            streamline[2 * (pm - i) + 0] = z[0] - deltam * tz[0]; 
            streamline[2 * (pm - i) + 1] = z[1] - deltam * tz[1];
        };
        // For cx,
        float H = 0.0f;
        for (int i = 0; i < p; i++)
        {

            printf("");
            printf("");
            // Calculate q values
            ls[2 * qm] = streamline[2 * i];
            ls[2 * qm + 1] = streamline[2 * i + 1];
            float gx = etfData[index + width * height]; //grad direction
            float gy = -etfData[index];
            for (int j = 1; j <= qm; j++)
            {
                ls[2 * (qm + j) + 0] = streamline[2 * i + 0] + j * deltan * gx;
                ls[2 * (qm + j) + 1] = streamline[2 * i + 1] + j * deltan * gy;
                ls[2 * (qm - j) + 0] = streamline[2 * i + 0] - j * deltan * gx;
                ls[2 * (qm - j) + 1] = streamline[2 * i + 1] - j * deltan * gy;
            }
            float Fs = 0.0f;
            for (int j = 0; j < q; j++)
            {
                //Integral, Equation 6
                float a = interpolateLinCUDA(lumiData, ls[2 * j + 0], ls[2 * j + 1], 0, width, height, true);
                float b = f[j];
                Fs += interpolateLinCUDA(lumiData, ls[2 * j + 0], ls[2 * j + 1], 0, width, height, true) * f[j] * deltan;
            }
            //cout << Fs << endl;
            //Integral, Equation 9
            H += gaussM[i] * Fs * deltam;
        }
     
        //Calculate Integral
        // printf("x, y, H: %d, %d, %f \n", x, y, H);
        outputData[index] = (float)!((H < 0) && (1 + tanh(H) < tau));

    
    }
    free(ls);
    free(streamline);
}

__device__ float interpolateLinCUDA(float* im, float x, float y, int z, int width, int height, bool clamp) {

   // get the neighboring points
    int xf = floor(x); // floor
    int yf = floor(y);
    int xc = xf + 1; //and ceil
    int yc = yf + 1;

    // compute the distances of the point to the floor-extreme point
    float yalpha = y - yf;
    float xalpha = x - xf;

    float tl = smartAccessorCUDA(im, xf, yf, z, width, height, clamp); // top-left
    float tr = smartAccessorCUDA(im, xc, yf, z, width, height, clamp); // ...
    float bl = smartAccessorCUDA(im, xf, yc, z, width, height, clamp);
    float br = smartAccessorCUDA(im, xc, yc, z, width, height, clamp);

    float topL = tr * xalpha + tl * (1.0f - xalpha);
    float botL = br * xalpha + bl * (1.0f - xalpha);

    float retv = botL * yalpha + topL * (1.0f - yalpha);

    return retv;
}

__device__ float smartAccessorCUDA(float* im, int x, int y, int z, int width, int height, bool clamp) {
    // // --------- HANDOUT  PS02 ------------------------------
    // return 0.0f; // change this

    // --------- SOLUTION PS02 ------------------------------
    float black = 0.0f;
    int x0 = x;
    int y0 = y;

    if (y >= height) {
        if (clamp) {
            y0 = height - 1;
        }
        else {
            return black;
        }
    }

    if (y < 0) {
        if (clamp) {
            y0 = 0;
        }
        else {
            return black;
        }
    }

    if (x >= width) {
        if (clamp) {
            x0 = width - 1;
        }
        else {
            return black;
        }
    }

    if (x < 0) {
        if (clamp) {
            x0 = 0;
        }
        else {
            return black;
        }
    }

    return im[x0 + width * y0 + width * height * z];
}

Image lineConstructionCUDA(const Image &im, float sigmam, float sigmac,float rho, float tau, float radius, int n, float eta) {
    //Draws the lines given t from ETF
    int width = im.width();
    int height = im.height();
    
    Image lumi = lumiChromi(im)[0];
    cout << "Starting ETF" << endl;
    Image etf = ETFCUDA(im, radius, n, eta);
    //etf = Image("Input/circleETF.png");
    etf = etf * 2 - 1;
    Image output = Image(width, height, 1);
    cout << "Starting Line Construction" << endl;
    float deltam = 1.0f;
    float deltan = 1.0f;
    float sigmas = 1.6 * sigmac;
    int pm = ceil(sigmam * 3);
    int qm = ceil(sigmas * 3);
    int p = 2 * pm + 1;
    int q = 2 * qm + 1;

    // Precompute gauss, f values
    float* gaussC = (float*) malloc(q * sizeof(float)); 
    memcpy(gaussC, calculateGaussValuesCUDA(sigmac, q).data(), q * sizeof(float));
    float* gaussS = (float*)malloc(q * sizeof(float));
    memcpy(gaussS, calculateGaussValuesCUDA(sigmas, q).data(), q * sizeof(float));
    float* gaussM = (float*)malloc(p * sizeof(float));
    memcpy(gaussM, calculateGaussValuesCUDA(sigmam, p).data(), p * sizeof(float));
    vector<float> f(q, 0.0f);
    for (int i = 0; i < q; i++)
    {
        f[i] = (gaussC[i] - rho * gaussS[i]); //DoG, Equation 7

    }

    float* fdata = (float*) malloc(q * sizeof(float));
    memcpy(fdata, f.data(), q * sizeof(float));

    cudaError_t cudaStatus;
    cudaDeviceReset();

    float* deviceLumiImageData;
    float* deviceETFImageData;
    float* deviceOutputImageData;
    float* hostLumiImageData = (float*)malloc(width * height * sizeof(float));
    float* hostETFImageData = (float*)malloc(width * height * 3 * sizeof(float));
    float* hostOutputImageData = (float*)malloc(width * height * sizeof(float));

    for (size_t i = 0; i < lumi.number_of_elements(); i++)
    {
        hostLumiImageData[i] = lumi(i);
    }
    for (size_t i = 0; i < etf.number_of_elements(); i++)
    {
        hostETFImageData[i] = etf(i);
    }

    cudaMalloc((void**) &deviceLumiImageData, width * height * sizeof(float));
    cudaMalloc((void**) &deviceETFImageData, width * height * 3 * sizeof(float));
    cudaMalloc((void**) &deviceOutputImageData, width * height * sizeof(float));
    cudaMemcpy(deviceLumiImageData, hostLumiImageData,
        width * height  * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(deviceETFImageData, hostETFImageData,
        width * height * 3 * sizeof(float),
        cudaMemcpyHostToDevice);

    float* deviceGaussC;
    float* deviceGaussS;
    float* deviceGaussM;
    float* deviceF;

    cudaMalloc((void**) &deviceGaussC, q * sizeof(float)); 
    cudaMalloc((void**) &deviceGaussS, q * sizeof(float));
    cudaMalloc((void**) &deviceGaussM, p * sizeof(float));
    cudaMalloc((void**) &deviceF, q * sizeof(float));

    cudaMemcpy(deviceGaussC, gaussC, q * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceGaussS, gaussS, q * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceGaussM, gaussM, p * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(deviceF, fdata, q * sizeof(float), cudaMemcpyHostToDevice);

    //between 1504 and 1520
    int renderSize = 16;
    dim3 dimGrid(renderSize, renderSize);
    dim3 dimBlock(ceil((float)width / renderSize), ceil((float)height / renderSize));
    cout << "Before kernel" << endl;
    lineKernel<<<dimBlock, dimGrid>>>(deviceLumiImageData, deviceOutputImageData, deviceETFImageData, 
        width, height, sigmam, sigmac, tau, 
        deviceGaussC, deviceGaussS, deviceGaussM, deviceF);
    cout << "After kernel" << endl;
    
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    cout << "afterCudaLastError" << endl;

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    cout << "after synchronize" << endl;

    cudaMemcpy(hostOutputImageData, deviceOutputImageData, 
        width * height * sizeof(float), cudaMemcpyDeviceToHost);

    cout << "after copy" << endl;
    for (size_t i = 0; i < output.number_of_elements(); i++)
    {
        output(i) = hostOutputImageData[i];

    }
    cout << "before median Filter" << endl;
    output = medianFilter(output, 1);
    cout << "after median Filter" << endl;


Error:
    cudaMemset(deviceLumiImageData, 0, width * height *
        sizeof(float));
    cudaMemset(deviceETFImageData, 0, width * height *
        3 * sizeof(float));
    cudaMemset(deviceOutputImageData, 0, width * height *
        sizeof(float));

    cudaFree(deviceLumiImageData);
    cudaFree(deviceETFImageData);
    cudaFree(deviceOutputImageData);
    cudaFree(deviceGaussC);
    cudaFree(deviceGaussS);
    cudaFree(deviceGaussM);
    cudaFree(deviceF);

    free(hostLumiImageData);
    free(hostETFImageData);
    free(hostOutputImageData);
    free(gaussC);
    free(gaussS);
    free(gaussM);
    free(fdata);

    return output;
}

vector<vector<float>> getStreamlineCUDA(Image tangent, float x, float y, int n, float stepsize) {
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
        tz = {interpolateLin(tangent, z[0], z[1], 0), interpolateLin(tangent, z[0], z[1], 1)}; 
        streamline[nm + i] = {z[0] + stepsize * tz[0], z[1] + stepsize * tz[1]};
        //backward
        z = streamline[nm - i + 1];
        tz = {interpolateLin(tangent, z[0], z[1], 0), interpolateLin(tangent, z[0], z[1], 1)};
        streamline[nm - i] = {z[0] - stepsize * tz[0], z[1] - stepsize * tz[1]};
    }
    return streamline;
}

vector<float> calculateGaussValuesCUDA(float sigma, float n) {
    //assume n is odd
    int offset = (n - 1) / 2;
    int filterSize = n;
    vector <float> fData(filterSize, 0.0f);
    for (int i = 0; i < filterSize; i++) {
        fData[i] = (1.0f / sigma / sqrt(2 * M_PI)) * exp(-(i - offset) * (i - offset) / (2.0f * sigma * sigma));
    }
    return fData;
}

Image medianFilter(const Image &im, int radius) {
    //For binary images only
    Image output = Image(im.width(), im.height(), 1);
    for (int x = 0; x < im.width(); x++)
    {
        for (int y = 0; y < im.height(); y++)
        {
            int numWhite = 0;
            for (int i = -radius; i <= radius; i++)
            {
                for (int j = -radius; j <= radius; j++)
                {
                    numWhite += (im.smartAccessor(x + i, y + j, 0, true) > 0.5f);
                    
                }
            }
            output(x, y, 0) = (float)(numWhite >= 5);
        }
    }
    return output;
}