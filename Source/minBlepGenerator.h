//
//  wavegen.h
//  Stagecraft Software
//
//  Created by Aaron Leese on 9/22/11.
//  Copyright 2011 Stagecraft Software. All rights reserved.
//

#pragma once
 
#include "JuceHeader.h"  

#if JUCE_INTEL
#define JUCE_SNAP_TO_ZERO(n)    if (! (n < -1.0e-8 || n > 1.0e-8)) n = 0;
#else
#define JUCE_SNAP_TO_ZERO(n)
#endif

class MinBlepGenerator {
    
    // SEE ....
    // http://www.kvraudio.com/forum/viewtopic.php?t=364256
    // http://www.cs.cmu.edu/~eli/papers/icmc01-hardsync.pdf
    // http://stackoverflow.com/questions/175312/bandlimited-waveform-generation
    
    // Basically, we need an oversampled, filtered, nonlinearity .... 1->0 ...
    // This will be added any time the waveform jumps ....
    // in order to to eliminate aliasing (basically, build a bandlimited wave)
    
    // ANTIALIASING FILTER ::::
    // since we are downsampling .... we can filter for better AA
    double coefficients[6];
    struct FilterState
    {
        double x1, x2, y1, y2;
    };
    int numChannels = 2;
    HeapBlock<FilterState> filterStates;
    double ratio, lastRatio;
    
public:
    double overSamplingRatio;
    int zeroCrossings;
    
    float lastValue;
    float lastDelta; // previous derivative ...
    
    // Tweaking the Blep F
    double proportionalBlepFreq;
    bool returnDerivative; // set this to return the FIRST DERIVATIVE of the blep (for first der. discontinuities)
    
    struct BlepOffset
    {
        double offset = 0;
        double freqMultiple = 0;
        double pos_change_magnitude = 0;
        double vel_change_magnitude = 0;
    };
    
    Array<BlepOffset, CriticalSection> currentActiveBlepOffsets;

public:
    
    MinBlepGenerator();
    ~MinBlepGenerator();
    
    Array<float> getMinBlepArray();
    Array<float> getMinBlepDerivArray();
    
    void setToReturnDerivative(bool derivative) {
        
        returnDerivative = derivative;
    }
    
    
    // TEMP 
    void setTest(double newTest);
    
    
    // Utility ....
    
    // SINC Function
    inline double SINC(double x)
    {
        double pix;
        
        if (x == 0.0f)
            return 1.0f;
        else
        {
            pix = double_Pi* x;
            return sin(pix) / pix;
        }
    }
    
    // Generate Blackman Window
    inline void BlackmanWindow(int n, double *w)
    {
        int m = n - 1;
        int i;
        double f1, f2, fm;
        
        fm = (double) m;
        for (i = 0; i <= m; i++) {
            f1 = (2.0f * double_Pi * (double) i) / fm;
            f2 = 2.0f * f1;
            w[i] = 0.42f - (0.5f * cos(f1)) + (0.08f * cos(f2));
        }
    }
    
    // Discrete Fourier Transform
    void DFT(int n, double *realTime, double *imagTime, double *realFreq, double *imagFreq) {
        int k, i;
        double sr, si, p;
        
        for (k = 0; k < n; k++) {
            realFreq[k] = 0.0f;
            imagFreq[k] = 0.0f;
        }
        
        for (k = 0; k < n; k++)
            for (i = 0; i < n; i++) {
                p = (2.0f * double_Pi * (double) (k * i)) / n;
                sr = cos(p);
                si = -sin(p);
                realFreq[k] += (realTime[i] * sr) - (imagTime[i] * si);
                imagFreq[k] += (realTime[i] * si) + (imagTime[i] * sr);
            }
    }
    
    // Inverse Discrete Fourier Transform
    void InverseDFT(int n, double *realTime, double *imagTime, double *realFreq, double *imagFreq) {
        int k, i;
        double sr, si, p;
        
        for (k = 0; k < n; k++) {
            realTime[k] = 0.0f;
            imagTime[k] = 0.0f;
        }
        
        for (k = 0; k < n; k++) {
            for (i = 0; i < n; i++) {
                p = (2.0f * double_Pi * (double) (k * i)) / n;
                sr = cos(p);
                si = -sin(p);
                realTime[k] += (realFreq[i] * sr) + (imagFreq[i] * si);
                imagTime[k] += (realFreq[i] * si) - (imagFreq[i] * sr);
            }
            realTime[k] /= n;
            imagTime[k] /= n;
        }
    }
    
    // Complex Absolute Value
    inline double cabs(double x, double y) {
        return sqrt((x * x) + (y * y));
    }
    
    // Complex Exponential
    inline void cexp(double x, double y, double *zx, double *zy) {
        float expx;
        
        expx = exp(x);
        *zx = expx * cos(y);
        *zy = expx * sin(y);
    }
    
    // Compute Real Cepstrum Of Signal
    void RealCepstrum(int n, double *signal, double *realCepstrum) {
        double *realTime, *imagTime, *realFreq, *imagFreq;
        int i;
        
        realTime = new double[n];
        imagTime = new double[n];
        realFreq = new double[n];
        imagFreq = new double[n];
        
        // Compose Complex FFT Input
        for (i = 0; i < n; i++) {
            realTime[i] = signal[i];
            imagTime[i] = 0.0f;
        }
        
        // Perform DFT
        
        DFT(n, realTime, imagTime, realFreq, imagFreq);
        
        // Calculate Log Of Absolute Value
        for (i = 0; i < n; i++) {
            realFreq[i] = log(cabs(realFreq[i], imagFreq[i]));
            imagFreq[i] = 0.0f;
        }
        
        // Perform Inverse FFT
        InverseDFT(n, realTime, imagTime, realFreq, imagFreq);
        
        // Output Real Part Of FFT
        for (i = 0; i < n; i++)
            realCepstrum[i] = realTime[i];
        
        delete realTime;
        delete imagTime;
        delete realFreq;
        delete imagFreq;
    }
    
    // Compute Minimum Phase Reconstruction Of Signal
    void MinimumPhase(int n, double *realCepstrum, double *minimumPhase) {
        int i, nd2;
        double *realTime, *imagTime, *realFreq, *imagFreq;
        
        nd2 = n / 2;
        realTime = new double[n];
        imagTime = new double[n];
        realFreq = new double[n];
        imagFreq = new double[n];
        
        if ((n % 2) == 1) {
            realTime[0] = realCepstrum[0];
            for (i = 1; i < nd2; i++)
                realTime[i] = 2.0f * realCepstrum[i];
            for (i = nd2; i < n; i++)
                realTime[i] = 0.0f;
        }
        else {
            realTime[0] = realCepstrum[0];
            for (i = 1; i < nd2; i++)
                realTime[i] = 2.0f * realCepstrum[i];
            realTime[nd2] = realCepstrum[nd2];
            for (i = nd2 + 1; i < n; i++)
                realTime[i] = 0.0f;
        }
        
        for (i = 0; i < n; i++)
            imagTime[i] = 0.0f;
        
        DFT(n, realTime, imagTime, realFreq, imagFreq);
        
        for (i = 0; i < n; i++)
            cexp(realFreq[i], imagFreq[i], &realFreq[i], &imagFreq[i]);
        
        InverseDFT(n, realTime, imagTime, realFreq, imagFreq);
        
        for (i = 0; i < n; i++)
            minimumPhase[i] = realTime[i];
        
        delete realTime;
        delete imagTime;
        delete realFreq;
        delete imagFreq;
    }
    
    
    // FILTER ::::::
    void createLowPass (const double frequencyRatio) {
        
        const double proportionalRate = (frequencyRatio > 1.0) ? 0.5 / frequencyRatio
        : 0.5 * frequencyRatio;
        
        const double n = 1.0 / std::tan (double_Pi * jmax (0.001, proportionalRate));
        const double nSquared = n * n;
        const double c1 = 1.0 / (1.0 + std::sqrt (2.0) * n + nSquared);
        
        setFilterCoefficients (c1,
                               c1 * 2.0f,
                               c1,
                               1.0,
                               c1 * 2.0 * (1.0 - nSquared),
                               c1 * (1.0 - std::sqrt (2.0) * n + nSquared));
    }
    void setFilterCoefficients (double c1, double c2, double c3, double c4, double c5, double c6) {
        const double a = 1.0 / c4;
        
        c1 *= a;
        c2 *= a;
        c3 *= a;
        c5 *= a;
        c6 *= a;
        
        coefficients[0] = c1;
        coefficients[1] = c2;
        coefficients[2] = c3;
        coefficients[3] = c4;
        coefficients[4] = c5;
        coefficients[5] = c6;
    }
    void resetFilters() {
        
        filterStates.clear ((size_t) numChannels);
    }
    void applyFilter (float* samples, int num, FilterState& fs) {
        
        while (--num >= 0)
        {
            const double in = *samples;
            
            double out = coefficients[0] * in
            + coefficients[1] * fs.x1
            + coefficients[2] * fs.x2
            - coefficients[4] * fs.y1
            - coefficients[5] * fs.y2;
            
#if JUCE_INTEL
            if (! (out < -1.0e-8 || out > 1.0e-8))
                out = 0;
#endif
            
            fs.x2 = fs.x1;
            fs.x1 = in;
            fs.y2 = fs.y1;
            fs.y1 = out;
            
            *samples++ = (float) out;
        }
    }
    float filterSample (float sample, FilterState& fs) {
        
        
        const double in = sample;
        
        double out = coefficients[0] * in
        + coefficients[1] * fs.x1
        + coefficients[2] * fs.x2
        - coefficients[4] * fs.y1
        - coefficients[5] * fs.y2;
        
#if JUCE_INTEL
        if (! (out < -1.0e-8 || out > 1.0e-8))
            out = 0;
#endif
        
        fs.x2 = fs.x1;
        fs.x1 = in;
        fs.y2 = fs.y1;
        fs.y1 = out;
        
        return (float) out;
        
    }
    
    
    // CUSTOM ::::
    void setLimitingFreq(float proportionOfSamplingRate);
    
    void buildBlep();
    void addBlep(BlepOffset newBlep);
    void addBlepArray(Array<BlepOffset> newBleps);
    
    Array<BlepOffset> getNextBleps();
    
    void processBlock(float* buffer, int numSamples);
    void rescale_bleps_to_buffer(float* buffer, int numSamples, float shiftBlepsBy=0);
    void process_currentBleps(float* buffer, int numSamples);
    

};



