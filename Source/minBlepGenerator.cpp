//
//  wavegen.cpp
//  Stagecraft Software
//
//  Created by Aaron Leese on 9/22/11.
//  Copyright 2011 Stagecraft Software. All rights reserved.
//


#include "JuceHeader.h" 
#include "WaveGenerator.h"

// STATIC ARRAYS - to house the minBlep and integral of the minBlep ...
static Array<float> minBlepArray;
static Array<float> minBlepDerivArray;

MinBlepGenerator::MinBlepGenerator() {
    
    overSamplingRatio = 16;
    zeroCrossings = 16;
    returnDerivative = false;
    proportionalBlepFreq = 0.5; // defaults to NyQuist ....
    
    lastValue = 0;
    lastDelta = 0;
    
    //currentOffsetInCircular = 0;
    
    // AA FILTER 
    zeromem (coefficients, sizeof (coefficients));
    
    numChannels = 2;
    filterStates.calloc (numChannels);
    
    ratio = 1;
    lastRatio = 1;
    
    createLowPass (ratio);
    resetFilters();

    
    buildBlep();
    
}
MinBlepGenerator::~MinBlepGenerator() {

    //
    
}

void MinBlepGenerator::setLimitingFreq(float proportionOfSamplingRate) {
    
    //
    // Instead of limiting to the sampling F,
    // We bring the maximum allowable F down to some known quantity
    // Doing this we can "tune" the blep to some desired F
    // So .... making wave-generators better ....
    
    // SINCE the buffer is only resized to 8x, we can only use blep adjustments down to 0.125
    proportionOfSamplingRate = jlimit<float>(0.0001, 1.0f, proportionOfSamplingRate);
    proportionalBlepFreq = proportionOfSamplingRate;
    
}

Array<float> MinBlepGenerator::getMinBlepArray() {
    
    return minBlepArray;
    
}
Array<float> MinBlepGenerator::getMinBlepDerivArray() {
    
    return minBlepDerivArray;
    
}


void MinBlepGenerator::clear() {
    
    jassert(currentActiveBlepOffsets.size() == 0);
    
    currentActiveBlepOffsets.clear();
    
}
bool MinBlepGenerator::isClear() {
    
    return (currentActiveBlepOffsets.size() > 0);
}

// MIN BLEP - freq domain calc
//
//  wavegen.cpp
//  Stagecraft Software
//
//  Created by Aaron Leese on 9/22/11.
//  Copyright 2011 Stagecraft Software. All rights reserved.
//


#include "JuceHeader.h"
#include "WaveGenerator.h"

// STATIC ARRAYS - to house the minBlep and integral of the minBlep ...
static Array<float> minBlepArray;
static Array<float> minBlepDerivArray;

MinBlepGenerator::MinBlepGenerator() {
    
    overSamplingRatio = 16;
    zeroCrossings = 16;
    returnDerivative = false;
    proportionalBlepFreq = 0.5; // defaults to NyQuist ....
    
    lastValue = 0;
    lastDelta = 0;
    
    // AA FILTER
    zeromem (coefficients, sizeof (coefficients));
    
    numChannels = 2;
    filterStates.calloc (numChannels);
    
    ratio = 1;
    lastRatio = 1;
    
    createLowPass (ratio);
    resetFilters();
    
    
    buildBlep();
    
}
MinBlepGenerator::~MinBlepGenerator() {
    
    //
    
}

void MinBlepGenerator::setLimitingFreq(float proportionOfSamplingRate) {
    
    //
    // Instead of limiting to the sampling F,
    // We bring the maximum allowable F down to some known quantity
    // Doing this we can "tune" the blep to some desired F
    // So .... making wave-generators better ....
    
    // SINCE the buffer is only resized to 8x, we can only use blep adjustments down to 0.125
    proportionOfSamplingRate = jlimit<float>(0.0001, 1.0f, proportionOfSamplingRate);
    proportionalBlepFreq = proportionOfSamplingRate;
    
}

Array<float> MinBlepGenerator::getMinBlepArray() {
    
    return minBlepArray;
    
}
Array<float> MinBlepGenerator::getMinBlepDerivArray() {
    
    return minBlepDerivArray;
    
}


void MinBlepGenerator::clear() {
    
    jassert(currentActiveBlepOffsets.size() == 0);
    
    currentActiveBlepOffsets.clear();
    
}
bool MinBlepGenerator::isClear() {
    
    return (currentActiveBlepOffsets.size() > 0);
}


// MIN BLEP - freq domain calc
void MinBlepGenerator::buildBlep() {
    
    // ALREADY built - so return ...
    if (minBlepArray.size() > 0) return;
    
    // BUILD the BLEP
    int i, n;
    double r, a, b;
    Array <double> buffer1;
    Array <double> buffer2;
    
    n = (zeroCrossings * 2 * overSamplingRatio) + 1;
    
    DBG("BUILD minBLEP - ratio " + String(overSamplingRatio) + " -> " + String(n));
    
    
    // Generate Sinc
    const float bandlimit = 0.9f;
    a = bandlimit * (double) -zeroCrossings;
    b = -a;
    for (i = 0; i < n; i++) {
        r = ((double) i) / ((double) (n - 1));
        
        buffer1.add( SINC(a + (r * (b - a))) );
        buffer2.add(0); // size ...
    }
    
    
    jassert(buffer1.size() == n);
    jassert(buffer2.size() == n);
    
    // Window Sinc
    BlackmanWindow(n, buffer2.getRawDataPointer());
    FloatVectorOperations::multiply(buffer1.getRawDataPointer(), buffer2.getRawDataPointer(), n);
    
    // Minimum Phase Reconstruction
    RealCepstrum(n, buffer1.getRawDataPointer(), buffer2.getRawDataPointer());
    MinimumPhase(n, buffer2.getRawDataPointer(), buffer1.getRawDataPointer());
    
    // Integrate Into MinBLEP
    minBlepArray.ensureStorageAllocated(n);
    minBlepDerivArray.ensureStorageAllocated(n);
    
    a = 0;
    double secondInt = 0;
    for (i = 0; i < n; i++) {
        
        a += buffer1[i]; // full integral ... so that we can normalize (make area=1)
        minBlepArray.add(a);
        
        // 2ND ORDER ::::
        secondInt += a;
        minBlepDerivArray.add(secondInt);
        
    }
    
    
    // Normalize
    double maxVal = minBlepArray.getUnchecked(n - 1);
    FloatVectorOperations::multiply(minBlepArray.getRawDataPointer(), 1.0/maxVal, n);
    
    // Normalize ...
    float max = FloatVectorOperations::findMaximum(minBlepDerivArray.getRawDataPointer(), n);
    jassert(max = minBlepDerivArray.getLast());
    FloatVectorOperations::multiply(minBlepDerivArray.getRawDataPointer(), 1.0/max, minBlepDerivArray.size());
    
    for (double ramp = 0; ramp < n; ramp++) {
        
        // 2ND ORDER ::::
        minBlepDerivArray.getRawDataPointer()[int(ramp)] -= double(ramp)/double(n-1);
    }
    
    
    jassert(fabsf(minBlepDerivArray[0]) < 0.01);
    jassert(fabsf(minBlepDerivArray[n-1]) < 0.01);
    
    
    // SUBTRACT 1 and invert so the signal (so it goes 1->0)
    FloatVectorOperations::add(minBlepArray.getRawDataPointer(), -1.f, minBlepArray.size());
    FloatVectorOperations::multiply(minBlepArray.getRawDataPointer(), -1.f, minBlepArray.size());
    
    jassert(fabsf(minBlepArray[n]) < 0.001);
    
}

void MinBlepGenerator::addBlep(BlepOffset newBlep) {
    
    //
    jassert(newBlep.offset <= 0);
    
    newBlep.freqMultiple = overSamplingRatio*proportionalBlepFreq;
    currentActiveBlepOffsets.add(newBlep);
    
}

void MinBlepGenerator::addBlepArray(Array<BlepOffset> newBleps) {
    
    //
    currentActiveBlepOffsets.addArray(newBleps, 0, newBleps.size());
    
}

Array<MinBlepGenerator::BlepOffset> MinBlepGenerator::getNextBleps() {
    
    
    Array<MinBlepGenerator::BlepOffset> newBleps;
    newBleps.addArray(currentActiveBlepOffsets, 0, currentActiveBlepOffsets.size());
    
    jassert(newBleps.size() == currentActiveBlepOffsets.size());
    
    // CLEAR the array ...
    currentActiveBlepOffsets.clearQuick();
    
    return newBleps;
    
}


// REAL TIME ::::: the core functions :::::
void MinBlepGenerator::processBlock(float* buffer, int numSamples) {
    
    
    // look for non-linearities ....
    jassert(numSamples > 0);
    
    // NON-LINEARITIES :::::
    // This is for processing detected nonlinearities about which we ONLY know the POSITION
    // process_nonlinearities(buffer, numSamples, nonlinearities);
    
    // GRAB the final value ....
    // just in case there is a nonlinearity at sample 0 of the next block ...
    // MUST be done BEFORE we ADD the bleps
    lastValue = buffer[numSamples-1];
    
    // Hmmm .... once in a while there is a nonlinearity at the edge ....
    // inwhich case, we probably shouldn't update the delta ...
    // jassert(lastDelta == 0);
    if ( int(currentActiveBlepOffsets.getLast().offset) != -(numSamples-1))
    {
        lastDelta = buffer[numSamples-1] - buffer[numSamples-2];
    }
    else // hacky .... hmmm ...
    {
        lastDelta = buffer[numSamples-2] - buffer[numSamples-3];
    }
    
    
    // PROCESS BLEPS :::::
    process_currentBleps(buffer, numSamples);
    
    
}

void MinBlepGenerator::rescale_bleps_to_buffer(float* buffer, int numSamples, float shiftBlepsBy) {
    
    
    // MUST be big enough to hold the entire wave after all ... and safety factor of 2
    jassert(currentActiveBlepOffsets.size() < 1000);
    
    for (int i = currentActiveBlepOffsets.size(); --i >= 0;)
    {
        
        // confusingly, the nonlinearities actually occued 1 sample BEFORE the number we get
        // this is because the detector detects that one JUST OCCURED, and then adds it with the fractional offset from the last sample
        BlepOffset blep = currentActiveBlepOffsets.getUnchecked(i); // ie 101.2
        
        
        // SCALE FREQ (to this bleps prop freq)
        blep.freqMultiple = overSamplingRatio*proportionalBlepFreq;
        
        // MODIFY :::: the exact offset ...
        // since this is an effect ... it manifests 1 sample later than the discontinuity
        float exactOffset = -blep.offset + shiftBlepsBy; // +1 here is NEEDED for flanger/chorus !!
        blep.offset = blep.offset - shiftBlepsBy; // starts compensating on the sample AFTER the blep ....
        
        // ACTIVE blep .... (not upcoming)
        if (exactOffset < 0) continue;
        
        // CHECK :::: further away than 1 buffer ... should never happen
        if (exactOffset > numSamples)
        {
            // LFOs have nonlinearities that affect the audio 1 sample later
            // ... so we can get edge cases here ....
            // simply roll it over to the next buffer ...
            DBG("OUT OF RANGE NONLINEARITY ??? " + String(exactOffset) );
            blep.offset = exactOffset-numSamples;
            currentActiveBlepOffsets.setUnchecked(i, blep);
            continue;
        }
        
        
        // CALCULATE the MAGNITUDE of the nonlinearity
        float magnitude_position = 0;
        float magnitude_velocity = 0;
        float currentDelta = 0;
        {
            
            // CALCULATE the integer (sample) offset, and the fractional (subsample) offset
            double intOffset = int(exactOffset);
            double fraction = modf(exactOffset, &intOffset);
            
            // UNLESS we're on the edge case, we get the most recent value from the buffer ...
            if (intOffset > 0)
                lastValue = buffer[int(intOffset) - 1];
            
            
            // 1st order (velocity)
            // MUST do this one first - since the 0th order may change the LastValue
            {
                
                // FIND the last and next deltas ... and compute the difference ....
                if (intOffset >= 2 ) lastDelta = buffer[int(intOffset)-1] - buffer[int(intOffset)-2];
                else if (intOffset >= 1 ) lastDelta = buffer[int(intOffset)-1] - lastValue;
                
                // DEFAULT :: assume flat ...
                currentDelta = 0;
                
                if (intOffset + 1 < numSamples)
                    currentDelta = buffer[int(intOffset)+1] - buffer[int(intOffset)];
                
                // CALCULATE change in velocity
                double change_in_delta = (currentDelta - lastDelta);
                double propDepth = proportionalBlepFreq;
                
                magnitude_velocity = -4*change_in_delta*(1/propDepth);
                
                //jassert(magnitude_velocity <= 1);
                
            }
            
            // 0th order (position)
            {
                // CALCULATE the magnitude of the 0 order nonlinearity *change in position*
                float extrapolated_last_pos = lastValue + lastDelta*(fraction);
                float extrapolated_jump_pos = buffer[int(intOffset)] - currentDelta*(1 - fraction);
                magnitude_position = extrapolated_last_pos - extrapolated_jump_pos;
            }
            
            
        }
        
        
        // TOO SMALL :::
        /// no need to compensating for tiny discontinuities
        if (fabsf(magnitude_position) < .001 && fabsf(magnitude_velocity) < .001)
        {
            currentActiveBlepOffsets.remove(i);
            continue;
        }
        
        // NEGLIGIBLE MAGNITUDES :::
        /// zero out any tiny effects here, so we don't waste time calculating them
        if (fabsf(magnitude_position) < .001) magnitude_position = 0;
        if (fabsf(magnitude_velocity) < .001) magnitude_velocity = 0;
        
        
        // ADD ::::
        // GAIN factors ... how big of a discontinutiy are we talking about ?
        blep.pos_change_magnitude = magnitude_position;
        blep.vel_change_magnitude = magnitude_velocity;
        
        // ALTER :::
        currentActiveBlepOffsets.setUnchecked(i, blep);
        
    }
    
    
}

void MinBlepGenerator::process_currentBleps(float* buffer, int numSamples) {
    
    
    
    // PROCESS ALL BLEPS -
    /// for each offset, copy a portion of the blep array to the output ....
    for (int i = currentActiveBlepOffsets.size(); --i >= 0;)
    {
        
        BlepOffset blep = currentActiveBlepOffsets[i];
        double adjusted_Freq = blep.freqMultiple;
        double exactPosition = blep.offset;
        
        
        // ADD the BLEP to the circular buffer ...
        for (float p = 0; p < numSamples; p++)
        {
            
            double blepPosExact = adjusted_Freq*(exactPosition + p + 1); // +1 because this needs to trigger on the LOW SAMPLE
            double blepPosSample = 0;
            double fraction = modf(blepPosExact, &blepPosSample);
            
            // LIMIT the scaling on the derivative array
            // otherwise, it can get TOO large
            double depthLimited = proportionalBlepFreq; //jlimit<double>(.1, 1, proportionalBlepFreq);
            double blepDeriv_PosExact = depthLimited*overSamplingRatio*(exactPosition + p + 1);
            double blepDeriv_Sample = 0;
            double fraction_Deriv = modf(blepDeriv_PosExact, &blepDeriv_Sample);
            
            
            // DONE ... we reached the end ...
            if (int(blepPosExact) > minBlepArray.size() && int(blepDeriv_PosExact) > minBlepArray.size())
                break;
            
            // BLEP has not yet occurred ...
            if (blepPosExact < 0)
                continue;
            
            
            // 0TH ORDER COMPENSATION ::::
            /// add the BLEP to compensate for discontinuties in the POSITION
            if ( fabs(blep.pos_change_magnitude) > 0 && blepPosSample < minBlepArray.size())
            {
                // LINEAR INTERPOLATION ::::
                float lowValue = minBlepArray.getRawDataPointer()[int(blepPosSample)];
                float hiValue = lowValue;
                
                if (int(blepPosSample) + 1 < minBlepArray.size())
                    hiValue = minBlepArray.getRawDataPointer()[int(blepPosSample) + 1];
                
                float delta = hiValue - lowValue;
                float exactValue = lowValue + fraction*delta;
                
                // SCALE by the discontinuity magnitude
                exactValue *= blep.pos_change_magnitude;
                
                // ADD to the thruput
                buffer[int(p)] += exactValue;
            }
            
            
            // 1ST ORDER COMPENSATION ::::
            /// add the BLEP DERIVATIVE to compensate for discontinuties in the VELOCITY
            if ( fabs(blep.vel_change_magnitude) > 0 && blepDeriv_PosExact < minBlepDerivArray.size())
            {
                
                // LINEAR INTERPOLATION ::::
                double lowValue = minBlepDerivArray.getRawDataPointer()[int(blepDeriv_PosExact)];
                double hiValue = lowValue;
                
                if (int(blepDeriv_PosExact) + 1 < minBlepDerivArray.size())
                    hiValue = minBlepDerivArray.getRawDataPointer()[int(blepDeriv_PosExact) + 1];
                
                double delta = hiValue - lowValue;
                double exactValue = lowValue + fraction_Deriv*delta;
                
                // SCALE by the discontinuity magnitude
                exactValue *= blep.vel_change_magnitude;
                
                // ADD to the thruput
                buffer[int(p)] += exactValue;
                
            }
            
        }
        
        // UPDATE ::::
        blep.offset = blep.offset + double(numSamples);
        if (blep.offset*adjusted_Freq > minBlepArray.size())
        {
            currentActiveBlepOffsets.remove(i);
        }
        else currentActiveBlepOffsets.setUnchecked(i, blep);
        
    }
    
    
    
}



void MinBlepGenerator::addBlep(BlepOffset newBlep) {
    
    //
    jassert(newBlep.offset <= 0);
    
    newBlep.freqMultiple = overSamplingRatio*proportionalBlepFreq;
    currentActiveBlepOffsets.add(newBlep);
    
}
void MinBlepGenerator::addBlepArray(Array<BlepOffset> newBleps) {
    
    //
    currentActiveBlepOffsets.addArray(newBleps, 0, newBleps.size());
    
}
Array<MinBlepGenerator::BlepOffset> MinBlepGenerator::getNextBleps() {
    

    Array<MinBlepGenerator::BlepOffset> newBleps;
    newBleps.addArray(currentActiveBlepOffsets, 0, currentActiveBlepOffsets.size());
    
    jassert(newBleps.size() == currentActiveBlepOffsets.size());
    
    // CLEAR the array ...
    currentActiveBlepOffsets.clearQuick();
    
    return newBleps;
    
}


// REAL TIME ::::: the core functions :::::
void MinBlepGenerator::processBlock(float* buffer, int numSamples) {
    
    
    // look for non-linearities ....
    jassert(numSamples > 0);
    
    // NON-LINEARITIES :::::
    // This is for processing detected nonlinearities about which we ONLY know the POSITION
    // process_nonlinearities(buffer, numSamples, nonlinearities);
   
    // GRAB the final value ....
    // just in case there is a nonlinearity at sample 0 of the next block ...
    // MUST be done BEFORE we ADD the bleps
    lastValue = buffer[numSamples-1];
    
    // Hmmm .... once in a while there is a nonlinearity at the edge ....
    // inwhich case, we probably shouldn't update the delta ...
    // jassert(lastDelta == 0);
    if ( int(currentActiveBlepOffsets.getLast().offset) != -(numSamples-1))
    {
        lastDelta = buffer[numSamples-1] - buffer[numSamples-2];
    }
    else // hacky .... hmmm ...
    {
        lastDelta = buffer[numSamples-2] - buffer[numSamples-3];
    }

    
    // PROCESS BLEPS :::::
    process_currentBleps(buffer, numSamples);
    
    
}

void MinBlepGenerator::rescale_bleps_to_buffer(float* buffer, int numSamples, float shiftBlepsBy) {
    
    
    // MUST be big enough to hold the entire wave after all ... and safety factor of 2
    jassert(currentActiveBlepOffsets.size() < 1000);
    
    for (int i = currentActiveBlepOffsets.size(); --i >= 0;)
    {
        
        // confusingly, the nonlinearities actually occued 1 sample BEFORE the number we get
        // this is because the detector detects that one JUST OCCURED, and then adds it with the fractional offset from the last sample
        BlepOffset blep = currentActiveBlepOffsets.getUnchecked(i); // ie 101.2
        
        
        // SCALE FREQ (to this bleps prop freq)
        blep.freqMultiple = overSamplingRatio*proportionalBlepFreq;
        
        // MODIFY :::: the exact offset ...
        // since this is an effect ... it manifests 1 sample later than the discontinuity
        float exactOffset = -blep.offset + shiftBlepsBy; // +1 here is NEEDED for flanger/chorus !!
        blep.offset = blep.offset - shiftBlepsBy; // starts compensating on the sample AFTER the blep ....
        
        // ACTIVE blep .... (not upcoming)
        if (exactOffset < 0) continue;
        
        // CHECK :::: further away than 1 buffer ... should never happen
        if (exactOffset > numSamples)
        {
            // LFOs have nonlinearities that affect the audio 1 sample later
            // ... so we can get edge cases here ....
            // simply roll it over to the next buffer ...
            DBG("OUT OF RANGE NONLINEARITY ??? " + String(exactOffset) );
            blep.offset = exactOffset-numSamples;
            currentActiveBlepOffsets.setUnchecked(i, blep);
            continue; 
        }
        
        
        // CALCULATE the MAGNITUDE of the nonlinearity
        float magnitude_position = 0;
        float magnitude_velocity = 0;
        float currentDelta = 0;
        {
            
            // CALCULATE the integer (sample) offset, and the fractional (subsample) offset
            double intOffset = int(exactOffset);
            double fraction = modf(exactOffset, &intOffset);
            
            // UNLESS we're on the edge case, we get the most recent value from the buffer ...
            if (intOffset > 0)
                lastValue = buffer[int(intOffset) - 1];
            
            
            // 1st order (velocity)
            // MUST do this one first - since the 0th order may change the LastValue
            {
                
                // FIND the last and next deltas ... and compute the difference ....
                if (intOffset >= 2 ) lastDelta = buffer[int(intOffset)-1] - buffer[int(intOffset)-2];
                else if (intOffset >= 1 ) lastDelta = buffer[int(intOffset)-1] - lastValue;
                
                // DEFAULT :: assume flat ...
                currentDelta = 0;
                
                if (intOffset + 1 < numSamples)
                    currentDelta = buffer[int(intOffset)+1] - buffer[int(intOffset)];
                
                // CALCULATE change in velocity
                double change_in_delta = (currentDelta - lastDelta);
                double propDepth = proportionalBlepFreq;
                
                magnitude_velocity = -4*change_in_delta*(1/propDepth);
                
                jassert(magnitude_velocity <= 1);
                
            }
            
            // 0th order (position)
            {
                // CALCULATE the magnitude of the 0 order nonlinearity *change in position*
                float extrapolated_last_pos = lastValue + lastDelta*(fraction);
                float extrapolated_jump_pos = buffer[int(intOffset)] - currentDelta*(1 - fraction);
                magnitude_position = extrapolated_last_pos - extrapolated_jump_pos;
            }
            
            
        }
        
        
        // TOO SMALL :::
        /// no need to compensating for tiny discontinuities
        if (fabsf(magnitude_position) < .001 && fabsf(magnitude_velocity) < .001)
        {
            currentActiveBlepOffsets.remove(i);
            continue;
        }
        
        // NEGLIGIBLE MAGNITUDES :::
        /// zero out any tiny effects here, so we don't waste time calculating them
        if (fabsf(magnitude_position) < .001) magnitude_position = 0;
        if (fabsf(magnitude_velocity) < .001) magnitude_velocity = 0;
        
        
        // ADD ::::
        // GAIN factors ... how big of a discontinutiy are we talking about ?
        blep.pos_change_magnitude = magnitude_position;
        blep.vel_change_magnitude = magnitude_velocity;
        
        // ALTER :::
        currentActiveBlepOffsets.setUnchecked(i, blep);
        
    }
    
    
}

void MinBlepGenerator::process_currentBleps(float* buffer, int numSamples) {
    
    
    
    // PROCESS ALL BLEPS -
    /// for each offset, copy a portion of the blep array to the output ....
    for (int i = currentActiveBlepOffsets.size(); --i >= 0;)
    {
        
        BlepOffset blep = currentActiveBlepOffsets[i];
        double adjusted_Freq = blep.freqMultiple;
        double exactPosition = blep.offset;
        
        
        // ADD the BLEP to the circular buffer ...
        for (float p = 0; p < numSamples; p++)
        {
            
            double blepPosExact = adjusted_Freq*(exactPosition + p + 1); // +1 because this needs to trigger on the LOW SAMPLE
            double blepPosSample = 0;
            double fraction = modf(blepPosExact, &blepPosSample);
            
            // LIMIT the scaling on the derivative array
            // otherwise, it can get TOO large
            double depthLimited = proportionalBlepFreq; //jlimit<double>(.1, 1, proportionalBlepFreq);
            double blepDeriv_PosExact = depthLimited*overSamplingRatio*(exactPosition + p + 1);
            double blepDeriv_Sample = 0;
            double fraction_Deriv = modf(blepDeriv_PosExact, &blepDeriv_Sample);
            
    
            // DONE ... we reached the end ...
            if (int(blepPosExact) > minBlepArray.size() && int(blepDeriv_PosExact) > minBlepArray.size())
                break;
            
            // BLEP has not yet occurred ...
            if (blepPosExact < 0)
                continue;
            
            
            // 0TH ORDER COMPENSATION ::::
            /// add the BLEP to compensate for discontinuties in the POSITION
            if ( fabs(blep.pos_change_magnitude) > 0 && blepPosSample < minBlepArray.size())
            {
                // LINEAR INTERPOLATION ::::
                float lowValue = minBlepArray.getRawDataPointer()[int(blepPosSample)];
                float hiValue = lowValue;
                
                if (int(blepPosSample) + 1 < minBlepArray.size())
                    hiValue = minBlepArray.getRawDataPointer()[int(blepPosSample) + 1];
                
                float delta = hiValue - lowValue;
                float exactValue = lowValue + fraction*delta;
                
                // SCALE by the discontinuity magnitude
                exactValue *= blep.pos_change_magnitude;
                
                // ADD to the thruput
                buffer[int(p)] += exactValue;
            }
            
            
            // 1ST ORDER COMPENSATION ::::
            /// add the BLEP DERIVATIVE to compensate for discontinuties in the VELOCITY
            if ( fabs(blep.vel_change_magnitude) > 0 && blepDeriv_PosExact < minBlepDerivArray.size())
            {
                
                // LINEAR INTERPOLATION ::::
                double lowValue = minBlepDerivArray.getRawDataPointer()[int(blepDeriv_PosExact)];
                double hiValue = lowValue;
                
                if (int(blepDeriv_PosExact) + 1 < minBlepDerivArray.size())
                    hiValue = minBlepDerivArray.getRawDataPointer()[int(blepDeriv_PosExact) + 1];
                
                double delta = hiValue - lowValue;
                double exactValue = lowValue + fraction_Deriv*delta;
                
                // SCALE by the discontinuity magnitude
                exactValue *= blep.vel_change_magnitude;
                
                // ADD to the thruput
                buffer[int(p)] += exactValue;
                
            }
            
        }
        
        // UPDATE ::::
        blep.offset = blep.offset + double(numSamples);
        if (blep.offset*adjusted_Freq > minBlepArray.size())
        {
            currentActiveBlepOffsets.remove(i);
        }
        else currentActiveBlepOffsets.setUnchecked(i, blep);
        
    }

    
    
}




