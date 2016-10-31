//
//  wavegen.cpp
//  Stagecraft Software
//
//  Created by Aaron Leese on 9/22/11.
//  Copyright 2011 Stagecraft Software. All rights reserved.
//


#include "JuceHeader.h" 
#include "WaveGenerator.h"


// DEBUG UI :::
static double test = 1;
bool clearWave = false;
void WaveGenerator::setTest(double newTest) {
    test = newTest;
}
void WaveGenerator::setClear(bool clear) {
    clearWave = clear;
}



// WAVE GEN :::::
WaveGenerator::WaveGenerator() {
    
    historyLength = 500;
    sampleRate = 0;
    
    myMode = NO_ANTIALIAS;
    myWaveType = square;
    
    slaveDeltaBase = 0;
    actualCurrentAngleDelta = 0;
    
    currentAngle = 0.01; // needed
    currentAngleSkewed = lastAngleSkewed = 0;
    pitchBendTarget = pitchBendActual = 1.0f;
    
    masterPitchOffset = 1;
    slavePitchOffset = 1;
    
    volume = 1;
    pan = 0.5;
    gainLast[0] = gainLast[1] = 0;
    skew = 0;
    lastSample = 0;
    
    stereo = false;
    
    for (int i=0; i < historyLength; i++)
        history.add(0);
    
    phaseAngleTarget = phaseAngleActual = 0; // expressed 0 - 2*PI
    
    // configure the blep multiple
    blepOvertoneDepth = 128; // high default ... buzzy
    
    
}

void WaveGenerator::prepareToPlay (double newSampleRate) {
    
    sampleRate = newSampleRate;
    
    // BUILD the appropriate BLEP step ....
    myBlepGenerator.buildBlep();
    
}
void WaveGenerator::changeBlepSize(float newOverSample) {
    
    jassertfalse; // need to rebuild, no ?
    
    myBlepGenerator.overSamplingRatio = newOverSample;
    
}

void WaveGenerator::setBlepOvertoneDepth(double newMult) {

    // SET the allowed depth of overtones
    // ** NOTE **
    // this is relative to the fundamental F
    // so each factor of 2 is an octave
    blepOvertoneDepth = newMult;
    
}

void WaveGenerator::clear() {
    
    currentAngle = phaseAngleTarget; // + phase !!!
    currentMasterAngle = phaseAngleTarget;
    
    pitchBendTarget = pitchBendActual = 1.0;
    slaveDeltaBase = 0.0;
    gainLast[0] = gainLast[1] = 0;
}

// FAST RENDER (AP) :::::
void WaveGenerator::renderNextBlock (AudioSampleBuffer& outputBuffer, int numSamples) {
    
    jassert(sampleRate != 0);
 
    if (slaveDeltaBase == 0.0) return;
    
    if (volume == 0 && gainLast[0] == 0 && gainLast[1] == 0)
        return;
    
    buildWave(numSamples);
    

    // ADD BAND-LIMITED (minBLEP) transitions :::
    if ( myMode == ANTIALIAS )
    {
        // Since we KNOW the intended F ... relative to F(sampling)
        // We can tweak the minBLEP to limit any harmonic above 4*(desired F)
        
        // TUNE the blep ....
        double freq = getCurrentPitchHz(); // Current, playing, Freq
        double relativeFreq = 2*freq/sampleRate; // 2 for Nyquist ...
        relativeFreq *= blepOvertoneDepth; // ie - up to the 3nd harmonic (2*2*2 -> 8x fundamental)
    
        myBlepGenerator.setLimitingFreq(relativeFreq); // up to the 2nd harmonic ..
        myBlepGenerator.processBlock( wave.getRawDataPointer(), numSamples);
        
    }
    
    
    // BUILD ::::
    for (int i=0; i < numSamples; i = i+20)
    {
        // just adding a sample every 20 or so to the history
        history.add(wave.getUnchecked(i));
        
        if (history.size() > historyLength)
            history.remove(0);
    }
    
    
    // PAN :::
    // CALCULATE the pan factors of each input, for each output ....
    // calculate the pan pots .... using -3 center position
    double panGain[2];
    double panL_dB = -3*pow(double(2.0f - 2.0f*pan), double(3.0));    // 0 => 0dB, .5 => -3 dB, 1 => -24db
    panGain[0] = Decibels::decibelsToGain(panL_dB);
    double panR_dB = -3*pow(double(2.0f*pan), 3.0);    // 0 => 0dB, .5 => -3 dB, 1 => -24db
    panGain[1] = Decibels::decibelsToGain(panR_dB);
    
    // STEREO :::
    double invert = 1;
    if (stereo) invert = -1;
    
    
    // COPY it to the outputbuffer ....
    outputBuffer.addFromWithRamp(0, 0, wave.getRawDataPointer(), numSamples, gainLast[0], volume*panGain[0]);
    outputBuffer.addFromWithRamp(1, 0, wave.getRawDataPointer(), numSamples, invert*gainLast[1], invert*volume*panGain[1]);
    
    gainLast[0] = volume*panGain[0];
    gainLast[1] = volume*panGain[1];
    
    
#if JUCE_DEBUG
    
    // Using dB -12 to +3 ... normalized to 0..1 for indicator
    float RMS = outputBuffer.getRMSLevel(0, 0, outputBuffer.getNumSamples());
    
    // Aaron - wtf?  How do we get NaN .. but we do .... hmmmm
    if (RMS != RMS)
    {
        DBG("NaN " + String(outputBuffer.getNumSamples())
            + " " + String(*outputBuffer.getReadPointer(0, 0)) );
    }
    
    if (RMS > 5 || RMS < 0)
    {
        DBG("WOAH! " + String(RMS)
            + " " + String(*outputBuffer.getReadPointer(0, 0)) );
    }

#endif
    
    
    
}
inline void WaveGenerator::buildWave (int numSamples) {
    
    if (slaveDeltaBase == 0.0) return;
    
    if (numSamples == 0) return;
    jassert(numSamples > 0);
    
    if (wave.size() != numSamples)
        wave.resize(numSamples);
    
    float* waveData = wave.getRawDataPointer();
    
    // LINEAR CHANGE for now .... over the numsamples (see above)
    double freqDelta = (pitchBendTarget - pitchBendActual)/double(numSamples);
    JUCE_SNAP_TO_ZERO(freqDelta);
    
    // SKEWING :::::
    // change the freqDelta - up to 2x the speed, when at PI, and corresponding slowness at 0 ...
    double phaseShiftPerSample = 0;
    if (phaseAngleTarget != phaseAngleActual)
    {
        // LINEAR RAMP from current phase to the target
        phaseShiftPerSample = (phaseAngleTarget - phaseAngleActual)/float(2*numSamples); // 2 here is smoooothing ...
        phaseAngleActual = phaseAngleActual + phaseShiftPerSample*numSamples;
    }
   
   
    // BUILD ::::
    for (int i=0; i < numSamples; i++)
    {
 
        bool masterBlepOccurred = false;
        
        // CHANGE the PITCH BEND (linear ramping)
        pitchBendActual += freqDelta;
        if (fabs(pitchBendActual - 1) < .00001) pitchBendActual = 1;
        
        // FOR CALCULATIONS,
        // note the current, actual delta (base pitch modified by pitch bends and phase shifting)
        actualCurrentAngleDelta = slaveDeltaBase*pitchBendActual + phaseShiftPerSample;
        
        // HARD SYNCING ::::::
        if (hardSync && masterDeltaBase != slaveDeltaBase )
        {
            // MASTER OSC DOES use pitch bending and phase shifting ...
            double actualCurrentMasterDelta = masterDeltaBase*pitchBendActual + phaseShiftPerSample;
            currentMasterAngle += actualCurrentMasterDelta;
            
            // ADD A BLEP ::::
            currentMasterAngle = fmod(currentMasterAngle, double(2*double_Pi));
            
            // MASTER (unskewed) rollover
            if (currentMasterAngle < actualCurrentMasterDelta)
            {
                double percAfterRoll = currentMasterAngle/actualCurrentMasterDelta;
                double percBeforeRoll = 1 - percAfterRoll;
                
                // ADD the blep ...
                {
                    MinBlepGenerator::BlepOffset blep;
                    blep.offset = percAfterRoll - double(i + 1);
                    
                    
                    // CALCULATE the MAGNITUDE of ths 2nd ORDER (VEL) discontinuity
                    // TRIG :: calculate the angle (rise/run) before and after the rollover
                    double delta = .0000001; // MIN
                    
                    double angleAtRoll = fmod(currentAngle + percBeforeRoll*actualCurrentAngleDelta, 2*double_Pi);
                    double angleBeforeRoll = fmod(angleAtRoll - delta, 2*double_Pi);
                    
                    
                    // SKEW ALL ANGLES !
                    double valueBeforeRoll = getValueAt(skewAngle(angleBeforeRoll));
                    double valueAtRoll = getValueAt(skewAngle(angleAtRoll));
                    
                    double valueAtZero = getValueAt(skewAngle(0));
                    double valueAfterZero = getValueAt(skewAngle(delta));
                    
                    
                    // CALCULATE the MAGNITUDE of ths 1st ORDER (POS) discontinuity
                    blep.pos_change_magnitude = valueAtRoll - getValueAt(skewAngle(0));
                    

                    // CALCULATE the skewed angular change AFTER the rollover
                    double angle_delta_after_roll = valueAfterZero - valueAtZero; // MODs based on the PITCH BEND ...
                    double angle_delta_before_roll = valueAtRoll - valueBeforeRoll; // MODs based on the PITCH BEND ...
                    
                    double change_in_delta = (angle_delta_after_roll - angle_delta_before_roll)*(1/(2*delta));
                    double depthLimited = jlimit<double>(.1, .5, myBlepGenerator.proportionalBlepFreq);
                    
                    // 1.3 here was experimentally determined ...
                    blep.vel_change_magnitude = 1.3*test*change_in_delta*(1/depthLimited);

                    // ADD
                    myBlepGenerator.addBlep(blep);
                }
                
                
                // MOVE the UNSKEWED ANGLE
                // so that it will actually roll over at this sub-sample ...
                // ESTIMATE !!!!! ERROR - this should be better, but we can't unskew (x = x/cos(x) is not solvable)
                currentAngle = 2*double_Pi - percBeforeRoll*actualCurrentAngleDelta;
                
                //
                masterBlepOccurred = true;
                
            }
           
        }
         
        
        // MOVE the ANGLE
        currentAngle += actualCurrentAngleDelta; // MODs based on the PITCH BEND ...
        currentAngle = fmod(double(currentAngle), double(2*double_Pi)); // ROLLOVER :::
        
        jassert(currentAngle == currentAngle);
    
        // APPLY SKEWING :::::
        currentAngleSkewed = skewAngle(currentAngle);
        
        
        // BUILD the antialiasing ....
        if (myMode != NO_ANTIALIAS && masterBlepOccurred == false && myWaveType != sine)
        {
            
            double actualCurrentAngleDeltaSkewed = currentAngleSkewed - lastAngleSkewed;
            if (actualCurrentAngleDeltaSkewed < 0) actualCurrentAngleDeltaSkewed += 2*double_Pi;
            
           
            // ROLLED through 2*PI
            if (myWaveType == square)
            {
                // :: SQUARE rolls twice - at PI and 2*PI ::::
                if (fmod(currentAngleSkewed, double_Pi) < actualCurrentAngleDeltaSkewed)
                {
                    double percAfterRoll = fmod(currentAngleSkewed, double_Pi)/actualCurrentAngleDeltaSkewed; // LINEAR interpolation
                    
                    MinBlepGenerator::BlepOffset blep;
                    
                    // CALCULATE subsample offset :::
                    blep.offset = percAfterRoll - double(i + 1);
                    
                    // MAGNITUDE of 1st order nonlinearity is 2 or -2 :::
                    if (fmod(currentAngleSkewed, 2*double_Pi) < actualCurrentAngleDeltaSkewed)
                        blep.pos_change_magnitude = -2;
                    else blep.pos_change_magnitude = 2;
                    
                    // NO CHANGE to slope - 0
                    blep.vel_change_magnitude = 0;
                    
                    // ADD
                    myBlepGenerator.addBlep(blep);
                
                }
            }
            else if (myWaveType == sawRise || myWaveType == sawFall ) // SAW
            {
             
                // SAW ROLLs only at PI
                if (fmod(currentAngleSkewed, 2*double_Pi) > actualCurrentAngleDeltaSkewed
                    && fmod(currentAngleSkewed, double_Pi) < actualCurrentAngleDeltaSkewed
                    )
                {
                    
                    double percAfterRoll = fmod(currentAngleSkewed, double_Pi)/actualCurrentAngleDeltaSkewed; // LINEAR interpolation
                    
                    // CALCULATE the OFFSET
                    MinBlepGenerator::BlepOffset blep;
                    blep.offset = percAfterRoll - double(i + 1);
                    
                    // MAGNITUDE of 1st order nonlinearity is 2 or -2 :::
                    if (myWaveType == sawRise) blep.pos_change_magnitude = -2;
                    else blep.pos_change_magnitude = 2;
                  
                    // NO CHANGE to slope - 0
                    blep.vel_change_magnitude = 0;
                    
                    // ADD
                    myBlepGenerator.addBlep(blep);
                
                }
            }
            else if (myWaveType == triangle)
            {
                
                    
                if (fmod(currentAngleSkewed + double_Pi/2, 2*double_Pi) < actualCurrentAngleDeltaSkewed
                    || fmod(currentAngleSkewed + 3*double_Pi/2, 2*double_Pi) < actualCurrentAngleDeltaSkewed
                    )
                {
                    
                    
                    double aboveNonlinearity = 0;
                    double percAfterRoll = 0;
                    
                    if (fmod(currentAngleSkewed + 3*double_Pi/2, 2*double_Pi) < actualCurrentAngleDeltaSkewed)
                    {
                        aboveNonlinearity = fmod(currentAngleSkewed + 3*double_Pi/2, 2*double_Pi);
                        percAfterRoll = aboveNonlinearity/actualCurrentAngleDeltaSkewed;
                    }
                    else // 3*double_Pi/2
                    {
                        aboveNonlinearity = fmod(currentAngleSkewed + double_Pi/2, 2*double_Pi);
                        percAfterRoll = aboveNonlinearity/actualCurrentAngleDeltaSkewed;
                    }
                    
                    MinBlepGenerator::BlepOffset blep;
                    blep.offset = percAfterRoll - double(i + 1);
                    
                    // SYMETRY :::::
                    // since this is a triangle
                    // ... we can average the two values,
                    // get the abs distance from 1
                    // and scale that to find the angle ....
                    
                    double nextValue = getValueAt(currentAngleSkewed);
                    double averageValue = (lastSample + nextValue)/2;
                    double slope = 1 - fabs(averageValue);
                    jassert(slope < 1);
                    
                    double sign = 1;
                    if (averageValue > 0) sign = -1;
                    
                    blep.pos_change_magnitude = 0;
                    
                    // SCALE the vel magnitude inversely with play speed
                    double depthLimited = jlimit<double>(.1, .5, myBlepGenerator.proportionalBlepFreq);
                    
                    // Assume nominal delta for all waves ... so ...
                    blep.vel_change_magnitude = sign*test*130*slope*(1/depthLimited);
                    
                    // ADD
                    myBlepGenerator.addBlep(blep);
                
                }
            }
            
            
        }
        
        
        *waveData = (float)getValueAt(currentAngleSkewed);
        
        // UPDATE the tracking variables ...
        // Used or computing exact values at rolls, etc.
        lastAngleSkewed = currentAngleSkewed; // NOTE the previous angle, for calculations
        lastSampleDelta = *waveData - lastSample;
        lastSample = *waveData; // NOTE* the most recent sample, for computation purposes
        
        waveData++;
    }
    
     
    jassert(wave.size() >= numSamples);
    
}

double WaveGenerator::skewAngle(double angle) {
    
    // APPLY PWM to angle ...
    
    // SKEW ANGLE :::::
    double skewedAngle = angle;
    if (skew != 0)
    {
        // ADD the INTEGRAL of SIN .. which is cos ... to the angle ...
        skewedAngle = angle + skew*cos(angle);
    }
    
        
    skewedAngle = fmod(skewedAngle, 2*double_Pi);
    if (skewedAngle < 0.0) skewedAngle += 2*double_Pi;

    jassert(skewedAngle >= 0 && skewedAngle <= 2*double_Pi);
    
    return skewedAngle;
}
double WaveGenerator::getValueAt(double angle) {
    
    jassert(angle >= 0 && angle <= 2*double_Pi);
    
    double currentSample = 0;
    
    if (myWaveType == sine) currentSample = getSine(angle);
    else if (myWaveType == sawRise)  currentSample = getSawRise(angle);
    else if (myWaveType == sawFall)  currentSample = getSawFall(angle);
    else if (myWaveType == triangle)  currentSample = getTriangle(angle);
    else if (myWaveType == square)  currentSample = getSquare(angle);
    else if (myWaveType == random)  currentSample = getRandom(angle);
    
    return currentSample;
}


// SLOW RENDER (LFO) ::::::
void WaveGenerator::moveAngleForward(int numSamples) {
    
    // Not moving ...
    if (numSamples == 0 || slaveDeltaBase == 0) return;
    
    double angle = currentAngle;
    int historyPoints = numSamples/20;
    
    // calculate MOD in the angle delta ...
    double modAngleDelta = slaveDeltaBase;
    
    // DISREGARDING any PITCH SHIFTING (for now, it's just not done anywhere on LFOs)
    jassert (pitchBendTarget == pitchBendActual);
   
    // CALCULATE the PHASE CHANGE per sample ...
    if (phaseAngleTarget != phaseAngleActual)
    {
        // LINEAR RAMP from current phase to the target
        float phaseShiftPerSample = (phaseAngleTarget - phaseAngleActual)/float(100*numSamples);
        modAngleDelta += phaseShiftPerSample;
        phaseAngleActual = phaseAngleActual + phaseShiftPerSample*numSamples;
    }
    
    
    // UPDATE the history
    for (int i=0; i < historyPoints; i++)
    {
        angle = currentAngle + i*(20*modAngleDelta);
        angle = fmod(angle, 2*double_Pi); // modulus
        angle = skewAngle(angle);
        
        history.add(getValueAt(angle));
        
        if (history.size() > historyLength)
            history.remove(0);
    }
    
    currentAngle = fmod(currentAngle + numSamples*modAngleDelta, 2*double_Pi); // ROLL
    
    jassert(currentAngle == currentAngle);
    
}
void WaveGenerator::moveAngleForwardTo(double newAngle) {
    

    double delta = newAngle - currentAngle;
    if (delta < 0) delta = delta + 2*double_Pi;
    
    double numSamples = delta/slaveDeltaBase;
    jassert(numSamples == numSamples);
    
    moveAngleForward(numSamples);
    
}

// LOAD / SAVE
void WaveGenerator::loadScene(ValueTree node) {
    
    // wavetype ....
    int type = (int)node.getProperty("waveType", 0);
    setWaveType((WaveType)type);
    
    masterPitchOffset = (double)node.getProperty("masterPitchOffset", 1.0);
    
    volume = (double)node.getProperty("volume", volume);
    pan = (double)node.getProperty("pan", pan);
    skew = (double)node.getProperty("skew", skew);
    
    blepOvertoneDepth = (double)node.getProperty("blepOvertoneDepth", 64.f);
    phaseAngleTarget = (double)node.getProperty("phaseAngleTarget", phaseAngleTarget);
    
    stereo = (bool)node.getProperty("stereo", false);
    
    // ALWAYS HARD SYNC (for now)
    hardSync = true; //(bool)node.getProperty("hardSync", false);
    
    // TONE :::
    float tone = (double)node.getProperty("toneOffset", 0);
    setToneOffset(tone);
    
    
}
void WaveGenerator::saveScene(ValueTree node) {
    
    node.setProperty("waveType", (int)myWaveType, 0);
    node.setProperty("masterPitchOffset", masterPitchOffset, 0);
    
    node.setProperty("volume", volume, 0);
    node.setProperty("pan", pan, 0);
    node.setProperty("skew", skew, 0);
    
    node.setProperty("blepOvertoneDepth", blepOvertoneDepth, 0);
    node.setProperty("phaseAngleTarget", phaseAngleTarget, 0);
    
    node.setProperty("stereo", stereo, 0);
    node.setProperty("hardSync", hardSync, 0); // always true!
    
    // TONE :::
    node.setProperty("toneOffset", getToneOffsetInSemis(), 0);
    
    
}




