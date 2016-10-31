//
//  wavegen.h
//  Stagecraft Software
//
//  Created by Aaron Leese on 9/22/11.
//  Copyright 2011 Stagecraft Software. All rights reserved.
//

#pragma once
 
#include "JuceHeader.h"  
#include "minBlepGenerator.h"

// ==== SIMPLER LFO - non AP ==========================================================================
class WaveGenerator
{
     
    MinBlepGenerator myBlepGenerator;
    
    double blepOvertoneDepth = 1; // tweak the "resolution" of the blep ... (multiple of highest harmonics allowed)
    
    // MASTER - HARD SYNC ::::
    double masterDeltaBase = 0; // for hard syncing
    double masterPitchOffset = 1; // relative offset ....
    double currentMasterAngle = 0;
    
    // SLAVE - HARD SYNC ::::
    double slaveDeltaBase = 0;
    double slavePitchOffset = 1; // relative offset ....
    double currentAngle = 0;
    double currentAngleSkewed = 0;
    double lastAngleSkewed = 0;
    
    
    // ACTUAL OUTPUT
    double lastSample = 0;
    double lastSampleDelta = 0;
    
    // PITCH BEND
    double pitchBendTarget = 0;
    double pitchBendActual = 0;
    
    // ACTUAL delta
    // - takes into account pitch bend and phase shift ...
    double actualCurrentAngleDelta = 0; // ACTUAL CURRENT DELTA (pitch)
    
    double volume = 1;
    double pan = 0.5; // [0..1]
    double gainLast[2] = {0, 0}; // for ramping ..
    double skew = 0; // [-1, 1]
    double sampleRate = 0;
    
    Array<float> wave; // hmmm ... faster to make this a pointer I bet ....
    Array<float> history; // a running averaged wave, for rendering purposes
    int historyLength = 500;
    
    bool stereo = false; // invert the right channel ...
    double phaseAngleTarget = 0;
    double phaseAngleActual = 0; // the target angle to get to (used for phase shifting)
    
    bool hardSync = true;

public:
    
    void setTest(double newTest);
    void setClear(bool clear);
    
    enum WaveType {
        
        sine        = 0,
        sawRise     = 1,
        sawFall     = 2,
        triangle    = 3,
        square      = 4,
        random      = 5
    };
    WaveType myWaveType;
    
    // these can be used to generate LFO waves
    // or high speed (synth) waveforms ...
    enum WaveMode {
        
        ANTIALIAS,
        BUILD_AA,
        NO_ANTIALIAS
                
    };
    WaveMode myMode;
    
    WaveGenerator ();
    
    void prepareToPlay (double newSampleRate);
    
    void setWaveType(WaveType newWaveType) {
        
        myWaveType = newWaveType;
        
        if (myWaveType == triangle) myBlepGenerator.setToReturnDerivative(true);
        else if (myWaveType == sine) myBlepGenerator.setToReturnDerivative(true);
        else  myBlepGenerator.setToReturnDerivative(false);
        
    }
    WaveType getWaveType() {
        
        return myWaveType;
        
    }
    
    void setMode(WaveMode mode) {
        
        myMode = mode;
        
        // BUILD the appropriate BLEP step ....
        if (myMode == ANTIALIAS)
        {
            myBlepGenerator.buildBlep();
        }
        
    }
    WaveMode getAntialiasMode() {
        
        return myMode;
    }

    // minBLEP :::
    void setBlepOvertoneDepth(double newMult);
    MinBlepGenerator* getBlepGenerator() {
        return &myBlepGenerator;
    }
    void changeBlepSize(float sampleRateForBlep);
    
    //
    double getSlaveDeltaBase() {
        
        return slaveDeltaBase;
    }
    double getMasterDeltaBase() {
        
        return masterDeltaBase;
    }
    double getAngleDeltaActual() {
        
        return actualCurrentAngleDelta;
    }
    
    void setMasterDelta(double newAngleDelta) {
        
        // MASTER OSC ::::
        masterDeltaBase = masterPitchOffset*newAngleDelta;
        
        // SLAVE OSC ::::
        slaveDeltaBase = slavePitchOffset*newAngleDelta;
        
    }
	double getCurrentAngle() {

		return currentAngle;
	};

    void setPitchSemitone(int midiNoteValue, double sampleRate) {
        
        double centerF = MidiMessage::getMidiNoteInHertz ( midiNoteValue ); // ....
        double cyclesPerSample = centerF / sampleRate;
        float angleDelta = cyclesPerSample * 2.0 * double_Pi;
        
        setMasterDelta(angleDelta);
    }
    void setPitchHz(double newFreq) {
        
        double cyclesPerSample = newFreq / sampleRate;
        float angleDelta = cyclesPerSample * 2.0 * double_Pi;
        setMasterDelta(angleDelta);
    }

    double getCurrentPitchHz() {
        
        // float angleDelta = cyclesPerSample * 2.0 * double_Pi;
        double cyclesPerSample = actualCurrentAngleDelta/( 2.0*double_Pi);
        double Freq = cyclesPerSample*sampleRate;
        
        return Freq;
        
    }
    
    Array<float> getBLEPArray();
    
    float getSkew() { return skew; }
    void setSkew(double newSkew) {
        
        jassert(newSkew >= -1 && newSkew <= 1);
        skew = newSkew;
    }
    
    double getPhaseTarget() {
       
        return phaseAngleTarget;
    }
    void setPhaseTarget(double angleToGetTo) {
        
        jassert(angleToGetTo >= -double_Pi && angleToGetTo <= double_Pi);
        phaseAngleTarget = angleToGetTo;
    }
    
    // PITCH MOD ::: shifts the MASTER angleDelta up/down in semitones ...
    void setPitchOffset(double pitchOffsetInSemitones) {
        
        // TONE is ALWAYS higher than master center (or has no effect)
        double slaveToneOffset = getToneOffsetInSemis();
        
        // Convert from semitones to * factor
        masterPitchOffset = (pow(double(2), double(pitchOffsetInSemitones/12)) );
       
        // UPDATE the slave freq
        double newToneOffset = pitchOffsetInSemitones + slaveToneOffset;
        slavePitchOffset = (pow(double(2), double(newToneOffset/12)) );
        
    }    
	double getPitchOffsetInSemis() {
        
        // return the Log pitch offset ....
        double pitchOffsetInSemis = 12*log2(masterPitchOffset);
        return pitchOffsetInSemis;
    }
    
    void setToneOffset(double newToneOffsetInSemis) {
        
        // Convert from semitones to * factor
        double masterSemis = getPitchOffsetInSemis();
        
        // TONE is ALWAYS higher than master center (or has no effect)
        double newToneOffset = masterSemis + newToneOffsetInSemis;
        slavePitchOffset = (pow(double(2), double(newToneOffset/12)) );
        
    }
    double getToneOffsetInSemis() {
        
        // return the Log pitch offset ....
        double toneOffsetInSemis = 12*log2(slavePitchOffset);
        
        // RELATIVE to master ...
        double masterSemis = getPitchOffsetInSemis();
        
        return (toneOffsetInSemis - masterSemis);
        
    }
    
    void setPitchBend(double newBendInSemiTones) {
        
        jassert(newBendInSemiTones == newBendInSemiTones);
        
        // For EQUAL TEMPERMENT :::::
        pitchBendTarget = pow (2.0, (newBendInSemiTones / 12.0));
         
		//jassert(pitchBendTarget != std::isin);
        jassert(pitchBendTarget == pitchBendTarget);
        
    }
	double getPitchBendSemis() {

		return 12*log2(pitchBendActual);
	}
	
    void setVolume(double dBMult) {
        
        if (dBMult <= -80) volume = 0;
        else volume = Decibels::decibelsToGain(dBMult);
	}
    void setPan(double proportionalPan) {
        
		jassert(proportionalPan >= 0);
		jassert(proportionalPan <= 1);
		pan = proportionalPan;
	}
    void setStereo(bool isStereo) {
        stereo = isStereo;
    }
    void setHardsync(bool shouldHardSync) {
        
        hardSync = shouldHardSync;
    }
    
    bool isStereo() { return stereo; }
    bool isHardSync() { return hardSync; }
    
    Array<float> getHistory() {
        
        return history;
    }
    
    void clear();
    
    // FAST RENDER (AP) :::::
    void renderNextBlock (AudioSampleBuffer& outputBuffer, int numSamples);
    inline void buildWave (int numSamples);
     
    // SLOW RENDER (LFO) ::::::
    void moveAngleForward(int numSamples);
    void moveAngleForwardTo(double newAngle);
    
    double getAngleAfter(double samplesSinceRollover) {
		
		// CALCULATE the WAVEFORM'S ANGULAR OFFSET
        // GIVEN a certain number of samples (since rollover) ....
        
        // since this is only done in LFO mode ....
        // reset phase so there is no changing ...
        phaseAngleActual = phaseAngleTarget;
       
        return slaveDeltaBase*samplesSinceRollover + phaseAngleActual;

	}
    double skewAngle(double angle);
     
    // GET value
    double getCurrentValue() {
        
        // SKEW the "current angle"
        double skewedAngle = skewAngle(currentAngle);
        return getValueAt(skewedAngle);

    }
    double getValueAt(double angle);
    
    
    // Wave calculations ...
	inline double getSine(double angle) {
		 
		const double sample = sin(angle);
		return sample;
		 
	}
	inline double getSawRise(double angle) {
        
        // remainder ....
		double sample = getSawFall(angle);
        
        // JUST INVERT IT NOW ... to get rising ...
        sample = -sample;
        
		return sample;
        
	}
	inline double getSawFall(double angle) {
        
        angle = fmod(angle + double_Pi, 2*double_Pi); // shift x
        double sample = angle/double_Pi - double(1); // computer as remainder
    
		return sample;
        
	}
	inline double getTriangle(double angle) {
 
		double sample = 0;
		
        // using a simple offset, we can make this a ramp up, then down ....
        //angle += double_Pi/4;
        angle = fmod(angle + double_Pi/2, double(2*double_Pi)); // ROLL
        
        double frac = angle/(2*double_Pi);

		if (frac < .5) sample = 2*frac; // RAMPS up
		else sample = 1 - 2*(frac - .5); // RAMPS down

		// SCALE and Y-OFFSET
		sample = (sample - .5)*2; // so it goes from -1 .. 1

		jassert(sample <= 1);
		jassert(sample >= -1);

		return sample;	 

	}
	inline double getSquare(double angle) {
        
        if (angle >= double_Pi) return -1;
        return 1;
        
	}
    inline double getRandom(double angle) {
        
        double r = Random::getSystemRandom().nextFloat();
        
        r = 2*(r - 0.5f); // scale to -1 .. 1
        r = jlimit(-10*slaveDeltaBase, 10*slaveDeltaBase, r);
        
        lastSample += r;
    
        // limit it based on the F ...
        lastSample = jlimit(-1.0, 1.0, lastSample);
        
        return lastSample;
        
	}

    
    // LOAD / SAVE ::::
    void loadScene(ValueTree node);
    void saveScene(ValueTree node);
};


 



