/*
 *  accelerateFFT.h
 *  Livetronica Studio
 *
 *  Created by Aaron Leese on 7/15/10.
 *  Copyright 2010 StageCraft Software. All rights reserved.
 *
 */


#pragma once

#include "JuceHeader.h"
#include "kiss_fft130/kiss_fft.h"
#include "kiss_fft130/tools/kiss_fftr.h"


class SpectralData;

class FFTProcessor : public AudioProcessor
{
	
private:
    
    Array<float> hanningWindow;
    Array<float> magnitudes;
    Array<float> phaseAngles;
    
    kiss_fft_cpx* input; // complex in
    kiss_fft_cpx* output; // complex result / output for FFT (input for iFFT)
    
    kiss_fft_cfg kiss_FFT_config;
    kiss_fft_cfg kiss_iFFT_config;
    
    kiss_fftr_cfg kiss_FFT_real_config;
    kiss_fftr_cfg kiss_iFFT_real_config;
    
    
    // For real only :::::
    kiss_fft_scalar* timedata; // for REAL only input
    
public:
    
    Array<float> incomingData;
    int recordMarker;
    int lastBufWritePos;

    int FFT_Resolution; 
    bool realOnlyFFT; // set the FFT to use only REAL input (no imaginery numbers)
    
public:
    
    bool stopCalculations;
	bool active; // for controlling if it repaints ...
    uint32 lastMillSecCount;
    
    // Pitch Detection
    StringArray notes; // = {"C ","C#","D ","D#","E ","F ","F#","G ","G#","A ","A#","B "};
    const float base_a4; //=440.0; // set A4=440Hz
    
    int currentRoot; // root notes, C,C#,etc.
    bool scale[12]; 
    float scaleIntervals[12];
    bool continuous;    
     
    unsigned short currentNote;
    unsigned short lastNote;
    
    unsigned int cents; // for pitch bending from input
    unsigned int detune; // for pitch wheel
    double currentNoteStartTime;
    
    float peakIndex;
    float currentSum;
	uint32 lastBeat;

	FFTProcessor();
	~FFTProcessor();	
	 
    int getResolution() {
        
        return incomingData.size();
    }
    
    void setResolution(int resolution, bool realDataOnly = true);
    
	float getFreqPeakPercent() {
        
        //
        float size = magnitudes.size();
        return peakIndex/size; 
    }

	float getCurrentSum() {
		 
		return currentSum;
	}
    
    Array<float> getHanningWindow() {
        
        return hanningWindow;
        
    }
    
    kiss_fft_cpx* getOutput() {
    
        return output;
    }
    
    kiss_fft_cpx* getCurrentFullFFT(int numSamples);
    
	// :::::::::::::: AudioProcessor Methods
    //============================================================================== PROCESSOR CALLBACK METHODS
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
	
    //==============================================================================
    void processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages) override;
    
    void processBlock (AudioBuffer<double>& buffer, MidiBuffer& midiMessages) override
    {
        // NOT implemented for double precision currently ...
        jassertfalse;
        
    }
    
    
    // ACTUALLY CALCULATE the FFT
    void calculateFFT();

    void getFFTfor(float* time, kiss_fft_cpx* freq);
    void getInverseFFTfor(kiss_fft_cpx* freq, float* timeData);
    
    Array<float> getCurrentFFT(int resolution);
    Array<float> getNoteResponse(int loMidiNote, int hiMidiNote, int noteSpan = 1);
    Array<float> buildInverseFFT();
    Array<float> getCurrentPhaseBins() {
        
        return phaseAngles;
        
    }
    
    Path getFullFFTPath();
    
    // Calculations ::::
    void calculateHistogram(AudioSampleBuffer* bufferToAnalyze, SpectralData& mySpectralAnalysis);
    
    void calculateKey(AudioSampleBuffer* bufferToAnalyze, SpectralData& mySpectralAnalysis);
    
    // AP methods :::
    
    bool acceptsMidi() const override                                           { return true; }
    bool producesMidi() const override                                          { return true; }
    
    double getTailLengthSeconds() const override                                { return 0.0; }
    
    //==============================================================================
    int getNumPrograms() override                                               { return 0; }
    int getCurrentProgram() override                                            { return 0; }
    void setCurrentProgram (int /*index*/) override                             {}
    const String getProgramName (int /*index*/) override                        { return String(); }
    void changeProgramName (int /*index*/, const String& /*name*/) override     {}
    
    //==============================================================================
    void getStateInformation (MemoryBlock&) override {}
    void setStateInformation (const void* data, int sizeInBytes) override {}
    
    
    //============================================================================== UI
    bool hasEditor() const override {
        return false;
    }
    
    AudioProcessorEditor* createEditor() override
    {
        return nullptr;
    }
    
    Component* createFFTGraph();
	Component* createHistogramGraph();

    
    const String getName() const override {
        
        return "FFT Processor";
    }

    
};


/// UI ::::::::::
class FFTGraphView : public Component,
                    public Timer
{
    
    ScopedPointer<Image> history;
    ScopedPointer<Graphics> historyGraphic;
    Path m_path;
    FFTProcessor* myFFTProcessor;
    
    bool refreshing;
    float animateCount;
    float animationLength = 100;
    
    ColourGradient spectrum;
    
    float maxPos = 0;
    float fadeMultiple = 0.8; //0.9;
    
public:
    
    FFTGraphView(FFTProcessor* FFTCalc);
    ~FFTGraphView();
    
    void setFadeMultiple(float fade);
    
private:
    
    void paint(Graphics& g) override;
    
    void update();
    
    void timerCallback() override;
    
    int getXforF(float freq);
    
    JUCE_LEAK_DETECTOR (FFTGraphView)
};
