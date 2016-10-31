/*
 *  FFTKiss.cpp
 *  Livetronica Studio
 *
 *  Created by Aaron Leese on 7/15/10.
 *  Copyright 2010 StageCraft Software. All rights reserved.
 *
 */


#include "FFTKiss.h"

/* :::::::::::::: RANGE OF FFT :::::::::: 
 
 When using the DFT it is necessary to know how to control the range and resolution of the spectrum. 
 The range is determined by the sampling time and is equal to ws/2. ( Nyquist F )
 The resolution is determined by the length of the sampling interval ~ 2p/NTs = ws/N .
 
*/


static double scaleLogPercentToFreqPercent(double logProp) {
    
    double fHz = 2*powf(10, 3*logProp + 1);
    double fProp = (fHz - 20)/22040;
    
    return fProp;
    
}


FFTProcessor::FFTProcessor() : base_a4(440.0) { 
    
    
	active = false;
    stopCalculations = false;
    
    notes.add("C ");
    notes.add("C#");
    notes.add("D ");
    notes.add("D#");
    notes.add("E ");
    notes.add("F ");
    notes.add("F#");
    notes.add("G ");
    notes.add("G#");
    notes.add("A ");
    notes.add("A#");
    notes.add("B ");
    notes.add("C ");
    
    currentNote = lastNote = 0;
    
    input = nullptr;
    output = nullptr;
    timedata = nullptr;
    
    currentSum = 0;
	lastBeat = 0;

    lastMillSecCount = 0;
    
    kiss_FFT_config = nullptr;
    kiss_iFFT_config = nullptr;
    kiss_FFT_real_config = nullptr;
    kiss_iFFT_real_config = nullptr;
    
    
    realOnlyFFT = true;
    
    // DEFAULT
	float res = kiss_fft_next_fast_size(16000);
    setResolution(res);
    
    
    // declare a bit array and set the bits appropriate to the scale selected ...
    // using equal temparment
    //
    //		 C    C#    D    D#   E    F     F#    G    G#   A    A#    B   C
    //		 1  16/15  9/8  6/5  5/4  4/3  45/32  3/2  8/5  5/3  9/5  15/8  2 
    
    scaleIntervals[0] = 1.0f;
    scaleIntervals[1] = 16.0f/15;
    scaleIntervals[2] = 9.0f/8;
    scaleIntervals[3] = 6.0f/5;
    scaleIntervals[4] = 5.0f/4;
    scaleIntervals[5] = 4.0f/3;
    scaleIntervals[6] = 45.0f/32;
    scaleIntervals[7] = 3.0f/2;
    scaleIntervals[8] = 8.0f/5;
    scaleIntervals[9] = 5.0f/3;
    scaleIntervals[10] = 9.0f/5;
    scaleIntervals[11] = 15.0f/8; 
    
    recordMarker = 0;
    lastBufWritePos = 0;
    currentNote = lastNote = 0;
    
    for (int i=0; i<12; i++)
        scale[i] = 1;
    
    
    // ADD main buses
    // TEMP hack - since Jules has set the bus creation to be linked to the preproccessor directive .. hmmm
    busArrangement.inputBuses.clear();
    busArrangement.inputBuses.add (AudioProcessorBus ("Input",  AudioChannelSet::canonicalChannelSet (2)));
    

}
FFTProcessor::~FFTProcessor() {
     
    if (input) kiss_fft_free(input);
    if (output) kiss_fft_free(output);
    if (timedata) kiss_fft_free(timedata);

    if (kiss_FFT_config != nullptr) free (kiss_FFT_config);
    if (kiss_iFFT_config != nullptr) free (kiss_iFFT_config);
    
    if (kiss_FFT_real_config != nullptr) free (kiss_FFT_real_config);
    if (kiss_iFFT_real_config != nullptr) free (kiss_iFFT_real_config);
    
}

void FFTProcessor::setResolution(int resolution, bool realDataOnly) {
    
    // LOCK the audio thread !
    const ScopedLock myScopedLock ( getCallbackLock() );
    
    // SET whether this FFT will include IMAGINERY data!
    realOnlyFFT = realDataOnly;
    
    resolution = kiss_fft_next_fast_size(resolution);
  
    // PROBLEM with KISS-FFT - above method doesn't always work ...
    if ( resolution%2 != 0) resolution++;

    // CHECK - already at this resolution ...
    if (FFT_Resolution == resolution) return;
    
    FFT_Resolution = resolution;
    
    incomingData.clear();
    incomingData.ensureStorageAllocated(FFT_Resolution);
    
    for (int i=0; i < FFT_Resolution; i++)
        incomingData.set(i, 0);
    
    recordMarker = 0;
    
    
    if (input != nullptr)
    {
        delete input;
        delete output;
        delete timedata;
    }
    
    // ALLOCATE the in/out for kiss FFT
    input = new kiss_fft_cpx[ FFT_Resolution ];
    output = new kiss_fft_cpx[ FFT_Resolution ];
    
    // for the simpler (real-only) FFT
    timedata = new kiss_fft_scalar[ FFT_Resolution ];
    
    // BUILD the HANNING WINDOW ::::
    hanningWindow.clear();
    for (int i=0; i < FFT_Resolution; i++)
    {
        float val = 0.5*(1.0f - cos(2.0f*(3.1415926)*(float)(i)/(float)(FFT_Resolution)));
        hanningWindow.add(val);
        input[i].i = 0;
    }
    
    
    // CLEAR the configs ....
    if (kiss_FFT_config != nullptr) free (kiss_FFT_config);
    if (kiss_iFFT_config != nullptr) free (kiss_iFFT_config);
    if (kiss_FFT_real_config != nullptr) free (kiss_FFT_real_config);
    if (kiss_iFFT_real_config != nullptr) free (kiss_iFFT_real_config);
    
    
    if (realOnlyFFT)
    {
        // For REAL ONLY input ....
        kiss_FFT_real_config = kiss_fftr_alloc(FFT_Resolution, 0, 0, 0);
        kiss_iFFT_real_config = kiss_fftr_alloc(FFT_Resolution, 1, 0, 0);
        kiss_FFT_config = nullptr;
        kiss_iFFT_config = nullptr;
    }
    else
    {
        // KISS FFT CONFIG
        kiss_FFT_config = kiss_fft_alloc(FFT_Resolution, 0, 0, 0);
        kiss_iFFT_config = kiss_fft_alloc(FFT_Resolution, 1, 0, 0);
        kiss_FFT_real_config = nullptr;
        kiss_iFFT_real_config = nullptr;

    }
    
    // ZERO the in/out buffers
    magnitudes.ensureStorageAllocated(FFT_Resolution);
	phaseAngles.ensureStorageAllocated(FFT_Resolution);
    
    FloatVectorOperations::clear(magnitudes.begin(), FFT_Resolution);
    FloatVectorOperations::clear(phaseAngles.begin(), FFT_Resolution);
    
    for (int i=0; i < FFT_Resolution; i++)
    {
		//magnitudes.add(0);
		//phaseAngles.add(0);
	 	timedata[i] = 0;
	}

    
}

// AP ::::::
void FFTProcessor::prepareToPlay (double sampleRate, int samplesPerBlock) {
    
    setPlayConfigDetails (2, 2, sampleRate, samplesPerBlock);
    
}
void FFTProcessor::releaseResources() { 
    
    // 
}

void FFTProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages) {
    
    // LOCK the audio thread !
    const ScopedLock myScopedLock ( getCallbackLock() );

    int numSamples = buffer.getNumSamples();
   
    // JUST UPDATE THE WINDOW ::::::
    // DONT do the calculation
   
    
    // IF the incoming audio is very quiet ....
    // then disable this calculation to save CPU
    if (buffer.getRMSLevel(0, 0, numSamples) + buffer.getRMSLevel(1, 0, numSamples) < 0.00001) active = false;
	else active = true;
    
    //
    lastMillSecCount = Time::getMillisecondCounter();
    
    // RECORD the data to the real part of the input ...
    for (int i=0; i < numSamples; i++)
    {
        float newData = *buffer.getReadPointer(0, i)/2 + *buffer.getReadPointer(1, i)/2;
        incomingData.set((recordMarker + i)%FFT_Resolution, newData);
    }
    
    //
    recordMarker = (recordMarker + numSamples)%incomingData.size();

}

// ACTUALLY CALCULATE the FFT
void FFTProcessor::calculateFFT() {
    
    // LOCK the audio thread !
    const ScopedLock myScopedLock ( getCallbackLock() );
    
    if (!timedata) return;
    
    if (realOnlyFFT)
    {
        // CONDITION :::: HANNING WINDOW
        // COPY the data to the real part of the input ...
        for (int i=0; i < incomingData.size(); i++)
        {
            // since the Hanning window has 0s at it's min and max ...
            // we don't need to worry about the discontinuity at incomingData[recordMarker]
            timedata[i] = incomingData[(recordMarker + i)%FFT_Resolution]*hanningWindow[i];
            
        }
        
        // FFT...
        kiss_fftr( kiss_FFT_real_config, timedata, output);
    }
    else // not sure if I ever need the full fft (when would input have imaginery portion?)
    {
        // CONDITION :::: HANNING WINDOW
        // COPY the data to the real part of the input ...
        for (int i=0; i < incomingData.size(); i++)
        {
            // since the Hanning window has 0s at it's min and max ...
            // we don't need to worry about the discontinuity at incomingData[recordMarker]
            input[i].r = incomingData[(recordMarker + i)%FFT_Resolution]*hanningWindow[i];
        }
        
        // FFT...
        kiss_fft( kiss_FFT_config, input, output);
    }
}
void FFTProcessor::getFFTfor(float* time, kiss_fft_cpx* freq) {
    
    // LOCK the audio thread !
    const ScopedLock myScopedLock ( getCallbackLock() );
    
    if (realOnlyFFT)
    {
        // CONDITION :::: HANNING WINDOW
        // COPY the data to the real part of the input ...
        for (int i=0; i < getResolution(); i++)
        {
            jassert(getResolution() == hanningWindow.size());
            
            // since the Hanning window has 0s at it's min and max ...
            // we don't need to worry about the discontinuity at incomingData[recordMarker]
            
            timedata[i] = time[i]*hanningWindow[i];
            jassert(timedata[i]  == timedata[i]);

        }
        
        // FFT...
        kiss_fftr( kiss_FFT_real_config, time, freq);
        jassert(freq->r == freq->r);
        
    }
    else // FULL FFT
    {
        // CONDITION :::: HANNING WINDOW
        // COPY the data to the real part of the input ...
        for (int i=0; i < getResolution(); i++)
        {
            // since the Hanning window has 0s at it's min and max ...
            // we don't need to worry about the discontinuity at incomingData[recordMarker]
            input[i].r = time[i]*hanningWindow[i];
            input[i].i = 0;
        }
        
        // FFT...
        kiss_fft( kiss_FFT_config, input, freq);
    }
    
}
void FFTProcessor::getInverseFFTfor(kiss_fft_cpx* freq, float* timeData) {
    
    // LOCK the audio thread !
    const ScopedLock myScopedLock ( getCallbackLock() );
    
    if (realOnlyFFT)
    {
        // FFT...
        kiss_fftri( kiss_iFFT_real_config, freq, timeData);
    }
    else
    {
        
        // DOES NOT EXIST ....
        jassertfalse;
        
    }

}


// GET FFT ::::
Array<float> FFTProcessor::getCurrentFFT(int resolution) {

    // LOCK the audio thread !
    const ScopedLock myScopedLock ( getCallbackLock() );
    
    if (!output)
    {
        // does this happen ???
        DBG("GAHHHHH!!!!!");
        jassertfalse;
        Array<float> empty;
        return empty;
    }
    
    // Do the math ....
    calculateFFT();
    
    // COPY and SCALE the data to the desired number of bins
	// the FFT is symetrical ... and we only use the first half ....
    float posLast = 0;
	 
    for (int i=0; i < resolution; i++)
    {

		float percentage = float(i)/float(resolution);
		float position = percentage*float(FFT_Resolution)*(.5); // .5 here because it comes out as symetrical image
        
        float intPosition = 0;
        float fraction = modff(position, &intPosition);
        
        if (magnitudes.size() <= position) magnitudes.add(0);
    
        if (intPosition + 1 == FFT_Resolution)
            continue;
        
        // LINEAR INTERPOLATION :::::
        float lowValue = sqrt(powf(output[int(intPosition)].r,2) + powf(output[int(intPosition)].i,2));
        float hiValue = sqrt(powf(output[int(intPosition+1)].r,2) + powf(output[int(intPosition+1)].i,2));
        
        float delta = hiValue - lowValue;
        float exactValue = lowValue + fraction*delta;
        
        // SET :::
        magnitudes.set(i, exactValue);
        posLast = position;
    }
 
    // update the max ....
    peakIndex = FloatVectorOperations::findMaximum(magnitudes.getRawDataPointer(), magnitudes.size() );

	return magnitudes;

}
kiss_fft_cpx* FFTProcessor::getCurrentFullFFT(int numSamples) {
    
    // Do the math ....
    if (realOnlyFFT)
    {
        float resolution = getResolution();
        
        // SHORT RECT window .... seems ok
        if (false)
        {
            // CONDITION :::: HANNING WINDOW
            // COPY the data to the real part of the input ...
            for (int i=0; i < resolution; i++)
            {
                // ZERO padding ...
                if (i >= numSamples)
                {
                    timedata[i] = 0;
                    continue;
                }
                
                // since the Hanning window has 0s at it's min and max ...
                // we don't need to worry about the discontinuity at incomingData[recordMarker]
                timedata[i] = incomingData[(recordMarker + i)%FFT_Resolution];
            }
            
        }
        // HANNING :::: 2 blocks
        else if (true)
        {
            // CONDITION :::: HANNING WINDOW
            // COPY the data to the real part of the input ...
            for (int i=0; i < resolution; i++)
            {
                // WINDOW at 2*blocksize ...
                float max = 2*numSamples;
                
                // ZERO padding ...
                if (i >= max)
                {
                    timedata[i] = 0;
                    continue;
                }
                
                float perc = float(i)/max;
                
                // since the Hanning window has 0s at it's min and max ...
                // we don't need to worry about the discontinuity at incomingData[recordMarker]
                timedata[i] = incomingData[(recordMarker + i)%FFT_Resolution]*hanningWindow[perc*resolution];
                
            }
        }
        
        
        
        
        // FFT...
        kiss_fftr( kiss_FFT_real_config, timedata, output);
    }
    else // not sure if I ever need the full fft (when would input have imaginery portion?)
    {
        // CONDITION :::: HANNING WINDOW
        // COPY the data to the real part of the input ...
        for (int i=0; i < incomingData.size(); i++)
        {
            // since the Hanning window has 0s at it's min and max ...
            // we don't need to worry about the discontinuity at incomingData[recordMarker]
            input[i].r = incomingData[(recordMarker + i)%FFT_Resolution]*hanningWindow[i];
        }
        
        // FFT...
        kiss_fft( kiss_FFT_config, input, output);
    }
    
    
    return output;
    
}
Array<float> FFTProcessor::buildInverseFFT () {
    
    // Assume we have already gotten the FFT
    
    if (realOnlyFFT)
    {   
        //void kiss_fftri(kiss_fftr_cfg cfg,const kiss_fft_cpx *freqdata,kiss_fft_scalar *timedata);
        kiss_fftri( kiss_FFT_real_config, output, magnitudes.getRawDataPointer());
    }
    else
    {
        //kiss_fftri( kiss_FFT_real_config, output, magnitudes.getRawDataPointer());
    }
    
    return magnitudes;
    
}

Path FFTProcessor::getFullFFTPath() {
    
    
    // LOCK the audio thread !
    //const ScopedLock myScopedLock ( getCallbackLock() );
    
    Path fftPath;
    
    // Do the math ....
    calculateFFT();
    
    for (int i=0; i < FFT_Resolution; i++)
    {
        float xProp = float(i)/float(FFT_Resolution);
       
        // displaying logarithmically ....
        float f = 22.03*powf(10, xProp*3); // here f will go from 20Hz to 20K
        
        float val = sqrt(powf(output[int(i)].r, 2) + powf(output[int(i)].i, 2));
        
        if (i == 0) fftPath.startNewSubPath(f, 10000*val);
        else fftPath.lineTo(f, 10000*val);
        
      
    }
    

    return fftPath;
}


/// UI ::::::::::
FFTGraphView::FFTGraphView(FFTProcessor* FFTCalc) {
    
    myFFTProcessor = FFTCalc;
    
    maxPos = 0;
    
    setSize(300, 300);
    
    history = new Image(Image::ARGB, 500, 300, true);
    historyGraphic = new Graphics(*history);
    
    Colour base = Colours::purple.withSaturation(.7).withBrightness(.4);
    
    spectrum = ColourGradient(base.withRotatedHue(0).withAlpha(0.0f), 0, 0,
                              base.withRotatedHue(.1), history->getWidth(), 0, false);
    
    spectrum.addColour(.2, base.withRotatedHue(0));
    spectrum.addColour(.4, base.withRotatedHue(.7));
    spectrum.addColour(.6, base.withRotatedHue(.5));
    spectrum.addColour(.75, base.withRotatedHue(.3));
    
    historyGraphic->setGradientFill(spectrum);
    
    animationLength = 100;
    animateCount = 0;
    startTimer(1000/10);
    
    const MessageManagerLock mmLock;
    
    // basically and overlay - no interaction
    setWantsKeyboardFocus(false);
    setInterceptsMouseClicks(false, false);
    removeMouseListener(this);
    
    
}
FFTGraphView::~FFTGraphView() {
    
    DBG("FFT GraphView delete ");
    stopTimer();
    
}

void FFTGraphView::setFadeMultiple(float fade) {
    
    fadeMultiple = fade;
    
}
void FFTGraphView::paint(Graphics& g) {
    
    // FADE out .... when the level falls too low ...
    // we animate 100 frames, then stop animating ... to save resources
    
    history->multiplyAllAlphas(fadeMultiple);
    
    if (m_path.getBounds().getHeight() != m_path.getBounds().getHeight()) return;
    if (m_path.getBounds().getWidth() != m_path.getBounds().getWidth()) return;
    
    
    float opacity = animateCount/animationLength; // fade in
    if (opacity > 0)
    {
        historyGraphic->setOpacity(opacity/8);
        historyGraphic->fillPath (m_path);
        
        historyGraphic->setOpacity(opacity);
        historyGraphic->strokePath (m_path, PathStrokeType(1.2));
    }
    
    if (refreshing == false)
    {
        history->desaturate();
    }
    
    int xDelta = jmax(2, getWidth()/200);
    int yDelta = jmax(2, getHeight()/100);
    
    history->moveImageSection(xDelta, -1.5 -yDelta, 0, 0, history->getWidth(), history->getHeight());
    
    // draw the FFT in the bottom half of the graph ....
    g.drawImage(*history, 2, 2, getWidth()-4, getHeight()-4, 0, 0, history->getWidth(), history->getHeight());
    
    
}
void FFTGraphView::update () {
    
    float w = history->getWidth();
    float h = history->getHeight();
    
    if (h < 20 || w < 20) return; // should never happen
    
    float dbCutoff = 20;
    
    m_path.clear();
    m_path.startNewSubPath (0, h*dbCutoff);
    
    float FFT_resolution = myFFTProcessor->getResolution();
    float view_resolution = history->getWidth();
    
    if ( myFFTProcessor->active )
    {
        Array<float> FFT = myFFTProcessor->getCurrentFFT(FFT_resolution);
        
        int lastIndex = 0;
        
        for (int xi = 0; xi < view_resolution; ++xi )
        {
            float xProp = float(xi)/float(view_resolution); // [0..1)
            
            // displaying logarithmically ....
            float f = scaleLogPercentToFreqPercent(xProp);
            
            int currentIndex = f*FFT_resolution;
            if (currentIndex == lastIndex) continue; // so we only draw the line when a new data point has been reached
            
            // get MAX value in this bin
            float y = 0;
            int span = int(currentIndex) - int(lastIndex);
            if (span > 1) y = FloatVectorOperations::findMaximum(&FFT.begin()[lastIndex], span);
            else y = FFT[int(currentIndex)];
            
            //
            y = Decibels::gainToDecibels(y);
            
            y = jlimit<float>(-dbCutoff, 120, y); // 100 is the positive (dB) cutoff ...
            
            lastIndex = currentIndex; // memory ...
            
            // DC values, essentially ...
            // which we still get (in the synth for example), but we don't really want to see ...
            // if (f < 0.002) y = 0;
            
            // NaN
            if (y != y) return;
            
            m_path.lineTo (history->getWidth()*xProp, -h*y); // 200 FUDGE FACTOR !!!!! I dont get the Y !?
            
        }
        
        
        m_path.lineTo (view_resolution, h*dbCutoff);
        
        
    }
    
    
    if (m_path.getBounds().getHeight() > .00001)
    {
        // RESCALE ::::: to fill the space ...
        m_path.scaleToFit(0, 0, w, h/2, false);
        m_path.applyTransform(AffineTransform::translation(0, h/2 - 2)); // Translate - down ...
    }
    
    repaint();
}
void FFTGraphView::timerCallback() {
    
    if ( ! isShowing()) return;
    
    uint32 current = Time::getMillisecondCounter();
    int elapsed = current - myFFTProcessor->lastMillSecCount;
    
    if (elapsed > 100) // FADE OUT
    {
        refreshing = false; // clears ...
        
        // FADE ...
        if (animateCount > 0)
        {
            animateCount -= 5;
            animateCount = jlimit(0.0f, animationLength, animateCount);
        }
    }
    else // FADE IN
    {
        refreshing = true; // clears ...
        if (animateCount < animationLength) animateCount += 10;
        
        animateCount = jlimit(0.0f, animationLength, animateCount);
    }
    
    update();
    
}
int FFTGraphView::getXforF(float freq) {
    
    float x = (powf(10, (freq/myFFTProcessor->getResolution())*3 + 1) - 10)/10000;
    return x*getWidth();
    
}


class FFTHistogram : public Component,
                    public Timer
{
    
    ScopedPointer<Image> history;
    ScopedPointer<Graphics> historyGraphic;
    Path m_path;
    FFTProcessor* myFFTProcessor;
    
    bool refreshing;
    float animateCount;
    float animationLength;
    
    ColourGradient spectrum;
    
public:
    
    FFTHistogram(FFTProcessor* FFTCalc) {
        
        myFFTProcessor = FFTCalc;
        
        
        // basically and overlay - no interaction
        setWantsKeyboardFocus(false);
        setInterceptsMouseClicks(false, false);
        removeMouseListener(this);
        
        setSize(300, 300);
        
        history = new Image(Image::ARGB, 500, 300, true);
        historyGraphic = new Graphics(*history);
        
        Colour base = Colours::purple.withSaturation(.7).withBrightness(.4);
        
        spectrum = ColourGradient(base.withRotatedHue(0).withAlpha(0.0f), 0, 0,
                                  base.withRotatedHue(.1), history->getWidth(), 0, false);
        
        spectrum.addColour(.2, base.withRotatedHue(0));
        spectrum.addColour(.4, base.withRotatedHue(.7));
        spectrum.addColour(.6, base.withRotatedHue(.5));
        spectrum.addColour(.75, base.withRotatedHue(.3));
        
        historyGraphic->setGradientFill(spectrum);
        
        animationLength = 100;
        animateCount = 0;
        
        
        startTimer(1000/24);
        
    }
    ~FFTHistogram() {
        
        DBG("FFT GraphView delete ");
        stopTimer();
        
    }
    
private:
    
    void paint(Graphics& g) override {
        
        g.fillAll(Colours::black);
        
        // FADE out .... when the level falls too low ...
        // we animate 100 frames, then stop animating ... to save resources
        
        
        int xDelta = 0; //jmax(2, getWidth()/70);
        int yDelta = 1; //jmax(3, getHeight()/50);
        
        history->moveImageSection(xDelta, -yDelta, 0, 0, history->getWidth(), history->getHeight());
        
        // draw the FFT in the bottom half of the graph ....
        g.drawImage(*history, 0, 0, getWidth(), getHeight(), 0, 0, history->getWidth(), history->getHeight());
        
    }
    
    void update () {
        
        
        float resolution = jmin(history->getWidth(), getWidth());
        
        if (true) // myFFTProcessor->active )
        {
            Array<float> FFT = myFFTProcessor->getCurrentFFT(resolution);
            
            int lastIndex = 0;
            float x = 0;
            
            juce::Rectangle<int> r(0, history->getHeight() - 1, history->getWidth(), 1);
            history->clear(r);
            
            for (int xi = 0; xi < resolution; ++xi )
            {
                float percentage = float(xi)/float(resolution); // [0..1)
                
                // displaying logarithmically ....
                float f = 22.03*powf(10, percentage*3); // here f will go from 20Hz to 20K
                f = (f)/22050; // scale 0->1
                
                int currentIndex = f*resolution;
                if (currentIndex == lastIndex) continue; // so we only draw the line when a new data point has been reached
                
                float y = FFT[int(currentIndex)];
                
                float saturation = FFT[int(currentIndex)]/30;
                float brightness = FFT[int(currentIndex)]/30;
                
                lastIndex = currentIndex; // memory ...
                
                // NaN
                if (y != y) return;
                
                Colour c = Colours::purple.withRotatedHue(float(xi)/float(resolution));
                c = c.withSaturation(saturation).withBrightness(brightness);
                
                float nextX = percentage*history->getWidth();
                
                while(x <= nextX)
                {
                    // ADDITIVE :::::
                    Colour current = history->getPixelAt(x, history->getHeight());
                    
                    float red = jmax(current.getRed(), c.getRed());
                    float green = jmax(current.getGreen(), c.getGreen());
                    float blue = jmax(current.getBlue(), c.getBlue());
                    
                    Colour newC(red, green, blue);
                    
                    history->setPixelAt(x, history->getHeight() - 1, newC);
                    x++;
                }
            }
            
            
            // SUM
            float sum = myFFTProcessor->getCurrentSum(); 
            
            for (int x = history->getWidth() - sum/5; x < history->getWidth(); x++)
                history->setPixelAt(x, history->getHeight() - 1, Colours::white);
            
            // TEST :::
            // if the sum is above a threshold ....
            // call this a "hit"
            // then look for the average time between hits ....
            // and guess at the BPM from that ....
            //uint32 current = Time::getMillisecondCounter();
            
        } 
        
        repaint();
    }
    
    void timerCallback() override {
        
        if ( ! isShowing()) return;
        
        uint32 current = Time::getMillisecondCounter();
        int elapsed = current - myFFTProcessor->lastMillSecCount;
        
        if (elapsed > 100) // FADE OUT
        {
            refreshing = false; // clears ...
            
            // FADE ...
            if (animateCount > 0)
            {
                animateCount -= 5;
                animateCount = jlimit(0.0f, animationLength, animateCount);
            }
            
            update();
            
        }
        else // FADE IN
        {
            refreshing = true; // clears ...
            if (animateCount < animationLength) animateCount += 10;
            
            animateCount = jlimit(0.0f, animationLength, animateCount);
            update();
        }
        
        repaint();
        
    }
    
    int getXforF(float freq) {
        
        float x = (powf(10, (freq/myFFTProcessor->getResolution())*3 + 1) - 10)/10000;
        return x*getWidth();
        
    }
    
    JUCE_LEAK_DETECTOR (FFTHistogram)
};


Component* FFTProcessor::createFFTGraph() {
    
    return new FFTGraphView(this);
    
}

Component* FFTProcessor::createHistogramGraph() {
    
    return new FFTHistogram(this);
    
}








