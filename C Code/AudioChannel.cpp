/*
 Copyright 2009 Music and Entertainment Technology Laboratory - Drexel University
 
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 
 http://www.apache.org/licenses/LICENSE-2.0
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

#include <stdlib.h>;
#include <stdio.h>;
#include <cstdio>;
#include <string.h>;
#include <math.h>;
#include <sys/time.h>;
#include "AudioChannel.h";
#include "CircularBuffer.h";
#include "AS3.h";


char out[200];
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging

//int minSampleSize 2048;
//int maxSampleSize 4096;

/*
 Class: AudioChannel.cpp
 
 This class manages the buffers and other data needed in audio storage, playback, and processing. There are circular buffers for the
 input and output audio streams. Data that is reused from frame to frame (e.g. filter coefficients) is stored in this class and is
 accessible in <ALFPackage.cpp> for use with the DSP functions contained in <DSPFuncs.c>.
*/	

/*
 Group: Constructor

 Constructor: AudioChannel
 
 Creates an AudioChannel object with the default parameters.
 
 Default parameters:
 
 * hopSize - the size, in samples, of the 'inAudioFrame' buffer. Default: 4096.
 * fs - the sampling rate of the audio under analysis. Set by 'initChannel' method. Default: 44100.
 * fftFlag - Boolean indicating whether or not an FFT has been performed on the current frame of audio. Default value: false.
 * magFlag - Boolean indicating if the magnitude spectrum has been computed on the fftFrame. Default: false.
 * centroidFlag - Boolean indicating if the spectral centroid has been performed on the current frame of data. Default: false.
 * centroidVal - Stores the value of the spectral centroid for feature computatons. Default: 0.
 * filterLen - WILL BECOME OBSOLETE. NEED TO REMOVE FROM AUDIOCHANNEL.CPP & .H AND ALFPACKAGE.CPP
 * filterProcessing - Boolean value indicating if processing functionality (aka filter) is currently in use. Default: 0;
 */
AudioChannel::AudioChannel(){
	
	hopSize = 4096;
	fs = 44100;
	fftFlag = false;
	magFlag = false;
	centroidFlag = false;
	centroidVal = 0;
	filterLen = 0;  //make obsolete
	filterProcessing = 0;
	
}

AudioChannel::~AudioChannel() {

	sprintf(out, "destroying AudioChannel \n"); DISP(out);
	free(inAudioFrame);
	free(outAudioFrame);
	free(fftFrame);
	free(freqData);
	free(magSpectrum);
	free(spectrumPrev);
	free(hannCoefficients);
	free(filter);
	
	delete [] roomSize;
	delete [] sourcePosition;
	delete [] micPosition;
	delete [] inBuffer;
	delete [] outBuffer;
}


/*
 Group: Utilities
 
 Function: clearFFTFrame()
 
 Clears the fftFrame buffer, which has the FFT dta.
 
 */
void AudioChannel::clearFFTFrame(){
	
	//not sure why memset isn't working here, so we use a loop
	int i;
	for(i = 0; i < fftSize; i++){ fftFrame[i] = 0;}
}
/*
 Function: clearInAudioFrame()
 
 Clears the inAudioFrame buffer that is the only buffer that can be read and written to by C++ and Actionscript. This function
 should be called when not writing a full frame as samples from the previous frame may remain and introduce error into 
 calculations that are performed.
 
*/
void AudioChannel::clearInAudioFrame(){
	memset(inAudioFrame, 0, hopSize*sizeof(float));
}
/* 
 Function: initChannel()
 
 This function is called with the allocation of every DATF (or ALF) object from ActionScript for a desired audio channel
 (left or right). This function initializes the circular buffers needed for input and output. It also initializes memory
 required for the feature analysis functions. Certain flags are initialized as well. See Section <Buffers>
 
 Parameters:
 
	hopSize - The desired size of the analysis frame, as selected by the user from ALF or DATF.
	fftSize - the size of the FFT to be computed on each frame when using spectral analysis fucntions.
	fs - The audio sampling rate of the user-selected audio
 
 */
void AudioChannel::initChannel(int _hopSize, int _fftSize, int _fs) {

	fftSize = _fftSize;
	hopSize = _hopSize;
	fs = _fs;
	if(hopSize < 2048){ audioOutputSize = 2*hopSize;} // = 2048 b/c this is the least # of samples flash can play at once
	else{ audioOutputSize = hopSize; }	
	inAudioFrame = (float*) malloc(sizeof(float)*(fftSize + 2));		// We initialize this to be 2 larger for the weird fftPacking thing
	outAudioFrame = (float*) malloc(sizeof(float)*4096);				// Maxsamplesize is the max number of samples flash can playback
	fftFrame = (float*) malloc(sizeof(float)*(fftSize + 2));
	freqData = (float*) malloc(sizeof(float)*(fftSize/2 + 1));			// Number of unique frequencies in the specific fftLength
	magSpectrum = (float*) malloc(sizeof(float)*(fftSize/2 + 1));		// Number of unique magnitudes ""			""
	spectrumPrev = (float *) calloc((fftSize/2 + 1),sizeof(float));
	hannCoefficients = (float *) malloc(sizeof(float)*fftSize*2);
	
	// These are the pointers that ActionScript will see. We use these instead of inAudioFrame and outAudioFrame
	// since a call to reInitAudioChannel may change inAudioFrame and outAudioFrame but we then just copy the 
	// value to inRefPtr and outRefPtr to prevent giving Actionscript new pointers to keep track of.
	//*inRefPtr = inAudioFrame;
	//outRefPtr = outAudioFrame;
		
	roomSize = new float[3];
	sourcePosition = new float[3];
	micPosition = new float[3];
	
	// Allocates our input circular buffer 
	int buffSize = 30000;
	inBuffer = new CircularBuffer[1];
	inBuffer->initBuffer(buffSize);
	
	// Allocates our output circular buffer 
	buffSize = 30000;
	outBuffer = new CircularBuffer[1];
	outBuffer->initBuffer(buffSize);
	
//	sprintf(out, "===== AudioChannel::initChannel ====="); DISP(out);
//	sprintf(out, "inBuffer: %i", inBuffer); DISP(out);
//	sprintf(out, "fftSize = %i hopSize = %i fs = %i audioOutputSize = %i", fftSize, hopSize, fs, audioOutputSize); DISP(out);
	
	inAudioFrameSamples = 0;
	outAudioFrameSamples = 0;
	circularBufferFlag = false;
	bufferReady = false;
	firstFrame = true;
	stereo = false;
}
/* 
 Function: reInitChannel()
 
 This function is called with the allocation of every DATF (or ALF) object from ActionScript for a desired audio channel
 (left or right). This function initializes the circular buffers needed for input and output. It also initializes memory
 required for the feature analysis functions. Certain flags are initialized as well. See Section <Buffers>
 
 Parameters:
 
	hopSize - The desired size of the analysis frame, as selected by the user from ALF or DATF.
	fftSize - the size of the FFT to be computed on each frame when using spectral analysis fucntions.
	fs - The audio sampling rate of the user-selected audio
 
 */
void AudioChannel::reInitChannel(int _hopSize, int _fftSize, int _fs, int numCh) {

	sprintf(out, "reInitializing Channel"); DISP(out);

	fftSize = _fftSize;
	hopSize = _hopSize;
	fs = _fs;
	
//	sprintf(out, "in reInitChannel: fs = %i, hopSize = %i, fftSize = %i", fs, hopSize, fftSize); DISP(out);	
//	sprintf(out, "inAudioFramePtr: %i\n fftFramePtr: %i\n freqDataPtr: %i\n magSpectrumPtr: %i\n spectrumPrevPtr: %i\n hannCoefPtr: %i\n",
//			inAudioFrame, fftFrame, freqData, magSpectrum, spectrumPrev, hannCoefficients); DISP(out);
	if(hopSize < 2048){ audioOutputSize = 2*hopSize;} // = 2048 b/c this is the least # of samples flash can play at once
	else{ audioOutputSize = hopSize; }	
	//inAudioFrame = (float*) realloc(sizeof(float)*(fftSize + 2));		// We initialize this to be 2 larger for the weird fftPacking thing
	free(inAudioFrame);
//	inAudioFrame = (float *)realloc(inAudioFrame, sizeof(float)*(fftSize + 2));
	inAudioFrame = (float *) calloc(2*hopSize, sizeof(float));
	fftFrame = (float *)realloc(fftFrame, sizeof(float)*(fftSize + 2));
	freqData = (float *)realloc(freqData, sizeof(float)*(fftSize/2 + 1));
	magSpectrum = (float *)realloc(magSpectrum, sizeof(float)*(fftSize/2 + 1))	;
	spectrumPrev = (float *)realloc(spectrumPrev, (fftSize/2 + 1)*sizeof(float));	
	hannCoefficients = (float *)realloc(hannCoefficients, sizeof(float)*fftSize*2);	
	
//	sprintf(out, "inAudioFramePtr: %i\n fftFramePtr: %i\n freqDataPtr: %i\n magSpectrumPtr: %i\n spectrumPrevPtr: %i\n hannCoefPtr: %i\n",
//			&inAudioFrame[0], fftFrame, freqData, magSpectrum, spectrumPrev, hannCoefficients); DISP(out);
	
	int i;
	for(i = 0; i < (fftSize/2 + 1); i++){
		spectrumPrev[i] = 0;
	}
	
	if(numCh == 1){stereo = false;}
	else if(numCh == 2){stereo = true;}
	else{sprintf(out, "Only Mono and Stereo files are supported!"); DISP(out);}

//for(i = 0; i < 10; i++){
//	sprintf(out, "inAudioFrame[%i] = %f", i, inAudioFrame[i]); DISP(out);
//}
//for(i = hopSize*2 - 10; i < hopSize*2; i++){
//	sprintf(out, "inAudioFrame[%i] = %f", i, inAudioFrame[i]); DISP(out);
//}
//	memset(spectrumPrev, 0, (fftSize/2 + 1));
//	fftFrame = (float*) realloc(sizeof(float)*(fftSize + 2));
//	freqData = (float*) realloc(sizeof(float)*(fftSize/2 + 1));			// Number of unique frequencies in the specific fftLength
//	magSpectrum = (float*) realloc(sizeof(float)*(fftSize/2 + 1));		// Number of unique magnitudes ""			""
//	spectrumPrev = (float *) calloc((fftSize/2 + 1),sizeof(float));
//	hannCoefficients = (float *) remalloc(sizeof(float)*fftSize*2);
	
	clearInAudioFrame();
	// inRefPtr = inAudioFrame;
}
/*
 Function: getCentroidVal()
 
 Returns the centroid of the audio spectrum from the current frame of data
 
 Returns:
 
 * centroidVal - the current frame's spectral centroid
 */
float AudioChannel::getCentroidVal() {return centroidVal;}
/*
 Function: getChannelName()
 
 Retrieves the name given to the AudioChannel object.
 
 Returns:
 
 * a pointer for the character array
 */
char *AudioChannel::getChannelName() {return channelName;}
/*
 Function: getFFTSize()
 
 Returns the size of the FFT used on the current audio frame
 
 Returns:
 
 * fftSize: the fftSize used for the data
 */
int AudioChannel::getFFTSize() { return  fftSize;}
/*
 Function: getFrameSize()
 
 Returns the user-specificed hopSize object
 
 Returns:
 
 * hopSize: the user's desired hopSize
 */
int AudioChannel::getHopSize() { return  hopSize;}
/*
 Function: getSampleRate()
 
 Returns the sampling rate of the audio under analysis
 
 Returns:
 
 * fs - the audio sampling rate
 */
int AudioChannel::getSampleRate() { return fs;}
/*
 Function: getOutAudioFrameSamples()
 
 Returns the current number of samples contained in outAudioFrame
 
 Returns:
 
 * outAudioFrameSamples - the current number of samples contained in outAudioFrame
 */
int AudioChannel::getOutAudioFrameSamples(){return outAudioFrameSamples;}
/*
	Funciton: resetChannel
	
	Parameters:
*/
void AudioChannel::resetChannel(){

	sprintf(out, "resetChannel in AudioChannel"); DISP(out);
	inBuffer->resetBuffer();
	outBuffer->resetBuffer();
	clearInAudioFrame();
	clearFFTFrame();
	clearFlags();
	outAudioFrameSamples = 0;
	firstFrame = true;
	stereo = false;
	bufferReady = false;
}
/*
 Function: setCentroidVal()
 
 Sets the value of the spectral centroid for the current frame
 
 Parameters:
 
 *val: a float specifying the centroid value of the current frame
 */
void AudioChannel::setCentroidVal(float val) {centroidVal = val;}
/*
 Function: setChannelName()
 
 Sets the name of the AudioChannel
 
 Parameters:
 
 *name: a pointer for the character array containing the desired name of the channel
 */
void AudioChannel::setChannelName(char name[]) {
	int i;
	for(i = 0; i < 4; i++){ channelName[i] = name[i];}
}
/*
 Function: setHopSize()
 
 Sets the value of the hopSize
 
 Parameters:
 
 *hop: An int which is the new hopSize
 */
void AudioChannel::setHopSize(int hop){
	hopSize = hop;
}
/*
 Function: setSampleRate()
 
 Sets the sampling rate of the audio under analysis
 
 Returns:
 
 * _fs - the audio sampling rate
 */
int AudioChannel::setSampleRate(int _fs) { fs = _fs;}
/*
 Function: setAudioFrameSamples()
 
 Sets the number of samples currently contained inthe InAudioFrame buffer
 
 Parameters:
 
 * numSamples: an int specifying the number of samples in the buffer
 
 */
void AudioChannel::setAudioFrameSamples(int numSamples) {inAudioFrameSamples = numSamples;}
/*
 Function: setOutAudioFrameSamples()
 
 Sets the number of samples currently contained inthe OutAudioFrame buffer
 
 Parameters:
 
 * numSamples: an int specifying the number of samples in the buffer
 
 */
void AudioChannel::setOutAudioFrameSamples(int numSamples) {outAudioFrameSamples = numSamples;}




void AudioChannel::setRoom(	float newRoomLength, float newRoomWidth, float newRoomHeight, 
							float newSourceLength, float newSourceWidth, float newSourceHeight, 
							float newMicLength, float newMicWidth, float newMicHeight, double newEchoStrength){
	//sprintf(out, "\nSet Room"); DISP(out);
	//sprintf(out, "sourceWidth = %4.15f, newSourceWidth = %4.15f", sourceWidth, newSourceWidth); DISP(out); 
	roomLength = newRoomLength;
	roomWidth = newRoomWidth;
	roomHeight = newRoomHeight;
	sourceLength = newSourceLength;
	sourceWidth = newSourceWidth;
	sourceHeight = newSourceHeight; 
	micLength = newMicLength;
	micWidth = newMicWidth;
	micHeight = newMicHeight;
	echoStrength = newEchoStrength;
	//sprintf(out, "sourceWidth = %4.15f, newSourceWidth = %4.15f", sourceWidth, newSourceWidth); DISP(out);
}
bool AudioChannel::checkRoom(	float newRoomLength, float newRoomWidth, float newRoomHeight, 
								float newSourceLength, float newSourceWidth, float newSourceHeight, 
								float newMicLength, float newMicWidth, float newMicHeight, double newEchoStrength){
	// See if parameters are the same
	//sprintf(out, "\nCheck Room"); DISP(out);
	//sprintf(out, "sourceWidth = %4.15f, newSourceWidth = %4.15f", sourceWidth, newSourceWidth); DISP(out);		
	if(	(roomLength == newRoomLength) &&
		(roomWidth == newRoomWidth) &&
		(roomHeight == newRoomHeight) &&
		(sourceLength == newSourceLength) &&
		(sourceWidth == newSourceWidth) &&
		(sourceHeight == newSourceHeight) &&
		(micLength == newMicLength) &&
		(micWidth == newMicWidth) &&
		(micHeight == newMicHeight) &&
		(echoStrength == newEchoStrength) )
	{
		newRoom = false;
		//sprintf(out, "false"); DISP(out);		
	} else{
		newRoom = true;
		//sprintf(out, "true"); DISP(out);		
	}
	//sprintf(out, "still in check room sourceWidth = %4.15f, newSourceWidth = %4.15f", sourceWidth, newSourceWidth); DISP(out);
//	if(roomLength == newRoomLength){sprintf(out, "roomLengths are equal"); DISP(out);	}
//	if(roomWidth == newRoomWidth){sprintf(out, "roomWidths are equal"); DISP(out);	}
//	if(roomHeight == newRoomHeight){sprintf(out, "roomHeight are equal"); DISP(out);	}
//	if(sourceLength == newSourceLength){sprintf(out, "sourceLengths are equal"); DISP(out);	}
//	if(sourceWidth == newSourceWidth){sprintf(out, "sourceWidths are equal"); DISP(out);	}
//	if(sourceHeight == newSourceHeight){sprintf(out, "sourceHeights are equal"); DISP(out);	}
//	if(micLength == newMicLength){sprintf(out, "micLengths are equal"); DISP(out);	}
//	if(micWidth == newMicWidth){sprintf(out, "micWidths are equal"); DISP(out);	}
//	if(micHeight == newMicHeight){sprintf(out, "micHeights are equal"); DISP(out);	}
//	if(echoStrength == newEchoStrength){sprintf(out, "echoStrengts are equal"); DISP(out);	}	
	//sprintf(out, "still in check room sourceWidth = %4.15f, newSourceWidth = %4.15f", sourceWidth, newSourceWidth); DISP(out);	
//	sprintf(out, "bools = %i %i %i", sourceLength == newSourceLength, sourceWidth == newSourceWidth, sourceHeight == newSourceHeight); DISP(out);
	return newRoom;
}
/*
	Group: Flag Management
	
	Function: clearFlags()
	
	This function clears all of the flags that are used to specify whether a given value
	has been calculated on the current frame. See the <Flags> section at the bottom of this page.

*/
void AudioChannel::clearFlags() {
	fftFlag = false;
	centroidFlag = false;
	magFlag = false;

}
/*
 Function: getCentroidFlag()
 
 Returns the status off the centroidFlag to determine if the centroid has been calculated
 
 Returns:
 
 * centroidFlag - the centroid flag for the current frame
 */
bool AudioChannel::getCentroidFlag() {return centroidFlag;}
/*
 Function: getCircularBufferFlag()
 
 Returns the status fo the circularBufferFlag to determine if audio output should be routed through the ouput circular buffer
 
 Returns:
 
 * circularBufferFlag - the flag indicating if the circular buffer is active
 */
bool AudioChannel::getCircularBufferFlag() {return circularBufferFlag;}
/*
 Function: getFFTFlag()
 
 Returns the status off the fftFlag to determine if the FFT has been performed on the current frame
 
 Returns:
 
 * fftFlag - the FFT flag for the current frame
 */
bool AudioChannel::getFFTFlag() {return fftFlag;}
/*
 Function: getMagFlag()
 
 Returns the status off the magflag to determine if the magnitude has been calculated on the current frame
 
 Returns:
 
 * magFlag - the magnitude flag for the current frame
 */
bool AudioChannel::getMagFlag() {return magFlag;}
/*
 Function: setCentroidFlag()
 
 Sets the status of the centroid flag
 
 Parameters:
 
 *flagVal: a boolean indicating whether or not the spectral centroid for the current frame has been performed
 */
void AudioChannel::setCentroidFlag(bool flagVal) {centroidFlag = flagVal;}
/*
 Function: setCircularBufferFlag()
 
 Sets the status of the circularBuffer flag
 
 Parameters:
 
 *flagVal: a boolean indicating the status of the inputCircular buffer
 */
void AudioChannel::setCircularBufferFlag(bool flagVal) {circularBufferFlag = flagVal;}
/*
 Function: setFFTFlag()
 
 Sets the status of the FFT flag
 
 Parameters:
 
 *flagVal: a boolean indicating whether or not the fft for the current frame has been computed
 */
void AudioChannel::setFFTFlag(bool flagVal) {fftFlag = flagVal;}
/*
 Function: setMagFlag()
 
 Sets the status of the magnitude flag
 
 Parameters:
 
 *flagVal: a boolean indicating whether or not the magnitude has been computed for the current frame's spectrum
 */
void AudioChannel::setMagFlag(bool flagVal) {magFlag = flagVal;}

/*
	Section: Public Members
	
	The AudioChannel class contains many buffers required for storing various aspects of the current audio frame
	as well as the input/output buffer access so that flash can read/write directly from/to the C namespace. The
	following description includes memory and flags.
	
	Audio Buffers:
	
	inAudioFrame - Before a frame can be processed, The DATF class writes samples directly to this buffer via a
		pointer. The buffer is then used to fill an input Circular Buffer for overlap processing.
	outAudioFrame - When the ALF requests audioplayback, outAudioFrame is the buffer that provides the output
		audio samples. The ALF class reads directly from this buffer via a pointer.			
	inBuffer - an instance of the CircularBuffer class. After inAudioFrame has been initialized by DATF. Methods in
		ALFPackage copy the data to inBuffer. This is required since features require overlap, so past samples are needed
		for future frames.		
	outBuffer - an instance of the CircularBuffer class. When the circularBufferFlag is set, a method the checkOutputBuffer
		method in ALF package loads outAudioFrame with samples from this buffer. This is necessary since filtering operations
		produce sequences that are larger than what can be stored in outAudioFrame (and played back at once by Flash). 
	
	Audio Framing:
	
	inAudioFrameSamples - the number of smples currently contained in inAudioFrame buffer.	
	outAudioFrameSamples - the number of samples currently contained in outAudioFrame.
	
			
	Spectral Data:
	
	fftFrame - after inAudioFrame is initialized, fftFrame is copied with the same data in order to allow for FFT
		computation with an in place algorithm. In future revisions, an FFT using separte input output memory will
		be implemented to avoid the copy.		
	freqData -  an array containing the allowable frequencies corresponding to the specified FFT size. This is needed
		for certain feature computation algorithms		
	magSpectrum - an array that contains the magnitude spectrum computed from fftFrame. This is required for certain
		feature computation algorithms.		
	spectrumPrev - an array containing the previous frame's magnitude spectrum. This is required for the spectral
		flux algorithm.		
	hannCoefficients - an array containing the the Hann Window coefficients required for a tapered window. The values
		in this array are multiplied after fftFrame has been initialized and before the FFT computation.
	
	Flags:
	
	firstFrame - a boolean indicating if the frame being processed is the first frame. Referenced by functions
			in ALF package.			
	bufferReady - a bool indicating if the buffer is ready (may be deprecated, need to check).	
	stereo -  a bool indicating if the AudioChannel is one of of a stereo pair.
	channelName - a string indicating the specified name of the AudioChannel.
	
*/

/*
	Section: Private Members
	
	Includes variables and flags.
			
	Audio Framing:

	fs - The sampling rate of the audio under analysis.	
	hopSize - The size of the current audio frame set by the user's desired sampling rate. Note that spectral
			  analysis is computed on twice the amount of data, that is a frame is twice the size of the hopSize. 
			  Only hopSize unique samples are read in on each frame, except the first where twice as many are read in.
	fftSize - The size of the fft used for computation. In general, this is not the same as hopSize
		because it requires powers of 2. fftSize = nextpow2(fftSize).	
	audioOutputSize - not sure what this one is doing here. Deprecate??
	
	Flags:
	
	fftFlag - This flag is true when the spectrum has been calculated on the current frame.	
	centroidFlag - This flag is true when the centroid has been calculated on the current frame.	
	magFlag - This flag is true when the magnitude spectrum has been calculated on the current frame.		
	circularBufferFlag - This flag is true when the circular output buffer is being used (i.e. audio 
		is being synthesized with ALF.
	
	Other:
	
	float centroidVal - Holds the value of the spectral centroid for the current frame;
*/
