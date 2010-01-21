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
#include <string.h>;
#include <sys/time.h>;
#include <math.h>;
#include "AS3.h";
#include "AudioChannel.h";
#include "DSPFunctions.h";

/*
	Class: ALFPackage.cpp
	
	This is the wrapper class that manages communication and interfacing between Actionscript and C++. Information
	on using alchemy can be found at http://forums.adobe.com/community/labs/alchemy?view=discussions and
	http://labs.adobe.com/technologies/alchemy/
*/

//some vars for timing functions
struct timeval tv;
time_t startTime, endTime;
double timeDiff;

int DEBUG = 0;
char outStr[200];

//initialize the left and right audio classes
AudioChannel *leftCh = NULL;
AudioChannel *rightCh = NULL;

// constants
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging

/*macro functions for timing functions, need to run START_TIME before a function, END_TIME and 
TIME_DIFF after a function is called*/

#define START_TIME { \
gettimeofday(&tv, NULL);\
startTime = tv.tv_usec;	\
}
#define END_TIME {\
gettimeofday(&tv, NULL);\
endTime = tv.tv_usec;\
}
#define TIME_DIFF(funcName) {\
sprintf(outStr, "it took %i msecs to execute %s \n", ((endTime - startTime)/1000), funcName);\
DISP(outStr);\
}

//float *rirFilter, *dataArray, *filterArray;
int prevFirLength = 0;
int count = 0;

//function prototypes
AS3_Val clearAudioFrame(void *self, AS3_Val args);
AS3_Val destroyAudioChannel(void *self, AS3_Val args);
AS3_Val initAudioChannel(void *self, AS3_Val args);
AS3_Val performFFT(void* self, AS3_Val args);
AS3_Val performIFFT(void* self, AS3_Val args);
AS3_Val getMagSpec(void *self, AS3_Val args);
AS3_Val getFlux(void *self, AS3_Val args);
AS3_Val getCentroid(void *self, AS3_Val args);
AS3_Val getIntensity(void *self, AS3_Val args);
AS3_Val getRolloff(void *self, AS3_Val args);
AS3_Val getBandwidth(void *self, AS3_Val args);
AS3_Val getLPC(void *self, AS3_Val args);
AS3_Val getHarmonics(void *self, AS3_Val args);
AS3_Val addReverb(void *self, AS3_Val args);
AS3_Val checkOutputBuffer(void *self, AS3_Val args);
AS3_Val resetFlags(void *self, AS3_Val args);
AS3_Val clearAudioBuffer(void *self, AS3_Val args);
AS3_Val setInputBuffer(void *self, AS3_Val args);
AS3_Val setHopSize(void* self, AS3_Val args);
AS3_Val checkInputBuffer(void *self, AS3_Val args);
AS3_Val reInitializeChannel(void* self, AS3_Val args);
AS3_Val getInAudioPtr(void *self, AS3_Val args);


//internal methods, not exposed to ActionScript
void filter(AudioChannel *ch, float *fir, int firLength);
void FFT(AudioChannel *ch);
void IFFT(AudioChannel *ch);
void computeMagnitudeSpectrum(AudioChannel *ch);
void computeCentroid(AudioChannel *ch);

/**************** Actionscript *******************/

/*
	Group: Alchemy Interface Utilities

	Function: clearAudioBuffer
	
		This clears the outBuffer (circular) of the audioChannel
	
	Parameters:
	
		*chPtr - a pointer to the audio channel containing the buffer to be cleared.

	See Also:
	
	<AudioChannel.cpp>
*/
AS3_Val clearAudioBuffer(void *self, AS3_Val args){
	
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	//Clearing audio buffer
	chPtr->outBuffer->clearBuffer();
	
	return 0;
}
/*	
	Function: clearAudioFrame
	
		This function clears the buffer inAudioFrame in the AudioChannel class. This buffer is shared between Actionscript 
		and C++. Samples are written to this buffer from Actionscript in order to be read and processed by C/C++ functions.
	
	Parameters: 
	
		*chPtr - a pointer to the audio channel containing the buffer to be cleared.
	
	See Also:
	
		- <AudioChannel.cpp>, <AudioChannel.cpp.initChannel()>

*/
AS3_Val clearAudioFrame(void *self, AS3_Val args){

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	//Clearing inAudioFrame
	chPtr->clearInAudioFrame();
		
	return 0;
}
/*
	Function: checkOutputBuffer
	
	This function is called prior to playback in <ALF>. See the ALF documentation for the proper usage.
	
	If using a mono track, checkOutputBuffer will see if there are enough samples to play in the audio channel.
	The buffer being used for playback will automatically be detected and the read and write pointers are 
	compared. If there are enough samples (greater than 2048), then the function returns that the buffer is ready
	for playback and also specifies the number of samples to be played.
	
	Parameters:
	
		*chPtr - a pointer to the audio channel containing the buffer to be cleared.
	
	Returns:
	
		An Array with the first element the status of the buffer (1 for ready, 0 for not ready) and the second element 
		isthe number of samples to be played.
*/
AS3_Val checkOutputBuffer(void *self, AS3_Val args) {
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	/*This function:
		1) Checks to see if the useCircularBuffer Flag is set in AudioChannel. This
			routes the output Buffer into the audio Chain 
			a) if so determines if there is at least a hopSize's worth of samples
			 that can be read
				i) if so, copy into chPtr->inAudioFrame and return a 'true'
				ii) if not, don't copy anything over, return a 'false'
		2) If the flag is not set, the input Buffer is routed into the audio chain 
			since no processing is requried.
	*/
	//sprintf(outStr, "=========== checkOutputBuffer ==========="); DISP(outStr);	
	int maxNumSamples = 4096; //the maximum number of samples that flash can playback
	int readPtr, writePtr, maxWritePtr, bufferSize, i;
	int numSamplesToPlay = 0, numSamplesToCopy = 0, diffSamples = 0, maxSamples = 0;
	int playReady = 0, bufferReady = 0;
	bool dumpSamples = false;

	if(chPtr->getCircularBufferFlag()) {
		//This is the case where the circular buffer flag is in use
/*
		// below lines are commented for debugging purposes
		readPtr = chPtr->outBuffer->getReadPtr();
		writePtr = chPtr->outBuffer->getWritePtr();
		maxWritePtr = chPtr->outBuffer->getMaxWritePtr();
		if(chPtr == leftCh){
			sprintf(outStr, "Enter COB: readPtr: %i, writePtr: %i, maxWritePtr: %i", readPtr, writePtr, maxWritePtr); DISP(outStr);
		}
*/		
		
		// diffSamples tells us how many samples are available based on the read and write pointers
		// of the output buffer
		diffSamples = chPtr->outBuffer->getPtrDiff();
		
		if(diffSamples > 0){
			if(diffSamples >= 2048 && diffSamples <= 4096) {
				numSamplesToCopy = diffSamples;
			} else if (diffSamples > 4096){
				numSamplesToCopy = 4096;
			}
			else {
				if(chPtr->filterProcessing) {
					numSamplesToCopy = 0;
					dumpSamples = false;
				}
				else{
					numSamplesToCopy = 0;
					dumpSamples = true;
				}
			}
		}
		
		if(diffSamples == 0) {
			if(chPtr->filterProcessing){
				numSamplesToCopy = 0;
				dumpSamples = false;
			}else{
				maxSamples = chPtr->outBuffer->getMaxPtrDiff();
				if(maxSamples >= 2048 && maxSamples <= 4096) {
					numSamplesToCopy = maxSamples;
				} else if (maxSamples > 4096){
					numSamplesToCopy = 4096;
				} else {
					numSamplesToCopy = maxSamples;
					dumpSamples = true;
				}
				chPtr->outBuffer->setWritePtr(maxSamples); //we need to re-sync the maxWrite and writePtrs
			}
		}
		
		//now we know how many samples we can play
		numSamplesToPlay = numSamplesToCopy;
		
		//determine if output is ready
		if(numSamplesToPlay >= 2048 || dumpSamples){
			playReady = true;
			for(i = 0; i < numSamplesToCopy; i++){
				chPtr->outAudioFrame[i] = chPtr->outBuffer->readBuffer();
			}
		} else{ playReady = false; }
		
		// lines below commented out for debugging purposes
/*		
		readPtr = chPtr->outBuffer->getReadPtr();
		writePtr = chPtr->outBuffer->getWritePtr();
		maxWritePtr = chPtr->outBuffer->getMaxWritePtr();
		if(chPtr == leftCh){
			sprintf(outStr, "End COB: readPtr: %i, writePtr: %i, maxWritePtr: %i", readPtr, writePtr, maxWritePtr); DISP(outStr);
			sprintf(outStr, "playReady %i, dumpSamples %i, numSamplesToPlay %i \n", playReady, dumpSamples, numSamplesToPlay); DISP(outStr);
		}
*/		
		
		
		//The flag below is very important. Any processing call should set this to '1' so that checkOutputBuffer
		//method knows that it is waiting on samples. checkOutputBuffer will always reset it to 0 since new input
		//audio could end at any frame. It is up the processing function in use (i.e. addReverb) to reset this to 1.
		chPtr->filterProcessing = 0;

	}else {
		
/*		
		readPtr = chPtr->inBuffer->getReadPtr();
		writePtr = chPtr->inBuffer->getWritePtr();
		maxWritePtr = chPtr->inBuffer->getMaxWritePtr();
		if(chPtr == leftCh){
			sprintf(outStr, "Enter COB  LEFT: readPtr: %i, writePtr: %i, maxWritePtr: %i", readPtr, writePtr, maxWritePtr); DISP(outStr);			
		}
		if(chPtr == rightCh){
			sprintf(outStr, "Enter COB RIGHT: readPtr: %i, writePtr: %i, maxWritePtr: %i", readPtr, writePtr, maxWritePtr); DISP(outStr);			
		}		
*/		
		diffSamples = chPtr->inBuffer->getPtrDiff();
		
		if(diffSamples == 0){ //case where there are no more samples tocopy
			numSamplesToCopy = 0;
			dumpSamples = true;
		} else {
			numSamplesToCopy = diffSamples;
		}
		
		// Now we need to make sure we aren't copying more than flash can possibly play
		if(numSamplesToCopy > maxNumSamples) { numSamplesToCopy = maxNumSamples; }
		
		// Get the number of audioSamples currently stored in outAudioFrame
		int outAS = chPtr->getOutAudioFrameSamples();
		int totalSamples = 0;
		
		// Figure out if the numberSamplesToCopy + current contents is too big
		if(outAS + numSamplesToCopy > maxNumSamples){ numSamplesToCopy = maxNumSamples - outAS; }
		totalSamples = outAS + numSamplesToCopy;
		//sprintf(outStr, "outAS = %i", outAS); DISP(outStr);

		// Copy samples into outAudioFrame		
		//sprintf(outStr, "numSamplesToCopy: %i totalSamples = %i", numSamplesToCopy, totalSamples); DISP(outStr);			
		for(i = 0; i < numSamplesToCopy; i++){
			chPtr->outAudioFrame[outAS + i] = chPtr->inBuffer->readBuffer();
		}
		
		if(totalSamples >= 2048){	//this is the min. num of samples we can play
			playReady = 1;
			numSamplesToPlay = totalSamples;
			chPtr->setOutAudioFrameSamples(0);
		} else if (dumpSamples && numSamplesToPlay != 0){	//if we've reached the end of the buffer, play whatever's left
			playReady = 1;
			numSamplesToPlay = totalSamples;
			chPtr->setOutAudioFrameSamples(0);
		} else if(dumpSamples && numSamplesToPlay == 0){
			playReady = 1;
			numSamplesToPlay = 0;
		}
		else {
			playReady = 0;			//if we don't have enough, we're not ready to play yet
			numSamplesToPlay = 0;
			chPtr->setOutAudioFrameSamples(totalSamples);
		}
		
		//now we set AudioFrameSamples to 0 we are removing all the samples
		//(this line might need to be deleted)
		chPtr->setAudioFrameSamples(0);	
		
//		 TODO: 
/*
		readPtr = chPtr->inBuffer->getReadPtr();
		writePtr = chPtr->inBuffer->getWritePtr();
		maxWritePtr = chPtr->inBuffer->getMaxWritePtr();
		if(chPtr == leftCh){
			sprintf(outStr, "End COB  LEFT: readPtr: %i, writePtr: %i, maxWritePtr: %i\n", readPtr, writePtr, maxWritePtr); DISP(outStr);			
		}
		if(chPtr == leftCh){
			sprintf(outStr, "End COB RIGHT: readPtr: %i, writePtr: %i, maxWritePtr: %i", readPtr, writePtr, maxWritePtr); DISP(outStr);			
		}		
*/		
				
	}
	
	// This flag is set to false here since CheckOutputBuffer should not be called until both channels (for stereo) have
	// been copied into C/C++ memory. If the flag is set to false before the right channel is read, the playback between the
	// left and right channels will not be synchronized.
	if(chPtr->firstFrame) chPtr->firstFrame = false;
	
	// These are the values we return
	AS3_Val checkBufferRes = AS3_Array("IntType, IntType", playReady, numSamplesToPlay);
	
	return checkBufferRes;
	AS3_Release(checkBufferRes);
}
/*
	Function: initAudioChannel
	
	This funciton initializes the AudioChannel (either left or right). If using stereo audio, both left and right channels have to be 
	initialized separately, but they can be initialized to the same parameters. The channel will allocate memory for the input and
	output buffers as well as arrays for performing DFT/IDFT and calculating spectral features. It sets flags that are used to keep
	track of what has been calculated on the current frame as well as what buffers are in use to default values. DO NOT change these
	values unless you know what you are doing.
	
	Parameters:
	
		channelType - A String, "left" or "right" that indicates which channel is being initialized
		sampleRate - The sample rate of the audio input.
		hopSize - The size in samples of each frame. If using ALF, this value will be automatically calculated
					when the frame rate (fps) is specified. 
	
	Returns:
	
		*chPtr - A pointer to the AudioChannel initialized in C.
	
	See Also:

		<AudioChannel.cpp>
*/
AS3_Val initAudioChannel(void *self, AS3_Val args) {
		
	int sampleRate, hopSize, frameSize;
	char *channelType;
	AudioChannel *chPtr;			// A pointer to return to Flash so it knows which channel to work with
	float *inAudPtr, *outAudPtr;	// A pointer for a float audioArray that will be initialized here
	int *samplesPtr;
	
	AS3_ArrayValue(args, "StrType, IntType, IntType", &channelType, &sampleRate, &hopSize);
	
	// Set frame size 
	frameSize = 2*hopSize;
	
	if(strcmp(channelType,"leftCh") == 0) {
	
		// Create the left AudioChannel
		leftCh = new AudioChannel[1];
		chPtr = leftCh;
		
		// Initialize audio channel
		leftCh->initChannel(hopSize, nextPowerOf2(frameSize), sampleRate);
		getFreq(leftCh->freqData, leftCh->getFFTSize(), leftCh->getSampleRate());	// Compute frequency array for later usage
		hannWindow(leftCh->hannCoefficients, leftCh->getFFTSize());					// Compute hann window for later usage
		inAudPtr = &(chPtr->inAudioFrame[0]);										// Pointer to where Flash will write to
		outAudPtr = &(chPtr->outAudioFrame[0]);										// Pointer to where Flash will read from
		samplesPtr = &(chPtr->inAudioFrameSamples);									// Pointer for flash to r/w the number of samples
		char name[] = "left";
		chPtr->setChannelName(name);
		
		// TODO: delete comment
//		sprintf(outStr, "fs = %i, hopSize = %i, frameSize = %i, fftSize = %i", chPtr->getSampleRate(), chPtr->getHopSize(), 2*chPtr->getHopSize(), chPtr->getFFTSize()); DISP(outStr);
		
	} else if(strcmp(channelType,"rightCh") == 0) {
						
		// Initialize the right AudioChannel
		rightCh = new AudioChannel[1];
		chPtr = rightCh;

		// Initialize audio channel
		rightCh->initChannel(hopSize, nextPowerOf2(frameSize), sampleRate);		
		getFreq(rightCh->freqData, rightCh->getFFTSize(), rightCh->getSampleRate());		// Compute frequency array for later usage
		hannWindow(rightCh->hannCoefficients, rightCh->getFFTSize());						// Compute hann window for later usage
		inAudPtr = &(chPtr->inAudioFrame[0]);												// Pointer to where Flash will play from
		outAudPtr = &(chPtr->outAudioFrame[0]);												// Pointer to where Flash will read from
		samplesPtr = &(chPtr->inAudioFrameSamples);											// Pointer for flash to r/w the number of samples
		char name[] = "right";
		chPtr->setChannelName(name);
		chPtr->stereo = true;
		leftCh->stereo = true;

	} else {sprintf(outStr, "INVALID AUDIO CHANNEL SPECIFIED \n"); DISP(outStr);}
	
//	sprintf(outStr, "chPtr = %d, inAudPtr = %d, outAudPtr = %d, samplesPtr = %d",
//			leftCh, inAudPtr, outAudPtr, samplesPtr); DISP(outStr);
	
	// Return the requried pointers in an array
	AS3_Val ptrArray = AS3_Array("PtrType, PtrType, PtrType, PtrType", chPtr, inAudPtr, outAudPtr, samplesPtr);
	
	return ptrArray;
	AS3_Release(ptrArray);
}
/*
	Function: reInitializeChannel
			
	Re-initializes the channel after a new song of different sample rate has been loaded
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to find the FFT size of.
		hop - The new hopSize value
		sampleRate - The new sample rate
		
	Returns:
	
*/
AS3_Val reInitializeChannel(void* self, AS3_Val args) {

	int hop, fs, numCh;
	
	AudioChannel *chPtr;
	//sprintf(outStr, "reInitChannel chPtr = %i", chPtr);
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType", &chPtr, &hop, &fs, &numCh);
	
	sprintf(outStr, "hop = %i, fs = %i, numCh = %i", hop, fs, numCh); DISP(outStr);
	chPtr->resetChannel();
	chPtr->setHopSize(hop);
	chPtr->setSampleRate(fs);
	chPtr->reInitChannel(hop, nextPowerOf2(2*hop), fs, numCh);
	getFreq(chPtr->freqData, chPtr->getFFTSize(), chPtr->getSampleRate());		// Compute frequency array for later usage
	hannWindow(chPtr->hannCoefficients, chPtr->getFFTSize());					// Compute hann window for later usage
	
	
	sprintf(outStr, "after reInitChannel: fs = %i, hopSize = %i, frameSize = %i, fftSize = %i, stereo = %i firstFrame = %i \n", 
					chPtr->getSampleRate(), chPtr->getHopSize(), 2*chPtr->getHopSize(), chPtr->getFFTSize(), chPtr->stereo, chPtr->firstFrame); DISP(outStr);	

	int readPtr, writePtr, maxWritePtr;
	// below lines are commented for debugging purposes
	readPtr = chPtr->outBuffer->getReadPtr();
	writePtr = chPtr->outBuffer->getWritePtr();
	maxWritePtr = chPtr->outBuffer->getMaxWritePtr();
	if(chPtr == leftCh){
		sprintf(outStr, "after reInitChannel OutBuff: readPtr: %i, writePtr: %i, maxWritePtr: %i", readPtr, writePtr, maxWritePtr); DISP(outStr);
	}
	if(chPtr == rightCh){
		sprintf(outStr, "after reInitChannel OutBuff: readPtr: %i, writePtr: %i, maxWritePtr: %i", readPtr, writePtr, maxWritePtr); DISP(outStr);
	}	
		
//
//	sprintf(outStr, "chPtr = %d, inAudPtr = %d, outAudPtr = %d, samplesPtr = %d",
//			leftCh, &chPtr->inAudioFrame[0], &(chPtr->outAudioFrame[0]), &(chPtr->inAudioFrameSamples)); DISP(outStr);
	return 0;
}
/*
	Function: destroyAudioChannel
	
	This function deletes the AudioChannel objects associated with an old ALF object. Should be called by
	DATF before initializing a new AudioChannels.
*/
AS3_Val destroyAudioChannel(void *self, AS3_Val args){
	
	delete [] leftCh;
	delete [] rightCh;
	
	return 0;
}

AS3_Val getInAudioPtr(void *self, AS3_Val args){
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType, IntType, IntType", &chPtr);
	
	AS3_Val newPtr = AS3_Array("PtrType", &(chPtr->inAudioFrame[0]));
	
//	sprintf(outStr, "setting new inaudptr: %i \n", &(chPtr->inAudioFrame[0])); DISP(outStr)	
	
	return newPtr;
	AS3_Release(newPtr);
}
/*
	Function: setInputBuffer
	
	This function copies the data from inAudioFrame in <AudioChannel.cpp> to to the inputBuffer in <AudioChannel.cpp>. It also
	copies the data into the fftFrame in <AudioChannel.cpp> for in-place computation of the frequency spectrum. The fftFrame
	is cleared prior to copying which zero pads the values at the end of the array. The values are hann windowed as they 
	are copied into this array. The read pointer for the inputBuffer is set accordingly so that it is the same upon exiting 
	the function as it was entering the function. 
	
	Parameters:
	
		*chPtr - A pointer to an AudioChannel.

*/
AS3_Val setInputBuffer(void *self, AS3_Val args){

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	unsigned int i = 0;
	int samp = chPtr->inAudioFrameSamples;
	
	// uncomment below for debugging
	
	//TODO: Delete comment
/*
	sprintf(outStr, "========= setInputBuffer==========");DISP(outStr);
	int readPtr = chPtr->inBuffer->getReadPtr();
	int writePtr = chPtr->inBuffer->getWritePtr();
	if(chPtr == leftCh){
		sprintf(outStr, "LEFT  readPtr is %i and writePtr is %i", readPtr, writePtr);DISP(outStr);
	}else{
		sprintf(outStr, "RIGHT readPtr is %i and writePtr is %i \n", readPtr, writePtr);DISP(outStr);
	}
	sprintf(outStr, "samp =  %i", samp);DISP(outStr);
*/
	
	// Write inAudioFrame into the input buffer
	for(i = 0; i < samp; i++){chPtr->inBuffer->writeBuffer(chPtr->inAudioFrame[i]);}
	
	
//	readPtr = chPtr->inBuffer->getReadPtr();
//	writePtr = chPtr->inBuffer->getWritePtr();
//	if(chPtr == leftCh){
//		sprintf(outStr, "left readPtr is %i and writePtr is %i", readPtr, writePtr);DISP(outStr);
//	}else{
//		sprintf(outStr, "right readPtr is %i and writePtr is %i \n", readPtr, writePtr);DISP(outStr);
//	}	
	
//	sprintf(outStr, "\nfirstFrame = %i", chPtr->firstFrame);DISP(outStr);
	// If it is past the first frame, kick the read pointer back a frame
	if(!chPtr->firstFrame){ chPtr->inBuffer->setReadPtr(-chPtr->getHopSize());}
	

	
//	readPtr = chPtr->inBuffer->getReadPtr();
//	writePtr = chPtr->inBuffer->getWritePtr();
//	if(chPtr == leftCh){
//		sprintf(outStr, "left readPtr is %i and writePtr is %i", readPtr, writePtr);DISP(outStr);
//	}else{
//		sprintf(outStr, "right readPtr is %i and writePtr is %i \n", readPtr, writePtr);DISP(outStr);
//	}	
	// Clear the FFT frame (zero pad)
	chPtr->clearFFTFrame();	
	
	// Hann window frame and copy to fftFrame
	if(!chPtr->stereo){
		
		// Mono case
		if(chPtr->firstFrame){
			for(i = 0; i < samp; i++){
				chPtr->fftFrame[i] = chPtr->inBuffer->readBuffer()*chPtr->hannCoefficients[i];
			}
		}else{
			for(i = 0; i < chPtr->getHopSize() + samp; i++){
				chPtr->fftFrame[i] = chPtr->inBuffer->readBuffer()*chPtr->hannCoefficients[i];
			}		
		}

		// Kick input buffer read pointer back a frame again so it is in it's original position	
		chPtr->inBuffer->setReadPtr(-chPtr->getHopSize());
			
	} else if (chPtr == leftCh){ 
		// We are in the left channel but we have two channels so we need to wait for the right channel
		// to be written before copying the data to fftFrame since we will perform spectral analysis only
		// on mono data to avoid excess computation.
	} else{	
//		sprintf(outStr, "firstFrame = %i", chPtr->firstFrame);DISP(outStr);
		// Convert to mono for spectral analysis
		if(chPtr->firstFrame){	
			for(i = 0; i < samp; i++){
				leftCh->fftFrame[i] = ((leftCh->inBuffer->readBuffer() + rightCh->inBuffer->readBuffer())/2) * leftCh->hannCoefficients[i];
//				if(i >= 400 && i <= 500){
//					sprintf(outStr, "x[%i] = %f", i, leftCh->fftFrame[i]);DISP(outStr);
//				}
			}
		}else{
			for(i = 0; i < chPtr->getHopSize() + samp; i++){
				leftCh->fftFrame[i] = ((leftCh->inBuffer->readBuffer() + rightCh->inBuffer->readBuffer())/2) * leftCh->hannCoefficients[i];
			}		
		}
//	readPtr = chPtr->inBuffer->getReadPtr();
//	writePtr = chPtr->inBuffer->getWritePtr();
//	if(chPtr == leftCh){
//		sprintf(outStr, "left readPtr is %i and writePtr is %i", readPtr, writePtr);DISP(outStr);
//	}else{
//		sprintf(outStr, "right readPtr is %i and writePtr is %i \n", readPtr, writePtr);DISP(outStr);
//	}
		// Since we read from the left and right channels we must set both pointers back
		leftCh->inBuffer->setReadPtr(-leftCh->getHopSize());		
		rightCh->inBuffer->setReadPtr(-rightCh->getHopSize());					
	}

	// uncomment below for debugging
// TODO: Delete comment
/*
	readPtr = chPtr->inBuffer->getReadPtr();
	writePtr = chPtr->inBuffer->getWritePtr();
	if(chPtr == leftCh){		
		sprintf(outStr, "End of SetInputBuffer: LEFT  readPtr is %i and writePtr is %i\n", readPtr, writePtr); DISP(outStr);
	}else{
		sprintf(outStr, "End of SetInputBuffer: RIGHT readPtr is %i and writePtr is %i\n", readPtr, writePtr); DISP(outStr);
	}
*/	
	
	

	return 0;
}
/*
 Function: checkInputBuffer
 
 Before copying new data to inAudioFrame, this function determines if the new data can be added. If there is insufficient 
 space in the audio channel's input buffer, adding more samples can eventually lead to the buffer over-running
 itself. This function determines if more data should be added based on the positions of the read and write pointers.
 The setFrame function in <DATF.as> uses the flag returned by this function to determine if a new frame should be set.
 
 Parameters:
 
	*chPtr - A pointer to an AudioChannel.
 
 Returns:
 
	returnVal: an array that contains a boolean flag (1 or 0) indicating if the channel's input buffer can be written to
 
 */
AS3_Val checkInputBuffer(void *self, AS3_Val args) {

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	int diffSamples = 0;
	int readSamples = 0;
	
	// uncomment for debugging
/*	
	int readPtr = chPtr->inBuffer->getReadPtr();
	int writePtr = chPtr->inBuffer->getWritePtr();
	int bufferSize = chPtr->inBuffer->getBufferSize();
	if(chPtr == leftCh){
		sprintf(outStr, "===== checkInputBuffer ====="); DISP(outStr);
		sprintf(outStr, "leftCh: readPtr %i writePTr %i bufferSize %i \n", readPtr, writePtr, bufferSize); DISP(outStr);
	}else {
		sprintf(outStr, "rightCh: readPtr %i writePTr %i bufferSize %i \n", readPtr, writePtr, bufferSize); DISP(outStr);
	}
*/
	
	
	diffSamples = chPtr->inBuffer->getPtrDiff();
	
	if(diffSamples > 2048){ readSamples = 0; }
	else{ readSamples = 1; }
	
	// uncomment for debugging
/*	
	if(chPtr == leftCh) {
		sprintf(outStr, "leftCh readSamples: %i \n", readSamples); DISP(outStr);
	} else{
		sprintf(outStr, "rightCh readSamples: %i \n", readSamples); DISP(outStr);
	}
*/	
	
	//return the requried pointers in an array
	AS3_Val returnVal = AS3_Array("IntType", readSamples);
	
	return returnVal;
	AS3_Release(returnVal);
}
/*
	Function: reset
	
		This funciton clears the flags that keep track of what features/spectrum etc. have been
		calculated in the current frame. This should be reset either at the end of a frame when
		no more data will be asked for or at the beginning of a frame PRIOR to calculating anything.

	See Also:
	
		<DATF.setFrame>

*/
AS3_Val resetFlags(void *self, AS3_Val args){
	
	leftCh->clearFlags();
	if(rightCh != NULL)	rightCh->clearFlags();
	
	return 0;
}

/*
	Function:
	
	Parameters:
	
	Returns:
	
	See Also:

*/
AS3_Val addReverb(void *self, AS3_Val args){

	// Invokes the RIR function and the filter method, in order to perform fast convolution based filtering
	
	float *rirFilter;
	char *active;
	AudioChannel *chPtr;			// AudioChannel which we will be doing the processing on
	
	// Default arguments, but make these arguments user specified from the DATF layer
	// default to a room size of 6m x 6m x 6m;
	double roomLen = 6, roomWid = 6, roomHt = 6, srcLen = 3, srcWid = 3, srcHt = 3, micLen = 3, micWid = 3.5, micHt = 3;
	double echoStrength = .5;
	
	// Read in user selected parameters
	AS3_ArrayValue(args, "StrType, IntType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType, DoubleType", 
				   &active,
				   &chPtr,
				   &roomLen, &roomWid, &roomHt, 
				   &srcLen, &srcWid, &srcHt, 
				   &micLen, &micWid, &micHt,
				   &echoStrength);
				   					  			   	
	if(strcmp(active,"on") == 0) {
		
		// We are brining reverb into the audio processing chain, so we need to indicate that the channel's output circular buffer is engaged
		chPtr->setCircularBufferFlag(true);
		
		int rirLen[1];
		int newFilter = 0;
		
		chPtr->roomSize[0] = (float)roomLen;	chPtr->sourcePosition[0] = (float)srcLen;	chPtr->micPosition[0] = (float)micLen;
		chPtr->roomSize[1] = (float)roomWid;	chPtr->sourcePosition[1] = (float)srcWid;	chPtr->micPosition[1] = (float)micLen;
		chPtr->roomSize[2] = (float)roomHt;		chPtr->sourcePosition[2] = (float)srcHt;	chPtr->micPosition[2] = (float)micLen;
		
		// Check to see if there are new room parameters
		bool needRoom = chPtr->checkRoom(chPtr->roomSize[0], chPtr->roomSize[1], chPtr->roomSize[2],
 									 	 chPtr->sourcePosition[0], chPtr->sourcePosition[1], chPtr->sourcePosition[2],
										 chPtr->micPosition[0], chPtr->micPosition[1], chPtr->micPosition[2], echoStrength);

		// Generate filter if needed
		if(needRoom){
			chPtr->filter = rir(chPtr->getSampleRate(), echoStrength, chPtr->micPosition, chPtr->roomSize, chPtr->sourcePosition, rirLen);
			chPtr->filterLen = rirLen[0];	
			chPtr->setRoom(chPtr->roomSize[0], chPtr->roomSize[1], chPtr->roomSize[2],
 						   chPtr->sourcePosition[0], chPtr->sourcePosition[1], chPtr->sourcePosition[2],
						   chPtr->micPosition[0], chPtr->micPosition[1], chPtr->micPosition[2], echoStrength);
		}
		
		// Tell the class that filter processing is in effect
		chPtr->filterProcessing = 1;
		filter(chPtr, chPtr->filter, chPtr->filterLen);
		
	} else if(strcmp(active,"off") == 0) {
		// No filtering needed, reverb is off so we remove it from the audio Processing chain
		chPtr->setCircularBufferFlag(false);
	} else {
		sprintf(outStr, "invalid 'active' flag for addReverb function!!! \n"); DISP(outStr);
		chPtr->setCircularBufferFlag(false);
	}

	return 0;
}
/*
	Group: DSP Functions

	Function: getBandwidth
	
		A function to calculate the spectral bandwidth of the current frame. If not already computed for the current frame, the 
		magnitude spectrum and spectral centroid will be calculated.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel on which to compute the bandwidth of.
	
	Returns:
	
		Bandwidth value, a float returned to AS as a Number.
	
	See Also:
	
		<DSPFunctions.c.bandwidth()>

*/
AS3_Val getBandwidth(void *self, AS3_Val args) {

	//Method has a dependency on FFT, MagSectrum, Centroid being computed first....
	AudioChannel *chPtr;
	float bandwidthVal;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	// Check for centroidFlag
	if(chPtr->getCentroidFlag() == 0) {
		computeCentroid(chPtr);		
	}
	
	// Compute bandwidth
	bandwidthVal= bandwidth(chPtr->magSpectrum, chPtr->freqData, chPtr->getCentroidVal(), chPtr->getFFTSize(), chPtr->getSampleRate());	
	if(bandwidthVal != bandwidthVal) bandwidthVal = 0;	// Handle NaN

	
	return AS3_Number(bandwidthVal);
}
/*
	Function: getCentroid
	
		A function to calculate the spectral centroid of the current frame. If not already computed for the current frame, the 
		magnitude spectrum will be calculated.	
	
	Parameters:

		*chPtr - A pointer to the audio channel on which to compute the bandwidth of.
	
	Returns:
	
		The centroid value as a float.
		
	See Also:
	
		<DSPFunctions.c.centroid()>

*/
AS3_Val getCentroid(void *self, AS3_Val args) {

	//Method has dependency on FFT, MagSpectrum being comptued first...
	AudioChannel *chPtr;
	float centroidVal;
	AS3_ArrayValue(args, "IntType", &chPtr);

	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr);
	}

	centroidVal = centroid(chPtr->magSpectrum, chPtr->freqData, chPtr->getFFTSize(), chPtr->getSampleRate());
	if(centroidVal != centroidVal) centroidVal = 0;	// Handles NaN

	// Store value in audio channel and set flag
	chPtr->setCentroidVal(centroidVal);
	chPtr->setCentroidFlag(true);	
	
	return AS3_Number(centroidVal);
}
/*
	Function: getFlux
	
		A function to calculate the spectral flux of the current frame. The current flux value is stored for use in
		calculation of the flux for the next frame.
		
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the flux of.
			
	Returns:
	
		The flux value as a float.
	
	See Also:
	
		<DSPFunctions.c.flux()>

*/
AS3_Val getFlux(void *self, AS3_Val args) {

	//Method has dependency on FFT, MagSpectrum being comptued first...
	AudioChannel *chPtr;
	float fluxVal;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr);
	}

	// Compute flux 
	fluxVal = flux(chPtr->magSpectrum, chPtr->spectrumPrev, chPtr->getFFTSize());
	if(fluxVal != fluxVal) fluxVal = 0;
	
	// Copy current spectrum to compute flux for next frame
	memcpy(chPtr->spectrumPrev, chPtr->magSpectrum, ( (chPtr->getFFTSize()/2) + 1)*sizeof(float) );
	
	return AS3_Number(fluxVal);
}
/*
	Function: 
		
	Parameters:
	
	Returns:
	
	See Also:

*/
AS3_Val getHarmonics(void *self, AS3_Val args) {

	// Dependencies - FFT is needed
	AudioChannel *chPtr;
	int error = 0;
	
	// Using fixed order now
	int order = 10;
	
	// Initialize arrays that will store the data
	float lpCE[2][order + 1];
	float resp[(chPtr->getFFTSize())/2 + 1];
	float harmPeaks[(chPtr->getFFTSize())/2 + 1];
	float harmFreqs[(chPtr->getFFTSize())/2 + 1];

	// Defaults to 10 requested harmonics unless the user specifies something else
	int i, reqHarms = 10, dbAdj = 3, foundHarms = 0;
	
	AS3_ArrayValue(args, "IntType, IntType", &chPtr, &reqHarms);
	if(DEBUG) {sprintf(outStr, "Reqesting LPC for data with length: %i and order: %i", chPtr->getFFTSize(), order); DISP(outStr);}
	if(DEBUG) {sprintf(outStr, "requesting %i harmonics for the the audio data \n", reqHarms); DISP(outStr);}
	
	// Some input checking to make sure the values are legit
	if(reqHarms <= 0 ) { reqHarms = 10;} 
	else if (reqHarms > (chPtr->getFFTSize())/2 + 1) { reqHarms = (chPtr->getFFTSize())/2 + 1;}
	else {}
	
	int numRows = 2;
	int numCols = order + 1;
	
	// First thing we do is check to see if we need to take an FFT or IFFT to get LPC data
	// FFT frame must contain the time domain samples for LPC
	if(chPtr->getFFTFlag()){IFFT(chPtr);}
	
	error = LPC(chPtr->fftFrame, chPtr->getFFTSize(), order, *lpCE, numRows, numCols);
	
	// Check for error in LPC formulation
	if(error == 1) { sprintf(outStr, "error in getHarmonics: LPC returned an error \n"); DISP(outStr);}
	else {
		
		// Initialize arrays for the impulse response
		float input[chPtr->getFFTSize()];			// This array has the impulse
		float output[chPtr->getFFTSize() + 2];		// This array has the output, which we will take FFT of
		for(i = 0; i < chPtr->getFFTSize(); i++) {	// Preallocate these to zero
			input[i] = 0; output[i] = 0;
		}
		input[0] = 1;								// Setting first element to 1 for unit impulse
		
		float numCoeffs[1]; numCoeffs[0] = 1;		// Init the numerator coeffs. its an all pole filter, so no zeros in the numerator
		float denomCoeffs[order + 1];				// Init the array for denom coeffs
		float gain = lpCE[1][0];					// Specify the gain explicitly
		for(i = 0; i < (order + 1); i++) {			// Initialize the denominator coefficients from the result of LPC
			denomCoeffs[i] = lpCE[0][i];
		}
		
		// Gen impulse response.... note we specify the length at 25% of fftSize to get a good ringing
		iirFilter(input, output, (chPtr->getFFTSize())/4, gain, numCoeffs, denomCoeffs, 1, (order + 1));
		
		// Take the FFT of the impulse response
		realft(output, chPtr->getFFTSize(), 1);	
		unpackFFT(output, chPtr->getFFTSize());
		
		// Get the magnitude response of IR
		magSpectrum(output, resp, chPtr->getFFTSize(), 1);
		
		// Now need to go back to frequency domain if we're not already in it
		if(!chPtr->getFFTFlag()){
			// Compute FFT here because we want frequency domain data
			FFT(chPtr);
		}
		
		// Obtain the magnitude spectrumm with dB's here
		magSpectrum(chPtr->fftFrame, chPtr->magSpectrum, chPtr->getFFTSize(), 1);
		
		// Vars for finding the max peak
		float alpha, beta, gamma, p;
		
		// Peak finding algorithm
		for(i = 1; i < (chPtr->getFFTSize())/2; i++) {
			// Detect region with candidate peak (s)
			if(chPtr->magSpectrum[i] > (resp[i] - dbAdj)) {
				// Look for a local maximum........
				if(chPtr->magSpectrum[i] > chPtr->magSpectrum[i-1] && chPtr->magSpectrum[i] > chPtr->magSpectrum[i+1]) {
					// We only add a candidate peak if it doesn't exceed reqHarms
					if(foundHarms < reqHarms) {	
								
						// Interpolate the maximum peak here
						alpha = chPtr->magSpectrum[i-1];
						beta = chPtr->magSpectrum[i];
						gamma =chPtr->magSpectrum[i+1];
						p = 0.5 * (alpha - gamma)/(alpha - 2*beta + gamma);
						
						// Corrected harmonic magnitude and frequency after correction
						harmPeaks[foundHarms] = beta - .25*(alpha - gamma)*p;
						harmFreqs[foundHarms] = ((float)chPtr->getSampleRate()/(float)chPtr->getFFTSize())*((float)i+p);
						foundHarms++;
						// Incrememt the foundHarms counter
					}
				}
			}
		}
	}
	
	// Outputs: number of harmonics found, ptr for harmoinc frequencies and amplitudes
	AS3_Val getHarmRes = AS3_Array("IntType, IntType, PtrType, PtrType", error, foundHarms, &harmPeaks, &harmFreqs);
	
	return getHarmRes;
	AS3_Release(getHarmRes);
}
/*
	Function: getIntensity
	
		This funciton calculates the spectral intensity and returns the value to Actionscript.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the intensity of.
	
	Returns:
	
		The intensity as a float.
	
	See Also:
	
		<DSPFunctions.c.intensity()>

*/
AS3_Val getIntensity(void *self, AS3_Val args) {

	// Method has dependency on MagSpectrum being comptued first...
	AudioChannel *chPtr;
	float intensityVal;
	AS3_ArrayValue(args, "IntType", &chPtr);

	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr);
	}
	
	// Compute intensity
	intensityVal = intensity(chPtr->magSpectrum, chPtr->getFFTSize());
	if(intensityVal != intensityVal) intensityVal = 0;					//Handle NaN
	
	return AS3_Number(intensityVal);
}
/*
	Function: 
		
	Parameters:
	
	Returns:
	
	See Also:

*/
AS3_Val getLPC(void *self, AS3_Val args) {

	int order, error = 0;
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType, IntType", &chPtr, &order);
	
	int numRows = 2;
	int numCols = order + 1;
	float lpCE[2][order + 1];
	
	// Call the LPC method...
	error = LPC(chPtr->inAudioFrame, chPtr->getHopSize(), order, *lpCE, numRows, numCols);
	if(error == 1) { sprintf(outStr, "error in getLPC: LPC returned an error \n"); DISP(outStr);}
	// We return the error to notify of an error in the computation
	// and we return the pointer to the coefficient array so we can read them out, if desired
	AS3_Val lpcRes = AS3_Array("IntType, PtrType", error, &lpCE);
	
	return lpcRes;
	AS3_Release(lpcRes);
}
/*
	Function: getComplexSpectrum
	
	This funciton calculates the complex frequency spectrum of the AudioChannel and returns a 
	pointer to the array containing the spectrum to C. 
		
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the intensity of.
	
	Returns:
	
		*fftFrame - A pointer to the magSpectrum array in the specified AudioChannel
	
	See Also:
	
		<AudioChannel.cpp>, <DSPFunctions.c.magSpectrum()>

*/
AS3_Val getComplexSpectrum(void *self, AS3_Val args){

	// Method has dependency on FFT being comptued first...
	AudioChannel *chPtr;
	float *fftPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
		
	// Check for fftFlag.....
	if(chPtr->getFFTFlag() == 0) {
		// Compute FFT if needed
		FFT(chPtr);
	}	
	
	fftPtr = &(chPtr->fftFrame[0]);

	// Return a magPtr so we can read the values in AS if we wish to access them
	return AS3_Ptr(fftPtr);
}
/*
	Function: getMagSpec
	
	This funciton calculates the magnitude spectrum of the AudioChannel and returns a pointer to 
	the array containing the spectrum to C. This is half of the spectrum since we are assuming
	a real input signal.
		
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the intensity of.
	
	Returns:
	
		*magSpectrum - A pointer to the magSpectrum array in the specified AudioChannel
	
	See Also:
	
		<AudioChannel.cpp>, <DSPFunctions.c.magSpectrum()>

*/
AS3_Val getMagSpec(void *self, AS3_Val args){

	//Method has dependency on FFT being comptued first...
	AudioChannel *chPtr;
	float *magPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);

	// Check for fftFlag.....
	if(chPtr->getFFTFlag() == 0) {
		// Compute FFT if needed
		FFT(chPtr);
	}

	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr);
	}	
	
	magPtr = &(chPtr->magSpectrum[0]);
	
	// Return a magPtr so we can read the values in AS if we wish to access them
	return AS3_Ptr(magPtr);
}
/*
	Function: getRolloff
	
		This funciton calculates the spectral rolloff and returns the value to Actionscript.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the rolloff for.
	
	Returns:
	
		The rolloff as a float.
	
	See Also:
	
		<DSPFuncs.rolloff()>

*/
AS3_Val getRolloff(void *self, AS3_Val args) {

	//Method has a dependency on FFT, Magnitude Spectrum being computed first
	AudioChannel *chPtr;
	float rolloffVal;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	if(chPtr->getMagFlag() == 0) {
		computeMagnitudeSpectrum(chPtr);
	}
	
	// Compute rolloff
	rolloffVal = rolloff(chPtr->magSpectrum, chPtr->getFFTSize(), chPtr->getSampleRate());
	if(rolloffVal != rolloffVal) rolloffVal = 0;
	
	return AS3_Number(rolloffVal);
}
/*
	Function: performIFFT
	
	This function calculates the IDFT of the audio frame for the specified AudioChannel.
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to compute the inverse DFT of.
	
	Returns:
	
		*ifftPtr - A pointer to an array containing the recovered signal.
		
	See Also:

	<DSPFuncs.realft()>, <DSPFuncs.repackFFT()>
	
	Notes:
	
		In general the funciton is used to reconstruct a signal that was already transformed using <performFFT>
*/
AS3_Val performIFFT(void* self, AS3_Val args) {

	//Only perform this method if FFT was computed first.......
	AudioChannel *chPtr;
	float *ifftPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	

	// Check to see if the FFT was performed
	if(chPtr->getFFTFlag() == 0) {
		sprintf(outStr, "No FFT was not computed for chPtr: %i !!!!!!!!\n", chPtr); DISP(outStr);
	} else {
		
		// An FFT is available to compute the inverse of
		repackFFT(chPtr->fftFrame, chPtr->getHopSize());
		realft(chPtr->fftFrame, chPtr->getHopSize(), -1);
		
		// Rescale for proper magnitude, numerical recipes does this...
		float scaleFactor = 2.0/(chPtr->getHopSize());	
		int i;
		for(i = 0; i < chPtr->getHopSize(); i++) {
			chPtr->fftFrame[i] *= scaleFactor;
		}
		
		chPtr->setFFTFlag(false);
	}
	
	ifftPtr = &(chPtr->fftFrame[0]);
	
	//return an ifftPtr so we can read the values in AS if we wish to access them
	return AS3_Ptr(ifftPtr);
}
/*
	Function: getFFTSize
	
	This function returns the the size of the FFT calculated on each frame
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to find the FFT size of.
	
	Returns:
	
		*fftSize - An int specifying the FFT size.
				
	See Also:

	<DSPFuncs.magSpectrum()>, <DATF.magSpectrum>
*/
AS3_Val getFFTSize(void* self, AS3_Val args) {

	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr);
	
	int fftSize = chPtr->getFFTSize();
	
	// Return the FFT Size
	return AS3_Number(chPtr->getFFTSize());
}
/*
	Function: setHopSize
		
	This function the hopSize
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to find the FFT size of.
		hop - The new hopSize value
	Returns:
	
*/
AS3_Val setHopSize(void* self, AS3_Val args) {

	int hop;
	
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr, &hop);
	
	//chPtr->setHopSize(hop);
	//chPtr->reInitChannel(hop, nextPowerOf2(2*hop));
	return 0;
}
/*
	Function: setSampleRate
		
	This function the SampleRate
	
	Parameters:
	
		*chPtr - A pointer to the audio channel to find the FFT size of.
		fs - The new SampleRate value
	Returns:
	
*/
AS3_Val setSampleRate(void* self, AS3_Val args) {

	int fs;
	
	AudioChannel *chPtr;
	AS3_ArrayValue(args, "IntType", &chPtr, &fs);
	
	chPtr->setSampleRate(fs);
	return 0;
}


//***************INTERNAL METHODS***************************************************

void FFT(AudioChannel *ch){
		
	int i;
	realft(ch->fftFrame, ch->getFFTSize(), 1);
	unpackFFT(ch->fftFrame, ch->getFFTSize());
	//realFFT(ch->fftFrame, ch->getFFTSize());
	ch->setFFTFlag(true);
}

void IFFT(AudioChannel *ch) {

	repackFFT(ch->fftFrame, ch->getFFTSize());
	realft(ch->fftFrame, ch->getFFTSize(), -1);
	
	float scaleFactor = 2.0/(ch->getFFTSize());	// Rescale for proper magnitude, numerical recipes does this...
	int i;
	for(i = 0; i < ch->getFFTSize(); i++) {
		ch->fftFrame[i] *= scaleFactor;
	}
	
	ch->setFFTFlag(false);
	
}

void filter(AudioChannel *ch, float *fir, int firLength) {
	//sprintf(outStr, "in filter firLength = %i", firLength); DISP(outStr);
	/*************************************************************
	 FILTER
	 This is an internal method that is used to implement fast convolution via frequency domain
		multiplication. it is called by the reverb method
	 INPUTS;
		*ch: a pointer for an AudioChannel object that will contain (among other things), the audio data
			and output buffers that the resulting, filtered data is placed i nto
		*fir: a pointer for an fir filter generated by some external means (i.e. filter method)
		firLength: an int specifying the length of the fir filter
	 OUTPUTS:
		none
	 */
	float *dataArray, *filterArray;
	float numOverlap;
	int convLen, i;
	bool processFilter = false;
	int minInputLen = 2048;
	
	// Filter will not execute unless there are at least 2048 samples ready from the input circular buffer
	int readPtr = ch->inBuffer->getReadPtr();
	int writePtr = ch->inBuffer->getWritePtr();
	int bufferSize = ch->inBuffer->getBufferSize();
	int diffSamples;
	
	// Calculate difference between read and write ptrs to determine if more samples should be read in
	if(writePtr - readPtr < 0){ 
		// Case where the writePtr wraps around
		diffSamples = writePtr + bufferSize - readPtr;
	}else {
		diffSamples = writePtr - readPtr;
	}
	
	if(diffSamples >= minInputLen){processFilter = true;}
	else{processFilter = false;}

	if(processFilter){
		convLen = diffSamples + firLength - 1;		// Size of convolution = legnth(seq1) + length(seq2) - 1
		
		// Now we need to figure out how large to allocate our arrays based on convLen
		// We need powers of 2 for efficiency and need to avoid circular convolution
		int pow2 = nextPowerOf2(convLen);
		dataArray = (float *) calloc((pow2 + 2), sizeof(float));
		filterArray = (float *) calloc((pow2 + 2), sizeof(float));
		
		// Need to copy the values in to the arrays we just allocated
		for(i = 0; i < diffSamples; i++){
			dataArray[i] = ch->inBuffer->readBuffer();
		}
		for(i = 0; i < firLength; i++){
			filterArray[i] = fir[i];
		}
		
		// Take the FFTs and unpack
		realft(dataArray, pow2, 1);
		unpackFFT(dataArray, pow2);
		realft(filterArray, pow2, 1);
		unpackFFT(filterArray, pow2);
		
		// Perform a frequency domain multiplication on data and filter arrays to perform convolution
		float realData, imagData, realFilter, imagFilter;
		for(i = 0; i < pow2; i = i+2) {
			realData = dataArray[i];
			imagData = dataArray[i+1];
			realFilter = filterArray[i];
			imagFilter = filterArray[i+1];
			filterArray[i] = realData*realFilter - imagData*imagFilter;
			filterArray[i+1] = imagData*realFilter + realData*imagFilter;
		}
		
		repackFFT(filterArray, pow2);		// Repack and take IFFT
		realft(filterArray, pow2, -1);
		
		float maxVal= 0.0;
		int clipping = 0;
		
		float scaleFactor = 2.0/(pow2);		// Need to re-scale the IFFT because numerical recipes weights it
		for (i = 0; i < pow2; i++) {
			filterArray[i] *= scaleFactor/2;
		}
		
		// Now need to write the values to the AudioChannel's output circular buffer
		int overlapSamples = convLen - diffSamples;
		
		// Figure out the amount overlap so we don't have sources clipping when they add together
		numOverlap = (float)convLen/(float)(diffSamples);
		numOverlap = ceil(numOverlap);
		
		int overlap = 2;
		
		for(i = 0; i< convLen; i++) {
			if(i < overlapSamples) {	
				// This means we add samples because they overlap with left over from previous frame
				ch->outBuffer->addBuffer(filterArray[i], 1);
			}else{						
				// If they are not overlapping, we just write to the buffer directly
				ch->outBuffer->writeBuffer(filterArray[i]);
			}
		}
		// Knock the write pointer back by the number of overlapping samples for the next frame
		ch->outBuffer->setWritePtr(-overlapSamples);
		
		// Free the memory we allocated
		free(dataArray);
		free(filterArray);
		
	}else{
		// If we don't have enough samples, we do nothing until there are more samples available
	}
	
}

void computeMagnitudeSpectrum(AudioChannel *ch) {

	// Check for fftFlag.....
	if(ch->getFFTFlag() == 0) {
		// Compute FFT if necessary
		FFT(ch);
	}
		
	// Default of '0' indicates that dB's will not be computed
	magSpectrum(ch->fftFrame, ch->magSpectrum, ch->getFFTSize(), 0);
	ch->setMagFlag(true);
}

void computeCentroid(AudioChannel *ch) {
	
	// Check for fftFlag
	if(ch->getMagFlag() == 0) {
		computeMagnitudeSpectrum(ch);
	}
	
	float centroidVal;
	centroidVal = centroid(ch->magSpectrum, ch->freqData, ch->getFFTSize(), ch->getSampleRate());
	ch->setCentroidVal(centroidVal);
	ch->setCentroidFlag(true);
}

//***************END INTERNAL METHODS***************************************************

//entry point for code
int main() {
	//define the methods exposed to ActionScript
	//typed as an ActionScript Function instance
	AS3_Val clearAudioBufferMethod = AS3_Function(NULL, clearAudioBuffer);
	AS3_Val clearAudioFrameMethod = AS3_Function(NULL, clearAudioFrame);
	AS3_Val initAudioChannelMethod = AS3_Function(NULL, initAudioChannel);
	AS3_Val destroyAudioChannelMethod = AS3_Function(NULL, destroyAudioChannel);
	AS3_Val performIFFTMethod = AS3_Function(NULL, performIFFT);
	AS3_Val getComplexSpectrumMethod = AS3_Function(NULL, getComplexSpectrum);
	AS3_Val getMagSpecMethod = AS3_Function(NULL, getMagSpec);
	AS3_Val getCentroidMethod = AS3_Function(NULL, getCentroid);
	AS3_Val getFluxMethod = AS3_Function(NULL, getFlux);
	AS3_Val getFFTSizeMethod = AS3_Function(NULL, getFFTSize);
	AS3_Val getIntensityMethod = AS3_Function(NULL, getIntensity);
	AS3_Val getRolloffMethod = AS3_Function(NULL, getRolloff);
	AS3_Val getBandwidthMethod = AS3_Function(NULL, getBandwidth);
	AS3_Val getLPCMethod = AS3_Function(NULL, getLPC);
	AS3_Val getHarmonicsMethod = AS3_Function(NULL, getHarmonics);
	AS3_Val addReverbMethod = AS3_Function(NULL, addReverb);
	AS3_Val checkOutputBufferMethod = AS3_Function(NULL, checkOutputBuffer);
	AS3_Val resetFlagsMethod = AS3_Function(NULL, resetFlags);
	AS3_Val setInputBufferMethod = AS3_Function(NULL, setInputBuffer);
	AS3_Val checkInputBufferMethod = AS3_Function(NULL, checkInputBuffer);
	AS3_Val setHopSizeMethod = AS3_Function(NULL, setHopSize);
	AS3_Val setSampleRateMethod = AS3_Function(NULL, setSampleRate);	
	AS3_Val reInitializeMethod = AS3_Function(NULL, reInitializeChannel);	
	AS3_Val getInAudioPtrMethod = AS3_Function(NULL, getInAudioPtr);	
	
	// construct an object that holds references to the functions
	AS3_Val result = AS3_Object("initAudioChannelC:AS3ValType", initAudioChannelMethod);
	AS3_SetS(result, "destroyAudioChannelC", destroyAudioChannelMethod);
	AS3_SetS(result, "setInputBufferC", setInputBufferMethod);
	AS3_SetS(result, "clearAudioFrameC", clearAudioFrameMethod);
	AS3_SetS(result, "clearAudioBufferC", clearAudioBufferMethod);
	AS3_SetS(result, "performIFFTC", performIFFTMethod);
	AS3_SetS(result, "getMagSpectrumC", getMagSpecMethod);
	AS3_SetS(result, "getComplexSpectrumC", getComplexSpectrumMethod);
	AS3_SetS(result, "getCentroidC", getCentroidMethod);
	AS3_SetS(result, "getFluxC", getFluxMethod);
	AS3_SetS(result, "getFFTSizeC", getFFTSizeMethod);
	AS3_SetS(result, "getIntensityC", getIntensityMethod);
	AS3_SetS(result, "getRolloffC", getRolloffMethod);
	AS3_SetS(result, "getBandwidthC", getBandwidthMethod);
	AS3_SetS(result, "getLPCC", getLPCMethod);
	AS3_SetS(result, "getHarmonicsC", getHarmonicsMethod);
	AS3_SetS(result, "addReverbC", addReverbMethod);
	AS3_SetS(result, "checkOutputBufferC", checkOutputBufferMethod);
	AS3_SetS(result, "checkInputBufferC", checkInputBufferMethod);
	AS3_SetS(result, "resetFlagsC", resetFlagsMethod);
	AS3_SetS(result, "setHopSizeC", setHopSizeMethod);
	AS3_SetS(result, "setSampleRateC", setSampleRateMethod);
	AS3_SetS(result, "reInitializeChannelC", reInitializeMethod);
	AS3_SetS(result, "getInAudioPtrC", getInAudioPtrMethod);
	
	// Release
	AS3_Release(setInputBufferMethod);
	AS3_Release(destroyAudioChannelMethod);
	AS3_Release(initAudioChannelMethod);
	AS3_Release(clearAudioFrameMethod);	
	AS3_Release(clearAudioBufferMethod);
	AS3_Release(performIFFTMethod);
	AS3_Release(getMagSpecMethod);
	AS3_Release(getComplexSpectrumMethod);	
	AS3_Release(getCentroidMethod);
	AS3_Release(getFluxMethod);
	AS3_Release(getIntensityMethod);
	AS3_Release(getRolloffMethod);
	AS3_Release(getBandwidthMethod);
	AS3_Release(getLPCMethod);
	AS3_Release(getHarmonicsMethod);
	AS3_Release(addReverbMethod);
	AS3_Release(checkOutputBufferMethod);
	AS3_Release(checkInputBufferMethod);
	AS3_Release(resetFlagsMethod);
	AS3_Release(setHopSizeMethod);
	AS3_Release(setSampleRateMethod);	
	AS3_Release(reInitializeMethod);
	AS3_Release(getInAudioPtrMethod);	
	
	// notify that we initialized -- THIS DOES NOT RETURN!
	AS3_LibInit(result);
	
	// should never get here!
	return 0;
}
