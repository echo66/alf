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
#include <math.h>;
#include <sys/time.h>;
#include "MathFloatFuncs.h";
#include "MathDoubleFuncs.h";
#include "DSPFunctions.h";
#include "AS3.h";

struct timeval tv1;
time_t startTime1, endTime1;
double timeDiff1;

const float PI = 3.141592653589793;
char out1[200];

//RIR parameters
#define N 3										// related to the number of virtual sources
#define NN (N * 2 + 1)							// number of virtual sources

#define SWAP(a,b)tempr=(a);(a)=(b);(b)=tempr;
#define DISP(lit) AS3_Trace(AS3_String(lit));	//macro function for printing data for debugging
#define START_TIME {		\
gettimeofday(&tv1, NULL);\
startTime = tv1.tv_usec;	\
}
#define END_TIME {			\
gettimeofday(&tv1, NULL);\
endTime = tv1.tv_usec;	\
}
#define TIME_DIFF(funcName) {	\
sprintf(out1, "it took %i msecs to execute %s \n", ((endTime1 - startTime1)/1000), funcName);\
DISP(out1);\
}

/*
	Class: DSPFunctions.c
*/

/*
	Group: Fourier Analysis
	
	Function: getFreq
	
	Calculates the frequency associated with each bin of the discrete fourier transform.
	
	Parameters:
	
		freq - Pointer to an array. The array will be filled with the frequency values.
		frameSize - The number of points in the discrete fourier transform.
		fs - The sample rate of the signal to take the transform of.
*/
void getFreq(float freq[], int frameSize, int fs){

	//Create frequency array
	int n;
	float fnyq = fs/2;								//Nyquist freq
	float deltaF =  fnyq/(frameSize/2);				//Distance between the center frequency of each bin
	for (n = 0; n < (frameSize/2) + 1; n++){
		freq[n] = deltaF*n;
	}
}
/*****************************************************
 Function: magSpec
 Returns the magnitude spectrum of FFT data
 
 Parameters:
 
	fft[] - a pointer to an fft array obtained using realFT and unpacked
	FFT[] - a pointer to an array, allocated outside, that holds the magnitude
	fftLength - an int specifying th elength of the FFT
	useDB - a boolean var indicating to return decibels (1) or no decibels (0)
 
 Returns:
	None, arrays are passed by reference
 
 ****************************************************/
void magSpectrum(float fft[], float FFT[], int fftLength, int useDB){
	
	unsigned int i,j = 0;
		
	if(useDB) { 
		for(i = 0; i <= fftLength; i = i + 2){
			FFT[j] = 20*log10(fabs(sqrt(pow(fft[i], 2) + pow(fft[i + 1], 2)))); 
			j++;
		}
	}else{
		for(i = 0; i <= fftLength; i = i + 2){
			FFT[j] = sqrt(pow(fft[i], 2) + pow(fft[i + 1], 2));		
			j++;
		}
	}	
	
}
/************************************************************
 Function: freqResp 
 Calculates the frequency response from an ALL POLE transfer function
 
 Parameters:
 lpCE - a pointer to an array holding the linear predictioin coefficients
 resp - a pointer holding an array of bin frequencies
 fftSize - the size (frameSize) of the FFT
 numRows - the number of rows in the lpCE array
 numCols - the number of cols inthe lpCE array
 fs - the sampling frequency
 useDB - a flag that indicates to return in decibels (1) or not (0)
 
 Returns:
 None, values are passed by reference
 *************************************************************/
void freqResp(float *lpCE, float *resp, int fftSize, int numRows, int numCols, int useDB) {
	
	float gain = *(lpCE + numCols);						//assign the gain value for access
	float freqInc = PI/(fftSize/2 + 1);
	float rePart, imPart, denom;
	int i, c;
	
	for(i =0; i < (fftSize/2 + 1); i++) { resp[i] = 0; }
	
	for(i = 0; i < (fftSize/2 + 1); i++) {
		rePart = 0;
		imPart = 0;
		
		for(c = 1; c < numCols; c++) {
			rePart += (*(lpCE + c))*cos((float)(-c*i)*freqInc);
			imPart += (*(lpCE + c))*sin((float)(-c*i)*freqInc);
		}
		
		denom = sqrt(pow((1 + rePart),2) + pow((imPart), 2));
		resp[i] += gain/denom;									//!!!important! notice the += sign to accumulate values from each coefficient
		if(useDB) {
			resp[i] = 20*log10(fabs(resp[i]));
		}
	}
}
/**********************************************************
 Group: DSP Algorithms
 
 Function: iirFilter  
 Performs filtering with a provided transfer function based on a direct form II -transpose structure
 
 Parameters:
	input - the input sequence that will be used to filter the audio signal
	ouput - the output sequence where the audio will be stored
	seqLen - the length of the input and output sequence (they must be the same)
	gain - the gain of the filter if any
	numCoeffs - an array specifying the numerator coefficients
	denomCoeffs - an array specifying the denominator coefficients
	numOrder - the number of numCoefficients
	denomOrder - the number of denomCoefficients
 
 Format:
 - denomCoeffs: (1 a1  a2  a3 ...... aM), order of denom = M
 - numCoeffs: (1  b1  b2 ....... bN), order of num = N 
 - for proper tf, should have M >= N
 
 Returns:
		None, arrays are passed by reference
 
 ********************************************************/
void iirFilter(float *input, float *output, int seqLen, float gain, float *numCoeffs, float *denomCoeffs, int numOrder, int denomOrder) {
	
	int i, n, d, e;
	float v[denomOrder];						//filter memory for delays
	for(i = 0; i < denomOrder; i++) v[i] = 0;	//init to zeros...
	
	//peform the filtering..........
	for(i = 0; i < seqLen; i++){
		
		//calculate v[n] = input[n] - a1*v[n-1] - a2*v[n-2] - ...... - aM*v[n-M]
		v[0] = input[i];
		for(d = 1; d < denomOrder; d++){
			v[0] -= denomCoeffs[d]*v[d];
		}
		
		//now calculate y[n] = b0*v[n] + b1*v[n-1] + .......+ bN*v[n-N]
		output[i] = 0;
		for(n = 0; n < numOrder; n++){
			output[i] += numCoeffs[n]*v[n];
		}
		output[i] *= gain;
		
		//now, need to shift memory in v[n] = v[n-1], v[n-1] = v[n-2] ......
		for(e = denomOrder - 1; e > 0; e--){
			v[e] = v[e-1];
		}
	}
}
/*********************************************************************
 Function: LPC
 Performs linear predictive analysis on an audio segment for a desired order. The algorithm 
 works by computing the autocorrelation of the sequency followed by the Levinson Recursion to 
 computed the prediction coefficients.
 
 Parameters:
 audioSeg - a pointer to an array containing the frame of audio of interest
 audioLength - the length of audioSeg ...MUST BE A POWER OF 2!!!!!
 order - the desired order of LP analysis
 lpCE - a pointer for a two dimensional array containing gain and coefficients (Coefficients 
 in first row, gain in second)
 numRows - the number of rows in lpCE
 numCols - the number of cols in lpCE
 
 Returns:
 Returns an integer indicating whether or not an error ocurred in 
 the algorithm (1 = error, 0 = no error)
 **********************************************************************/
int LPC(float *audioSeg, int audioLength, int order, float *lpCE, int numRows, int numCols) {
	int error = 0;
	int c, p, i;
	
	if (order < 0)	error = 1;					//can't have negative order prediction coefficients
	else if (order > audioLength) error = 1;	//can't have more prediction coefficients than samples
	else {
		
		/**************************AUTOCORRELATION CODE HERE*****************************************
		 *	the following implementation utilizes the FFT algorithms to compute the autocorrelation of
		 *		a sequence. The data is copied into a new array which is padded to twice the length
		 *		with trailing zeros to avoid the circular convolution.
		 ********************************************************************************************/
		
		//needs to be twice as long to account for the circular convolutiion. The plus 2 is for repacking the FFT
		float *corrData = (float * ) calloc((audioLength*2), sizeof(float));
		float *corrDataOut = (float * ) calloc((audioLength*2), sizeof(float));
		float *twid = (float * ) calloc(audioLength*2, sizeof(float));
		float *invTwid = (float * ) calloc(audioLength*2, sizeof(float));
		//now need to copy the audio into the array
		for (i = 0; i < audioLength; i++) {
			corrData[i] = audioSeg[i];
		}
		
		computeTwiddleFactors(twid, audioLength*2, 1);
		realFFT(corrData, audioLength*2, twid, corrDataOut, 1);
		unpackFrequency(corrDataOut, audioLength*2);


		//now multiply the FFT by its conjugate....
		float RE, IM;
		for(i = 0; i < (audioLength*2); i = i+2) {
			RE = corrDataOut[i];
			IM = corrDataOut[i+1];
			corrDataOut[i] = RE * RE - (IM * -IM);
			corrDataOut[i+1] = IM * RE + (-IM * RE);
		}
		
		//repack the FFT and take the ifft
		computeTwiddleFactors(invTwid, audioLength*2, -1);
		pack(corrDataOut, (audioLength*2));
		realFFT(corrDataOut, (audioLength*2), invTwid, corrData, -1);			
	
		//rescale the FFT to compensate for frameSize weighting
		float scaleFactor = 2.0/(audioLength*2);
		for(i = 0; i < audioLength; i++) {
			corrData[i] = corrData[i] * scaleFactor;
		}
		
		//***********************************LEVINSON RECURSION FOLLOWS*********************************
		double lpcTemp[order + 1][order + 1];			//this array stores the partial correlaton coefficients
		double A = 0;									//this is the gain computed from the predicition error
		double E;		 								//this is the zeroth order predictor
		double ki;										//initial value of the weighting factor
		double sum;										//temp variable for storing values
		int j = 0;
		
		for(c = 0; c < order+1; c++) {					//initialize the lpc Temp array to zeros
			for(p = 0; p < order + 1; p++){				//for troubleshooting purposes
				lpcTemp[c][p] = 0;
			}
		}
		
		lpcTemp[0][0] = 0;									//initializations before recursion
		E = corrData[0];									//for the zeroth order predictor
		
		for(i = 1; i < order + 1; i++) {					//begin calculating the partial correlation coefficients
			sum = 0;										//temp sum for the calculation
			ki = 0;											//initial value of the weighting factor
			//compute the terms under the summation
			for(j = 1; j <= i - 1; j++){
				sum = sum + lpcTemp[j][i-1]*corrData[i-j];
			}
			ki = (corrData[i] - sum)/E;						//compute the weighting factor
			lpcTemp[i][i] = ki;								//store the weighting factor
			
			//recursively compute partio corr. coefficients
			for(j = 1;j <= i - 1; j++){
				lpcTemp[j][i] = lpcTemp[j][i-1] - (ki*lpcTemp[i-j][i-1]);
			}
			E = (1 - pow(ki,2))*E;							//updatethe prediction error
		}
		
		sum = 0;											//assign the pth order coefficients to an output vector and compute the gain A
		for(p = 1; p < (order + 1); p++){
			sum = sum + lpcTemp[p][order]*corrData[p];
		}
		A = corrData[0] - sum;
		A = sqrt(A);
		
		for(i = 0; i < numRows; i++) {
			for(j = 1; j < numCols; j++) {
				if(i == 0){
					*(lpCE + ((numCols*i) + j)) = -lpcTemp[j][order];
				}
			}
		}
		//need to manually set the gain and the leading coefficient of 1
		*lpCE = 1;
		*(lpCE + numCols) = A;
		
		free(corrData);
		free(corrDataOut);
		free(twid);
		free(invTwid);
	}
	return error;
}
/**********************************************************************************
 Function: rir 
 Generates a room impulse response for the specified room dimensions, speaker
 and microphone positions. An FIR represents the RIR
 
 Parameters:
 fs - the sample rae we wish to operate at
 refCo - the reflection coeffcients, a float between 0 and 1 (ecch strength)
 mic - the 3 dimensional positions of the microphone, in meters (LXWXH)
 room - the 3 dimensional room dimensions (L X W X H)
 src - the 3 dimensional position of the source (L X W X H)
 rirLen - the length of the resulting FIR filter
 
 Returns:
 Returns a float* for the resulting FIR filter
 **********************************************************************************/
float* rir(int fs, float refCo, float mic[], float room[], float src[], int rirLen[]){	
	int i, j, k;
	
	// Index for the sequence
	// nn=-n:1:n;
	float nn[NN];
	for(i=(-1 * N);i<=N;i++) {
		nn[i+N] = (float)i;		
	}
	
	
	// Part of equations 2, 3 & 4
	// rms = nn + 0.5 - 0.5*(-1).^nn;
	// srcs=(-1).^(nn);
	float rms[NN], srcs[NN];
	for(i=0;i<NN;i++) {
		rms[i] = nn[i] + 0.5 - (0.5 * ((float)pow(-1,(double)nn[i])));
		srcs[i] = (float)pow(-1,nn[i]);
	}	
	
	
	// Equation 2
	// xi=srcs*src(1)+rms*rm(1)-mic(1);
	// Equation 3
	// yj=srcs*src(2)+rms*rm(2)-mic(2);
	// Equation 4
	// zk=srcs*src(3)+rms*rm(3)-mic(3);
	float xi[NN], yj[NN], zk[NN];
	for(i=0;i<NN;i++) {
		xi[i] = srcs[i] * src[0] + rms[i] * room[0] - mic[0];
		yj[i] = srcs[i] * src[1] + rms[i] * room[1] - mic[1];
		zk[i] = srcs[i] * src[2] + rms[i] * room[2] - mic[2];
	}
	
	
	// Convert vectors to 3D matrices
	// [i,j,k]=meshgrid(xi,yj,zk);
	float meshOut[NN][NN][3*NN];
	meshgrid_float(xi,yj,zk,&meshOut[0][0][0],NN,NN,3*NN);
	
	
	// Equation 5
	// d=sqrt(i.^2+j.^2+k.^2);
	float d[NN][NN][NN];
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				d[i][j][k] = sqrt(pow(meshOut[i][j][k],2) + \
								  pow(meshOut[i][j][k + NN],2) + \
								  pow(meshOut[i][j][k + (2 * NN)],2));			
			}
		}
	}
	
	
	// Similar to Equation 6
	// time=round(fs*d/343)+1;
	float timeMat[NN][NN][NN];	
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				timeMat[i][j][k] = round_float(((fs * d[i][j][k] / 343) + 1));
			}
		}
	}
	
	
	// Convert vectors to 3D matrices
	// [e,f,g]=meshgrid(nn, nn, nn);
	float meshOutefg[NN][NN][3*NN];
	meshgrid_float(nn,nn,nn,&meshOutefg[0][0][0],NN,NN,3*NN);
	
	
	// Equation 9
	// c=r.^(abs(e)+abs(f)+abs(g));
	float c[NN][NN][NN];
	double constSum;
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				constSum = abs_float(meshOutefg[i][j][k]) + abs_float(meshOutefg[i][j][k + NN]) + abs_float(meshOutefg[i][j][k + 2 * NN]);
				c[i][j][k] = (float)(pow((double)refCo,constSum));
			}
		}
	}
	
	
	// Equation 10
	// e=c./d;
	float e[NN][NN][NN];
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				e[i][j][k] = c[i][j][k] / d[i][j][k];
			}
		}
	}
	
	
	// Equation 11
	// h=full(sparse(time(:),1,e(:)));
	int len = (int)pow(NN,3);
	
	// Left channel
	float* retVal = (float *)malloc(sizeof(float)*4);
	maxabs3D_float(&timeMat[0][0][0],NN,NN,NN,retVal);
	rirLen[0] = (int)retVal[3];
	float* rirArr = (float *)calloc(rirLen[0], sizeof(float));
	for(k=0;k<NN;k++) {
		for(j=0;j<NN;j++) {
			for(i=0;i<NN;i++) {
				// TOOK MINUS ONE AWAY FROM BELOW TO OFFSET SO THAT RIR[0] = 0
				rirArr[(int)timeMat[i][j][k] - 1] = rirArr[(int)timeMat[i][j][k] - 1] + e[i][j][k];
			}
		}
	}
	free(retVal);
	
	
	// Setting the time domain representation of the rir for the specified channel
	float* retVal2 = (float *)malloc(sizeof(float)*2);
	maxabs1D_float(rirArr,rirLen[0],retVal2);
	float sum = 0;
	for(i = 0;i < rirLen[0];i++) {
		rirArr[i] = rirArr[i] / retVal2[1];
		sum+=rirArr[i];
	}
	free(retVal2);
	
	return rirArr;	
}
/*********************************************************************
 Group: Spectral Features
 Function: bandwidth
 Computes the centroid on a frame-by-frame basis for a vector of sample data

 Parameters:
	x[] - array of FFT magnitude values
	fs - sample frequency	
	winLength - window length

 Returns:
	Returns a float that is the spectral bandwidth of the given audio frame
 
*********************************************************************/
float bandwidth(float spectrum[], float freq[], float centroid, int winLength, int fs){
	
	float *diff = (float *) malloc(sizeof(float)*(floor(winLength/2) + 1));
	int i;
	float band = 0;
	
	//Create frequency array
	float fnyq = fs/2;									//Nyquist freq
	float deltaF =  fnyq/(winLength/2);				//Distance between the center frequency of each bin
	for (i = 0; i < floor(winLength/2) + 1; i++){
		freq[i] = deltaF*(i);
	}
	//Find the distance of each frequency from the centroid
	for (i = 0; i < floor(winLength/2)+1; i++){
		diff[i] = fabs(centroid - freq[i]);	
		
	}
	
	//Weight the differences by the magnitude
	for (i = 0; i < floor(winLength/2)+1; i++){
		band = band + diff[i]*spectrum[i]/(winLength/2);
	}
	
	free(diff);
	return band;
}
/*********************************************************************

 Function: centroid
 Calculates the spectral centroid 

 Parameters:
	spectrum[] - the MAGNITUDE spectrum of the data to compute the centroid of
	fs - the sample frequency
	winLength - the number of points of the FFT taken to compute the associated spectrum

 Returns:
	Returns a float that is the centroid for the given frame

*********************************************************************/
float centroid(float spectrum[], float freq[], int winLength, int fs){
	
	int i;
	float centVal;
	float sumNum = 0;
	float sumDen = 0;
		
	//Calculate Centroid - sum of the frequencies weighted by the magnitude spectrum dided by 
	//the sum of the magnitude spectrum

	for (i = 0; i < (winLength/2) + 1; i++){
		sumNum = spectrum[i]*freq[i] + sumNum;
		sumDen = spectrum[i] + sumDen;
	}
	
	centVal = sumNum/sumDen;

	return centVal;
}
/*
	Function: flux
	Calculates the spectral flux.
	
	Parameters:
	
		spectrum - Pointer to the current spectrum
		spectrumPrev - Pointer to the spectrum from the previous frame
		winLength - The length of the DFT
*/
float flux(float spectrum[], float spectrumPrev[], int winLength){
	
	int i;
	
	//Calculate Flux
	float fluxVal = 0;
	for (i = 0; i < (winLength/2) + 1; i++){
		fluxVal = pow((spectrum[i] - spectrumPrev[i]),2) + fluxVal;
	}
	
	return fluxVal;
}
/*********************************************************************
 Function: intensity
 Calculates the spectral energy

 Parameters:
	spectrum[] - the MAGNITUDE spectrum of the data to 
	winLength - the window length
 
 Returns:
	Returns a float that is the energy for the given frame

*********************************************************************/
float intensity(float spectrum[], int winLength){

	//Find the total energy of the magnitude spectrum
	float totalEnergy = 0;
	int n;
	for (n = 0; n < (winLength/2) + 1; n++){
		totalEnergy = totalEnergy + spectrum[n];
	}

	return totalEnergy;
}
/*********************************************************************
 Function: rolloff
 Calculates the spectral centroid 

 Parameters:
	spectrum[] - the MAGNITUDE spectrum of the data to compute the centroid of
	fs - the sample frequency
	winLength - the window lenghth specified earlier

 Returns:
	Returns a float that is the centroid for the given frame
 
*********************************************************************/
float rolloff(float spectrum[], int winLength, int fs){
	
	float rollPercent = 0.85;
	float *freq = (float *) malloc(sizeof(float)*((winLength/2) + 1));	
	
	//Create frequency array
	float fnyq = fs/2;								//Nyquist freq
	float deltaF =  fnyq/(winLength/2);			//Distance between the center frequency of each bin
	int n;
	for (n = 0; n < (winLength/2) + 1; n++){
		freq[n] = deltaF*(n);
	}
	
	/*
	* Calculate Rolloff
	*/
	
	//Find the total energy of the magnitude spectrum
	float totalEnergy = 0;
	for (n = 0; n < (winLength/2) + 1; n++){
		totalEnergy = totalEnergy + spectrum[n];
	}
	
	//Find the index of the rollof frequency
	float currentEnergy = 0;
	int k = 0;
	while(currentEnergy <= totalEnergy*rollPercent && k <= winLength/2){
		currentEnergy = currentEnergy + spectrum[k];
		k++;
	
	}
		
	//Output the rollof frequency	
	float rollFreq = freq[k-1];
	free(freq);
	return rollFreq;
}
/************************************************************************
*	Function:  hannWindow
*
*	Parameters:		hann[] - An array that will contain the Hann coefficients.
*					winLength - The number of coefficients to be calculated
*
*	Returns:		Replaces the values in hann[] with the windowed values
*
*************************************************************************/
void hannWindow(float hann[], int winLength){

	int n;
	for (n = 0; n < winLength; n++){
		hann[n] = 0.5*(1 - cos(PI*2*(n)/(winLength - 1)));
	}

}
/************************************************************************
*	Function:  nextPowerOf2
*
*	Parameters:		number - The number to find the next highest power of two for.
*
*	Returns:		An integer which is the next highest power of two above the argument.
*
*************************************************************************/
int nextPowerOf2(int number){

	unsigned int count = 0;
	number--;
	while(pow(2,count) < sizeof(int)*8){
		
		number = number | number >> (int) pow(2,count);	
		count++;
	}
	number++;
	
    return number;
}



/*******************************************************************************
 Function: polarToComplex
 Converts polar numbers to complex numbers
 
 Parameters:
 mag - magnitude
 phase - phase
 ans - output array ans[0] = real, ans[1] = imag
 ******************************************************************************/
void polarToComplex(float mag, float phase, float* ans) {
    ans[0] = mag * cos(phase);
    ans[1] = mag * sin(phase);
}
/*******************************************************************************
 Function: computeTwiddleFactors
 Pre computes the twiddle factors needed in the FFTfunction
 
 Parameters:
 twiddle - array of size 2*fftLength
 fftLength - fftLength
 ******************************************************************************/
void computeTwiddleFactors(float* twiddle, int fftLength, float sign) {
    int k;
    float temp[2];

    for (k = 0; k < fftLength / 2; k++) {
        polarToComplex(1, sign * 2 * PI * (k) / fftLength, temp);
        twiddle[(2 * k)] = temp[0];
        twiddle[(2 * k) + 1] = temp[1];
    }
}
/*******************************************************************************
 Function: FFT
 Function for the fast fourier transform of an array with 
 length fftLength. fftLength must be a power of two.
 
 Parameters:
 x - input vector
 fftLength - length of fft (varies with each iteration
 ******************************************************************************/
void FFT(float* x, int fftLength, float* twiddles, float* output, int sign) {
    float* scratch = (float*) malloc(sizeof (float) * (2 * fftLength));
    FFTHelper(x, fftLength, output, scratch, twiddles, fftLength);

    int i = 0;
    if (sign == -1) {
        for (i = 0; i < fftLength; i++) {
            output[i] /= fftLength;
            output[i + fftLength] /= fftLength;
        }
    }

    free(scratch);
}
/*******************************************************************************
 Function: realFFT
 Function for the fast fourier transform of a real valued array
 with length fftLength. fftLength must be a power of two.
 
 Parameters:
 x - input vector
 fftLength - length of fft (varies with each iteration
 ******************************************************************************/
void realFFT(float* x, int fftLength, float* twiddles, float *out, int sign) {

    float xNew[fftLength];
    float scratch[2 * fftLength];
    float halfTwiddles[fftLength];

    int imagStart = fftLength / 2;
    int half = fftLength / 2;
    int quarter = half / 2;
    int i;

    /* Rearrange the original array. Even indexes have become the real part and
     * the odd indicies have become the imaginary parts */
    for (i = 0; i < half; i++) {
        xNew[i] = x[2 * i];
        xNew[i + (imagStart)] = x[2 * i + 1];
    }


    /* If we are taking the FFT */
    if (sign == 1) {

        computeTwiddleFactors(halfTwiddles, half, sign);
        /* FFT of new array */
        FFTHelper(xNew, half, out, scratch, halfTwiddles, imagStart);


        /* Manipulate tempOut for correct FFT */
        float temp1[2];
        float temp2[2];

        for (i = 1; i < quarter; i++) {
            temp1[0] = out[i];
            temp1[1] = out[i + (imagStart)];

            temp2[0] = out[half - i];
            temp2[1] = out[half - i + (imagStart)];

            out[i] = (0.5)*(temp1[0] + temp2[0]
                    + sign * twiddles[2 * i]*(temp1[1] + temp2[1])
                    + twiddles[(2 * i) + 1]*(temp1[0] - temp2[0]));

            out[i + (imagStart)] = (0.5)*(temp1[1] - temp2[1]
                    - sign * twiddles[2 * i]*(temp1[0] - temp2[0])
                    + twiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));

            out[half - i] = (0.5)*(temp1[0] + temp2[0]
                    - sign * twiddles[2 * i]*(temp1[1] + temp2[1])
                    - twiddles[(2 * i) + 1]*(temp1[0] - temp2[0]));

            out[half - i + (imagStart)] = (0.5)*(-temp1[1] + temp2[1]
                    - sign * twiddles[2 * i]*(temp1[0] - temp2[0])
                    + twiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));
        }

        temp1[0] = out[0];
        temp1[1] = out[(imagStart)];

        out[0] = temp1[0] + temp1[1];
        out[(imagStart)] = temp1[0] - temp1[1];
    }

    /* Inverse FFT of real signal */
    if (sign == -1) {

        /* Manipulate tempOutput for correct FFT */
        float temp1[2];
        float temp2[2];

        for (i = 1; i < quarter; i++) {
            temp1[0] = xNew[i];
            temp1[1] = xNew[i + (imagStart)];

            temp2[0] = xNew[half - i];
            temp2[1] = xNew[half - i + (imagStart)];

            xNew[i] = (0.5)*(temp1[0] + temp2[0]
                    + sign * twiddles[2 * i] * (temp1[1] + temp2[1])
                    - twiddles[(2 * i) + 1] * (temp1[0] - temp2[0]));

            xNew[i + (imagStart)] = (0.5)*(temp1[1] - temp2[1]
                    - sign * twiddles[2 * i]*(temp1[0] - temp2[0])
                    - twiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));

            xNew[half - i] = (0.5)*(temp1[0] + temp2[0]
                    - sign * twiddles[2 * i]*(temp1[1] + temp2[1])
                    + twiddles[(2 * i) + 1]*(temp1[0] - temp2[0]));

            xNew[half - i + (imagStart)] = (0.5)*(-temp1[1] + temp2[1]
                    - sign * twiddles[2 * i]*(temp1[0] - temp2[0])
                    - twiddles[(2 * i) + 1]*(temp1[1] + temp2[1]));
        }

        temp1[0] = xNew[0];
        temp1[1] = xNew[(imagStart)];

        xNew[0] = temp1[0] + temp1[1];
        xNew[(imagStart)] = temp1[0] - temp1[1];

        xNew[0] *= (0.5);
        xNew[imagStart] *= (0.5);

        computeTwiddleFactors(halfTwiddles, half, sign);
        FFTHelper(xNew, half, out, scratch, halfTwiddles, imagStart);
//		for(i = 0; i < fftLength; i++){
//			out[i] = out[i]/(fftLength/2);
//		}	
    }
}
/*******************************************************************************
 Function: FFTHelper
 Calcualtes the fast fourier transform of length N, N must be a power of two
 
 Parameters
 x - input vector
 fftLength - length of fft (varies with each iteration
 X - output vector
 scratch - empty array of size fftLength for use in the function
 twiddle - array of twiddle values
 twiddleLength - original fft length (never changes)
 ******************************************************************************/
void FFTHelper(float* x, int fftLength, float* X, float* scratch,
        float* twiddle, int imagStart) {
    int k, m, n;
    int skip;
    /* int imagStart = fftLength; */
    int evenItr = fftLength & 0x55555555;

    float* E, *D;
    float* Xp, *Xp2, *XStart;
    float temp[2], temp2[2];

    /* Special Case */
    if (fftLength == 1) {
        X[0] = x[0];
        X[1] = x[imagStart];
        return;
    }

    E = x;

    for (n = 1; n < fftLength; n *= 2) {
        XStart = evenItr ? scratch : X;
        skip = (fftLength) / (2 * n);
        Xp = XStart;
        Xp2 = XStart + (fftLength / 2);
        for (k = 0; k != n; k++) {

            temp[0] = twiddle[2 * (k * skip)];
            temp[1] = twiddle[2 * (k * skip) + 1];

            for (m = 0; m != skip; ++m) {
                D = E + (skip);

                temp2[0] = (*D * temp[0]) - (*(D + imagStart) * temp[1]);
                temp2[1] = (*D * temp[1]) + (*(D + imagStart) * temp[0]);

                *Xp = *E + temp2[0];
                *(Xp + imagStart) = *(E + imagStart) + temp2[1];

                *Xp2 = *E - temp2[0];
                *(Xp2 + imagStart) = *(E + imagStart) - temp2[1];

                Xp = Xp + 1;
                Xp2 = Xp2 + 1;
                E = E + 1;
            }
            E = E + skip;
        }
        E = XStart;
        evenItr = !evenItr;
    }
}

/*
	Function: pack
	
	Parameters:
	
		*in - A pointer to an array that contains the output data from <realFFT> after computing the *INVERSE*.
		fftLength - The number of FFT points used in calculating the spectrum.
*/
void pack(float* in, int fftLength) {
    int k;	
		
	in[1] = in[fftLength];		// Set second element to the real part at fs/2
	in[fftLength] = 0;			// Set imaginary component at fs/2 to zero
}
/*
	Function: unpackFrequency
	
	The output of RealFFT is in the form of Re[0], Re[1] ... Re[N/2], Im[0], Im[1], ..., Im[N/2]. This function
	changes the order to Re[0], Im[0], Re[1], Im[1], ... Re[N/2], Im[N/2]
	
	Parameters:
	
		*in - A pointer to an array that contains the output data from <realFFT>
		fftLength - The length of the FFT calculated.
*/
void unpackFrequency(float* in, int fftLength) {

	int k;
	float *temp = (float *) malloc(sizeof(float)*(fftLength*2));
    for (k = 0; k <= fftLength + 2; k++) {
        temp[k] = in[k];
    }
	
    for (k = 0; k <= fftLength/2; k++) {
        in[2*k] = temp[k];
		in[2*k + 1] = temp[k + fftLength / 2];
    }

	in[1] = 0;
	in[2*fftLength - 1] = temp[fftLength];
	free(temp);
}
/*
	Function: unpackTime
	
	This function reorders the values to be in the proper order after resynthesis.
	
	Parameters:
	
		*in - A pointer to an array that contains the output data from <realFFT>
		fftLength - The length of the FFT calculated.
*/
void unpackTime(float* in, int fftLength) {

	int k;
	float *temp = (float *) malloc(sizeof(float)*(fftLength*2));
    for (k = 0; k <= fftLength + 2; k++) {
        temp[k] = in[k];
    }

    for (k = 0; k <= fftLength/2; k++) {
        in[2*k] = temp[k];
		in[2*k + 1] = temp[k + fftLength / 2];
    }
	
	in[2*fftLength - 1] = temp[fftLength];
	free(temp);
}