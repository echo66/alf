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

#ifndef DSPFUNCTIONS_H
#define DSPFUNCTIONS_H

extern void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned long n);
extern void four1(float data[], unsigned long nn, int isign);
void realft(float data[], unsigned long n, int isign);
void unpackFFT(float data[], long dataLength);
void repackFFT(float data[], long dataLength);
void magSpectrum(float fft[], float FFT[], int fftLength, int useDB);
float centroid(float spectrum[], float freq[], int frameSize, int fs);
float intensity(float spectrum[], int winLength);
void hannWindow(float hann[], int winLength);
float rolloff(float spectrum[], int winLength, int fs);
float bandwidth(float spectrum[], float freq[], float centroid, int winLength, int fs);
void getFreq(float freq[], int frameSize, int fs);
int LPC(float *audioSeg, int audioLength, int order, float *lpCE, int numRows, int numCols);
void freqResp(float *lpCE, float *resp, int fftSize, int numRows, int numCols, int useDB);
float flux(float spectrum[], float spectrumPrev[], int winLength);
void iirFilter(float *input, float *output, int seqLen, float gain, float *numCoeffs, float *denomCoeffs, int numOrder, int denomOrder);
float* rir(int fs, float refCo, float mic[], float room[], float src[], int rirLen[]);
int nextPowerOf2(int number);
void FFTHelper(float* x, int fftLength, float* X, float* scratch, float* twiddle, int twiddleLength);
float* realFFT(float *x, int fftLength);
float* FFT(float* x, int fftLength);
void computeTwiddleFactors(float* twiddle, int N);
void polarToComplex(float mag, float phase, float* ans);

#endif