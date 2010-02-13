#!/bin/bash
echo

which g++
echo

echo Compiling Source files.......

#Compile files
g++ ALFpackage.cpp -c -o ALFpackage.o;
g++ AudioChannel.cpp -c -o AudioChannel.o;
g++ CircularBuffer.cpp -c -o CircularBuffer.o
g++ DSPFunctions.c -c -o DSPFunctions.o;
g++ MathDoubleFuncs.c -c -o MathDoubleFuncs.o;
g++ MathFloatFuncs.c -c -o MathFloatFuncs.o;
g++ AudioChannel.o CircularBuffer.o ALFpackage.o DSPFunctions.o MathDoubleFuncs.o MathFloatFuncs.o -O3 -Wall -swc -o ALFPackage.swc;

echo Done!
echo

#Move SWC file
echo Moving SWC to Parent Directory...
mv ./ALFPackage.swc ../;
echo Done!
echo
echo Removing Binary Files...

#Cleaning up binaries - you messy betty
rm ALFpackage.o
rm AudioChannel.o
rm CircularBuffer.o
rm DSPFunctions.o
rm MathDoubleFuncs.o
rm MathFloatFuncs.o

echo Done!
echo
echo Process Complete 

say -v Zarvox finished