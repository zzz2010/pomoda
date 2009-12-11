/* MiscUtilities.h		Miscellaneous Utilities

   Copyright 2004, Wong Chi Kwong, all rights reserved.

   This module contains miscellaneous utility functions.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/
#pragma once
#ifndef __MISC_UTILITIES_H__
#define __MISC_UTILITIES_H__

#include "TypeNLimit.h"
#include "stdio.h"

#define truncateRight(value, offset)			( (value) >> (offset) << (offset) )
#define truncateLeft(value, offset)				( (value) << (offset) >> (offset) )
// alignBoundary must be power of 2
#define downwardAlign(offset, alignBoundary)	( ((offset) + (alignBoundary) - 1) & (- (alignBoundary)) )
#define upwardAlign(offset, alignBoundary)		( (offset) & (- (alignBoundary)) )
#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )
#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define med3Index(a, b, c, ia, ib, ic)			( a<b ? (b<c ? ib : a<c ? ic : ia) : (b>c ? ib : a>c ? ic : ia))
#define swap(a, b, t);							t = a; a = b; b = t;

void QuickSortSingle(VAL*  key, const VAL numItem);
VAL checkDuplicate(int *input, const VAL numItem, const int minValue, const int maxValue, char* text);
VAL leadingZero(const VAL input);
VAL ceilLog2(const VAL input);
VAL floorLog2(const VAL input);
VAL power(const VAL base, const VAL power);
void formatVALAsBinary(const VAL input, char* output, VAL bitGroup);
VAL getRandomSeed();
VAL setStartTime();
double getElapsedTime(VAL startTime);
void printElapsedTime(FILE *file, const VAL printHour, const VAL printMin, const VAL printSec, 
					  const VAL secNumberOfDecimal, const double seconds);

VAL reverseBit(VAL x);
void initializeVAL(VAL *startAddr, const VAL length, const VAL initValue);
void initializeCHAR(CHAR *startAddr, const VAL length, const CHAR initValue);
VAL numberOfMatchInVAL(VAL *startAddr, const VAL length, const VAL searchValue);
VAL numberOfMatchInCHAR(CHAR *startAddr, const VAL length, const CHAR searchValue);

void bitCopyNoDestOffset(VAL *destinationAddress, const VAL *sourceAddress,
							VAL sourceBitOffset, VAL copyLengthInBit);
void bitCopyNoDestBitOffset(VAL *destinationAddress, VAL destinationWordOffset,
							const VAL *sourceAddress, VAL sourceWordOffset,
							VAL sourceBitOffset, VAL copyLengthInBit);
VAL bitCopy(VAL *destinationAddress, VAL destinationWordOffset, VAL destinationBitOffset,
			 const VAL *sourceAddress, VAL sourceBitOffset, VAL copyLengthInBit);

#endif

