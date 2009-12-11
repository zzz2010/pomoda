/* MiscUtilities.c		Miscellaneous Utilities

   Copyright 2004, Wong Chi Kwong, all rights reserved.

   This module contains miscellaneous utility functions.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MiscUtilities.h"

void QuickSortSingle(VAL*  key, const VAL numItem) {

	VAL s;
	VAL i, j;
	VAL m, l, n;
	VAL b, c;
	VAL tempKey;
	VAL partitionKey;

	if (numItem <= 7) {	 // Insertion sort on smallest arrays
		for (i = 1; i < numItem; i++) {
			tempKey = key[i];
			for (j = i; j > 0 && key[j-1] > tempKey; j--) {
				key[j] = key[j-1];
			}
			if (j != i) {
				key[j] = tempKey;
			}
		}
		return;
	}

	if (numItem > 40) {    // Big arrays, pseudomedian of 9
		s = numItem / 8;
		l = med3(key[0], key[s], key[2 * s]);
		m = med3(key[numItem / 2 - s], key[numItem / 2], key[numItem / 2 + s]);
		n = med3(key[numItem - 1 - 2 * s], key[numItem - 1 - s], key[numItem - 1]);
	} else {	// Mid-size, med of 3
		l = key[0];
		m = key[numItem / 2];
		n = key[numItem - 1];
	}
	partitionKey = med3(l, m, n);

	b = 0;
	c = numItem - 1;

	for (;;) {
		while (b <= c && key[b] <= partitionKey) {
			b++;
		}
		while (b < c && key[c] > partitionKey) {
			c--;
		}
		if (b >= c) {
			break;
		}
		swap(key[b], key[c], tempKey);
		b++;
		c--;
	}

	c = numItem - b;

	if (b > c) {
		QuickSortSingle(key + b, c);
	}
	QuickSortSingle(key, b);
	if (b <= c) {
		QuickSortSingle(key + b, c);
	}

}

VAL checkDuplicate(int *input, const VAL numItem, const int minValue, const int maxValue, char* text) {

	VAL *present;
	VAL i;
	char defaultText[17] = "checkDuplicate()";

	if (text == NULL) {
		text = defaultText;
	}

	present = (VAL*)malloc((maxValue - minValue + 1) * sizeof(VAL));
	initializeVAL(present, maxValue - minValue + 1, 0);

	for (i=0; i<numItem; i++) {
		if (input[i] >= minValue && input[i] <= maxValue) {
			if (present[input[i] - minValue] > 0) {
				fprintf(stderr, "%s : Item %u and %u contains duplicate value of %d\n", 
							    text, present[input[i] - minValue], i, input[i]);
				free(present);
				return FALSE;
			}
			present[input[i] - minValue] = i;
		}
	}

	free(present);
	return TRUE;

}


VAL leadingZero(const VAL input) {

	VAL l;
	const static VAL leadingZero8bit[256] = {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
											 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
											 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
											 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
											 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	if (input & 0xFFFF0000) {
		if (input & 0xFF000000) {
			l = leadingZero8bit[input >> 24];
		} else {
			l = 8 + leadingZero8bit[input >> 16];
		}
	} else {
		if (input & 0x0000FF00) {
			l = 16 + leadingZero8bit[input >> 8];
		} else {
			l = 24 + leadingZero8bit[input];
		}
	}
	return l;

}

VAL ceilLog2(const VAL input) {

	if (input <= 1) {
		return 0;
	}

	return BITS_IN_WORD - leadingZero(input - 1);

}

VAL floorLog2(const VAL input) {

	if (input <= 1) {
		return 0;
	}

	return BITS_IN_WORD - leadingZero(input) - 1;
}

VAL power(const VAL base, const VAL power) {

	VAL i;
	VAL result = 1;

	for (i=0; i<power; i++)	{
		result *= base;
	}

	return result;

}

void formatVALAsBinary(const VAL input, char* output, VAL bitGroup) {

	int i, j=0;

	for (i=0; i<BITS_IN_WORD; i++) {
		if ((input & (input << i >> i)) >> (BITS_IN_WORD - i - 1)) {
			output[j] = '1';
		} else {
			output[j] = '0';
		}
		j++;
		if (bitGroup > 0 && bitGroup < BITS_IN_WORD) {
			if ((i+1) % bitGroup == 0) {
				output[j] = ' ';
				j++;
			}
		}
	}
	output[j] = '\0';

}

VAL getRandomSeed() {

	time_t timer;
	VAL t;
	VAL byteToCopy;
	void *startAddr;

	time(&timer);

	byteToCopy = min(sizeof(time_t), sizeof(VAL));
	startAddr = &timer + sizeof(time_t) - byteToCopy;

	memcpy(&t, startAddr, byteToCopy);

	return t;
}

VAL setStartTime() {

	return (VAL)clock();
	
}

double getElapsedTime(VAL startTime) {

	return  (double)((VAL)clock() - startTime) / (double)CLOCKS_PER_SEC;

}

void printElapsedTime(FILE *file, const VAL printHour, const VAL printMin, const VAL printSec, 
					  const VAL secNumberOfDecimal, const double seconds) {

	VAL hour, min;
	double sec;
	char secondDisplay[7] = "%.0f s";
	
	#ifdef DEBUG
	if (printHour == TRUE && printMin == FALSE && printSec == TRUE) {
		fprintf(stderr, "printElapsedTime(): Cannot skip minute only!\n");
		exit(1);
	}
	if (secNumberOfDecimal > 9) {
		fprintf(stderr, "printElapsedTime(): secNumberOfDecimal > 9!\n");
		exit(1);
	}
	#endif

	secondDisplay[2] += (char)secNumberOfDecimal;

	sec = seconds;
	min = (VAL)(seconds / 60);
	if (printSec == FALSE && printMin == TRUE) {
		if (seconds - min * 60 >= 30) {
			min++;
		}
	}
	if (printMin == TRUE) {
		sec -= min * 60;
	}

	hour = min / 60;
	if (printMin == FALSE) {
		if (min = hour * 60 >= 30) {
			hour++;
		}
	}
	if (printHour == TRUE) {
		min -= hour * 60;
	}

	if (printHour == TRUE) {
        fprintf(file, "%d h ", hour);
	}
	if (printMin == TRUE) {
		fprintf(file, "%d m ", min);
	}
	if (printSec == TRUE) {
		fprintf(file, secondDisplay, sec);
	}
	
	fprintf(file, "\n");

}

VAL reverseBit(VAL x)
{
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return((x >> 16) | (x << 16));

}

void initializeVAL(VAL *startAddr, const VAL length, const VAL initValue) {

	VAL i;

	for (i=0; i<length; i++) {
		startAddr[i] = initValue;
	}

}

void initializeCHAR(CHAR *startAddr, const VAL length, const CHAR initValue) {

	VAL i;

	for (i=0; i<length; i++) {
		startAddr[i] = initValue;
	}

}

VAL numberOfMatchInVAL(VAL *startAddr, const VAL length, const VAL searchValue) {

	VAL i;
	VAL numberOfMatch = 0;

	for (i=0; i<length; i++) {
		if (startAddr[i] == searchValue) {
			numberOfMatch++;
		}
	}

	return numberOfMatch;

}

VAL numberOfMatchInCHAR(CHAR *startAddr, const VAL length, const CHAR searchValue) {

	VAL i;
	VAL numberOfMatch = 0;

	for (i=0; i<length; i++) {
		if (startAddr[i] == searchValue) {
			numberOfMatch++;
		}
	}

	return numberOfMatch;

}


// destinationAddress + up to the next 4 words boundary (depending on copy length) will be overridden
// sourceAddress + multiple of 4 words (depending on copy length) will be accessed
// calling program must ensure that those address, although not directly useful/used as it seems,
// must be safely overridden/accessed
// The remaining bits in the resulting ending word, if the resulting bit offset > 0, 
// are guaranteed to be cleared as 0
// The bits in the resulting ending word are undefined if the resulting bit offset = 0
// The remaining words (to make up the last 4 word multiple) are undefined

void bitCopyNoDestOffset(VAL *destinationAddress, const VAL *sourceAddress,
							VAL sourceBitOffset, VAL copyLengthInBit) {

	VAL i;
	VAL rightShift;
	VAL copyLeftBuffer[4], copyRightBuffer[4];
	VAL copyWordLength, copyWordLengthRoundTo4;

	#ifdef DEBUG
	if (copyLengthInBit == 0) {
		fprintf(stderr, "bitCopyNoDestOffset() : copyLengthInBit = 0!\n");
		exit(1);
	}
	#endif

	copyWordLength = (copyLengthInBit + BITS_IN_WORD_MINUS_1) / BITS_IN_WORD;

    if (sourceBitOffset == 0) {
		memcpy(destinationAddress, sourceAddress, copyWordLength * 4); 
	} else {
		rightShift = BITS_IN_WORD - sourceBitOffset;
		copyWordLengthRoundTo4 = (copyWordLength + 3) & FOUR_MULTIPLE_MASK;
		for (i=0; i<copyWordLengthRoundTo4; i+=4) {
			// This is supposed to generate SSE2 codes
			// Need to check rewrite using intrinics is necessary
			copyLeftBuffer[0] = sourceAddress[i + 0] << sourceBitOffset;
			copyLeftBuffer[1] = sourceAddress[i + 1] << sourceBitOffset;
			copyLeftBuffer[2] = sourceAddress[i + 2] << sourceBitOffset;
			copyLeftBuffer[3] = sourceAddress[i + 3] << sourceBitOffset;
			copyRightBuffer[0] = sourceAddress[i + 1] >> rightShift;
			copyRightBuffer[1] = sourceAddress[i + 2] >> rightShift;
			copyRightBuffer[2] = sourceAddress[i + 3] >> rightShift;
			copyRightBuffer[3] = sourceAddress[i + 4] >> rightShift;
			destinationAddress[i + 0] = copyLeftBuffer[0] | copyRightBuffer[0];
			destinationAddress[i + 1] = copyLeftBuffer[1] | copyRightBuffer[1];
			destinationAddress[i + 2] = copyLeftBuffer[2] | copyRightBuffer[2];
			destinationAddress[i + 3] = copyLeftBuffer[3] | copyRightBuffer[3];
		}
	}

	if (copyLengthInBit % BITS_IN_WORD > 0) {
		destinationAddress[copyWordLength - 1] = truncateRight(destinationAddress[copyWordLength - 1], 
													BITS_IN_WORD - (copyLengthInBit % BITS_IN_WORD));
	}

}

void bitCopyDestWordOffsetOnly(VAL *destinationAddress, VAL destinationWordOffset,
							const VAL *sourceAddress, VAL sourceBitOffset, VAL copyLengthInBit) {

	VAL i;
	VAL rightShift;
	VAL copyLeftBuffer[4], copyRightBuffer[4];
	VAL copyWordLength, copyWordLengthRoundTo4, wordToNext4WordBoundary;
	VAL *destAddr;
	const VAL *srcAddr;

	#ifdef DEBUG
	if (copyLengthInBit == 0) {
		fprintf(stderr, "bitCopyDestWordOffsetOnly() : copyLengthInBit = 0!\n");
		exit(1);
	}
	#endif

	copyWordLength = (copyLengthInBit + BITS_IN_WORD_MINUS_1) / BITS_IN_WORD;
	destAddr = destinationAddress + destinationWordOffset;
	srcAddr = sourceAddress;

	wordToNext4WordBoundary = (FOUR_MULTIPLE_MASK - destinationWordOffset) % 4;

    if (sourceBitOffset == 0) {
		memcpy(destAddr, srcAddr, copyWordLength * 4); 
	} else {
		rightShift = BITS_IN_WORD - sourceBitOffset;
		for (i=0; i<wordToNext4WordBoundary; i++) {
			destAddr[i] = (srcAddr[i] << sourceBitOffset) | 
						  (srcAddr[i+1] >> rightShift);
		}
		destAddr += wordToNext4WordBoundary;
		srcAddr += wordToNext4WordBoundary;
		copyWordLengthRoundTo4 = (copyWordLength - wordToNext4WordBoundary + 3) & FOUR_MULTIPLE_MASK;
		for (i=0; i<copyWordLengthRoundTo4; i+=4) {
			// This is supposed to generate SSE2 codes
			// Need to check rewrite using intrinics is necessary
			copyLeftBuffer[0] = srcAddr[i + 0] << sourceBitOffset;
			copyLeftBuffer[1] = srcAddr[i + 1] << sourceBitOffset;
			copyLeftBuffer[2] = srcAddr[i + 2] << sourceBitOffset;
			copyLeftBuffer[3] = srcAddr[i + 3] << sourceBitOffset;
			copyRightBuffer[0] = srcAddr[i + 1] >> rightShift;
			copyRightBuffer[1] = srcAddr[i + 2] >> rightShift;
			copyRightBuffer[2] = srcAddr[i + 3] >> rightShift;
			copyRightBuffer[3] = srcAddr[i + 4] >> rightShift;
			destAddr[i + 0] = copyLeftBuffer[0] | copyRightBuffer[0];
			destAddr[i + 1] = copyLeftBuffer[1] | copyRightBuffer[1];
			destAddr[i + 2] = copyLeftBuffer[2] | copyRightBuffer[2];
			destAddr[i + 3] = copyLeftBuffer[3] | copyRightBuffer[3];
		}
	}

	if (copyLengthInBit % BITS_IN_WORD > 0) {
		destinationAddress[copyWordLength - 1] = truncateRight(destinationAddress[copyWordLength - 1], 
													BITS_IN_WORD - (copyLengthInBit % BITS_IN_WORD));
	}

}

// The remaining bits in destinationAddress, if destinationBitOffset > 0, must be cleared as 0

VAL bitCopy(VAL *destinationAddress, VAL destinationWordOffset, VAL destinationBitOffset,
			 const VAL *sourceAddress, VAL sourceBitOffset, VAL copyLengthInBit) {

	VAL i;
	VAL rightShift;
	VAL copyLeftBuffer[4], copyRightBuffer[4];
	VAL copyWordLength, copyWordLengthRoundTo4, wordToNext4WordBoundary;
	VAL *destAddr;
	const VAL *srcAddr;

	#ifdef DEBUG
	if (copyLengthInBit == 0) {
		fprintf(stderr, "bitCopy() : copyLengthInBit = 0!\n");
		exit(1);
	}
	#endif

	destAddr = destinationAddress + destinationWordOffset;
	srcAddr = sourceAddress;

	if (destinationBitOffset > 0) {
		destAddr[0] = destAddr[0] | (srcAddr[0] << sourceBitOffset >> destinationBitOffset);
		if (destinationBitOffset < sourceBitOffset) {
			destAddr[0] = destAddr[0] |
						(srcAddr[1] >> destinationBitOffset >> (BITS_IN_WORD - sourceBitOffset));
		}
		if (copyLengthInBit > BITS_IN_WORD - destinationBitOffset) {
			destAddr++;
			srcAddr += (sourceBitOffset + BITS_IN_WORD - destinationBitOffset) / BITS_IN_WORD;
			sourceBitOffset = (sourceBitOffset + BITS_IN_WORD - destinationBitOffset) % BITS_IN_WORD;
			copyLengthInBit -= BITS_IN_WORD - destinationBitOffset;
			destinationWordOffset++;
		} else {
			if ((destinationBitOffset + copyLengthInBit) % BITS_IN_WORD > 0) {
				destAddr[0] = truncateRight(destAddr[0], BITS_IN_WORD - destinationBitOffset - copyLengthInBit);
			}
			return 0;
		}
	}

	copyWordLength = (copyLengthInBit + BITS_IN_WORD_MINUS_1) / BITS_IN_WORD;

	if (sourceBitOffset == 0) {
		memcpy(destAddr, srcAddr, copyWordLength * 4); 
	} else {
		wordToNext4WordBoundary = (FOUR_MULTIPLE_MASK - destinationWordOffset) % 4;
		rightShift = BITS_IN_WORD - sourceBitOffset;
		for (i=0; i<wordToNext4WordBoundary; i++) {
			destAddr[i] = (srcAddr[i] << sourceBitOffset) | 
						(srcAddr[i+1] >> rightShift);
		}
		if (wordToNext4WordBoundary >= copyWordLength) {
			if (copyLengthInBit % BITS_IN_WORD > 0) {
				destAddr[copyWordLength - 1] = truncateRight(destAddr[copyWordLength - 1], 
								BITS_IN_WORD - (copyLengthInBit % BITS_IN_WORD));
			}
			return 0;
		}
		destAddr += wordToNext4WordBoundary;
		srcAddr += wordToNext4WordBoundary;
		copyWordLength -= wordToNext4WordBoundary;
		copyWordLengthRoundTo4 = (copyWordLength + 3) & FOUR_MULTIPLE_MASK;
		for (i=0; i<copyWordLengthRoundTo4; i+=4) {
			// This is supposed to generate SSE2 codes
			// Need to check rewrite using intrinics is necessary
			copyLeftBuffer[0] = srcAddr[i + 0] << sourceBitOffset;
			copyLeftBuffer[1] = srcAddr[i + 1] << sourceBitOffset;
			copyLeftBuffer[2] = srcAddr[i + 2] << sourceBitOffset;
			copyLeftBuffer[3] = srcAddr[i + 3] << sourceBitOffset;
			copyRightBuffer[0] = srcAddr[i + 1] >> rightShift;
			copyRightBuffer[1] = srcAddr[i + 2] >> rightShift;
			copyRightBuffer[2] = srcAddr[i + 3] >> rightShift;
			copyRightBuffer[3] = srcAddr[i + 4] >> rightShift;
			destAddr[i + 0] = copyLeftBuffer[0] | copyRightBuffer[0];
			destAddr[i + 1] = copyLeftBuffer[1] | copyRightBuffer[1];
			destAddr[i + 2] = copyLeftBuffer[2] | copyRightBuffer[2];
			destAddr[i + 3] = copyLeftBuffer[3] | copyRightBuffer[3];
		}
	}

	if (copyLengthInBit % BITS_IN_WORD > 0) {
		destAddr[copyWordLength - 1] = truncateRight(destAddr[copyWordLength - 1], 
							BITS_IN_WORD - (copyLengthInBit % BITS_IN_WORD));
	}

	return 0;

}
