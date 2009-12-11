/* TextConverter.c		Text Converter

   Copyright 2004, Wong Chi Kwong, all rights reserved.

   This module contains miscellaneous text conversion functions.

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

#include "TextConverter.h"
#include "MiscUtilities.h"

/*******************************************************************************
 * Generate the character map for dna sequence {a,c,g,t,n}
 * -> remove the need to read in the external character map file
 * -> Added by Ooi Hong Sain 12/01/2005
 *
 * Changed 27/02/2005
 * Important: changed the initialize charMap[i] to INVALID_CHAR. Previously, it
 * was initialized to 0, hence during the FMIDecompressText, the string will 
 * be terminated when character {a} was found.
 ******************************************************************************/
void GenerateCharMapACGTN(CHAR *charMap) {
	int i;

	for (i=0; i<CHAR_MAP_SIZE; i++) {
		charMap[i] = INVALID_CHAR;
	}

	charMap['a'] = (char)0;
	charMap['c'] = (char)1;
	charMap['g'] = (char)2;
	charMap['t'] = (char)3;
	charMap['n'] = (char)4;
	charMap['A'] = (char)0;
	charMap['C'] = (char)1;
	charMap['G'] = (char)2;
	charMap['T'] = (char)3;
	charMap['N'] = (char)4;
}

VAL ReadCharMap(CHAR *charMap, const char *inputFileName, const CHAR defaultMapping) {

	FILE *inputFile;
	char c;
	VAL v, alphabetSize;

	inputFile = fopen(inputFileName, "r");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadCharMap() : Cannot open character map!\n");
		exit(1);
	}

	for (v=0; v<CHAR_MAP_SIZE; v++) {
		charMap[v] = defaultMapping;
	}

	alphabetSize = 0;

	while (!feof(inputFile)) {
		fscanf(inputFile, " %c %u \n", &c, &v);
		if (v > CHAR_MAP_SIZE) {
			fprintf(stderr, "ReadCharMap() : Invalid charMap!\n");
			return 0;
		}
		charMap[c] = (char)v;
		if (v > alphabetSize) {
			alphabetSize = v;
		}
	}

	fclose(inputFile);

	alphabetSize++;

	return alphabetSize;

}

void GenerateReverseCharMap(const CHAR *charMap, CHAR *reverseCharMap) {

	VAL i, j;

	for (i=0; i<CHAR_MAP_SIZE; i++) {
		reverseCharMap[i] = INVALID_CHAR;
		for (j=0; j<CHAR_MAP_SIZE; j++) {
			if (charMap[j] == i) {
				reverseCharMap[i] = (CHAR)j;
				break;
			}
		}
	}

}

VAL BitPerWordPackedChar(const VAL alphabetSize) {

	VAL bitPerChar;

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerWordPackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif

	bitPerChar = ceilLog2(alphabetSize);

	#ifdef DEBUG
	if (bitPerChar > BITS_IN_WORD) {
		fprintf(stderr, "BitPerWordPackedChar() : bitPerChar > BITS_IN_WORD!\n");
		exit(1);
	}
	#endif

	// Return the largest number of bit that does not affect packing efficiency
	if (BITS_IN_WORD / (BITS_IN_WORD / bitPerChar) > bitPerChar) {
		bitPerChar = BITS_IN_WORD / (BITS_IN_WORD / bitPerChar);
	}
	return bitPerChar;

}

VAL TextLengthFromWordPacked(VAL wordPackedLength, VAL bitPerChar, VAL lastWordLength) {

	return (wordPackedLength - 1) * (BITS_IN_WORD / bitPerChar) + lastWordLength;

}

VAL WordPackedLengthFromText(VAL textLength, VAL bitPerChar) {

	return (textLength + (BITS_IN_WORD / bitPerChar) - 1) / (BITS_IN_WORD / bitPerChar);

}

VAL LastWordLength(VAL textLength, VAL bitPerChar) {

	return textLength % (BITS_IN_WORD / bitPerChar);

}

VAL BitPerBytePackedChar(const VAL alphabetSize) {

	VAL bitPerChar;

	#ifdef DEBUG
	if (alphabetSize < 2) {
		fprintf(stderr, "BitPerBytePackedChar() : alphabetSize < 2!\n");
		exit(1);
	}
	#endif

	bitPerChar = ceilLog2(alphabetSize);

	#ifdef DEBUG
	if (bitPerChar > BITS_IN_BYTE) {
		fprintf(stderr, "BitPerBytePackedChar() : bitPerChar > BITS_IN_BYTE!\n");
		exit(1);
	}
	#endif

	// Return the largest number of bit that does not affect packing efficiency
	if (BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar) > bitPerChar) {
		bitPerChar = BITS_IN_BYTE / (BITS_IN_BYTE / bitPerChar);
	}
	return bitPerChar;
}

VAL TextLengthFromBytePacked(VAL bytePackedLength, VAL bitPerChar, VAL lastByteLength) {

	return (bytePackedLength - 1) * (BITS_IN_BYTE / bitPerChar) + lastByteLength;

}

VAL BytePackedLengthFromText(VAL textLength, VAL bitPerChar) {

	return (textLength + (BITS_IN_BYTE / bitPerChar) - 1) / (BITS_IN_BYTE / bitPerChar);

}

CHAR LastByteLength(VAL textLength, VAL bitPerChar) {

	return (CHAR)(textLength % (BITS_IN_BYTE / bitPerChar));

}

void ConvertTextToWordPacked(const CHAR *input, VAL *output, const CHAR *charMap, const VAL alphabetSize, const VAL textLength) {

	VAL bitPerChar, charPerWord;
	VAL i, j, k;
	VAL c;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = 0;
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			c = c | (charMap[input[j+k]] << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerWord < textLength) {
		c = 0;
		j = i * charPerWord;
		for (k=0; j+k < textLength; k++) {
			c = c | (charMap[input[j+k]] << (BITS_IN_WORD - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertTextToBytePacked(const CHAR *input, CHAR *output, const CHAR *charMap, const VAL alphabetSize, const VAL textLength) {

	VAL bitPerChar, charPerByte;
	VAL i, j, k;
	CHAR c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = 0;
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			c = c | (charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}
	if (i * charPerByte < textLength) {
		c = 0;
		j = i * charPerByte;
		for (k=0; j+k < textLength; k++) {
			c = c | (charMap[input[j+k]] << (BITS_IN_BYTE - (k+1) * bitPerChar));
		}
		output[i] = c;
	}

}

void ConvertWordPackedToText(const VAL *input, CHAR *output, const CHAR *reverseCharMap, const VAL alphabetSize, const VAL textLength) {

	VAL bitPerChar, charPerWord;
	VAL i, j, k;
	VAL c;

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	for (i=0; i<textLength/charPerWord; i++) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; k<charPerWord; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerWord < textLength) {
		c = input[i];
		j = i * charPerWord;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_WORD - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertBytePackedToText(const CHAR *input, CHAR *output, const CHAR *reverseCharMap, const VAL alphabetSize, const VAL textLength) {

	VAL bitPerChar, charPerByte;
	VAL i, j, k;
	CHAR c;

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	for (i=0; i<textLength/charPerByte; i++) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; k<charPerByte; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}
	if (i * charPerByte < textLength) {
		c = input[i];
		j = i * charPerByte;
		for (k=0; j+k<textLength; k++) {
			output[j+k] = reverseCharMap[c >> (BITS_IN_BYTE - bitPerChar)];
			c <<= bitPerChar;
		}
	}

}

void ConvertWordPackedToBytePacked(const VAL *input, CHAR *output, const VAL alphabetSize, const VAL textLength) {

	VAL i, j, k;
	VAL c;
	VAL bitPerBytePackedChar;
	VAL bitPerWordPackedChar;
	VAL charPerWord;
	VAL charPerByte;
	VAL bytePerIteration;
	VAL byteProcessed = 0;
	VAL wordProcessed = 0;
	VAL mask, shift;
	
	VAL buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerBytePackedChar;
	charPerByte = BITS_IN_BYTE / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		c = input[wordProcessed];
		for (i=0; i<charPerWord; i++) {
			buffer[i] = c >> shift;
			c <<= bitPerWordPackedChar;
		}
		wordProcessed++;

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = 0;
			for (j=0; j<charPerByte; j++) {
				c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
				k++;
			}
			output[byteProcessed] = (CHAR)c;
			byteProcessed++;
		}

	}

	c = input[wordProcessed];
	for (i=0; i < textLength - wordProcessed * charPerWord; i++) {
		buffer[i] = c >> shift;
		c <<= bitPerWordPackedChar;
	}

	k = 0;
	while (byteProcessed * charPerByte < textLength) {
		c = 0;
		for (j=0; j < textLength - wordProcessed * charPerWord; j++) {
			c |= buffer[k] << (BITS_IN_BYTE - (j+1) * bitPerBytePackedChar);
			k++;
		}
		output[byteProcessed] = (CHAR)c;
		byteProcessed++;
	}

}

void ConvertBytePackedToWordPacked(const CHAR *input, VAL *output, const VAL alphabetSize, const VAL textLength) {

	VAL i, j, k;
	VAL c;
	VAL bitPerBytePackedChar;
	VAL bitPerWordPackedChar;
	VAL charPerWord;
	VAL charPerByte;
	VAL bytePerIteration;
	VAL byteProcessed = 0;
	VAL wordProcessed = 0;
	VAL mask, shift;
	
	VAL buffer[BITS_IN_WORD];

	bitPerBytePackedChar = BitPerBytePackedChar(alphabetSize);
	bitPerWordPackedChar = BitPerWordPackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerBytePackedChar;
	charPerWord = BITS_IN_WORD / bitPerWordPackedChar;

	bytePerIteration = charPerWord / charPerByte;
	mask = truncateRight(ALL_ONE_MASK, BITS_IN_WORD - bitPerWordPackedChar);
	shift = BITS_IN_WORD - BITS_IN_BYTE + bitPerBytePackedChar - bitPerWordPackedChar;

	while ((wordProcessed + 1) * charPerWord < textLength) {

		k = 0;
		for (i=0; i<bytePerIteration; i++) {
			c = (VAL)input[byteProcessed] << shift;
			for (j=0; j<charPerByte; j++) {
				buffer[k] = c & mask;
				c <<= bitPerBytePackedChar;
				k++;
			}
			byteProcessed++;
		}

		c = 0;
		for (i=0; i<charPerWord; i++) {
			c |= buffer[i] >> bitPerWordPackedChar * i;
		}
		output[wordProcessed] = c;
		wordProcessed++;

	}

	k = 0;
	for (i=0; i < (textLength - wordProcessed * charPerWord - 1) / charPerByte + 1; i++) {
		c = (VAL)input[byteProcessed] << shift;
		for (j=0; j<charPerByte; j++) {
			buffer[k] = c & mask;
			c <<= bitPerBytePackedChar;
			k++;
		}
		byteProcessed++;
	}

	c = 0;
	for (i=0; i<textLength - wordProcessed * charPerWord; i++) {
		c |= buffer[i] >> bitPerWordPackedChar * i;
	}
	output[wordProcessed] = c;

}

void ConvertTextToCode(const CHAR *input, CHAR *output, const CHAR *charMap, const VAL textLength) {

	VAL i;

	for (i=0; i< textLength; i++) {
		output[i] = charMap[input[i]];
	}

}

void ConvertCodeToText(const CHAR *input, CHAR *output, const CHAR *reverseCharMap, const VAL textLength) {

	VAL i;

	for (i=0; i< textLength; i++) {
		output[i] = reverseCharMap[input[i]];
	}

}

VAL ReadTextAsWordPacked(const char *inputFileName, const CHAR *charMap, const VAL alphabetSize, VAL *targetAddress, const VAL maxTextLength) {

	FILE *inputFile;
	CHAR *buffer;
	VAL charPerWord;
	VAL charRead;
	VAL charProcessed = 0, wordProcessed = 0;
	VAL charPerBuffer;

	inputFile = fopen(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadTextAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	charPerBuffer = PACKED_BUFFER_SIZE / charPerWord * charPerWord;

	buffer = (unsigned char*)MMUnitAllocate(charPerBuffer);

	charRead = (VAL)fread(buffer, 1, charPerBuffer, inputFile);
	while (charRead > 0 && charProcessed + charRead < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, charRead);
		wordProcessed += charRead / charPerWord;
		charProcessed += charRead;
		charRead = (VAL)fread(buffer, 1, charPerBuffer, inputFile);
	}

	if (charRead > 0 && charProcessed < maxTextLength) {
		ConvertTextToWordPacked(buffer, targetAddress + wordProcessed, charMap, alphabetSize, min(charRead, maxTextLength - charProcessed));
		charProcessed += charRead;
	}

	MMUnitFree(buffer, charPerBuffer);

	fclose(inputFile);

	return charProcessed;

}

VAL ReadBytePackedAsWordPacked(const char *inputFileName, const VAL alphabetSize, VAL *targetAddress, const VAL maxTextLength) {

	FILE *inputFile;
	CHAR *buffer1, *buffer2;
	VAL charPerByte, charPerWord;
	VAL charPerBuffer, wordPerBuffer;
	VAL charProcessed = 0, wordProcessed = 0;
	VAL byteRead, tempByteRead;
	VAL charInLastBuffer;
	VAL bufferSize;

	inputFile = fopen(inputFileName, "rb");

	if (inputFile == NULL) {
		fprintf(stderr, "ReadBytePackedAsWordPacked() : Cannot open inputFileName!\n");
		exit(1);
	}

	charPerByte = BITS_IN_BYTE / BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / BitPerWordPackedChar(alphabetSize);
	bufferSize = PACKED_BUFFER_SIZE / charPerByte / charPerWord * charPerByte * charPerWord;

	charPerBuffer = bufferSize * charPerByte;
	wordPerBuffer = charPerBuffer / charPerWord;

	buffer1 = (unsigned char*)MMUnitAllocate(bufferSize);
	buffer2 = (unsigned char*)MMUnitAllocate(bufferSize);

	byteRead = (VAL)fread(buffer1, 1, bufferSize, inputFile);
	tempByteRead = (VAL)fread(buffer2, 1, bufferSize, inputFile);

	while (tempByteRead > 1 && charProcessed + charPerBuffer < maxTextLength) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, charPerBuffer);
		charProcessed += charPerBuffer;
		wordProcessed += wordPerBuffer;
		memcpy(buffer1, buffer2, bufferSize);
		byteRead = tempByteRead;
		tempByteRead = (VAL)fread(buffer2, 1, bufferSize, inputFile);
	}

	if (tempByteRead > 1) {
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, maxTextLength - charProcessed);
		charProcessed += charPerBuffer;
	} else {
		if (tempByteRead == 1) {
			charInLastBuffer = charPerBuffer - charPerByte + buffer2[0];
		} else {
			charInLastBuffer = (byteRead - 2) * charPerByte + buffer1[byteRead - 1];
		}
		ConvertBytePackedToWordPacked(buffer1, targetAddress + wordProcessed, alphabetSize, min(maxTextLength - charProcessed, charInLastBuffer));
		charProcessed += charInLastBuffer;
	}

	MMUnitFree(buffer1, bufferSize);
	MMUnitFree(buffer2, bufferSize);

	fclose(inputFile);

	return charProcessed;

}

void SaveText(const char *outputFileName, const CHAR *text, const VAL textLength) {

	FILE *outputFile;

	outputFile = fopen(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveText() : Cannot open output file!\n");
		exit(1);
	}

	fwrite(text, sizeof(CHAR), textLength, outputFile);
	fclose(outputFile);

}

void SaveBytePacked(const char *outputFileName, const CHAR *bytePacked, const VAL textLength, const VAL alphabetSize) {

	FILE *outputFile;
	VAL bitPerChar, charPerByte, bytePackedLen;
	CHAR lastByteLen;
	CHAR zero = 0;

	outputFile = fopen(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveBytePacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerByte = BITS_IN_BYTE / bitPerChar;

	bytePackedLen = BytePackedLengthFromText(textLength, bitPerChar);
	lastByteLen = LastByteLength(textLength, bitPerChar);

	fwrite(bytePacked, sizeof(CHAR), bytePackedLen, outputFile);
	if (lastByteLen == 0) {
		fwrite(&zero, sizeof(CHAR), 1, outputFile);
	}
	fwrite(&lastByteLen, sizeof(CHAR), 1, outputFile);
	fclose(outputFile);

}

void SaveWordPacked(const char *outputFileName, const VAL *wordPacked, const VAL textLength, const VAL alphabetSize) {

	FILE *outputFile;
	VAL bitPerChar, charPerWord, wordPackedLen;
	VAL lastWordLen;
	VAL zero = 0;

	outputFile = fopen(outputFileName, "wb");

	if (outputFile == NULL) {
		fprintf(stderr, "SaveWordPacked() : Cannot open output file!\n");
		exit(1);
	}

	bitPerChar = BitPerWordPackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	wordPackedLen = WordPackedLengthFromText(textLength, bitPerChar);
	lastWordLen = LastWordLength(textLength, bitPerChar);

	fwrite(wordPacked, sizeof(VAL), wordPackedLen, outputFile);
	if (lastWordLen == 0) {
		fwrite(&zero, sizeof(VAL), 1, outputFile);
	}
	fwrite(&lastWordLen, sizeof(VAL), 1, outputFile);
	fclose(outputFile);

}

FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, CHAR *packedOutput, const VAL alphabetSize, 
								  const VAL packedLengthPerLoad, VAL *textLength, VAL *textLengthForThisLoad) {

	FILE *packedFile;
	VAL len, packedFileLenForThisLoad, packedFileLen;
	CHAR lastByteLength;
	VAL bitPerChar, charPerWord;

	packedFile = fopen(inputFileName, "rb");

	if (packedFile == NULL) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	bitPerChar = BitPerBytePackedChar(alphabetSize);
	charPerWord = BITS_IN_WORD / bitPerChar;

	fseek(packedFile, -1, SEEK_END);
	packedFileLen = ftell(packedFile);
	if ((int)packedFileLen < 0) {
		fprintf(stderr, "InitialLoadPackedIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}
	fread(&lastByteLength, sizeof(CHAR), 1, packedFile);

	len = TextLengthFromBytePacked(packedFileLen, bitPerChar, lastByteLength);

	if (lastByteLength == 0 && (packedFileLen - 1) % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = 0;
		fseek(packedFile, -((int)(2+packedLengthPerLoad)), SEEK_END);
		*textLength = len;
		*textLengthForThisLoad = 0;
		return packedFile;
	}

	if (packedFileLen % packedLengthPerLoad == 0) {
		packedFileLenForThisLoad = packedLengthPerLoad;
	} else {
		packedFileLenForThisLoad = packedFileLen % packedLengthPerLoad;
	}
	fseek(packedFile, -1, SEEK_END);

	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	fread(packedOutput, sizeof(CHAR), packedFileLenForThisLoad, packedFile);
	fseek(packedFile, -((int)packedFileLenForThisLoad), SEEK_CUR);
	if (packedFileLen > packedFileLenForThisLoad) {
		fseek(packedFile, -((int)packedLengthPerLoad), SEEK_CUR);
	}

	*textLength = len;
	*textLengthForThisLoad = TextLengthFromBytePacked(packedFileLenForThisLoad, bitPerChar, lastByteLength);

	return packedFile;

}

void LoadPackedIncFromEnd(FILE *packedFile, CHAR *packedOutput, const VAL packedLengthPerLoad) {
	
	fread(packedOutput, sizeof(CHAR), packedLengthPerLoad, packedFile);
	fseek(packedFile, -(2*(int)packedLengthPerLoad), SEEK_CUR);

}


FILE *InitialLoadTextIncFromEnd(const char* inputFileName, CHAR *textOutput, const VAL textLengthPerLoad, VAL *textLength, VAL *textLengthForThisLoad) {

	FILE *textFile;
	VAL len, textLenForThisLoad;

	textFile = fopen(inputFileName, "rb");

	if (textFile == NULL) {
		fprintf(stderr, "InitialLoadTextIncFromEnd() : Cannot open inputFileName!\n");
		exit(1);
	}

	fseek(textFile, 0, SEEK_END);
	len = ftell(textFile);
	if ((int)len < 0) {
		fprintf(stderr, "InitialLoadTextIncFromEnd(): Cannot determine file length!\n");
		exit(1);
	}

	textLenForThisLoad = len % textLengthPerLoad;

	if (textLenForThisLoad > 0) {
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
		fread(textOutput, sizeof(CHAR), textLenForThisLoad, textFile);
		fseek(textFile, -((int)textLenForThisLoad), SEEK_END);
	}

	*textLength = len;
	*textLengthForThisLoad = textLenForThisLoad;

	return textFile;
}

void LoadTextIncFromEnd(FILE *textFile, CHAR *textOutput, const VAL textLengthPerLoad) {

	if (ftell(textFile) < (int)textLengthPerLoad) {
		fprintf(stderr, "LoadTextIncFromEnd(): file pointer is not correctly placed!\n");
		exit(1);
	}

	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);
	fread(textOutput, sizeof(CHAR), textLengthPerLoad, textFile);
	fseek(textFile, -((int)textLengthPerLoad), SEEK_CUR);

}
