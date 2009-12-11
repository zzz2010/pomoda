/* TextConverter.h		Text Converter

   Copyright 2004, Wong Chi Kwong, all rights reserved.

   This module contains miscellaneous text conversion functions.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/

#ifndef __TEXTCONVERTOR_H__
#define __TEXTCONVERTOR_H__

#include "TypeNLimit.h"
#include "MemManager.h"

#define INVALID_CHAR 0xFF
#define CHAR_MAP_SIZE 256
#define PACKED_BUFFER_SIZE			(PACKED_BUFFER_SIZE_IN_WORD * BYTES_IN_WORD)
#define PACKED_BUFFER_SIZE_IN_WORD	65536

// charMap is a char array of size 256. The index of the array is the input text value
// and the content of the array is the output text value. e.g. A -> 0, C -> 1
// If the value of an entry = INVALID_CHAR, the indexed text value is an invalid input

// Character map functions
VAL ReadCharMap(CHAR *charMap, const char *inputFileName, const CHAR defaultMapping);
void GenerateReverseCharMap(const CHAR *charMap, CHAR *reverseCharMap);
void GenerateCharMapACGTN(CHAR *charMap);

// Word packed text functions
VAL BitPerWordPackedChar(const VAL alphabetSize);
VAL TextLengthFromWordPacked(VAL wordPackedLength, VAL bitPerChar, VAL lastWordLength);
VAL WordPackedLengthFromText(VAL textLength, VAL bitPerChar);
VAL LastWordLength(VAL textLength, VAL bitPerChar);

// Byte packed text functions
VAL BitPerBytePackedChar(const VAL alphabetSize);
VAL TextLengthFromBytePacked(VAL bytePackedLength, VAL bitPerChar, VAL lastByteLength);
VAL BytePackedLengthFromText(VAL textLength, VAL bitPerChar);
CHAR LastByteLength(VAL textLength, VAL bitPerChar);

// Conversion functions
void ConvertTextToWordPacked(const CHAR *input, VAL *output, const CHAR *charMap, const VAL alphabetSize, const VAL textLength);
void ConvertTextToBytePacked(const CHAR *input, CHAR *output, const CHAR *charMap, const VAL alphabetSize, const VAL textLength);
void ConvertWordPackedToText(const VAL *input, CHAR *output, const CHAR *reverseCharMap, const VAL alphabetSize, const VAL textLength);
void ConvertBytePackedToText(const CHAR *input, CHAR *output, const CHAR *reverseCharMap, const VAL alphabetSize, const VAL textLength);
void ConvertWordPackedToBytePacked(const VAL *input, CHAR *output, const VAL alphabetSize, const VAL textLength);
void ConvertBytePackedToWordPacked(const CHAR *input, VAL *output, const VAL alphabetSize, const VAL textLength);
void ConvertTextToCode(const CHAR *input, CHAR *output, const CHAR *charMap, const VAL textLength);
void ConvertCodeToText(const CHAR *input, CHAR *output, const CHAR *reverseCharMap, const VAL textLength);

// Full load function
VAL ReadTextAsWordPacked(const char *inputFileName, const CHAR *charMap, const VAL alphabetSize, VAL *targetAddress, const VAL maxTextLength);
VAL ReadBytePackedAsWordPacked(const char *inputFileName, const VAL alphabetSize, VAL *targetAddress, const VAL maxTextLength);

// Save functions
void SaveText(const char *outputFileName, const CHAR *text, const VAL textLength);
void SaveBytePacked(const char *outputFileName, const CHAR *wordPacked, const VAL textLength, const VAL alphabetSize);
void SaveWordPacked(const char *outputFileName, const VAL *wordPacked, const VAL textLength, const VAL alphabetSize);

// Incremental load functions (start from end of text)
FILE *InitialLoadPackedIncFromEnd(const char* inputFileName, CHAR *packedOutput, const VAL alphabetSize, const VAL packedLengthPerLoad, VAL *textLength, VAL *textLengthForThisLoad);
void LoadPackedIncFromEnd(FILE *packedFile, CHAR *packedOutput, const VAL packedLengthPerLoad);
FILE *InitialLoadTextIncFromEnd(const char* inputFileName, CHAR *textOutput, const VAL textLengthPerLoad, VAL *textLength, VAL *textLengthForThisLoad);
void LoadTextIncFromEnd(FILE *textFile, CHAR *textOutput, const VAL textLengthPerLoad);

#endif
