#pragma once
#ifndef __TYPENLIMIT_H__
#define __TYPENLIMIT_H__

#include <limits.h>

#define BITS_IN_WORD 32
#define BITS_IN_WORD_MINUS_1 31
#define BITS_IN_WORD_MASK 0x0000001F
#define BITS_IN_WORD_SHIFT 5
#define BITS_IN_HALF_WORD 16
#define BITS_IN_4_WORD 128
#define BITS_IN_4_WORD_MINUS_1 127
#define BITS_IN_4_WORD_SHIFT 7
#define FIRST_BIT_MASK 0x80000000
#define ALL_BUT_FIRST_BIT_MASK 0x7FFFFFFF
#define ALL_ONE_MASK 0xFFFFFFFF
#define FOUR_MULTIPLE_MASK 0xFFFFFFFC
#define BITS_IN_BYTE 8
#define BITS_IN_BYTE_SHIFT 3
#define BYTES_IN_WORD 4
//#define VAL     VAL
#define CHAR	unsigned char
#define BOOL    VAL
#define TRUE    1
#define FALSE   0
#define __restrict__				// For compatibility of C compilers without support of __restrict__ keyword

#define MAX_FILENAME_LEN 256

#endif
