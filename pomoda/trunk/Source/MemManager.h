/* MemManager.h		Memory Manager

   Copyright 2004, Wong Chi Kwong, all rights reserved.

   This module provides memory management functions.

   Four types of memory allocation types are provided:
   unit, pool, temporary, and bulk. See Memory type remarks below.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/

#ifndef __MEM_MANAGER_H__
#define __MEM_MANAGER_H__

#include "TypeNLimit.h"

#define MAX_ALIGN	16
#define MIN_ALIGN	1

#define RECORD_GRAND_TOTAL

//	Memory type:
//
//		unit memory:	allocation managed by malloc() individually;
//						to be used for large and less frequently accessed items
//						allocation can be freed individually at any time
//		pool memory:	pre-allocated memory pool for items with varying sizes
//						allocation cannot be freed by individually
//						to be used for small and frequently accessed items
//		temp memory:	temporary use granted from pool memory
//						allocation is allocated and freed like the items in a stack
//						pool memory allocation is disabled while temporary memory is in use
//		bulk memory:	pre-allocated memory pool for items with the same size
//						to be used for massively numbered items
//						memory address of dispatched items can be calculated by dispatch index


#ifdef DEBUG
#define Mem(mmBulk, index)  MMBulkAddress(mmBulk, index)
#else
#define Mem(mmBulk, index) 	(void*)&(mmBulk->directory[index >> mmBulk->itemPerAllocationInPowerOf2][(index & mmBulk->indexMask) * mmBulk->itemSize])
#endif

typedef struct MMPool {
	VAL poolSize;
	VAL poolByteDispatched;				// Includes any spillover and memory skipped for align
	VAL poolByteSkippedForAlign;
	VAL poolByteSpillover;				// Exclude spillover pointers
	void *firstSpillOverAddress;		// if pool is freed, = address of mmPool
	VAL currentTempByteDispatched;		// Includes any spillover
	VAL currentTempByteSpillover;		// Exclude spillover pointers
	VAL maxTotalByteDispatched;			// The max of pool memory + temp memory dispatched
} MMPool;


typedef struct MMBulk {
	VAL itemSize;
	VAL itemPerAllocationInPowerOf2;
	VAL boundaryCushionSize;			// boundary cushion is a piece of memory allocated so that the memory around items can be safely referenced
	VAL indexMask;
	VAL currentDirectoryEntry;
	VAL nextUnusedItem;
	VAL directorySize;
	CHAR **directory;			// if bulk is freed, = NULL
} MMBulk;

typedef struct MMMaster {
	VAL currentUnitByteAllocated;
	VAL maxUnitByteAllocated;
	VAL maxNumberOfPools;
	MMPool **mmPool;
	VAL maxNumberOfBulks;
	MMBulk **mmBulk;
	VAL maxTotalByteAllocated;
	VAL maxTotalByteDispatched;
} MMMaster;

void *MMMalloc(const VAL memSize);
void MMFree(void *address);
void MMMasterInitialize(const VAL maxNumberOfPools, const VAL maxNumberOfBulks);
void MMMasterFreeAll();
VAL MMMasterCurrentTotalByteAllocated();
VAL MMMasterCurrentTotalByteDispatched();
VAL MMMasterMaxTotalByteAllocated();
VAL MMMasterMaxTotalByteDispatched();
void MMMasterSetMaxTotalByteAllocated();
void MMMasterSetMaxTotalByteDispatched();
void MMMasterPrintReport(FILE *output, const VAL withUnitDetails, const VAL withPoolDetails, const VAL withBulkDetails);

void *MMUnitAllocate(const VAL memSize);
void MMUnitFree(void *address, const VAL memSize);
VAL MMUnitCurrentByteAllocated();
VAL MMUnitMaxByteAllocated();
void MMUnitPrintReport(FILE *output);

MMPool *MMPoolCreate(const VAL poolSize);
VAL MMPoolIsActive(const MMPool *mmPool);
void MMPoolSetInactive(MMPool *mmPool);
VAL MMPoolCurrentTotalByteAllocated(const MMPool *mmPool);
VAL MMPoolCurrentTotalByteDispatched(const MMPool *mmPool);
VAL MMPoolMaxTotalByteDispatched(const MMPool *mmPool);
MMPool *MMPoolFree(MMPool *mmPool);
void MMPoolDestory(MMPool *mmPool);
void *MMPoolDispatch(MMPool *mmPool, const VAL memSize);
void MMPoolReturn(MMPool *mmPool, void *address, const VAL memSize);		// Dummy function
void MMPoolPrintReport(MMPool *mmPool, FILE *output);

void *MMTempDispatch(MMPool *mmPool, const VAL memsize);
void MMTempReturn(MMPool *mmPool, void *address, const VAL memSize);
void MMTempPrintReport(MMPool *mmPool, FILE *output);

MMBulk *MMBulkCreate(MMPool *mmPool, const VAL itemSize, const VAL itemPerAllocationInPowerOf2, 
					 VAL const boundaryCushionSize, VAL const directorySize);
VAL MMBulkIsActive(const MMBulk *mmBulk);
void MMBulkSetInactive(MMBulk *mmBulk);
VAL MMBulkByteAllocated(const MMBulk *mmBulk);
VAL MMBulkByteDispatched(const MMBulk *mmBulk);
VAL MMBulkUnitDispatched(const MMBulk *mmBulk);
void MMBulkFree(MMBulk *mmBulk);
void MMBulkDestory(MMBulk *mmBulk);
VAL MMBulkDispatch(MMBulk *mmBulk);
void *MMBulkAddress(const MMBulk *mmBulk, const VAL index);
MMPool *MMBulkFindPoolUsed(const MMBulk *mmBulk);
void MMBulkPrintReport(MMBulk *mmBulk, FILE *output);

void MMBulkSave(MMBulk *mmBulk, FILE *output);
MMBulk *MMBulkLoad(MMPool *mmPool, FILE *input);


#endif
