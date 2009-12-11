/* MemManager.c		Memory Manager

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stdafx.h"
#include "MiscUtilities.h"
#include "MemManager.h"

MMMaster mmMaster;

void *MMMalloc(const VAL memSize) {

	void *address;

	address = malloc(memSize);
	if (address == NULL) {
		fprintf(stderr, "MMMalloc() : cannot allocate memory!\n");
		exit(1);
	}
	return address;

}

void MMFree(void *address) {

	free(address);

}

void MMMasterInitialize(const VAL maxNumberOfPools, const VAL maxNumberOfBulks) {

	VAL i;

	mmMaster.maxTotalByteAllocated = 0;
	mmMaster.maxTotalByteDispatched = 0;
	mmMaster.currentUnitByteAllocated = 0;
	mmMaster.maxUnitByteAllocated = 0;

	mmMaster.maxNumberOfBulks = maxNumberOfBulks;
	mmMaster.maxNumberOfPools = maxNumberOfPools;
	if (maxNumberOfBulks > 0) {
		mmMaster.mmBulk = (MMBulk**)malloc(sizeof(MMBulk*) * maxNumberOfBulks);
		for (i=0; i<maxNumberOfBulks; i++) {
			mmMaster.mmBulk[i] = NULL;
		}
	} else {
		mmMaster.mmBulk = NULL;
	}
	if (maxNumberOfPools > 0) {
		mmMaster.mmPool = (MMPool**)malloc(sizeof(MMPool*) * maxNumberOfPools);
		for (i=0; i<maxNumberOfPools; i++) {
			mmMaster.mmPool[i] = NULL;
		}
	} else {
		mmMaster.mmPool = NULL;
	}

}

void MMMasterFreeAll() {

	VAL i;

	for (i=0; i < mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] != NULL) {
			if (MMBulkIsActive(mmMaster.mmBulk[i])) {
				MMBulkFree(mmMaster.mmBulk[i]);
			}
			if (MMBulkFindPoolUsed(mmMaster.mmBulk[i]) == NULL) {
				MMUnitFree(mmMaster.mmBulk[i], sizeof(MMBulk));
			}
			mmMaster.mmBulk[i] = NULL;
		}
	}
	free(mmMaster.mmBulk);

	for (i=0; i < mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] != NULL) {
			if (MMPoolIsActive(mmMaster.mmPool[i])) {
				MMPoolFree(mmMaster.mmPool[i]);
			}
			free(mmMaster.mmPool[i]);
			mmMaster.mmPool[i] = NULL;
		}
	}
	free(mmMaster.mmPool);

}

VAL MMMasterCurrentTotalByteAllocated() {

	VAL i;
	VAL currentTotalByteAllocated;

	// unit memory allocated
	currentTotalByteAllocated = mmMaster.currentUnitByteAllocated;

	// pool and temp memory allocated
	for (i=0; i < mmMaster.maxNumberOfPools; i++) {
        if (mmMaster.mmPool[i] != NULL && MMPoolIsActive(mmMaster.mmPool[i])) {
			currentTotalByteAllocated += MMPoolCurrentTotalByteAllocated(mmMaster.mmPool[i]);
		}
	}

	// bulk memory allocated
	for (i=0; i < mmMaster.maxNumberOfBulks; i++) {
        if (mmMaster.mmBulk[i] != NULL && MMBulkIsActive(mmMaster.mmBulk[i])) {
			currentTotalByteAllocated += MMBulkByteAllocated(mmMaster.mmBulk[i]);
		}
	}

	return currentTotalByteAllocated;

}

VAL MMMasterCurrentTotalByteDispatched() {

	VAL i;
	VAL currentTotalByteDispatched;

	// unit memory dispatched
	currentTotalByteDispatched = mmMaster.currentUnitByteAllocated;

	// pool and temp memory dispatched
	for (i=0; i < mmMaster.maxNumberOfPools; i++) {
        if (mmMaster.mmPool[i] != NULL && MMPoolIsActive(mmMaster.mmPool[i])) {
			currentTotalByteDispatched += MMPoolCurrentTotalByteDispatched(mmMaster.mmPool[i]);
		}
	}

	// bulk memory dispatched
	for (i=0; i < mmMaster.maxNumberOfBulks; i++) {
        if (mmMaster.mmBulk[i] != NULL && MMBulkIsActive(mmMaster.mmBulk[i])) {
			currentTotalByteDispatched += MMBulkByteDispatched(mmMaster.mmBulk[i]);
		}
	}

	return currentTotalByteDispatched;

}

VAL MMMasterMaxTotalByteAllocated() {

	VAL currentTotalByteAllocated;

	currentTotalByteAllocated = MMMasterCurrentTotalByteAllocated();

	if (currentTotalByteAllocated > mmMaster.maxTotalByteAllocated) {
		return currentTotalByteAllocated;
	} else {
		return mmMaster.maxTotalByteAllocated;
	}

}

VAL MMMasterMaxTotalByteDispatched() {

	VAL currentTotalByteDispatched ;

	currentTotalByteDispatched = MMMasterCurrentTotalByteDispatched();

	if (currentTotalByteDispatched > mmMaster.maxTotalByteDispatched) {
		return currentTotalByteDispatched;
	} else {
		return mmMaster.maxTotalByteDispatched;
	}

}

void MMMasterSetMaxTotalByteAllocated() {

	VAL currentTotalByteAllocated;
	
	currentTotalByteAllocated = MMMasterCurrentTotalByteAllocated();

	if (currentTotalByteAllocated > mmMaster.maxTotalByteAllocated) {
		mmMaster.maxTotalByteAllocated = currentTotalByteAllocated;
	}

}

void MMMasterSetMaxTotalByteDispatched() {

	VAL currentTotalByteDispatched;
	
	currentTotalByteDispatched = MMMasterCurrentTotalByteDispatched();

	if (currentTotalByteDispatched > mmMaster.maxTotalByteDispatched) {
		mmMaster.maxTotalByteDispatched = currentTotalByteDispatched;
	}

}

void MMMasterPrintReport(FILE *output, const VAL withUnitDetails, const VAL withPoolDetails, const VAL withBulkDetails) {

	VAL i;

	fprintf(output, "Maximum amount of memory allocated: %u\n", MMMasterMaxTotalByteAllocated());
	fprintf(output, "Maximum amount of memory dispatch:  %u\n", MMMasterMaxTotalByteDispatched());

	if (withUnitDetails) {
		fprintf(output, "\n");
		MMUnitPrintReport(output);
	}

	if (withPoolDetails) {
		for (i=0; i<mmMaster.maxNumberOfPools; i++) {
			if (mmMaster.mmPool[i] != NULL) {
				fprintf(output, "\nPool number %u\n", i);
				MMPoolPrintReport(mmMaster.mmPool[i], output);
			}
		}
	}

	if (withBulkDetails) {
		for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
			if (mmMaster.mmBulk[i] != NULL) {
				fprintf(output, "\nBulk number %u\n", i);
				MMBulkPrintReport(mmMaster.mmBulk[i], output);
			}
		}
	}

}

void *MMUnitAllocate(const VAL memSize) {

	void *temp;

	#ifdef DEBUG
	if (memSize == 0) {
		fprintf(stderr, "MMUnitAllocate() : memSize = 0!\n");
		exit(1);
	}
	#endif

	temp = malloc(memSize);
	if (temp == NULL) {
		fprintf(stderr, "MMUnitAllocate() : cannot allocate memory!\n");
		exit(1);
	}

	mmMaster.currentUnitByteAllocated += memSize;

	return temp;

}

void MMUnitFree(void *address, const VAL memSize) {

	#ifdef DEBUG
	if (address == NULL) {
		fprintf(stderr, "MMUnitFree() : address = NULL!\n");
		exit(1);
	}
	if (mmMaster.currentUnitByteAllocated < memSize) {
		fprintf(stderr, "MMUnitFree() : currentUnitByteAllocated < memSize!\n");
		exit(1);
	}
	#endif

	free(address);

	#ifdef RECORD_GRAND_TOTAL
	MMMasterSetMaxTotalByteAllocated();
	MMMasterSetMaxTotalByteDispatched();
	#endif

	if (mmMaster.currentUnitByteAllocated > mmMaster.maxUnitByteAllocated) {
		mmMaster.maxUnitByteAllocated = mmMaster.currentUnitByteAllocated;
	}
	mmMaster.currentUnitByteAllocated -= memSize;

}

VAL MMUnitCurrentByteAllocated() {

	return mmMaster.currentUnitByteAllocated;

}

VAL MMUnitMaxByteAllocated() {

	if (mmMaster.currentUnitByteAllocated > mmMaster.maxUnitByteAllocated) {
		return mmMaster.currentUnitByteAllocated;
	} else {
		return mmMaster.maxUnitByteAllocated;
	}

}

void MMUnitPrintReport(FILE *output) {

	fprintf(output, "Maximum amount of unit memory allocated: %u\n", MMUnitMaxByteAllocated());

}

MMPool *MMPoolCreate(const VAL poolSize) {

	MMPool *mmPool;
	VAL i;

	#ifdef DEBUG
	if (poolSize < sizeof(MMPool)) {
		fprintf(stderr, "MMPoolCreate() : poolSize < MMPool!\n");
		exit(1);
	}
	#endif

	mmPool = (MMPool*)malloc(poolSize);
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolCreate() : cannot allocate memory!\n");
		exit(1);
	}

	if (poolSize / MAX_ALIGN * MAX_ALIGN != poolSize) {
		fprintf(stderr, "MMPoolCreate() : poolSize must be multiple of MAX_ALIGN (%u)!\n", MAX_ALIGN);
		exit(1);
	}

	mmPool->poolSize = poolSize;
	mmPool->poolByteDispatched = sizeof(MMPool);
	mmPool->poolByteSkippedForAlign = 0;
	mmPool->poolByteSpillover = 0;
	mmPool->firstSpillOverAddress = NULL;
	mmPool->currentTempByteDispatched = 0;
	mmPool->currentTempByteSpillover = 0;
	mmPool->maxTotalByteDispatched = 0;

	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] == NULL) {
			mmMaster.mmPool[i] = mmPool;
			return mmPool;
		}
	}

	fprintf(stderr, "MMPoolCreate() : number of pools > maxNumberOfPools!\n");
	exit(1);

}

VAL MMPoolIsActive(const MMPool *mmPool) {

	return ((mmPool->firstSpillOverAddress) != (void*)mmPool);

}
void MMPoolSetInactive(MMPool *mmPool) {

	if (mmPool->firstSpillOverAddress != NULL) {
		fprintf(stderr, "MMPoolSetInactive() : spillover memory not freed yet!\n");
		exit(1);
	}

	mmPool->firstSpillOverAddress = (void*)mmPool;
}

VAL MMPoolCurrentTotalByteAllocated(const MMPool *mmPool) {

	return mmPool->poolSize + mmPool->poolByteSpillover + mmPool->currentTempByteSpillover;

}

VAL MMPoolCurrentTotalByteDispatched(const MMPool *mmPool) {

	return mmPool->poolByteDispatched + mmPool->currentTempByteDispatched;

}

VAL MMPoolMaxTotalByteDispatched(const MMPool *mmPool) {

	VAL currentTotalByteDispatched;

	currentTotalByteDispatched = MMPoolCurrentTotalByteDispatched(mmPool);
	
	if (currentTotalByteDispatched > mmPool->maxTotalByteDispatched) {
		return currentTotalByteDispatched;
	} else {
		return mmPool->maxTotalByteDispatched;
	}

}

MMPool *MMPoolFree(MMPool *mmPool) {

	MMPool *dummyMMPool;
	VAL i;
	void *temp1, *temp2;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolFree(): mmPool = NULL!\n");
		exit(1);
	}
	#endif

	#ifdef RECORD_GRAND_TOTAL
	MMMasterSetMaxTotalByteAllocated();
	MMMasterSetMaxTotalByteDispatched();
	#endif

	dummyMMPool = (MMPool*)malloc(sizeof(MMPool));
	if (dummyMMPool == NULL) {
		fprintf(stderr, "MMPoolFree() : cannot allocate memory!\n");
		exit(1);
	}

	// Fill spillover memory
	temp1 = mmPool->firstSpillOverAddress;
	while (temp1 != NULL) {
		temp2 = *((void**)temp1);
		free(temp1);
		temp1 = temp2;
	}
	mmPool->firstSpillOverAddress = NULL;

	dummyMMPool->poolByteDispatched = mmPool->poolByteDispatched;
	dummyMMPool->poolByteSpillover = mmPool->poolByteSpillover;
	dummyMMPool->currentTempByteDispatched = mmPool->currentTempByteDispatched;
	dummyMMPool->currentTempByteSpillover = mmPool->currentTempByteSpillover;
	dummyMMPool->firstSpillOverAddress = mmPool->firstSpillOverAddress;
	dummyMMPool->maxTotalByteDispatched = mmPool->maxTotalByteDispatched;
	dummyMMPool->poolByteSkippedForAlign = mmPool->poolByteSkippedForAlign;
	dummyMMPool->poolSize = mmPool->poolSize;

	MMPoolSetInactive(dummyMMPool);

	// Update master directory
	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] == mmPool) {
			mmMaster.mmPool[i] = dummyMMPool;
			free(mmPool);
			return dummyMMPool;
		}
	}

	fprintf(stderr, "MMPoolFree() : cannot locate pool in master!\n");
	exit(1);

}

void MMPoolDestory(MMPool *mmPool) {

	VAL i;
	MMPool *temp;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolDestory(): mmPool = NULL!\n");
		exit(1);
	}
	#endif

	if (MMPoolIsActive(mmPool)) {
		temp = MMPoolFree(mmPool);
	} else {
		temp = mmPool;
	}

	// Update master directory
	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] == temp) {
			mmMaster.mmPool[i] = NULL;
			free(temp);
			temp = NULL;
		}
	}

	if (temp != NULL) {
		fprintf(stderr, "MMPoolDestory() : cannot locate pool in master!\n");
		exit(1);
	}

}

void *MMPoolDispatch(MMPool *mmPool, const VAL memSize) {

	void **temp;
	VAL totalPoolMemoryUsed, nextPoolMemoryOffset;
	VAL align, skipForAlign;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolDispatch(): mmPool = NULL!\n");
		exit(1);
	}
	if (memSize == 0) {
		fprintf(stderr, "MMPoolDispatch(): memSize = 0!\n");
		exit(1);
	}
	#endif

	totalPoolMemoryUsed = mmPool->poolByteDispatched - mmPool->poolByteSpillover +
						  mmPool->currentTempByteDispatched - mmPool->currentTempByteSpillover;
	nextPoolMemoryOffset = mmPool->poolByteDispatched - mmPool->poolByteSpillover;

	// Calculate the number of byte to skip in order to align the memory dispatched
	align = 1 << (BITS_IN_WORD - leadingZero(memSize - 1));
	if (align > MAX_ALIGN) {
		align = MAX_ALIGN;
	}
	if (align < MIN_ALIGN) {
		align = MIN_ALIGN;
	}
	skipForAlign = downwardAlign(nextPoolMemoryOffset, align) - nextPoolMemoryOffset;

	if (totalPoolMemoryUsed + memSize + skipForAlign <= mmPool->poolSize) {
		temp = (void**)(((char*)mmPool) + nextPoolMemoryOffset + skipForAlign);
		mmPool->poolByteSkippedForAlign += skipForAlign;
		mmPool->poolByteDispatched += memSize + skipForAlign;
		return temp;
	} else {
		// Spillover
		// Allocate for linked list pointer as well
		temp = (void**)malloc(memSize + sizeof(void*));
		if (temp == NULL) {
			fprintf(stderr, "MMPoolDispatch(): cannot allocate memory!\n");
			exit(1);
		}
		// Add spillover memory to linked list
		*temp = mmPool->firstSpillOverAddress;
		mmPool->firstSpillOverAddress = temp;
		mmPool->poolByteSpillover += memSize;
		mmPool->poolByteDispatched += memSize;
		return (char*)temp + sizeof(void*);
	}
		
}

void MMPoolReturn(MMPool *mmPool, void *address, const VAL memSize) {
	
	// Dummy function

}

void MMPoolPrintReport(MMPool *mmPool, FILE *output) {

	fprintf(output, "Pool Size     : %u\n", mmPool->poolSize);
	fprintf(output, "   Dispatched : %u\n", mmPool->poolByteDispatched);
	fprintf(output, "     - Skipped for alignment : %u\n", mmPool->poolByteSkippedForAlign);
	fprintf(output, "     - Spillover             : %u\n", mmPool->poolByteSpillover);
	fprintf(output, "Maximum amount of memory dispatched including temp memory : %u\n", 
			MMPoolMaxTotalByteDispatched(mmPool));

}

void *MMTempDispatch(MMPool *mmPool, const VAL memSize) {

	void **temp;
	VAL totalPoolMemoryUsed, nextTempMemoryOffset;
	VAL alignedMemSize;
	void **pointerToLastSpilloverAddress;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMTempDispatch(): mmPool = NULL!\n");
		exit(1);
	}
	if (memSize == 0) {
		fprintf(stderr, "MMTempDispatch(): memSize = 0!\n");
		exit(1);
	}
	#endif

	alignedMemSize = downwardAlign(memSize, MAX_ALIGN);

	totalPoolMemoryUsed = mmPool->poolByteDispatched - mmPool->poolByteSpillover +
						  mmPool->currentTempByteDispatched - mmPool->currentTempByteSpillover;
	nextTempMemoryOffset = mmPool->currentTempByteDispatched - mmPool->currentTempByteSpillover;

	if (totalPoolMemoryUsed + alignedMemSize <= mmPool->poolSize) {
		temp = (void**)(((char*)mmPool) + mmPool->poolSize - nextTempMemoryOffset - alignedMemSize);
		mmPool->currentTempByteDispatched += alignedMemSize;
		return temp;
	} else {
		// Spillover
		// Locate the last spillover memory
		pointerToLastSpilloverAddress = &(mmPool->firstSpillOverAddress);
		temp = (void**)(*pointerToLastSpilloverAddress);
		while (temp != NULL) {
			pointerToLastSpilloverAddress = temp;
			temp = (void**)*pointerToLastSpilloverAddress;
		}
		// Allocate for linked list pointer as well
		temp = (void**)malloc(memSize + sizeof(void*));
		if (temp == NULL) {
			fprintf(stderr, "MMTempDispatch(): cannot allocate memory!\n");
			exit(1);
		}
		*pointerToLastSpilloverAddress = temp;
		*temp = NULL;
		mmPool->currentTempByteDispatched += memSize;
		mmPool->currentTempByteSpillover += memSize;
		return (char*)temp + sizeof(void*);
	}
		
}

void MMTempReturn(MMPool *mmPool, void *address, const VAL memSize) {

	void **temp;
	VAL alignedMemSize;
	void **pointerToLastButOneSpillover;
	void *spilloverPointerAddress;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMTempReturn(): mmPool = NULL!\n");
		exit(1);
	}
	if (memSize == 0) {
		fprintf(stderr, "MMTempReturn(): memSize = 0!\n");
		exit(1);
	}
	#endif

	alignedMemSize = downwardAlign(memSize, MAX_ALIGN);

	if (address >= (void*)mmPool && address <= (void*)((char*)mmPool + mmPool->poolSize)) {
		// No need to record the global level max memory dispatched/allocated
		// because memory pool is allocated as a whole and fluctuation across pools should not be counted
		if (mmPool->poolByteDispatched + mmPool->currentTempByteDispatched > mmPool->maxTotalByteDispatched) {
			mmPool->maxTotalByteDispatched = mmPool->poolByteDispatched + mmPool->currentTempByteDispatched;
		}
		mmPool->currentTempByteDispatched -= alignedMemSize;
	} else {
		#ifdef RECORD_GRAND_TOTAL
		MMMasterSetMaxTotalByteAllocated();
		MMMasterSetMaxTotalByteDispatched();
		#endif
		// Spillover
		spilloverPointerAddress = (void*)((char*)address - sizeof(void*));
		// Locate the last spillover memory
		pointerToLastButOneSpillover = &(mmPool->firstSpillOverAddress);
		temp = (void**)(*pointerToLastButOneSpillover);
		while (*temp != NULL) {
			pointerToLastButOneSpillover = temp;
			temp = (void**)*pointerToLastButOneSpillover;
		}
		if (*pointerToLastButOneSpillover != spilloverPointerAddress) {
			fprintf(stderr, "MMTempReturn(): address != lastSpilloverAddress!\n");
			exit(1);
		}
		free(spilloverPointerAddress);
		*pointerToLastButOneSpillover = NULL;

		if (mmPool->poolByteDispatched + mmPool->currentTempByteDispatched > mmPool->maxTotalByteDispatched) {
			mmPool->maxTotalByteDispatched = mmPool->poolByteDispatched + mmPool->currentTempByteDispatched;
		}
		mmPool->currentTempByteDispatched -= memSize;
		mmPool->currentTempByteSpillover -= memSize;
	}

}

void MMTempPrintReport(MMPool *mmPool, FILE *output) {

	MMPoolPrintReport(mmPool, output);

}

MMBulk *MMBulkCreate(MMPool *mmPool, const VAL itemSize, const VAL itemPerAllocationInPowerOf2, 
					 const VAL boundaryCushionSize, const VAL directorySize) {

	VAL i;
	MMBulk *mmBulk;

	#ifdef DEBUG
	if (itemSize == 0) {
		fprintf(stderr, "MMBulkCreate() : itemSize = 0!\n");
		exit(1);
	}
	if (itemPerAllocationInPowerOf2 >= BITS_IN_WORD) {
		fprintf(stderr, "MMBulkCreate() : itemPerAllocationInPowerOf2 >= BITS_IN_WORD!\n");
		exit(1);
	}
	#endif

	if (mmPool == NULL) {
		mmBulk = (MMBulk*)MMUnitAllocate(sizeof(MMBulk));
	} else {
		mmBulk = (MMBulk*)MMPoolDispatch(mmPool, sizeof(MMBulk));
	}

	mmBulk->itemSize = itemSize;
	mmBulk->itemPerAllocationInPowerOf2 = itemPerAllocationInPowerOf2;
	mmBulk->boundaryCushionSize = boundaryCushionSize;
	mmBulk->indexMask = truncateLeft(ALL_ONE_MASK,  BITS_IN_WORD - itemPerAllocationInPowerOf2);
	mmBulk->currentDirectoryEntry = 0;
	mmBulk->nextUnusedItem = 0;
	mmBulk->directorySize = directorySize;

	if (mmPool == NULL) {
		mmBulk->directory = (unsigned char**)MMUnitAllocate(sizeof(unsigned char*) * directorySize);
	} else {
		mmBulk->directory = (unsigned char**)MMPoolDispatch(mmPool, sizeof(unsigned char*) * directorySize);
	}

	//Allocate memory for the first directory entry
	mmBulk->directory[0] = (unsigned char*)malloc(boundaryCushionSize * 2 + (itemSize << itemPerAllocationInPowerOf2));
	if (mmBulk->directory[0] == NULL) {
		fprintf(stderr, "MMBulkCreate() : cannot allocate memory!\n");
		exit(1);
	}

	//Advance the address by boundaryCushionSize
	mmBulk->directory[0] += boundaryCushionSize;

	for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] == NULL) {
			mmMaster.mmBulk[i] = mmBulk;
			return mmBulk;
		}
	}

	fprintf(stderr, "MMBulkCreate() : number of bulks > maxNumberOfBulk!\n");
	exit(1);

}

VAL MMBulkIsActive(const MMBulk *mmBulk) {

	return (mmBulk->directory != (void*)mmBulk);

}

void MMBulkSetInactive(MMBulk *mmBulk) {

	if (mmBulk->directory != NULL) {
	}
	mmBulk->directory = (unsigned char**)mmBulk;

}

VAL MMBulkByteAllocated(const MMBulk *mmBulk) {

	return (mmBulk->currentDirectoryEntry + 1) *
			(mmBulk->boundaryCushionSize * 2 + (mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2));

}

VAL MMBulkByteDispatched(const MMBulk *mmBulk) {

	return (mmBulk->currentDirectoryEntry) *
			(mmBulk->boundaryCushionSize * 2 + (mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2)) +
			mmBulk->boundaryCushionSize * 2 +
			mmBulk->itemSize * mmBulk->nextUnusedItem;

}

VAL MMBulkUnitDispatched(const MMBulk *mmBulk) {

	return mmBulk->currentDirectoryEntry * (1 << mmBulk->itemPerAllocationInPowerOf2) + mmBulk->nextUnusedItem;

}

void MMBulkFree(MMBulk *mmBulk) {

	VAL i;

	#ifdef RECORD_GRAND_TOTAL
	MMMasterSetMaxTotalByteAllocated();
	MMMasterSetMaxTotalByteDispatched();
	#endif

	for (i=0; i<=mmBulk->currentDirectoryEntry; i++) {
		free(mmBulk->directory[i] - mmBulk->boundaryCushionSize);
	}

	if (MMBulkFindPoolUsed(mmBulk) == NULL) {
        MMUnitFree(mmBulk->directory, sizeof(unsigned char*) * mmBulk->directorySize);
	}

	mmBulk->directory = NULL;

	MMBulkSetInactive(mmBulk);

}

void MMBulkDestory(MMBulk *mmBulk) {

	VAL i;
	MMBulk *temp;

	#ifdef DEBUG
	if (mmBulk == NULL) {
		fprintf(stderr, "MMBulkDestory(): mmBulk = NULL!\n");
		exit(1);
	}
	#endif

	if (MMBulkIsActive(mmBulk)) {
		MMBulkFree(mmBulk);
	}

	temp = mmBulk;

	// Update master directory
	for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] == temp) {
			mmMaster.mmBulk[i] = NULL;
			if (MMBulkFindPoolUsed(temp) == NULL) {
				MMUnitFree(temp, sizeof(MMBulk));
			}
			temp = NULL;
		}
	}

	if (temp != NULL) {
		fprintf(stderr, "MMBulkDestory() : cannot locate bulk in master!\n");
		exit(1);
	}

}
VAL MMBulkDispatch(MMBulk *mmBulk) {

	if (mmBulk->nextUnusedItem >> mmBulk->itemPerAllocationInPowerOf2) {
		mmBulk->currentDirectoryEntry++;
		if (mmBulk->currentDirectoryEntry >= mmBulk->directorySize) {
			fprintf(stderr, "MMBulkDispatch() : memory directory size overflow!\n");
			exit(1);
		}
		//Allocate memory for the next directory entry
		mmBulk->directory[mmBulk->currentDirectoryEntry] = (unsigned char*)malloc(mmBulk->boundaryCushionSize * 2 + (mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2));
		if (mmBulk->directory[mmBulk->currentDirectoryEntry] == NULL) {
			fprintf(stderr, "MMBulkDispatch() : cannot allocate memory!\n");
			exit(1);
		}
		//Advance the address by boundaryCushionSize
		mmBulk->directory[mmBulk->currentDirectoryEntry] += mmBulk->boundaryCushionSize;
		mmBulk->nextUnusedItem = 0;
	}
	return ((mmBulk->currentDirectoryEntry << mmBulk->itemPerAllocationInPowerOf2) | mmBulk->nextUnusedItem++);

}

void *MMBulkAddress(const MMBulk *mmBulk, const VAL index) {

	#ifdef DEBUG
	if (index >= (((mmBulk->currentDirectoryEntry+1) << mmBulk->itemPerAllocationInPowerOf2) | mmBulk->nextUnusedItem)) {
		fprintf(stderr, "MMBulkAddress() : index out of range!\n");
		exit(1);
	}
	#endif

	return &(mmBulk->directory[index >> mmBulk->itemPerAllocationInPowerOf2][(index & mmBulk->indexMask) * mmBulk->itemSize]);
}

MMPool *MMBulkFindPoolUsed(const MMBulk *mmBulk) {

	VAL i;
	void *temp;

	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] != NULL) {
			if ((void*)mmBulk >= (void*)mmMaster.mmPool[i] &&
				(void*)mmBulk <= (void*)((char*)mmMaster.mmPool[i] + mmMaster.mmPool[i]->poolSize)) {
				return mmMaster.mmPool[i];
			}
			temp = mmMaster.mmPool[i]->firstSpillOverAddress;
			while (temp != NULL) {
				if ((void*)((char*)temp + sizeof(void*)) == (void*)mmBulk) {
					return mmMaster.mmPool[i];
				}
				temp = *((void**)temp);
			}
		}
	}

	return NULL;

}

void MMBulkPrintReport(MMBulk *mmBulk, FILE *output){

	fprintf(output, "Memory allocated  : %u\n", MMBulkByteAllocated(mmBulk));
	fprintf(output, "Memory dispatched : %u\n", MMBulkByteDispatched(mmBulk));

}

void MMBulkSave(MMBulk *mmBulk, FILE *output) {

	VAL i;

	fwrite(&mmBulk->itemSize, sizeof(VAL), 1, output);
	fwrite(&mmBulk->itemPerAllocationInPowerOf2, sizeof(VAL), 1, output);
	fwrite(&mmBulk->boundaryCushionSize, sizeof(VAL), 1, output);
	fwrite(&mmBulk->currentDirectoryEntry, sizeof(VAL), 1, output);
	fwrite(&mmBulk->nextUnusedItem, sizeof(VAL), 1, output);
	fwrite(&mmBulk->directorySize, sizeof(VAL), 1, output);

	for (i=0; i<mmBulk->currentDirectoryEntry; i++) {
		fwrite(mmBulk->directory[i], mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2, 1, output);
	}

	if (mmBulk->nextUnusedItem > 0) {
		fwrite(mmBulk->directory[i], mmBulk->itemSize * mmBulk->nextUnusedItem, 1, output);
	}

}

MMBulk *MMBulkLoad(MMPool *mmPool, FILE *input) {

	VAL i;
	MMBulk *mmBulk;

	mmBulk = (MMBulk*)MMPoolDispatch(mmPool, sizeof(MMBulk));

	fread(&mmBulk->itemSize, sizeof(VAL), 1, input);
	fread(&mmBulk->itemPerAllocationInPowerOf2, sizeof(VAL), 1, input);
	fread(&mmBulk->boundaryCushionSize, sizeof(VAL), 1, input);
	fread(&mmBulk->currentDirectoryEntry, sizeof(VAL), 1, input);
	fread(&mmBulk->nextUnusedItem, sizeof(VAL), 1, input);
	fread(&mmBulk->directorySize, sizeof(VAL), 1, input);

	mmBulk->indexMask = truncateLeft(ALL_ONE_MASK,  BITS_IN_WORD - mmBulk->itemPerAllocationInPowerOf2);

	mmBulk->directory = (unsigned char**)MMPoolDispatch(mmPool, sizeof(CHAR*) * mmBulk->directorySize);

	for (i=0; i<mmBulk->currentDirectoryEntry; i++) {
		mmBulk->directory[i] = (unsigned char*)malloc(mmBulk->boundaryCushionSize * 2 + 
									(mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2));
		if (mmBulk->directory[i] == NULL) {
			fprintf(stderr, "MMBulkLoad() : cannot allocate memory!\n");
			exit(1);
		}

		//Advance the address by boundaryCushionSize
		mmBulk->directory[i] += mmBulk->boundaryCushionSize;
		
		fread(mmBulk->directory[i], mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2, 1, input);
	}

	mmBulk->directory[i] = (unsigned char*)malloc(mmBulk->boundaryCushionSize * 2 + 
								(mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2));
	if (mmBulk->directory[i] == NULL) {
		fprintf(stderr, "MMBulkLoad() : cannot allocate memory!\n");
		exit(1);
	}

	//Advance the address by boundaryCushionSize
	mmBulk->directory[i] += mmBulk->boundaryCushionSize;

	if (mmBulk->nextUnusedItem > 0) {
		fread(mmBulk->directory[i], mmBulk->itemSize * mmBulk->nextUnusedItem, 1, input);
	}


	for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] == NULL) {
			mmMaster.mmBulk[i] = mmBulk;
			return mmBulk;
		}
	}

	fprintf(stderr, "MMBulkLoad() : number of bulks > maxNumberOfBulk!\n");
	exit(1);

}
