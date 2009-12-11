/* MemManagerDummy.h		Memory Manager Dummy Header File

   Copyright 2004, Wong Chi Kwong, all rights reserved.

   This header file translates interfaces to Memory Manager to ordinary memory
   management functions provided by C (malloc and free).

   This is to be used when standard C memory management functions are prefered.

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.

*/

#ifndef __MEM_MANAGER_DUMMY_H__
#define __MEM_MANAGER_DUMMY_H__

#define MMUnitAllocate(memSize)					MMMalloc(memSize)
#define MMUnitFree(address, memSize)			MMFree(address)
#define MMTempDispatch(mmPool, memSize)			MMMalloc(memSize)
#define MMTempReturn(mmPool, address, memSize)	MMFree(address)
#define MMPoolDispatch(mmPool, memSize)			MMMalloc(memSize)
#define MMPoolReturn(mmPool, address, memSize)	MMFree(address)

#define MMPool	void

#endif