#ifndef MDC_H
#define MDC_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// c interface to AliMDC                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifdef __cplusplus
extern "C" {
#endif

void* alimdcCreate(int compress, int filterMode, 
		   const char* localRunDB, int rdbmsRunDB,
		   const char* alienHostRunDB, const char* alienDirRunDB,
		   double maxSizeTagDB, const char* fileNameTagDB);
int   alimdcOpen(void* alimdc, int mode, const char* fileName);
int   alimdcProcessEvent(void* alimdc, void* event, int isIovecArray);
int   alimdcClose(void* alimdc);
void  alimdcDelete(void* alimdc);

#ifdef __cplusplus
}
#endif

#endif
