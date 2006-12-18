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
		   double maxSizeTagDB, const char* fileNameTagDB);
int   alimdcOpen(void* alimdc, int mode, const char* fileName);
int   alimdcProcessEvent(void* alimdc, void* event, int isIovecArray);
int   alimdcGetTotalFileSize(void* alimdc);
int   alimdcClose(void* alimdc);
void  alimdcDelete(void* alimdc);
void  alimdcEnableDebug();

#ifdef __cplusplus
}
#endif

#endif
