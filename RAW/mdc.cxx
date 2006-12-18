/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// c interface to AliMDC                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "mdc.h"
#include "AliMDC.h"
#include "AliLog.h"

void* alimdcCreate(int compress, int filterMode, 
		   double maxSizeTagDB, const char* fileNameTagDB)
{
// create an AliMDC object

  return new AliMDC(compress, kFALSE, AliMDC::EFilterMode(filterMode), 
		    maxSizeTagDB, fileNameTagDB);

}

int alimdcOpen(void* alimdc, int mode, const char* fileName)
{
// open a new raw DB

  return ((AliMDC*)alimdc)->Open(AliMDC::EWriteMode(mode), fileName);
}

int alimdcProcessEvent(void* alimdc, void* event, int isIovecArray)
{
// process one event

  return ((AliMDC*)alimdc)->ProcessEvent(event, isIovecArray);
}

int alimdcGetTotalFileSize(void* alimdc)
{
// return the total current file size

  return ((AliMDC*)alimdc)->GetTotalSize();
}

int alimdcClose(void* alimdc)
{
// close the raw DB

  return ((AliMDC*)alimdc)->Close();
}

void  alimdcDelete(void* alimdc)
{
// delete the AliMDC object

  delete (AliMDC*)alimdc;
}

void  alimdcEnableDebug()
{
// enable debug and log messages

  AliLog::EnableDebug(kTRUE);
  AliLog::SetGlobalLogLevel(AliLog::kMaxType);
  AliLog::SetGlobalDebugLevel(AliLog::kMaxType);
  AliLog::SetPrintRepetitions(kFALSE);
  AliLog::SetHandleRootMessages(kTRUE);
}
