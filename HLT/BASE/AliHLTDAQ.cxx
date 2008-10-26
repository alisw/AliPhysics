// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTDAQ.cxx
    @author Matthias Richter
    @date   24.10.2008
    @brief  Virtual Interface to the AliDAQ class.
*/

#include "AliHLTDAQ.h"
#include "AliHLTLogging.h"
#include "TClass.h"
#include "TSystem.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDAQ)

AliHLTDAQ::AliHLTDAQ()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTDAQ* AliHLTDAQ::fgpInstance=NULL;
const char* AliHLTDAQ::fgkImplName="AliHLTDAQInterfaceImplementation";
const char* AliHLTDAQ::fgkImplLibrary="libHLTrec.so";

AliHLTDAQ::~AliHLTDAQ()
{
  // see header file for class documentation
}

Int_t       AliHLTDAQ::DetectorID(const char *detectorName)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorID(detectorName);
  return -1;
}

const char *AliHLTDAQ::DetectorName(Int_t detectorID)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorName(detectorID);
  return NULL;
}

Int_t       AliHLTDAQ::DdlIDOffset(const char *detectorName)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlIDOffset(detectorName);
  return -1;
}

Int_t       AliHLTDAQ::DdlIDOffset(Int_t detectorID)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlIDOffset(detectorID);
  return -1;
}

const char *AliHLTDAQ::DetectorNameFromDdlID(Int_t ddlID, Int_t &ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorNameFromDdlID(ddlID, ddlIndex);
  return NULL;
}

Int_t       AliHLTDAQ::DetectorIDFromDdlID(Int_t ddlID, Int_t &ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDetectorIDFromDdlID(ddlID, ddlIndex);
  return -1;
}

Int_t       AliHLTDAQ::DdlID(const char *detectorName, Int_t ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlID(detectorName, ddlIndex);
  return -1;
}

Int_t       AliHLTDAQ::DdlID(Int_t detectorID, Int_t ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlID(detectorID, ddlIndex);
  return -1;
}

const char *AliHLTDAQ::DdlFileName(const char *detectorName, Int_t ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlFileName(detectorName, ddlIndex);
  return NULL;
}

const char *AliHLTDAQ::DdlFileName(Int_t detectorID, Int_t ddlIndex)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtDdlFileName(detectorID, ddlIndex);
  return NULL;
}

Int_t       AliHLTDAQ::NumberOfDdls(const char *detectorName)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtNumberOfDdls(detectorName);
  return -1;
}

Int_t       AliHLTDAQ::NumberOfDdls(Int_t detectorID)
{
  // see header file for class documentation
  if (!fgpInstance) GetInstance();
  if (fgpInstance) return fgpInstance->VirtNumberOfDdls(detectorID);
  return -1;
}

AliHLTDAQ* AliHLTDAQ::GetInstance()
{
  // see header file for class documentation
  int iLibResult=0;
  if (!fgpInstance) {
    AliHLTLogging log;
    TClass* pCl=NULL;
    ROOT::NewFunc_t pNewFunc=NULL;
    do {
      pCl=TClass::GetClass(fgkImplName);
    } while (!pCl && (iLibResult=gSystem->Load(fgkImplLibrary))==0);
    if (iLibResult>=0) {
      if (pCl && (pNewFunc=pCl->GetNew())!=NULL) {
	void* p=(*pNewFunc)(NULL);
	if (p) {
	  fgpInstance=reinterpret_cast<AliHLTDAQ*>(p);
	  if (!fgpInstance) {
	    log.Logging(kHLTLogError, "AliHLTDAQ::Instance", "HLT Analysis", "type cast to AliHLTDAQ instance failed");
	  }
	} else {
	  log.Logging(kHLTLogError, "AliHLTDAQ::Instance", "HLT Analysis", "can not create AliHLTDAQ instance from class descriptor");
	}
      } else {
	log.Logging(kHLTLogError, "AliHLTDAQ::Instance", "HLT Analysis", "can not find AliHLTDAQ class descriptor");
      }
    } else {
      log.Logging(kHLTLogError, "AliHLTDAQ::Instance", "HLT Analysis", "can not load libHLTrec library");
    }
  }
  return fgpInstance;
}
