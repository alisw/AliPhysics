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

/** @file   AliHLTEsdManager.cxx
    @author Matthias Richter
    @date   
    @brief  Manager for merging and writing of HLT ESDs
*/

#include "AliHLTEsdManager.h"
#include "TSystem.h"
#include "TClass.h"
#include "TROOT.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTEsdManager)

AliHLTEsdManager::AliHLTEsdManager()
  :
  AliHLTLogging()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTEsdManager::~AliHLTEsdManager()
{
  // see header file for class documentation
}

const char* AliHLTEsdManager::fgkImplName="AliHLTEsdManagerImplementation";
const char* AliHLTEsdManager::fgkImplLibrary="libHLTrec.so";


AliHLTEsdManager* AliHLTEsdManager::New()
{
  // see header file for class documentation
  int iLibResult=0;
  AliHLTEsdManager* instance=NULL;
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
	instance=reinterpret_cast<AliHLTEsdManager*>(p);
	if (!instance) {
	  log.Logging(kHLTLogError, "AliHLTEsdManager::New", "ESD handling", "type cast to AliHLTEsdManager instance failed");
	}
      } else {
	log.Logging(kHLTLogError, "AliHLTEsdManager::New", "ESD handling", "can not create AliHLTEsdManager instance from class descriptor");
      }
    } else {
      log.Logging(kHLTLogError, "AliHLTEsdManager::New", "ESD handling", "can not find AliHLTEsdManager class descriptor");
    }
  } else {
    log.Logging(kHLTLogError, "AliHLTEsdManager::New", "ESD handling", "can not load libHLTrec library");
  }
  return instance;
}

void AliHLTEsdManager::Delete(AliHLTEsdManager* instance)
{
  // see header file for class documentation
  if (!instance) return;

  // check if the library is still there in order to have the
  // destructor available
  TClass* pCl=TClass::GetClass(fgkImplName);
  if (!pCl) {
    AliHLTLogging log;
    log.Logging(kHLTLogError, "AliHLTEsdManager::Delete", "ESD handling", "potential memory leak: libHLTrec library not available, skipping destruction %p", instance);    
    return;
  }

  delete instance;
}
