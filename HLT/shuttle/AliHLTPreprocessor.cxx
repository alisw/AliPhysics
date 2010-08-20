// $Id: AliHLTPreprocessor.cxx 23039 2007-12-13 20:53:02Z richterm $

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


// @file   AliHLTPreprocessor.cxx
// @author Matthias Richter
// @brief  Container for HLT module preprocessors, acts to the outside as
//         HLT preprocessor used by the Offline Shuttle 
// 

#include "AliHLTPreprocessor.h"
#include "AliHLTModulePreprocessor.h"
#include "AliHLTSystem.h"
#include "AliHLTModuleAgent.h"
#include "TSystem.h"

ClassImp(AliHLTPreprocessor)

AliHLTPreprocessor::AliHLTPreprocessor(AliShuttleInterface* shuttle) 
  :
  AliPreprocessor(fgkHLTPreproc, shuttle),
  fProcessors(),
  fActiveDetectors(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  // run types according to 
  // http://alice-ecs.web.cern.ch/alice-ecs/runtypes_3.16.html

  // PHOS (retrieve Huffman tables)
  AddRunType("STANDALONE");

  // TPC (retrieve Huffman tables and temperature data)
  AddRunType("PHYSICS");
  AddRunType("COSMIC");
  AddRunType("LASER");
  AddRunType("PEDESTAL");
  AddRunType("PULSER");

  // TRD
  AddRunType("PEDESTAL");
  AddRunType("STANDALONE");
 
  fProcessors.SetOwner();

}

const char* AliHLTPreprocessor::fgkHLTPreproc = "HLT";

/** HLT default component libraries */
const char* AliHLTPreprocessor::fgkHLTDefaultShuttleLibs[]= {
  "libAliHLTUtil.so", 
  "libAliHLTRCU.so", 
  "libAliHLTTPC.so", 
  "libAliHLTComp.so", 
  "libAliHLTSample.so",
  //"libAliHLTPHOS.so",
  //"libAliHLTMUON.so",
  "libAliHLTTRD.so",
  "libAliHLTGlobal.so",
  "libAliHLTTrigger.so",
  NULL
};

AliHLTPreprocessor::~AliHLTPreprocessor()
{
  // see header file for function documentation
}

void AliHLTPreprocessor::Initialize(Int_t run, UInt_t startTime, 
			UInt_t endTime) 
{
  // see header file for function documentation
  fRun = run;
  fStartTime = startTime;
  fEndTime = endTime;

  // TODO: read a configuration object from OCDB
  // configure
  // - component libraries

  // retrieve list of active detectors from previous run.
  fActiveDetectors = atoi(AliPreprocessor::GetRunParameter("detectorMask"));

  //   TString msg("Preprocessor for HLT initialized for run: ");
  //   msg += run;
  //   Log(msg.Data());
  
  // load component libraries
  TString libs;
  const char** deflib=fgkHLTDefaultShuttleLibs;
  while (*deflib) {
    if (gSystem->Load(*deflib)==0) {
      Log(Form("HLT component library %s loaded", *deflib));
    }
    
    deflib++;
  }

  for (AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent();
       pAgent!=NULL;
       pAgent=AliHLTModuleAgent::GetNextAgent()) {
    AliHLTModulePreprocessor* pProc=pAgent->GetPreprocessor();
    if (pProc) 
      {

	// test if pProc is necessary, if not, take next one
	if((pProc->GetModuleNumber() & fActiveDetectors) == 0)
	  {
	    TString msg;
	    msg.Form("%s not needed", pProc->GetModuleID());
	    Log(msg.Data());
	    continue;
	  }
	
	pProc->SetShuttleInterface(this);
	pProc->Initialize(run, startTime, endTime);
	fProcessors.Add(pProc);
	TString msg;
	msg.Form("added preprocessor %p with ID %s for module %p", pProc, pProc->GetModuleID(), pAgent);
	Log(msg.Data());
      }
  }
}

UInt_t AliHLTPreprocessor::Process(TMap* dcsAliasMap)
{
  // see header file for function documentation
  UInt_t retVal = 0;

  if (!GetHLTStatus()) {
    return 0;
  }

  bool bAllFailed=true;
  TObjLink *lnk = NULL;
  lnk=fProcessors.FirstLink();
  while (lnk) {
    AliHLTModulePreprocessor* pProc=dynamic_cast<AliHLTModulePreprocessor*>(lnk->GetObject());
    if (pProc) {
      UInt_t result=pProc->Process(dcsAliasMap);
      if (result) {
	TString msg;
	msg.Form("preprocessor for module %s failed with error code %d", pProc->GetName(), result);
	Log(msg.Data());
      } else {
	bAllFailed=false;
      }
    }
    lnk = lnk->Next();
  }

  if (bAllFailed) return 1;
  return retVal;
}


Bool_t AliHLTPreprocessor::ProcessDCS()
{
  // see header file for function documentation
  return kFALSE;
}
