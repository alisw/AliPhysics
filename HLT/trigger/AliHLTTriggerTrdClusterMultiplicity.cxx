// $Id$
//-*- Mode: C++ -*-
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Svein Lindal <svein.lindal@fys.uio.no>                 *
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

/// @file   AliHLTTriggerTrdClusterMultiplicity.cxx
/// @author Svein Lindal
/// @date   2009-08-19
/// @brief  HLT Minimum Ionizing Particle (MIP) trigger for PHOS

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTTriggerTrdClusterMultiplicity.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDUtils.h"
#include "AliHLTTRDCluster.h"
#include "AliTRDgeometry.h"
#include "TObjString.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerTrdClusterMultiplicity)

AliHLTTriggerTrdClusterMultiplicity::AliHLTTriggerTrdClusterMultiplicity() 
  : AliHLTTrigger()
  , fClusterMult(10)
  , fClusterArray(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlts
}

AliHLTTriggerTrdClusterMultiplicity::~AliHLTTriggerTrdClusterMultiplicity()
{
  // see header file for class documentation
}

const char* AliHLTTriggerTrdClusterMultiplicity::GetTriggerName() const
{
  // see header file for class documentation
  return "TrdClusterMultiplicityTrigger";
}

AliHLTComponent* AliHLTTriggerTrdClusterMultiplicity::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerTrdClusterMultiplicity;
}

int AliHLTTriggerTrdClusterMultiplicity::DoTrigger()
{
  // see header file for class documentation

  TString description;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(AliHLTTRDDefinitions::fgkClusterDataType);
       pBlock!=NULL; pBlock=GetNextInputBlock()) {

    AliHLTTRDUtils::ReadClusters(fClusterArray, pBlock->fPtr, pBlock->fSize);

    AliTRDcluster* inCluster;
    Int_t stack[5] = {0,0,0,0,0};
    for(Int_t count=0; count<fClusterArray->GetEntriesFast(); count++)
      {
	inCluster=(AliTRDcluster*)fClusterArray->UncheckedAt(count);
	stack[AliTRDgeometry::GetStack(inCluster->GetDetector())]++;
      }

    Bool_t trigger=kFALSE;
    for(Int_t count=0; count<5; count++)
      {
	if(stack[count]>fClusterMult)trigger=kTRUE;
      }

    if(trigger){
      description.Form("Event contains at least one stack with at least %i TRD clusters", fClusterMult);
      SetDescription(description.Data());
	
      // Enable the detectors for readout.
      GetReadoutList().Enable( AliHLTReadoutList::kTRD );
	
      // Add the available HLT information for readout too.
      GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "TRD");
	
      //Set trigger decision
      TriggerEvent(true);

    }else{
      description.Form("No stack in this event contains enough clusters");
      SetDescription(description.Data());
      TriggerEvent(false);
    }
    
    HLTDebug(description.Data());
    fClusterArray->Delete();

    return 0;
  }
  return 0;
}

int AliHLTTriggerTrdClusterMultiplicity::DoInit(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;

  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }

  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  }

  fClusterArray = new TClonesArray("AliTRDcluster");
  return 0;
}

int AliHLTTriggerTrdClusterMultiplicity::DoDeinit()
{
  // see header file for class documentation

  fClusterArray->Delete();
  delete fClusterArray;

  return 0;
}

int AliHLTTriggerTrdClusterMultiplicity::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
 
  return 0;
}

int AliHLTTriggerTrdClusterMultiplicity::Configure(const char* arguments){
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;

      if (argument.CompareTo("-MultiplicityThresh")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Setting multiplicity threshold to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fClusterMult=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 

      else {
	HLTError("unknown argument: %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}
