// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk                                        *
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

/// @file   AliHLTTriggerITSMultiplicity.cxx
/// @author Gaute Ovrebekk
/// @date   2009-10-22
/// @brief  HLT trigger component for cluster multiplicity
///         in ITS.


// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerITSMultiplicity.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "TObjString.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTITSSpacePointData.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerITSMultiplicity)

AliHLTTriggerITSMultiplicity::AliHLTTriggerITSMultiplicity()
  : AliHLTTrigger()
  , fnClusters(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

const char* AliHLTTriggerITSMultiplicity::fgkOCDBEntry="HLT/ConfigHLT/ITSMultiplicityTrigger";

AliHLTTriggerITSMultiplicity::~AliHLTTriggerITSMultiplicity()
{
  // see header file for class documentation
}

const char* AliHLTTriggerITSMultiplicity::GetTriggerName() const
{
  // see header file for class documentation
  return "ITSMultiplicityTrigger";
}

AliHLTComponent* AliHLTTriggerITSMultiplicity::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerITSMultiplicity;
}

int AliHLTTriggerITSMultiplicity::DoTrigger()
{
  // see header file for class documentation
  int iResult=0;
  int numberOfClusters=-1;
  TString description;
 
  const AliHLTComponentBlockData* iter = NULL;

  if(!IsDataEvent()) return 0;

  for(iter = GetFirstInputBlock(kAliHLTDataTypeClusters); iter != NULL; iter = GetNextInputBlock()){

      if(iter->fDataType!=(kAliHLTAnyDataType|kAliHLTDataOriginITSSPD) &&
         iter->fDataType!=(kAliHLTAnyDataType|kAliHLTDataOriginITSSDD) &&
         iter->fDataType!=(kAliHLTAnyDataType|kAliHLTDataOriginITSSSD))
         continue;
      
      const AliHLTITSClusterData* clusterData = (const AliHLTITSClusterData*) iter->fPtr;
      Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;
      numberOfClusters += nSpacepoint;

  }

  if (iResult>=0 && numberOfClusters>=0) {
    
    if (numberOfClusters>=fnClusters) {
      description.Form("Event contains %d cluster(s)", numberOfClusters);
      SetDescription(description.Data());
      // Enable the central detectors for readout.
      GetReadoutList().Enable(
			      AliHLTReadoutList::kITSSPD |
			      AliHLTReadoutList::kITSSDD |
			      AliHLTReadoutList::kITSSSD |
			      AliHLTReadoutList::kTPC |
			      AliHLTReadoutList::kTRD |
			      AliHLTReadoutList::kTOF |
			      AliHLTReadoutList::kHMPID |
			      AliHLTReadoutList::kPHOS
			      );
      // Add the available HLT information for readout too.
      GetTriggerDomain().Add("CLUSTERS", "ITS ");
      TriggerEvent(true);
      return 0;
    }
    description.Form("No clusters matching the tresholds found (min clusters %d  ", fnClusters);
  } else {
    description.Form("No input blocks found");
  }
  SetDescription(description.Data());
  TriggerEvent(false);
  return iResult;
}

int AliHLTTriggerITSMultiplicity::DoInit(int argc, const char** argv)
{
  // see header file for class documentation

  // first configure the default
  int iResult=0;
  if (iResult>=0) iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  return iResult;
}

int AliHLTTriggerITSMultiplicity::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTriggerITSMultiplicity::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) {
    entry=fgkOCDBEntry;
  }

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerITSMultiplicity::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -nclusters
  if (argument.CompareTo("-nclusters")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fnClusters=argument.Atoi();
    return 2;
  }    
  
  // unknown argument
  return -EINVAL;
}
