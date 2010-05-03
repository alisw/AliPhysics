// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Svein Lindal <svein.lindal@gmail.com>                 *
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

/// @file   AliHLTTriggerCaloClusterEnergy.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @date   2009-08-17
/// @brief  BASE class for energy threshold trigger for Calorimeters
///      

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTTriggerCaloClusterEnergy.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTCaloClusterReader.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "TRefArray.h"
#include "TString.h"
#include "TMap.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerCaloClusterEnergy)

AliHLTTriggerCaloClusterEnergy::AliHLTTriggerCaloClusterEnergy(TString detector) : 
  AliHLTTrigger(),
  fEThreshold(0.0),
  fClustersRefs(NULL),
  fDetector(detector),
  fClusterReader(NULL),
  fgkOCDBEntry(""), 
  fgkInputDataType()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlts

  fClusterReader = new AliHLTCaloClusterReader();
  fClustersRefs = new TRefArray();

}

//const char* AliHLTTriggerCaloClusterEnergy::fgkOCDBEntry="HLT/ConfigHLT/CaloClusterEnergyTrigger";

AliHLTTriggerCaloClusterEnergy::~AliHLTTriggerCaloClusterEnergy() {
  // see header file for class documentation
  if (fClusterReader)
    delete fClusterReader;
  fClusterReader = NULL;

  if(fClustersRefs)
    delete fClustersRefs;
  fClustersRefs = NULL;
}

Int_t AliHLTTriggerCaloClusterEnergy::DoTrigger() {
  // see header file for class documentation
  
  Int_t iResult = 0;


  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  //Try the caloclusterstruct input

  
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(fgkInputDataType); pBlock!=NULL; pBlock=GetNextInputBlock()) {
    AliHLTCaloClusterHeaderStruct *caloClusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
    fClusterReader->SetMemory(caloClusterHeader);
    
    AliHLTCaloClusterDataStruct * caloClusterStruct;
    while( (caloClusterStruct = fClusterReader->NextCluster()) != 0) {
      if (TriggerOnCluster(caloClusterStruct)) {
	return iResult;
      }
    }
  }

  //Try the ESD input
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  
  if (esd != NULL) {
    esd->GetStdContent();

    Int_t ncc = GetClustersFromEsd(esd, fClustersRefs); 
    
    for (Int_t i = 0; i < ncc ; i++) {
      
      AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(fClustersRefs->At(i));
      if(TriggerOnCluster(cluster)) {
	return iResult;
      }
    }
  }

  // If we got to this point then we did not find any clusters with E > fEThreshold
  // generate negative trigger decision
  TString description;
  description.Form("No %s clusters containing energy > %.02f GeV found.", fDetector.Data(), fEThreshold);
  SetDescription(description.Data());
  TriggerEvent(false);
  return iResult;

}


template <class T>
Bool_t AliHLTTriggerCaloClusterEnergy::TriggerOnCluster(T* cluster) {
  
  if (cluster->E() > fEThreshold) {

    //We have a cluster satisfying trigger criteria
    TString description;
    description.Form("Event contains at least one %s cluster with energy > %.02f GeV.", fDetector.Data(), fEThreshold);
    SetDescription(description.Data());
    
    // Enable the detectors for readout.
    GetReadoutList().Enable( AliHLTReadoutList::kPHOS );
    
    // Add the available HLT information for readout too.
    GetTriggerDomain().Add(kAliHLTAnyDataTypeID, fDetector.Data());
    
    //Set trigger decision
    TriggerEvent(kTRUE);
    
    return kTRUE;
  } 


  return kFALSE;

}




int AliHLTTriggerCaloClusterEnergy::DoInit(int argc, const char** argv) {
  // see header file for class documentation

  // first configure the default
  int iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0) {
    iResult=ConfigureFromArgumentString(argc, argv);
    HLTImportant("Trigger threshold set from argument string:  %.02f GeV:", fEThreshold ); 
  } else if ( iResult >=0 ) {
    HLTImportant("Trigger threshold set from OCDB database entry:  %.02f GeV:", fEThreshold ); 
  }
  return iResult;
}

int AliHLTTriggerCaloClusterEnergy::DoDeinit() {

  // see header file for class documentation
 
  return 0;
}

int AliHLTTriggerCaloClusterEnergy::Reconfigure(const char* cdbEntry, const char* /*chainId*/) {
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) entry=fgkOCDBEntry;

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerCaloClusterEnergy::ScanConfigurationArgument(int argc, const char** argv) {
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -maxpt
  if (argument.CompareTo("-energy")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fEThreshold=argument.Atof();
    return 2;
  }    
  
  // unknown argument
  return -EINVAL;
}

void AliHLTTriggerCaloClusterEnergy::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) {
  // see header file for documentation
  constBase = sizeof(AliHLTTriggerDecision) + sizeof(AliHLTDomainEntry)*14;
  inputMultiplier = 1;
}


void AliHLTTriggerCaloClusterEnergy::GetOCDBObjectDescription( TMap* const targetMap) {
  // Get a list of OCDB object description.
  if (!targetMap) return;
  targetMap->Add(new TObjString(fgkOCDBEntry),
		 new TObjString(Form("%s threshold trigger OCDB object", fDetector.Data()) ) 
		 );
}
