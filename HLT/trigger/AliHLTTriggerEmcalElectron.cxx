// $Id: AliHLTTriggerEmcalElectron.cxx 50471 2011-07-07 09:50:47Z fronchet $
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: marcelfigueredo@gmail.com                             *
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

#include "AliHLTTriggerEmcalElectron.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
#include "AliHLTCaloClusterReader.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "TRefArray.h"
#include "TString.h"
#include "TMap.h"
#include "AliESDtrack.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerEmcalElectron)

AliHLTTriggerEmcalElectron::AliHLTTriggerEmcalElectron() : 
  AliHLTTrigger(),
  fEThreshold(0.0),
  fEoverPThreshold(0.8),
  fEoverPLimit(1.2),
  fClustersRefs(NULL),
  fDetector("EMCAL"),
  fClusterReader(NULL),
  fOCDBEntry(""), 
  fInputDataType()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlts

  fClusterReader = new AliHLTCaloClusterReader();
  fClustersRefs = new TRefArray();

}


AliHLTTriggerEmcalElectron::~AliHLTTriggerEmcalElectron() {
  // see header file for class documentation
  if (fClusterReader)
    delete fClusterReader;
  fClusterReader = NULL;

  if(fClustersRefs)
    delete fClustersRefs;
  fClustersRefs = NULL;
}

//Trigger name
const char* AliHLTTriggerEmcalElectron::GetTriggerName() const {
  // see header file for class documentation
  return "EmcalElectronTrigger";
}

AliHLTComponent* AliHLTTriggerEmcalElectron::Spawn() {
  // see header file for class documentation
  return new AliHLTTriggerEmcalElectron;
}

Int_t AliHLTTriggerEmcalElectron::DoTrigger() {
  // see header file for class documentation
  
  Int_t iResult = 0;


  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  //Try the caloclusterstruct input

/*  
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(fInputDataType); pBlock!=NULL; pBlock=GetNextInputBlock()) {
    AliHLTCaloClusterHeaderStruct *caloClusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
    fClusterReader->SetMemory(caloClusterHeader);
    
    AliHLTCaloClusterDataStruct * caloClusterStruct;
    while( (caloClusterStruct = fClusterReader->NextCluster()) != 0) {
      if (TriggerOnEoverP(caloClusterStruct)) {
	return iResult;
      }
    }
  }*/

  //Try the ESD input
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  
  if (esd != NULL) {
    esd->GetStdContent();

    Int_t ncc = GetClustersFromEsd(esd, fClustersRefs); 
    
    for (Int_t i = 0; i < ncc ; i++) {
      
      AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(fClustersRefs->At(i));
      if(TriggerOnEoverP(cluster,esd)) { 
	return iResult;
      }
    }
  }

  // If we got to this point then we did not find any clusters with E > fEThreshold
  // generate negative trigger decision
  TString description;
  description.Form("No %s clusters corresponding to Energy >  %.02f GeV and %.02f < E/P < %.02f", fDetector.Data(),fEThreshold,fEoverPThreshold,fEoverPLimit);
  SetDescription(description.Data());
  TriggerEvent(false);
  return iResult;

}


template <class T>
Bool_t AliHLTTriggerEmcalElectron::TriggerOnEoverP(T* cluster,AliESDEvent *esd) {
  
  if (cluster->E() > fEThreshold) {
    
      if(!cluster->IsEMCAL()) return kFALSE;
      //Double_t clusterEnergy=cluster->E();  
      //Int_t nmatched=cluster->GetNTracksMatched();    
      Int_t trackindex=cluster->GetTrackMatchedIndex(); 
      if(trackindex<0)return kFALSE;
      AliESDtrack* track = esd->GetTrack(trackindex);
      if(!track)return kFALSE;     
      Double_t EoverP=cluster->E()/track->P(); 
      if(EoverP<fEoverPThreshold)return kFALSE;
      if(EoverP>fEoverPLimit)return kFALSE;

   //We have a cluster satisfying trigger criteria
    TString description;
    description.Form("Event contains at least one %s cluster with energy greater than %.02f corresponding to  %.02f < E/P < %.02f", fDetector.Data(),fEThreshold,fEoverPThreshold,fEoverPLimit);
    SetDescription(description.Data());
    

    // Enable the detectors for readout.

    //GetReadoutList().Enable(AliHLTReadoutList::kPHOS);
//     SetCaloReadoutList("EMCAL");  //FR
    SetCaloReadoutList();  //FR
    
    
    
    // Add the available HLT information for readout too.
    GetTriggerDomain().Add(kAliHLTAnyDataTypeID, fDetector.Data());
    
    //Set trigger decision
    TriggerEvent(kTRUE);
    
    return kTRUE;
  } 


  return kFALSE;

}

int AliHLTTriggerEmcalElectron::DoInit(int argc, const char** argv) {
  // see header file for class documentation

  // first configure the default
  int iResult=ConfigureFromCDBTObjString(fOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0) {
    iResult=ConfigureFromArgumentString(argc, argv);
    HLTImportant("Trigger threshold set from argument string:  %.02f GeV: E/P:%.02f and %.02f", fEThreshold,fEoverPThreshold,fEoverPLimit);   
  } else if ( iResult >=0 ) {
    HLTImportant("Trigger threshold set from OCDB database entry:  %.02f GeV: E/P:%.02f and %.02f", fEThreshold,fEoverPThreshold,fEoverPLimit);  
  }
  return iResult;
}

int AliHLTTriggerEmcalElectron::DoDeinit() {

  // see header file for class documentation
 
  return 0;
}

int AliHLTTriggerEmcalElectron::Reconfigure(const char* cdbEntry, const char* /*chainId*/) {
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) entry=fOCDBEntry;

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerEmcalElectron::ScanConfigurationArgument(int argc, const char** argv) {
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=2;
  TString argument=argv[i];

  // -maxpt
  if (argument.CompareTo("-energy")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fEThreshold=argument.Atof(); // 
    return 2;
  }    

  if (argument.CompareTo("-minEoverP")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fEoverPThreshold=argument.Atof(); // 
    return 2;
  } 
  
    if (argument.CompareTo("-maxEoverP")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fEoverPLimit=argument.Atof(); // 
    return 2;
  } 
  // unknown argument
  return -EINVAL;
}

void AliHLTTriggerEmcalElectron::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier) {
  // see header file for documentation
  constBase = sizeof(AliHLTTriggerDecision) + sizeof(AliHLTDomainEntry)*14;
  inputMultiplier = 1;
}


void AliHLTTriggerEmcalElectron::GetOCDBObjectDescription( TMap* const targetMap) {
  // Get a list of OCDB object description.
  if (!targetMap) return;
  targetMap->Add(new TObjString(fOCDBEntry),
		 new TObjString(Form("%s threshold trigger OCDB object", fDetector.Data()) ) 
		 );
}


