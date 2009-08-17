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

/// @file   AliHLTTriggerPhosClusterEnergy.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @date   2009-08-17
/// @brief  HLT energy threshold trigger for PHOS
///      

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTTriggerPhosClusterEnergy.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerPhosClusterEnergy)

AliHLTTriggerPhosClusterEnergy::AliHLTTriggerPhosClusterEnergy() 
  : AliHLTTrigger()
  , fEThreshold(0.0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlts
}

const char* AliHLTTriggerPhosClusterEnergy::fgkOCDBEntry="HLT/ConfigHLT/PhosClusterEnergyTrigger";

AliHLTTriggerPhosClusterEnergy::~AliHLTTriggerPhosClusterEnergy()
{
  // see header file for class documentation
}

const char* AliHLTTriggerPhosClusterEnergy::GetTriggerName() const
{
  // see header file for class documentation
  return "PhosClusterEnergyTrigger";
}

AliHLTComponent* AliHLTTriggerPhosClusterEnergy::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerPhosClusterEnergy;
}

int AliHLTTriggerPhosClusterEnergy::DoTrigger()
{
  // see header file for class documentation

  TString description;

  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  
  if (esd != NULL) {
    esd->GetStdContent();
    
    //Loop over Calorimeter clusters
    Int_t ncc = esd->GetNumberOfCaloClusters();
    for (Int_t i = 0; i < ncc ; i++) {
      AliESDCaloCluster * cluster = esd->GetCaloCluster(i);
      
      // Trigger condition: PHOS clusters with energy > fEThreshold
      if (cluster->IsPHOS() && cluster->E() > fEThreshold) {
	description.Form("Event contains at least one PHOS cluster with energy > %.02f GeV.", fEThreshold);
	SetDescription(description.Data());

	// Enable the detectors for readout.
	GetReadoutList().Enable( AliHLTReadoutList::kPHOS );
	
	// Add the available HLT information for readout too.
	GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "PHOS");
	
	//Set trigger decision
	TriggerEvent(true);
	
	return 0;
      }
    }
  }
  
  // If we got to this point then we did not find any tracks with E > fEThreshold
  // generate negative trigger decision
  description.Form("No PHOS clusters containing energy > %.02f GeV found.", fEThreshold);
  SetDescription(description.Data());
  TriggerEvent(false);
  return 0;
}

int AliHLTTriggerPhosClusterEnergy::DoInit(int argc, const char** argv) {
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

int AliHLTTriggerPhosClusterEnergy::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTriggerPhosClusterEnergy::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) entry=fgkOCDBEntry;

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerPhosClusterEnergy::ScanConfigurationArgument(int argc, const char** argv)
{
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

void AliHLTTriggerPhosClusterEnergy::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase = sizeof(AliHLTTriggerDecision) + sizeof(AliHLTDomainEntry)*14;
  inputMultiplier = 1;
}
