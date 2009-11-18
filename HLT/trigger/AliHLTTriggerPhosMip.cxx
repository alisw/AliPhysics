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

/// @file   AliHLTTriggerPhosMip.cxx
/// @author Svein Lindal
/// @date   2009-08-19
/// @brief  HLT Minimum Ionizing Particle (MIP) trigger for PHOS

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTTriggerPhosMip.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerPhosMip)

AliHLTTriggerPhosMip::AliHLTTriggerPhosMip() 
  : AliHLTTrigger()
  , fEMin(0.0)
  , fEMax(0.0)
  , fNCellsMax(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlts
}

const char* AliHLTTriggerPhosMip::fgkOCDBEntry="HLT/ConfigHLT/PhosMipTrigger";

AliHLTTriggerPhosMip::~AliHLTTriggerPhosMip()
{
  // see header file for class documentation
}

const char* AliHLTTriggerPhosMip::GetTriggerName() const
{
  // see header file for class documentation
  return "PhosMipTrigger";
}

AliHLTComponent* AliHLTTriggerPhosMip::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerPhosMip;
}

int AliHLTTriggerPhosMip::DoTrigger()
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
      
      // Trigger condition: PHOS clusters within energy range, cells < fNCellsMax
      if ( cluster->IsPHOS() && 
	   cluster->E() > fEMin && 
	   cluster->E() < fEMax && 
	   cluster->GetNCells() <= fNCellsMax ) {
	  
	description.Form("Event contains at least one PHOS cluster satisfying MIP criteria");
	SetDescription(description.Data());
	
	// Enable the detectors for readout.
	GetReadoutList().Enable( AliHLTReadoutList::kPHOS );
	
	// Add the available HLT information for readout too.
	GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "PHOS");
	
	//Set trigger decision
	TriggerEvent(true);
	
	return 0;
      } /// if trigger criteria
    }  /// cluster loop
  }
  
  // If we got to this point then we did not find any good MIP candidates
  // generate negative trigger decision
  description.Form("No PHOS clusters satisfying MIP criteria found");
  SetDescription(description.Data());
  TriggerEvent(false);
  return 0;
}

int AliHLTTriggerPhosMip::DoInit(int argc, const char** argv) {
  // see header file for class documentation

  // first configure the default
  int iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0) {
    HLTInfo("Trigger configuration from OCDB database entry:  emin %.03f GeV \n emax %.03f, \n ncells: %i", fEMin, fEMax, fNCellsMax ); 
    iResult=ConfigureFromArgumentString(argc, argv);
    HLTInfo("Trigger configuration overwritten from command line:  emin %.03f GeV \n emax %.03f, \n ncells: %i", fEMin, fEMax, fNCellsMax ); 
   } else if ( iResult >=0 ) {
    HLTInfo("Trigger configuration from OCDB database entry:  emin %.03f GeV \n emax %.03f, \n ncells: %i", fEMin, fEMax, fNCellsMax ); 
  }
  return iResult;
}

int AliHLTTriggerPhosMip::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTriggerPhosMip::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) entry=fgkOCDBEntry;

  return ConfigureFromCDBTObjString(entry);
}

int AliHLTTriggerPhosMip::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -mine
  if (argument.CompareTo("-emin")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fEMin=argument.Atof();
    return 2;
  }    

  ///-maxe
  else if (argument.CompareTo("-emax")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fEMax=argument.Atof();
    return 2;
  }    

  ///-ncells
  else if (argument.CompareTo("-ncells")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fNCellsMax=argument.Atoi();
    return 2;
  }    
  
  // unknown argument
  return -EINVAL;
}

void AliHLTTriggerPhosMip::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase = sizeof(AliHLTTriggerDecision) + sizeof(AliHLTDomainEntry)*14;
  inputMultiplier = 1;
}
