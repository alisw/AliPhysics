// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kyrre Skjerdal                                        *
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

/// @file   AliHLTUpcTriggerComponent.cxx
/// @author Kyrre Skjerdal
/// @date   2010-04-16
/// @brief  HLT trigger component for Ultra-Peripheral Collisions

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTUpcTriggerComponent.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"
ClassImp(AliHLTUpcTriggerComponent)

const char* AliHLTUpcTriggerComponent::GetTriggerName() const
{
  //See header file for documentation
  return "UpcTrigger";
}

AliHLTComponent* AliHLTUpcTriggerComponent::Spawn()
{
  //See header file for documentation
  return new AliHLTUpcTriggerComponent;
}

int AliHLTUpcTriggerComponent::DoTrigger()
{
  //See header file for documentation
  HLTInfo("Entering DoTrigger()");
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  
  if (esd != NULL){
      esd->GetStdContent();     
      //We skip events with no primary vertex reconstructed
      if(!(PrimaryVertexReconstructed(esd))){
	HLTWarning("No primary vertex reconstructed");
	SetDescription("No reconstructed primary vertex. Not a candidate for an Ultra-peripheral Collison");
	TriggerEvent(false);
	return 0;
      }
      
      Int_t nGoodRec = 0;
      Int_t charge1 = 0;
      Int_t charge2 = 0;
      
      for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++){
	cout << "Track number: " << esd->GetNumberOfTracks() << endl;
	AliESDtrack* track = esd->GetTrack(i);

	Int_t nItsClusters = track->GetNcls(0);
	Int_t nTpcClusters = track->GetNcls(1);

	Float_t bxy = 0;
	Float_t bz = 0;
	track->GetImpactParametersTPC(bxy, bz);
	
	//Check if the track is comming from the primary vertex
	Bool_t isPrimary = kFALSE;
	if(fabs(bxy) < 1.5 && fabs(bz) < 1.5){
	  isPrimary = kTRUE;
	}
	
	Int_t charge = track->Charge();    
	
	if(nItsClusters > 3 && nTpcClusters > 50 && isPrimary){
	  nGoodRec++;
	  if(nGoodRec == 1){
	    charge1 = charge;
	  } else if(nGoodRec == 2){
	    charge2 = charge;
	  }
	}
	
      }
      
      //Demand two good tracks with opposite charge
      if(nGoodRec == 2){
	if(charge1 == -charge2){
	  SetDescription("Event is a candidate for an Ultra-peripheral Collision.");
	  // Enable the central detectors for readout.
	  GetReadoutList().Enable(

				  AliHLTReadoutList::kITSSPD |
				  AliHLTReadoutList::kITSSDD |
				  AliHLTReadoutList::kITSSSD |
				  AliHLTReadoutList::kTPC    |
				  AliHLTReadoutList::kTRD    |
				  AliHLTReadoutList::kTRD    |
				  AliHLTReadoutList::kTOF    |
				  AliHLTReadoutList::kHMPID  |
				  AliHLTReadoutList::kPHOS   
				  );
	  // Add the available HLT information for readout too.
	  GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "ISPD");
	  GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "ISDD");
	  GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "ISSD");
	  GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "ITPC");
	  GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "ITRD");
	  GetTriggerDomain().Add(kAliHLTAnyDataTypeID, "IPHOS");
	  TriggerEvent(true);
	  return 0;
	}
      }
      
  }else{
    HLTFatal("No ESD found");
    SetDescription("Not a candidate for an Ultra-peripheral Collision.");
    TriggerEvent(false);
    return 0;
  }
  SetDescription("Not a candidate for an Ultra-peripheral Collision.");
  TriggerEvent(false);
  return 0;
}

void AliHLTUpcTriggerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //See header file for documentation
  
  constBase = sizeof(AliHLTTriggerDecision) + sizeof(AliHLTDomainEntry)*14;
  inputMultiplier = 1;
}


Bool_t AliHLTUpcTriggerComponent::PrimaryVertexReconstructed(const AliESDEvent *event) const
{ 
  //See header file for documentation
  
  return (event->GetPrimaryVertex()->GetNContributors() > 0);
}
