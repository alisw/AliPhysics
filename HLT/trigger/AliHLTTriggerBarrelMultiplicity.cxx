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

/// @file   AliHLTTriggerBarrelMultiplicity.cxx
/// @author Matthias Richter
/// @date   2009-06-30
/// @brief  HLT trigger component for charged particle multiplicity in
///         the central barrel.

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerBarrelMultiplicity.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerBarrelMultiplicity)

AliHLTTriggerBarrelMultiplicity::AliHLTTriggerBarrelMultiplicity()
  : AliHLTTrigger()
  , fPtMin(0.0)
  , fPtMax(0.0)
  , fMinTracks(1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTriggerBarrelMultiplicity::~AliHLTTriggerBarrelMultiplicity()
{
  // see header file for class documentation
}

const char* AliHLTTriggerBarrelMultiplicity::GetTriggerName() const
{
  // see header file for class documentation
  return "BarrelMultiplicityTrigger";
}

AliHLTComponent* AliHLTTriggerBarrelMultiplicity::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerBarrelMultiplicity;
}

int AliHLTTriggerBarrelMultiplicity::DoTrigger()
{
  // see header file for class documentation
  const TObject* obj = GetFirstInputObject(kAliHLTAllDataTypes, "AliESDEvent");
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(obj));
  TString description;
  TString ptcut;
  if (esd != NULL) {
    esd->GetStdContent();
    
    unsigned int numberOfTracks=0;
    for (Int_t i = 0; i < esd->GetNumberOfTracks(); i++) {
      AliESDtrack* track = esd->GetTrack(i);
      if (track && track->Pt() >= fPtMin &&
	  (fPtMax<=fPtMin || track->Pt() < fPtMax)) {
	numberOfTracks++;
      }
    }

    if (fPtMax>fPtMin) {
      ptcut.Form(" %.02f GeV/c <= pt < %.02f GeV/c", fPtMin, fPtMax);
    } else {
      ptcut.Form(" pt >= %.02f GeV/c", fPtMin);
    }
    if (numberOfTracks>=fMinTracks) {
      description.Form("Event contains %d track(s) with ", numberOfTracks);
      description+=ptcut;
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
      GetTriggerDomain().Add("CLUSTERS", "TPC ");
      TriggerEvent(true);
      return 0;
    }
  }
  description.Form("No tracks matching the tresholds found in the central barrel (min tracks %d, %s)",
		   fMinTracks, ptcut.Data());
  SetDescription(description.Data());
  TriggerEvent(false);
  return 0;
}
