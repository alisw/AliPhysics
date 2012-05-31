
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Oystein Djuvsland                                     *
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

/// @file   AliHLTTriggerDetectorGeomRectangle.cxx
/// @author Oystein Djuvsland
/// @date   2009-10-08
/// @brief  HLT class describing simple rectangular geometry of (sub-)detectors.
///         Used for the AliHLTTriggerBarrelGeomMultiplicity classes

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerDetectorGeomRectangle.h"
#include "AliHLTTriggerDetectorGeom.h"
#include "TVector3.h"
#include "TMath.h"


AliHLTTriggerDetectorGeomRectangle::AliHLTTriggerDetectorGeomRectangle() : AliHLTTriggerDetectorGeom()
{
  // See header file for class documentation
}

AliHLTTriggerDetectorGeomRectangle::~AliHLTTriggerDetectorGeomRectangle()
{
  // See header file for class documentation
}

Bool_t AliHLTTriggerDetectorGeomRectangle::IsInDetector(Double_t point[3])
{
  // See header file for class documentation
  TVector3 trackPos(point);
  
  Double_t phi = 0;

  if(trackPos.Phi() < 0) phi = trackPos.Phi() + 2*TMath::Pi();
  else phi = trackPos.Phi();

  if(trackPos.Eta() >= fEtaMin && 
     trackPos.Eta() <= fEtaMax &&
     phi >= fPhiMin &&
     phi <= fPhiMax)
    {
      //      printf("Checking point: Eta: %f Phi: %f against cuts EtaMin: %f EtaMax: %f PhiMin: %f PhiMax: %f - true\n", trackPos.Eta(), phi, fEtaMin, fEtaMax, fPhiMin, fPhiMax);
      return true;
    }

  //  printf("Checking point: Eta: %f Phi: %f against cuts EtaMin: %f EtaMax: %f PhiMin: %f PhiMax: %f - false\n", trackPos.Eta(), phi, fEtaMin, fEtaMax, fPhiMin, fPhiMax);

  if(trackPos.Eta() >= fEtaMin && 
     trackPos.Eta() <= fEtaMax &&
     phi >= fPhiMin &&
     phi <= fPhiMax)
    {
      //      printf("Checking point: Eta: %f Phi: %f against cuts EtaMin: %f EtaMax: %f PhiMin: %f PhiMax: %f - true\n", trackPos.Eta(), phi, fEtaMin, fEtaMax, fPhiMin, fPhiMax);
      return true;
    }

  //  printf("Checking point: Eta: %f Phi: %f against cuts EtaMin: %f EtaMax: %f PhiMin: %f PhiMax: %f - false\n", trackPos.Eta(), phi, fEtaMin, fEtaMax, fPhiMin, fPhiMax);

  return false;
}

