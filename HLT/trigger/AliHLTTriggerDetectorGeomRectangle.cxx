
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
  
  if(trackPos.Eta() >= fEtaMin && 
     trackPos.Eta() <= fEtaMax &&
     trackPos.Phi() >= fPhiMin &&
     trackPos.Phi() <= fPhiMax)
    {
      return true;
    }
  return false;
}

