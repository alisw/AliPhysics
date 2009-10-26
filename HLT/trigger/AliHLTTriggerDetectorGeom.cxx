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

/// @file   AliHLTTriggerDetectorGeom.cxx
/// @author Oystein Djuvsland
/// @date   2009-10-08
/// @brief  HLT class describing simple geometry of (sub-)detectors.
///         Used for the AliHLTTriggerBarrelGeomMultiplicity class

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerDetectorGeom.h"
#include <ostream>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerDetectorGeom)

AliHLTTriggerDetectorGeom::AliHLTTriggerDetectorGeom()
: TObject(),
  fEtaMin(0),
  fEtaMax(0),
  fPhiMin(0),
  fPhiMax(0),
  fName('\0')
{
  // See header file for class documentation
  for(Int_t i = 0; i < 3; i++) 
    {
      fInitalPoint[i] = 0;
      fNormVector[i] = 0;
    }
}
  
AliHLTTriggerDetectorGeom::~AliHLTTriggerDetectorGeom()
{
  // See header file for class documentation
}

void AliHLTTriggerDetectorGeom::SetInitialPoint(Double_t *point)
{
  // See header file for class documentation
  for(int i = 0; i < 3; i++)
    {
      fInitalPoint[i] = point[i];
    }
}


void AliHLTTriggerDetectorGeom::SetNormVector(Double_t *nVector)
{
  // See header file for class documentation
  for(int i = 0; i < 3; i++)
    {
      fNormVector[i] = nVector[i];
    }
}

void AliHLTTriggerDetectorGeom::GetInitialPoint(Double_t *point)
{
  // See header file for class documentation
  for(int i = 0; i < 3; i++)
    {
      point[i] = fInitalPoint[i];
    }
}

void AliHLTTriggerDetectorGeom::GetNormVector(Double_t *vec)
{
  // See header file for class documentation
  for(int i = 0; i < 3; i++)
    {
      vec[i] = fNormVector[i];
    }
}

void AliHLTTriggerDetectorGeom::PrintDetectorGeom(std::ostream &out)
{
  // See header file for class documentation

  out << "Name: " << fName << std::endl;
  out << "Eta Min: " << fEtaMin << std::endl;
  out << "Eta Max: " << fEtaMax << std::endl;
  out << "Phi Min: " << fPhiMin << std::endl;
  out << "Phi Max: " << fPhiMax << std::endl;
  out << "Initial Point: {" << fInitalPoint[0] << ", " << fInitalPoint[1] << ", " << fInitalPoint[2] << "}" << std::endl; 
  out << "Normal Vector: {" << fNormVector[0] << ", " << fNormVector[1] << ", " << fNormVector[2] << "}" << std::endl; 
}

