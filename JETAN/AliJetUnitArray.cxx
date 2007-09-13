/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//------------------------------------------------------------------
//  Unit used by UA1 algorithm
//  Authors: Sarah Blyth (LBL/UCT)
//           Magali Estienne (IReS) (new version for JETAN)
//------------------------------------------------------------------

#include "AliJetUnitArray.h"


ClassImp(AliJetUnitArray)

AliJetUnitArray::AliJetUnitArray():
  fUnitEnergy(0.0),
  fUnitEta(0.0),
  fUnitPhi(0.0),
  fUnitDeta(0.),
  fUnitDphi(0.),
  fUnitID(0),
  fUnitTrackID(0),
  fUnitNum(0),
  fUnitClusterID(0),
  fUnitFlag(kOutJet),
  fUnitCutFlag(kPtSmaller),
  fUnitSignalFlag(kBad), 
  fUnitDetectorFlag(kTpc),
  fUnitPx(0.),
  fUnitPy(0.),
  fUnitPz(0.),
  fUnitMass(0.)
{
  // Default constructor
}  

AliJetUnitArray::AliJetUnitArray(Int_t absId, Int_t esdId, Float_t eta, Float_t phi, Float_t en, Float_t px, Float_t py, Float_t pz, Float_t Deta, Float_t Dphi, AliJetFinderUnitDetectorFlagType_t det, AliJetFinderUnitFlagType_t inout, AliJetFinderUnitCutFlagType_t cut, Float_t mass, Int_t clusId):
  fUnitEnergy(en),
  fUnitEta(eta),
  fUnitPhi(phi),
  fUnitDeta(Deta),
  fUnitDphi(Dphi),
  fUnitID(absId),
  fUnitTrackID(esdId),
  fUnitNum(0),
  fUnitClusterID(clusId),
  fUnitFlag(inout),
  fUnitCutFlag(cut),
  fUnitSignalFlag(kBad), 
  fUnitDetectorFlag(det),
  fUnitPx(px),
  fUnitPy(py),
  fUnitPz(pz),
  fUnitMass(mass)
{
  // Constructor 2
}

AliJetUnitArray::~AliJetUnitArray()
{
  // Destructor 
}
	
Bool_t AliJetUnitArray::operator>(AliJetUnitArray unit) const
{
  // Greater than operator used by sort
  if( fUnitEnergy > unit.GetUnitEnergy())
    return kTRUE;
  else 
    return kFALSE;
}

Bool_t AliJetUnitArray::operator<( AliJetUnitArray unit) const
{
  // Less than operator used by sort
  if( fUnitEnergy < unit.GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}

Bool_t AliJetUnitArray::operator==( AliJetUnitArray unit) const
{
  // equality operator used by sort
  if( fUnitEnergy == unit.GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}
