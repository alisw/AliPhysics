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


/* $Id$ */

//_________________________________________________________________________
//  Unit used by UA1 algorithm
// --
//*-- Author: Sarah Blyth (LBL/UCT)
// --
// Revised Version for JETAN 
// -- Magali Estienne (IReS)

//#include <vector>

#include <TVector3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>

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
  fUnitCutFlag2(kPtSmaller),
  fUnitSignalFlag(kBad), 
  fUnitDetectorFlag(kTpc),
  fUnitPx(0.),
  fUnitPy(0.),
  fUnitPz(0.),
  fUnitMass(0.),
  fV(0),
  fVc(0),
  fVn(0),
  fVet(0)
{
  // Default constructor
}  

AliJetUnitArray::AliJetUnitArray(Int_t absId, Int_t esdId, Float_t eta, Float_t phi, Float_t en, Float_t Deta, Float_t Dphi, AliJetFinderUnitDetectorFlagType_t det, AliJetFinderUnitFlagType_t inout, AliJetFinderUnitCutFlagType_t cut, AliJetFinderUnitCutFlagType_t cut2, AliJetFinderUnitSignalFlagType_t signal,Float_t mass, Int_t clusId):
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
  fUnitCutFlag2(cut2),
  fUnitSignalFlag(signal), 
  fUnitDetectorFlag(det),
  fUnitPx(0.),
  fUnitPy(0.),
  fUnitPz(0.),
  fUnitMass(mass),
  fV(0),
  fVc(0),
  fVn(0),
  fVet(0)
{
  //abs ID (in a eta,phi grid, track ID in ESD, eta, phi, energy, px, py, pz, Deta, Dphi, detector flag, in/out jet, mass

  // Constructor 2
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
  fUnitCutFlag2(kPtSmaller),
  fUnitSignalFlag(kBad),
  fUnitDetectorFlag(det),
  fUnitPx(px),
  fUnitPy(py),
  fUnitPz(pz),
  fUnitMass(mass),
  fV(0),
  fVc(0),
  fVn(0),
  fVet(0)
{
  // Constructor 2
}
 
//------------------------------------------------------------------------
AliJetUnitArray::~AliJetUnitArray()
{
  // Destructor 
}

//------------------------------------------------------------------------
void AliJetUnitArray::SetUnitSignalFlagC(Bool_t init, AliJetFinderUnitSignalFlagType_t flag)
{
  // Set signal flag of the charged particle
  if(init){
    if(!fVc.empty())
      fVc.clear();
  }
  else fVc.push_back(flag);
}

//------------------------------------------------------------------------
void AliJetUnitArray::SetUnitSignalFlagN(Bool_t init, AliJetFinderUnitSignalFlagType_t flag)
{
  // Set signal flag of the neutral cell
  if(init){
    if(!fVn.empty())
      fVn.clear();
  }
  else fVn.push_back(flag);
}

//------------------------------------------------------------------------
void AliJetUnitArray::SetUnitEtN(Bool_t init, Float_t et)
{
  // Set transverse energy of the neutral cell
  if(init){
    if(!fVet.empty())
      fVet.clear();
  }
  else fVet.push_back(et);
}


//------------------------------------------------------------------------
void AliJetUnitArray::SetUnitPxPyPz(Bool_t init, vector<Float_t> v3)
{
  // Set momentum components of the charged particle
  if(init)
    {
      if(!fV.empty()){
	fV.clear();
      }
    }
  else{
    fV.push_back(v3);
  }
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::GetUnitSignalFlagC(Int_t ind, AliJetFinderUnitSignalFlagType_t &flagc)
{
  // Get signal flag of the charged particle
  if(ind <= (Int_t)fVc.size())
    {
      flagc = (AliJetFinderUnitSignalFlagType_t)fVc[ind];
      return kTRUE;
    }
  else return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::GetUnitSignalFlagN(Int_t ind, AliJetFinderUnitSignalFlagType_t &flagn)
{
  // Get signal flag of the neutral cell
  if(ind <= (Int_t)fVn.size())
    {
      flagn = (AliJetFinderUnitSignalFlagType_t)fVn[ind];
      return kTRUE;
    }
  else return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::GetUnitEtN(Int_t ind, Float_t &et)
{
  // Get transverse energy of the neutral cell
  if(ind <= (Int_t)fVet.size())
    {
      et = (Float_t)fVet[ind];
      return kTRUE;
    }
  else return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::GetUnitPxPyPz(Int_t ind, Float_t &px, Float_t &py, Float_t &pz)
{
  // Get momentum components of the charged particle
  if(ind <= (Int_t)fV.size())
    {
      px = (Float_t)fV[ind][0];
      py = (Float_t)fV[ind][1];
      pz = (Float_t)fV[ind][2];
      return kTRUE;
    }
  else return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::GetUnitPxPyPzE(Int_t ind, Float_t &px, Float_t &py, Float_t &pz, Float_t &en)
{
// Get 4-momentum components of the charged particle
  if(ind <= (Int_t)fV.size())
    {
      px = (Float_t)fV[ind][0];
      py = (Float_t)fV[ind][1];
      pz = (Float_t)fV[ind][2];
      en = TMath::Sqrt(px*px+py*py+pz*pz);
      return kTRUE;
    }
  else return kFALSE;
}

//------------------------------------------------------------------------
Float_t  AliJetUnitArray::EtaToTheta(Float_t arg) const
{
  // Eta to theta transformation
  return 2.*atan(exp(-arg));
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::operator>(AliJetUnitArray unit) const
{
  // Greater than operator used by sort
  if( fUnitEnergy > unit.GetUnitEnergy())
    return kTRUE;
  else 
    return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::operator<( AliJetUnitArray unit) const
{
  // Less than operator used by sort
  if( fUnitEnergy < unit.GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::operator==( AliJetUnitArray unit) const
{
  // equality operator used by sort
  if( fUnitEnergy == unit.GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}
