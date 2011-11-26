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
// -- Magali Estienne (magali.estienne@subatech.in2p3.fr)

#include "AliJetUnitArray.h"

class TVector3;
class TLorentzVector;
class TClonesArray;

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
  fVc(0),
  fVn(0),
  fUnitTrackRef(new TRefArray),
  fUnitCellRef(new TRefArray),
  fUnitClusterRef(new TRefArray)
{
  // Default constructor
}  

AliJetUnitArray::AliJetUnitArray(Int_t absId, Int_t esdId, Float_t eta, Float_t phi, Float_t en, Float_t Deta, Float_t Dphi, AliJetFinderUnitDetectorFlagType_t det, AliJetFinderUnitFlagType_t inout, AliJetFinderUnitCutFlagType_t cut, AliJetFinderUnitCutFlagType_t cut2, AliJetFinderUnitSignalFlagType_t /*signal*/,Float_t mass, Int_t clusId):
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
  fUnitSignalFlag(kBad), 
  fUnitDetectorFlag(det),
  fUnitPx(0.),
  fUnitPy(0.),
  fUnitPz(0.),
  fUnitMass(mass),
  fVc(0),
  fVn(0),
  fUnitTrackRef(new TRefArray),
  fUnitCellRef(new TRefArray),
  fUnitClusterRef(new TRefArray)
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
  fVc(0),
  fVn(0),
  fUnitTrackRef(new TRefArray),
  fUnitCellRef(new TRefArray),
  fUnitClusterRef(new TRefArray)
{
  // Constructor 2
}
 
AliJetUnitArray::AliJetUnitArray(const AliJetUnitArray& rUnit):
  TObject(rUnit),
  fUnitEnergy(rUnit.fUnitEnergy),
  fUnitEta(rUnit.fUnitEta),
  fUnitPhi(rUnit.fUnitPhi),
  fUnitDeta(rUnit.fUnitDeta),
  fUnitDphi(rUnit.fUnitDphi),
  fUnitID(rUnit.fUnitID),
  fUnitTrackID(rUnit.fUnitTrackID),
  fUnitNum(rUnit.fUnitNum),
  fUnitClusterID(rUnit.fUnitClusterID),
  fUnitFlag(rUnit.fUnitFlag),
  fUnitCutFlag(rUnit.fUnitCutFlag),
  fUnitCutFlag2(rUnit.fUnitCutFlag2),
  fUnitSignalFlag(rUnit.fUnitSignalFlag),
  fUnitDetectorFlag(rUnit.fUnitDetectorFlag),
  fUnitPx(rUnit.fUnitPx),
  fUnitPy(rUnit.fUnitPy),
  fUnitPz(rUnit.fUnitPz),
  fUnitMass(rUnit.fUnitMass),
  fVc(rUnit.fVc),
  fVn(rUnit.fVn),
  fUnitTrackRef(rUnit.fUnitTrackRef),
  fUnitCellRef(rUnit.fUnitCellRef),
  fUnitClusterRef(rUnit.fUnitClusterRef)
{
  // Copy constructor
}

AliJetUnitArray& AliJetUnitArray::operator=(const AliJetUnitArray& rhs)
{
  // Assignment
    if (this != &rhs) {
	fUnitEnergy = rhs.fUnitEnergy;
	fUnitEta = rhs.fUnitEta;
	fUnitPhi = rhs.fUnitPhi;
	fUnitDeta = rhs.fUnitDeta;
	fUnitDphi = rhs.fUnitDphi;
	fUnitID = rhs.fUnitID;
	fUnitTrackID = rhs.fUnitTrackID;
	fUnitNum = rhs.fUnitNum;
	fUnitClusterID = rhs.fUnitClusterID;
	fUnitFlag = rhs.fUnitFlag;
	fUnitCutFlag = rhs.fUnitCutFlag;
	fUnitCutFlag2 = rhs.fUnitCutFlag2;
	fUnitSignalFlag = rhs.fUnitSignalFlag;
	fUnitDetectorFlag = rhs.fUnitDetectorFlag;
	fUnitPx = rhs.fUnitPx;
	fUnitPy = rhs.fUnitPy;
	fUnitPz = rhs.fUnitPz;
	fUnitMass = rhs.fUnitMass;
	fVc = rhs.fVc;
	fVn = rhs.fVn;
	fUnitTrackRef = rhs.fUnitTrackRef;
	fUnitCellRef = rhs.fUnitCellRef;
	fUnitClusterRef = rhs.fUnitClusterRef;
    }
  return *this;

}


//------------------------------------------------------------------------
AliJetUnitArray::~AliJetUnitArray()
{
  // Destructor 
  delete fUnitTrackRef;
  delete fUnitCellRef;
  delete fUnitClusterRef;

}

void AliJetUnitArray::ClearUnitTrackRef()
{
  fUnitTrackRef->Clear();
}

void AliJetUnitArray::ClearUnitCellRef()
{
  fUnitCellRef->Clear();
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
Float_t  AliJetUnitArray::EtaToTheta(Float_t arg) const
{
  // Eta to theta transformation
  return 2.*atan(exp(-arg));
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::operator>(const AliJetUnitArray* unit) const
{
  // Greater than operator used by sort
  if( fUnitEnergy > unit->GetUnitEnergy())
    return kTRUE;
  else 
    return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::operator<(const AliJetUnitArray* unit) const
{
  // Less than operator used by sort
  if( fUnitEnergy < unit->GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}

//------------------------------------------------------------------------
Bool_t AliJetUnitArray::operator==(const AliJetUnitArray* unit) const
{
  // equality operator used by sort
  if( fUnitEnergy == unit->GetUnitEnergy())
    return kTRUE;
  else
    return kFALSE;
}
