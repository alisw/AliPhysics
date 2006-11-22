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
 
//---------------------------------------------------------------------
// Jet class 
// Stores the output of a jet algorithm
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------
 
#include <Riostream.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "AliJet.h"
ClassImp(AliJet)
  
AliJet::AliJet():
  fNInput(0),
  fNJets(0),
  fEtAvg(0),
  fInJet(0),
  fMultiplicities(0),
  fNCells(0),
  fPtFromSignal(0),
  fJets(0),
  fEtaIn(0),
  fPhiIn(0),
  fPtIn(0)
{
  // Default constructor
  fJets = new TClonesArray("TLorentzVector",1000);
  fInJet = TArrayI();
  fPtIn = TArrayF();
  fEtaIn = TArrayF();
  fPhiIn = TArrayF();
  fPtFromSignal = TArrayF();
  fMultiplicities = TArrayI();
  fNCells = TArrayI();
} 

////////////////////////////////////////////////////////////////////////

AliJet::~AliJet()
{
  // destructor
  if (fJets) {
    fJets->Delete();
    delete fJets;
  }
}

////////////////////////////////////////////////////////////////////////

Bool_t AliJet::OutOfRange(Int_t i, const char *s) const
{
  // checks if i is a valid index. s = name of calling method
  if (i >= fNJets || i < 0) {
    cout << s << " Index " << i << " out of range" << endl;
    return kTRUE;
  }
  return kFALSE;
}

////////////////////////////////////////////////////////////////////////

TLorentzVector* AliJet::GetJet(Int_t i)
{
  // returns i-jet
  if (OutOfRange(i, "AliJet::GetJet:")) return 0;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv; 
}

////////////////////////////////////////////////////////////////////////

Int_t AliJet::GetMultiplicity(Int_t i) const
{
  // gets multiplicity of i-jet
  if (OutOfRange(i, "AliJet::GetMultiplicity:")) return 0;
  return fMultiplicities[i];
}

////////////////////////////////////////////////////////////////////////

Int_t AliJet::GetNCell(Int_t i) const
{
  // gets number of cell of i-jet
  if (OutOfRange(i, "AliJet::GetNCell:")) return 0;
  return fNCells[i];
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetPx(Int_t i)
{
// Get Px component of jet i
  if (OutOfRange(i, "AliJet::GetPx:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->Px();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetPy(Int_t i)
{
// Get Py component of jet i
  if (OutOfRange(i, "AliJet::GetPy:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->Py();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetPz(Int_t i)
{
// Get Pz component of jet i
  if (OutOfRange(i, "AliJet::GetPz:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->Pz();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetP(Int_t i)
{
// Get momentum of jet i
  if (OutOfRange(i, "AliJet::GetP:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->P();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetE(Int_t i)
{
// Get energy of jet i
  if (OutOfRange(i, "AliJet::GetE:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->E();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetPt(Int_t i)
{
// Get transverse momentum of jet i
  if (OutOfRange(i, "AliJet::GetPt:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->Pt();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetEta(Int_t i)
{
// Get eta of jet i
  if (OutOfRange(i, "AliJet::GetEta:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->Eta();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetPhi(Int_t i)
{
// Get phi of jet i
  if (OutOfRange(i, "AliJet::GetPhi:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return ( (lv->Phi() < 0) ? (lv->Phi()) + 2. * TMath::Pi() : lv->Phi());
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetTheta(Int_t i)
{
// Get theta of jet i
  if (OutOfRange(i, "AliJet::GetTheta:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->Theta();
}

////////////////////////////////////////////////////////////////////////

Double_t AliJet::GetMass(Int_t i)
{
// Get invariant mass of jet i
  if (OutOfRange(i, "AliJet::GetMass:")) return -1e30;
  TLorentzVector *lv = (TLorentzVector*) fJets->At(i);
  return lv->M();
}

////////////////////////////////////////////////////////////////////////
 

void AliJet::AddJet(Double_t px, Double_t py, Double_t pz, Double_t e)
{
// Add new jet to the list
  new ((*fJets)[fNJets++]) TLorentzVector(px,py,pz,e);
}

////////////////////////////////////////////////////////////////////////

void AliJet::SetInJet(Int_t* j)
{
  // set information of which input object belongs
  // to each jet. If zero, object was not assigned to
  // a jet, if n,positive, it was assiged to jet n
  // if n, negative, it is within cone of jet n, but
  // it did not passed the user cuts. filled in by AliJetFinder
  if (fNInput>0) fInJet.Set(fNInput, j);
}

////////////////////////////////////////////////////////////////////////

void AliJet::SetEtaIn(Float_t* r)
{
  if (fNInput>0) fEtaIn.Set(fNInput, r);
}

////////////////////////////////////////////////////////////////////////

void AliJet::SetPtIn(Float_t* pt)
{
  if (fNInput>0) fPtIn.Set(fNInput, pt);
}

////////////////////////////////////////////////////////////////////////

void AliJet::SetPhiIn(Float_t* x)
{
  if (fNInput>0) fPhiIn.Set(fNInput, x);
}

////////////////////////////////////////////////////////////////////////

void AliJet::SetPtFromSignal(Float_t* p)
{
  // set information of percentage of pt of jets
  // coming from signal (ie Pythia)
  if (fNJets>0) fPtFromSignal.Set(fNJets, p);
}

////////////////////////////////////////////////////////////////////////

void AliJet::SetMultiplicities(Int_t* m)
{
  // set information of jet multiplicities
  // filled in by AliJetFinder
  if (fNJets>0) fMultiplicities.Set(fNJets, m);
}

////////////////////////////////////////////////////////////////////////

void AliJet::SetNCells(Int_t* n)
{
  if (fNJets>0) fNCells.Set(fNJets, n);
}

////////////////////////////////////////////////////////////////////////

void AliJet::ClearJets(Option_t *option)
{
  // reset all values
  fJets->Clear(option);
  fNInput=0;
  fNJets=0;
  fMultiplicities.Set(0);
  fInJet.Set(0); 
  fPtFromSignal.Set(0);
  fPhiIn.Set(0);
  fEtaIn.Set(0);
  fPtIn.Set(0);
  fNCells.Set(0);
}

////////////////////////////////////////////////////////////////////////

void AliJet::PrintJets()
{
// Print jet information
  if (fNJets == 0) {
    cout << " AliJet::PrintJets: There are no jets in this event " << endl;
    return;
  }
  cout << " AliJet::PrintJets: There are " << fNJets
       << " jets in this event" << endl;
  for(Int_t i=0;i<fNJets;i++) {
    cout << "   Jet " << i << " (px,py,pz,en)=(" << GetPx(i)
	 << "," << GetPy(i)
	 << "," << GetPz(i)
	 << "," << GetE(i) 
	 << ")" << endl;
    cout << "         (pt,eta,phi)=(" << GetPt(i)
	 << "," << GetEta(i)
	 << "," << GetPhi(i) << ")" << endl;
    cout << "         # of tracks =" << GetMultiplicity(i) << endl;
  }
}
