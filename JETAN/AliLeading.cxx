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
// Class to find and store the leading particle in event and
// store its correlation to associated particles
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "AliLeading.h"
#include "AliJetReader.h"
#include "AliJetReaderHeader.h"

ClassImp(AliLeading)

////////////////////////////////////////////////////////////////////////

AliLeading::AliLeading():
  fNassoc(0),
  fLeading(0),
  fCorr(0),
  fnBin(45),
  fLow(-TMath::Pi()/2.0),
  fFound(kFALSE)
{
  // Constructor
  fLeading = new TLorentzVector(0.,0.,0.,0.);
  fCorr    = TArrayI(fnBin);
}

////////////////////////////////////////////////////////////////////////

AliLeading::~AliLeading()
{
  // Destructor
  delete fLeading;
}

////////////////////////////////////////////////////////////////////////

void AliLeading::FindLeading(AliJetReader *reader)
{
  // find leading particle in the array of lorentz vectors
  // lvArray and fill the correlation histogram
  
  AliJetReaderHeader* header = reader->GetReaderHeader();
  
  TClonesArray* lvArray = reader->GetMomentumArray();
  Int_t nIn = lvArray->GetEntries();
  
  // find max
  Double_t ptMax = 0.0;
  Int_t idxMax = -1;
  for (Int_t i = 0; i < nIn; i++){
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    if ((reader->GetCutFlag(i) == 1)             &&
	lv->Pt()   > ptMax                       && 
	lv->Eta()  > header->GetFiducialEtaMin() &&
	lv->Eta()  < header->GetFiducialEtaMax()){
      ptMax  = lv->Pt();
      idxMax = i;
    }
  }
  
  if (idxMax == -1) {
    fFound = kFALSE;
    Reset();
    return;
  }
  
  // fill correlation array
  *fLeading = *((TLorentzVector*) lvArray->At(idxMax));
  fFound = kTRUE;
  
  fNassoc = 0;  
  for (Int_t i = 0; i < nIn; i++) {
    if (i == idxMax) continue;
    TLorentzVector *lv = (TLorentzVector*) lvArray->At(i);
    if ( (reader->GetCutFlag(i) == 1) &&
	 lv->Eta()  > header->GetFiducialEtaMin() &&
	 lv->Eta()  < header->GetFiducialEtaMax()) {
      Double_t dphi = fLeading->DeltaPhi(*lv);
      if (dphi < fLow) dphi = 2.0 * TMath::Pi() + dphi;
      // find bin and fill array
      Int_t iBin = (Int_t) 
	TMath::Floor((dphi - fLow)
		     *((Double_t) fnBin) / (2.0 * TMath::Pi()));
      fCorr.AddAt(fCorr.At(iBin)+1,iBin);
      fNassoc++;
    }
  }
}

////////////////////////////////////////////////////////////////////////

void AliLeading::Reset()

{
// Reset leading particle information
  fLeading->SetPxPyPzE(0., 0., 0., 0.);
  fNassoc=0;
  fCorr.Reset();
}

////////////////////////////////////////////////////////////////////////

void AliLeading::PrintLeading()
{
// Print leading particle information
  if (fNassoc<0) {
    cout << " No leading particle in this event" << endl;
    return;
  }
  cout << " Leading particle: " << endl;
  cout << "    (px,py,pz,e) = (" << fLeading->Px() << ","
       << fLeading->Py() << "," << fLeading->Pz() << ","
       << fLeading->E() << ")" << endl;
  cout << "    (pt,eta,phi) = (" << fLeading->Pt() << ","
       << fLeading->Eta() << "," << fLeading->Phi() << ")" << endl;
  cout << "    " << fNassoc << " associated particles." << endl;
}

////////////////////////////////////////////////////////////////////////
 
Double_t AliLeading::GetE()
{
  return fLeading->E();
}
 
////////////////////////////////////////////////////////////////////////
 
Double_t AliLeading::GetPt()
{
  return fLeading->Pt();
}
 
////////////////////////////////////////////////////////////////////////
 
Double_t AliLeading::GetEta()
{
  return fLeading->Eta();
}
 
////////////////////////////////////////////////////////////////////////
 
Double_t AliLeading::GetPhi()
{
  // get phi of leading
  return ( (fLeading->Phi() < 0) ? 
	   (fLeading->Phi()) + 2. * TMath::Pi() : 
	   fLeading->Phi());
}
                                                                                
