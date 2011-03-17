/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Author: Svein Lindal <slindal@fys.uio.no>                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliAnaConvCorrBase.cxx
/// @author Svein Lindal
/// @brief  Base class for analysation of conversion particle - track correlations


#include "AliAnaConvCorrBase.h"
#include "AliAODTrack.h"

#include "TClonesArray.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TList.h"
#include "TNtuple.h"
#include "AliAODConversionParticle.h"


#include <iostream>

// Gamma - jet correlation analysis task
// Authors: Svein Lindal


using namespace std;
ClassImp(AliAnaConvCorrBase)

//________________________________________________________________________
AliAnaConvCorrBase::AliAnaConvCorrBase(TString name) : TObject(),
  fName(name),
  fHistograms(NULL),
  fTriggerPt(0), 
  fCorrelatedPt(0),
  fNPhiBins(32)
{
  //Constructor
  fTBins[0] = 7.0;
  fTBins[1] = 5.0;
  fTBins[2] = 3.0;

  fCBins[0] = 5.0;
  fCBins[1] = 3.0;
  fCBins[2] = 1.5;

  for(int iIso = 0; iIso < 2; iIso++) {
    fHEtaPhiPt[iIso] = NULL;
    fHdEta[iIso] = NULL;
    fHdPhi[iIso] = NULL;
    fHNTriggers[iIso] = NULL;
    for(Int_t tBin = 0; tBin < fNTBins; tBin++) {
      for(Int_t cBin = 0; cBin < fNCBins; cBin++) {
	fHdPhiBins[iIso][tBin][cBin] = NULL;
      }
    }
  }

 
}


//________________________________________________________________________________
AliAnaConvCorrBase::~AliAnaConvCorrBase() {
  
  ///destructor
  // delete[] fCBins;
  // delete[] fTBins;

}


//________________________________________________________________________
void AliAnaConvCorrBase::CreateHistograms() {
  //Create histograms add, to outputlis

  fHistograms = new TList();
  fHistograms->SetOwner(kFALSE);
  fHistograms->SetName(fName);



  for(int iIso = 0; iIso < 2; iIso++) {
    fHEtaPhiPt[iIso] = new TH2F(
				  Form("%s_deta_dphi_corr_%s", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
				  Form("%s_deta_dphi_corr_%s", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
				  200, -2., 2., fNPhiBins, -TMath::PiOver2(), 3*TMath::PiOver2());
      fHistograms->Add(fHEtaPhiPt[iIso]);
  }

  for(int iIso = 0; iIso < 2; iIso++) {
    fHdEta[iIso] = new TH1F(Form("%s_deta_corr_%s", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
			    Form("%s_deta_corr_%s", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
			    200, -2, 2);
    fHistograms->Add(fHdEta[iIso]);
  }

  for(int iIso = 0; iIso < 2; iIso++) {
    fHdPhi[iIso] = new TH1F(Form("%s_dphi_corr_%s", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
			    Form("%s_dphi_corr_%s", fName.Data(),  (iIso==0)?"nonIso":"isolated"),
			    fNPhiBins, -TMath::PiOver2(), 3*TMath::PiOver2());
    fHistograms->Add(fHdPhi[iIso]);
  }
  
  for(int iIso = 0; iIso < 2; iIso++) {
      fHNTriggers[iIso] = new TH1I(Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
				   Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
				   fNTBins, 0, fNTBins );
      fHistograms->Add(fHNTriggers[iIso]);
  }


  for(int iIso = 0; iIso < 2; iIso++) {
    for(Int_t tBin = 0; tBin < fNTBins; tBin++) {
      for(Int_t cBin = 0; cBin < fNCBins; cBin++) {
	fHdPhiBins[iIso][tBin][cBin] =  new TH1F(Form("%s_%s_phi_corr_tBin_%d_cBin_%d", fName.Data(), (iIso==0)?"nonIso":"isolated", tBin, cBin),
						 Form("%s particles dPhi %f<tPt< %f  %f<cPt<%f", (iIso==0)?"not isolated":"isolated", GetLowerBinLimit(tBin, fTBins), GetUpperBinLimit(tBin, fTBins), 
						      GetLowerBinLimit(cBin, fCBins), GetUpperBinLimit(cBin, fCBins)),
						 fNPhiBins, -TMath::PiOver2(), 3*TMath::PiOver2());
	fHistograms->Add(fHdPhiBins[iIso][tBin][cBin]);
      }
    }
  }
}


///____________________________________________________________________________
void AliAnaConvCorrBase::FillTriggerCounters(Float_t tPt, Bool_t isolated){ 
  //Fill histogram with trigger counters
  fHNTriggers[0]->Fill(GetTriggerBin(tPt));
  if(isolated)   fHNTriggers[isolated]->Fill(GetTriggerBin(tPt));


}

///_____________________________________________________________________________
void AliAnaConvCorrBase::FillHistograms(Float_t tPt, Float_t cPt, Float_t dPhi, Float_t dEta, Bool_t isolated) {
  //Fill histograms
  Float_t ptFrac = cPt/tPt;

  Int_t triggerBin = GetTriggerBin(tPt);
  Int_t corrBin = GetCorrBin(cPt);

  if(triggerBin < 0 ||corrBin < 0) return;

  fHEtaPhiPt[0]->Fill(dEta, dPhi, cPt);
  fHdEta[0]->Fill(dEta, cPt);
  fHdPhi[0]->Fill(dPhi, cPt);
  fHdPhiBins[0][triggerBin][corrBin]->Fill(dPhi);
  

  if(isolated) {
    fHEtaPhiPt[isolated]->Fill(dEta, dPhi, cPt);
    fHdEta[isolated]->Fill(dEta, ptFrac);
    fHdPhi[isolated]->Fill(dPhi, ptFrac);
    fHdPhiBins[isolated][triggerBin][corrBin]->Fill(dPhi);

  }
}

///____________________________________________________________________________
Int_t AliAnaConvCorrBase::GetTriggerBin(Float_t pt) const {
  //Get trigger bin
  for(Int_t i = 0; i < fNTBins; i++) {
    if(pt > fTBins[i]) return i;
  }

  return -1;
}

///____________________________________________________________________________
Int_t AliAnaConvCorrBase::GetCorrBin(Float_t pt) const {
  //Get correlation particle bin
  for(Int_t i = 0; i < fNCBins; i++) {
    if(pt > fCBins[i]) return i;
  }

  return -1;
}



///______________________________________________________________________________
Float_t AliAnaConvCorrBase::GetLowerBinLimit(const Int_t bin, const Float_t * const bins) const {
  //Get lower bin limit for bin 
  return bins[bin];
}

///______________________________________________________________________________
Float_t AliAnaConvCorrBase::GetUpperBinLimit(const Int_t bin, const Float_t * const bins) const {
  //Get upper bin limit for bin 

  Float_t limit = -999.;
  if(bin < 1)
    limit = 999.;
  else
    limit = bins[bin - 1];

  return limit;

}


//_______________________________________________________________________________

void AliAnaConvCorrBase::PrintStatistics()  { 
  //Print some statistics between each file
  
  for(Int_t i = 0; i < fNTBins; i++) {
    Int_t nTrig = (Int_t) fHNTriggers[0]->GetBinContent(i+1);
    cout << "triggers: " << nTrig << endl;
    for(int j = 0; j < fNCBins; j++) {
      cout << fHdPhiBins[0][i][j]->GetEntries() << "/" << ((nTrig>0)? fHdPhiBins[0][i][j]->GetEntries()/nTrig : 0) << ";  ";
    }
    cout << endl;
  }
}
