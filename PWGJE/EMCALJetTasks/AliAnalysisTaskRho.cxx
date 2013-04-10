// $Id$
//
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fRhoName".apppend("_Scaled").
//
// Authors: R.Reed, S.Aiola

#include "AliAnalysisTaskRho.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"

ClassImp(AliAnalysisTaskRho)

//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho() : 
  AliAnalysisTaskRhoBase("AliAnalysisTaskRho"),
  fNExclLeadJets(0)
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name, histo),
  fNExclLeadJets(0)
{
  // Constructor.
}


//________________________________________________________________________
Bool_t AliAnalysisTaskRho::Run() 
{
  // Run the analysis.

  fRho->SetVal(0);
  if (fRhoScaled)
    fRhoScaled->SetVal(0);

  if (!fJets)
    return kFALSE;

  const Int_t Njets   = fJets->GetEntries();

  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  if (fNExclLeadJets > 0) {
    for (Int_t ij = 0; ij < Njets; ++ij) {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(ij));
      if (!jet) {
	AliError(Form("%s: Could not receive jet %d", GetName(), ij));
	continue;
      } 

      if (!AcceptJet(jet))
        continue;

      if (jet->Pt() > maxJetPts[0]) {
	maxJetPts[1] = maxJetPts[0];
	maxJetIds[1] = maxJetIds[0];
	maxJetPts[0] = jet->Pt();
	maxJetIds[0] = ij;
      } else if (jet->Pt() > maxJetPts[1]) {
	maxJetPts[1] = jet->Pt();
	maxJetIds[1] = ij;
      }
    }
    if (fNExclLeadJets < 2) {
      maxJetIds[1] = -1;
      maxJetPts[1] = 0;
    }
  }

  static Double_t rhovec[999];
  Int_t NjetAcc = 0;

  // push all jets within selected acceptance into stack
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {

    // exlcuding lead jets
    if (iJets == maxJetIds[0] || iJets == maxJetIds[1])
      continue;

    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
    if (!jet) {
      AliError(Form("%s: Could not receive jet %d", GetName(), iJets));
      continue;
    } 

    if (!AcceptJet(jet))
      continue;

    rhovec[NjetAcc] = jet->Pt() / jet->Area();
    ++NjetAcc;
  }


  if (NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);
    fRho->SetVal(rho);

    if (fRhoScaled) {
      Double_t rhoScaled = rho * GetScaleFactor(fCent);
      fRhoScaled->SetVal(rhoScaled);
    }
  }

  return kTRUE;
} 
