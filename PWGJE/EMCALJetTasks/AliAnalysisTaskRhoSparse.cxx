// $Id: AliAnalysisTaskRho.cxx 58408 2012-09-03 07:00:58Z loizides $
//
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fRhoName".apppend("_Scaled").
//
// Authors: R.Reed, S.Aiola, M.Connors

#include "AliAnalysisTaskRhoSparse.h"

#include <TClonesArray.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"

ClassImp(AliAnalysisTaskRhoSparse)

//________________________________________________________________________
AliAnalysisTaskRhoSparse::AliAnalysisTaskRhoSparse() : 
  AliAnalysisTaskRhoBase("AliAnalysisTaskRhoSparse"),
  fHistOccCorrvsCent(0),
  fNExclLeadJets(0),
  fRhoCMS(0)
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskRhoSparse::AliAnalysisTaskRhoSparse(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name, histo),
  fHistOccCorrvsCent(0),
  fNExclLeadJets(0),
  fRhoCMS(0)
{
  // Constructor.
}

//________________________________________________________________________
void AliAnalysisTaskRhoSparse::UserCreateOutputObjects()
{
  if (!fCreateHisto)
    return;

  AliAnalysisTaskRhoBase::UserCreateOutputObjects();
  
  fHistOccCorrvsCent = new TH2F("OccCorrvsCent", "OccCorrvsCent", 101, -1, 100, 2000, 0 , 2);
  fOutput->Add(fHistOccCorrvsCent);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoSparse::Run() 
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
  Double_t TotaljetArea=0;
  Double_t TotaljetAreaPhys=0;

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

    //cout << "jetpt: " << jet->Pt() << " jetArea: " << jet->Area() <<endl;

    if (!AcceptJet(jet))
      continue;


    if(jet->Pt()>0.01) rhovec[NjetAcc] = jet->Pt() / jet->Area();
    //cout << "ACCEPTED: jetpt: " << jet->Pt() << " jetArea: " << jet->Area() << " jetRho: " << rhovec[NjetAcc] <<endl;
    if(jet->Pt()>0.01) TotaljetAreaPhys+=jet->Area();
    TotaljetArea+=jet->Area();
    ++NjetAcc;
  }

  const Double_t TpcMaxPhi = TMath::Pi()*2.;

  const Double_t TpcArea = TpcMaxPhi * 2.*(0.7);
  Double_t OccCorr=0.0;
  //cout << "Area Physical: " << TotaljetAreaPhys << " total: " << TotaljetArea <<endl;
  //if(TotaljetArea>0) OccCorr=TotaljetAreaPhys/TotaljetArea;
  if(TpcArea>0) OccCorr=TotaljetAreaPhys/TpcArea;
 
  fHistOccCorrvsCent->Fill(fCent, OccCorr);


  if (NjetAcc > 0) {
    //find median value
    Double_t rho = TMath::Median(NjetAcc, rhovec);

    if(fRhoCMS){
      rho = rho * OccCorr;
    }

    fRho->SetVal(rho);


    if (fRhoScaled) {
      Double_t rhoScaled = rho * GetScaleFactor(fCent);
      fRhoScaled->SetVal(rhoScaled);
    }
  }

  return kTRUE;
} 
