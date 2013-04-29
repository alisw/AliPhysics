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
  fRhoCMS(0),
  fSigJetsName("SJets")
{
  // Constructor.
}

//________________________________________________________________________
AliAnalysisTaskRhoSparse::AliAnalysisTaskRhoSparse(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name, histo),
  fHistOccCorrvsCent(0),
  fNExclLeadJets(0),
  fRhoCMS(0),
  fSigJetsName("SJets")
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
Bool_t AliAnalysisTaskRhoSparse::IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2)
{
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i)
  {
    Int_t jet1Track = jet1->TrackAt(i);
    for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j)
    {
      Int_t jet2Track = jet2->TrackAt(j);
      if (jet1Track == jet2Track)
        return kTRUE;
    }
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoSparse::IsJetSignal(AliEmcalJet* jet)
{
  if(jet->Pt()>5){
      return kTRUE;
  }else{
    return kFALSE;
  }
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

  TClonesArray *sigjets = 0;
  sigjets= dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fSigJetsName));

  Int_t NjetsSig = 0;
  if(sigjets) NjetsSig = sigjets->GetEntries();
 

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

    TotaljetArea+=jet->Area();
    
    if(jet->Pt()>0.1){
      TotaljetAreaPhys+=jet->Area();
    }

    if (!AcceptJet(jet))
      continue;

 

   // Search for overlap with signal jets
    Bool_t isOverlapping = kFALSE;
    if(sigjets){
      for(Int_t j=0;j<NjetsSig;j++)
	{
	  AliEmcalJet* signalJet = static_cast<AliEmcalJet*>(sigjets->At(j));
	  if(!AcceptJet(signalJet))
	    continue;
	  if(!IsJetSignal(signalJet))     
	    continue;
	  
	  if(IsJetOverlapping(signalJet, jet))
	    {
	      isOverlapping = kTRUE;
	      break;
	    }
	}
    }

    if(isOverlapping) 
      continue;

    if(jet->Pt()>0.1){
      rhovec[NjetAcc] = jet->Pt() / jet->Area();
      ++NjetAcc;
    }
    

  }

  Double_t OccCorr=0.0;
  if(TotaljetArea>0) OccCorr=TotaljetAreaPhys/TotaljetArea;
 
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
