//
// Jet mass response analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskEmcalJetMassResponse.h"

ClassImp(AliAnalysisTaskEmcalJetMassResponse)

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassResponse::AliAnalysisTaskEmcalJetMassResponse() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetMassResponse", kTRUE),
  fContainerBase(0),
  fMinFractionShared(0),
  f1JetMassAvg(0),
  fh2PtJet1DeltaMNoSub(0),
  fh2PtJet2DeltaMNoSub(0),
  fh3PtJet1DeltaPtDeltaMCheat(0),
  fh3PtJet2DeltaPtDeltaMCheat(0),
  fh3PtJet1DeltaPtDeltaM(0),
  fh3PtJet2DeltaPtDeltaM(0),
  fh3PtJet1MJet1MJet2(0),
  fh3PtJet2MJet1MJet2(0),
  fh2PtJet1DeltaPtVecSub(0)
{
  // Default constructor.

  fh2PtJet1DeltaMNoSub         = new TH2F*[fNcentBins];
  fh2PtJet2DeltaMNoSub         = new TH2F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet1MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet2MJet1MJet2          = new TH3F*[fNcentBins];
  fh2PtJet1DeltaPtVecSub       = new TH2F*[fNcentBins];
 
  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1DeltaMNoSub[i]        = 0;
    fh2PtJet2DeltaMNoSub[i]        = 0;
    fh3PtJet1DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet2DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet1DeltaPtDeltaM[i]      = 0; 
    fh3PtJet2DeltaPtDeltaM[i]      = 0;
    fh3PtJet1MJet1MJet2[i]         = 0;
    fh3PtJet2MJet1MJet2[i]         = 0;
    fh2PtJet1DeltaPtVecSub[i]      = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassResponse::AliAnalysisTaskEmcalJetMassResponse(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fMinFractionShared(0),
  f1JetMassAvg(0),
  fh2PtJet1DeltaMNoSub(0),
  fh2PtJet2DeltaMNoSub(0),
  fh3PtJet1DeltaPtDeltaMCheat(0),
  fh3PtJet2DeltaPtDeltaMCheat(0),
  fh3PtJet1DeltaPtDeltaM(0),
  fh3PtJet2DeltaPtDeltaM(0),
  fh3PtJet1MJet1MJet2(0),
  fh3PtJet2MJet1MJet2(0),
  fh2PtJet1DeltaPtVecSub(0)
{
  // Standard constructor.

  fh2PtJet1DeltaMNoSub         = new TH2F*[fNcentBins];
  fh2PtJet2DeltaMNoSub         = new TH2F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet1MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet2MJet1MJet2          = new TH3F*[fNcentBins];
  fh2PtJet1DeltaPtVecSub       = new TH2F*[fNcentBins];
 
  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1DeltaMNoSub[i]        = 0;
    fh2PtJet2DeltaMNoSub[i]        = 0;
    fh3PtJet1DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet2DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet1DeltaPtDeltaM[i]      = 0; 
    fh3PtJet2DeltaPtDeltaM[i]      = 0;
    fh3PtJet1MJet1MJet2[i]         = 0;
    fh3PtJet2MJet1MJet2[i]         = 0;
    fh2PtJet1DeltaPtVecSub[i]      = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassResponse::~AliAnalysisTaskEmcalJetMassResponse()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassResponse::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsPt  = 200;
  const Double_t minPt = -50.;
  const Double_t maxPt = 150.;

  const Int_t nBinsM  = 150;
  const Double_t minM = -50.;
  const Double_t maxM = 100.;

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = TString::Format("fh2PtJet1DeltaMNoSub_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{M}_{jet}",histName.Data());
    fh2PtJet1DeltaMNoSub[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtJet1DeltaMNoSub[i]);

    histName = TString::Format("fh2PtJet2DeltaMNoSub_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{M}_{jet}",histName.Data());
    fh2PtJet2DeltaMNoSub[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtJet2DeltaMNoSub[i]);

    histName = TString::Format("fh3PtJet1DeltaPtDeltaMCheat_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{T};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet1DeltaPtDeltaMCheat[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet1DeltaPtDeltaMCheat[i]);

    histName = TString::Format("fh3PtJet2DeltaPtDeltaMCheat_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{T};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet2DeltaPtDeltaMCheat[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet2DeltaPtDeltaMCheat[i]);

    histName = TString::Format("fh3PtJet1DeltaPtDeltaM_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{T};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet1DeltaPtDeltaM[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet1DeltaPtDeltaM[i]);

    histName = TString::Format("fh3PtJet2DeltaPtDeltaM_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet2};#delta#it{p}_{T};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet2DeltaPtDeltaM[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet2DeltaPtDeltaM[i]);

    histName = TString::Format("fh3PtJet1MJet1MJet2_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{M}_{jet2}",histName.Data());
    fh3PtJet1MJet1MJet2[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet1MJet1MJet2[i]);

    histName = TString::Format("fh3PtJet2MJet1MJet2_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet2};#it{M}_{jet1};#it{M}_{jet2}",histName.Data());
    fh3PtJet2MJet1MJet2[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet2MJet1MJet2[i]);

    histName = TString::Format("fh2PtJet1DeltaPtVecSub_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{T}",histName.Data());
    fh2PtJet1DeltaPtVecSub[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet1DeltaPtVecSub[i]);
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassResponse::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassResponse::FillHistograms()
{
  // Fill histograms.

  AliInfo(Form("%s",GetName()));

  AliEmcalJet* jet1 = NULL;

  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      AliEmcalJet *jet2 = jet1->ClosestJet();
      if(!jet2) return -1;

      Double_t fraction = jetCont->GetFractionSharedPt(jet1);
      if(fMinFractionShared>0. && fraction<fMinFractionShared) continue;

      Double_t ptJet1 = jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();
      Double_t ptJet2 = jet2->Pt();
      Double_t massJet1 = GetJetMass(jet1);//jet1->M();
      Double_t massJet2 = jet2->M();

      Double_t deltaPt = ptJet1 - ptJet2;
      Double_t deltaM  = massJet1 - massJet2;

      fh2PtJet1DeltaMNoSub[fCentBin]->Fill(ptJet1,jet1->M()-jet2->M());
      fh2PtJet2DeltaMNoSub[fCentBin]->Fill(ptJet2,jet1->M()-jet2->M());

      fh3PtJet1DeltaPtDeltaM[fCentBin]->Fill(ptJet1,deltaPt,deltaM);
      fh3PtJet2DeltaPtDeltaM[fCentBin]->Fill(ptJet2,deltaPt,deltaM);

      fh3PtJet1MJet1MJet2[fCentBin]->Fill(ptJet1,massJet1,massJet2);
      fh3PtJet2MJet1MJet2[fCentBin]->Fill(ptJet2,massJet1,massJet2);

      TLorentzVector vpC = GetSubtractedVectorCheat(jet1);
      fh3PtJet1DeltaPtDeltaMCheat[fCentBin]->Fill(vpC.Perp(),vpC.Perp()-jet2->Pt(),vpC.M()-jet2->M());
      fh3PtJet2DeltaPtDeltaMCheat[fCentBin]->Fill(ptJet2,vpC.Perp()-jet2->Pt(),vpC.M()-jet2->M());
    }
  }

  return kTRUE;
}

//________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassResponse::GetSubtractedVector(AliEmcalJet *jet) {
  //get subtracted vector
  TLorentzVector vpS;
  if(f1JetMassAvg) {
    Double_t pt = jet->Pt() - GetRhoVal(fContainerBase)*jet->Area();
    TLorentzVector vpB; vpB.SetPtEtaPhiE(GetRhoVal(fContainerBase)*jet->Area(),0.,0.,f1JetMassAvg->Eval(pt));
    TLorentzVector vpAAboost; vpAAboost.SetPtEtaPhiM(jet->Pt(),0.,0.,jet->M());
    TLorentzVector vpSboost = vpAAboost - vpB;
    vpS.SetPtEtaPhiM(vpSboost.Perp(),jet->Eta(),jet->Phi(),vpSboost.M());
  }
  return vpS;
}

//________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassResponse::GetSubtractedVectorCheat(AliEmcalJet *jet) {
  //get subtracted vector taking pT and mass difference from MC match
  TLorentzVector vpS;
  AliEmcalJet *jet2 = jet->ClosestJet();
  if(jet2) {
    TLorentzVector vpAAboost; vpAAboost.SetPtEtaPhiM(jet->Pt(),0.,0.,jet->M());
    TLorentzVector vpPPboost; vpPPboost.SetPtEtaPhiM(jet2->Pt(),0.,0.,jet2->M());
    Double_t dpt = vpAAboost.Perp()-vpPPboost.Perp();
    /*
      Double_t dm = vpAAboost.M() - vpPPboost.M();
      Double_t dE = TMath::Sqrt(dpt*dpt + dm*dm);
    */
    Double_t dE = vpAAboost.E()-vpPPboost.E();
    TLorentzVector vpB; vpB.SetPtEtaPhiE(dpt,0.,0.,dE);
    TLorentzVector vpSboost = vpAAboost - vpB;
    vpS.SetPtEtaPhiM(vpSboost.Perp(),jet->Eta(),jet->Phi(),vpSboost.M());
  }
  return vpS;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMassResponse::GetJetMass(AliEmcalJet *jet) {

  if(f1JetMassAvg) {
    TLorentzVector vpS = GetSubtractedVector(jet);
  
    AliEmcalJet *jet2 = jet->ClosestJet();
    if(jet2) fh2PtJet1DeltaPtVecSub[fCentBin]->Fill(vpS.Perp(),vpS.Perp()-jet2->Pt());

    return vpS.M();
  }
  else
    return jet->M();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassResponse::RetrieveEventObjects() {
  //
  // retrieve event objects
  //

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;

}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetMassResponse::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

