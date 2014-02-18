//
// Jet mass response analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
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
  fJetMassAvg(0),
  fh3PtJet1DeltaPtDeltaM(0),
  fh3PtJet2DeltaPtDeltaM(0),
  fh3PtJet1MJet1MJet2(0),
  fh3PtJet2MJet1MJet2(0)
{
  // Default constructor.

  fh3PtJet1DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet1MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet2MJet1MJet2          = new TH3F*[fNcentBins];
 
  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtJet1DeltaPtDeltaM[i]     = 0; 
    fh3PtJet2DeltaPtDeltaM[i]     = 0;
    fh3PtJet1MJet1MJet2[i]        = 0;
    fh3PtJet2MJet1MJet2[i]        = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassResponse::AliAnalysisTaskEmcalJetMassResponse(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fMinFractionShared(0),
  fJetMassAvg(0),
  fh3PtJet1DeltaPtDeltaM(0),
  fh3PtJet2DeltaPtDeltaM(0),
  fh3PtJet1MJet1MJet2(0),
  fh3PtJet2MJet1MJet2(0)
{
  // Standard constructor.

  fh3PtJet1DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet1MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet2MJet1MJet2          = new TH3F*[fNcentBins];
 
  for (Int_t i = 0; i < fNcentBins; i++) {
    fh3PtJet1DeltaPtDeltaM[i]     = 0; 
    fh3PtJet2DeltaPtDeltaM[i]     = 0; 
    fh3PtJet1MJet1MJet2[i]        = 0;
    fh3PtJet2MJet1MJet2[i]        = 0;
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

  const Int_t nBinsPt  = 250;
  const Double_t minPt = -50.;
  const Double_t maxPt = 200.;

  const Int_t nBinsM  = 150;
  const Double_t minM = -50.;
  const Double_t maxM = 100.;

  // const Int_t nBinsArea = 100;
  // const Double_t minArea = 0.;
  // const Double_t maxArea = 1.;

  // const Int_t nBinsNConst = 100;
  // const Double_t minNConst = 0.;
  // const Double_t maxNConst = 500.;

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
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

      fh3PtJet1DeltaPtDeltaM[fCentBin]->Fill(ptJet1,deltaPt,deltaM);
      fh3PtJet2DeltaPtDeltaM[fCentBin]->Fill(ptJet2,deltaPt,deltaM);

      fh3PtJet1MJet1MJet2[fCentBin]->Fill(ptJet1,massJet1,massJet2);
      fh3PtJet2MJet1MJet2[fCentBin]->Fill(ptJet2,massJet1,massJet2);
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMassResponse::GetJetMass(AliEmcalJet *jet) {

  Double_t subM = jet->M() - fJetMassAvg;
  // Double_t scale = deltaM / jet->M();

  // Double_t mt2 = jet->E()*jet->E() - jet->Pt()*jet->Pt();
  // Double_t et2 = jet->M()*jet->M() + jet->Pt()*jet->Pt();

  // Double_t pxScale = jet->Px()*scale;
  // Double_t pyScale = jet->Py()*scale;
  // Double_t pzScale = jet->Pz()*scale;
  
  // Printf("scaled jet 4-vector: %f-%f-%f-%f",pxScale,pyScale,pzScale);

  return subM;

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

