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
#include <TTree.h>

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
#include "AliParticleContainer.h"

#include "AliAODEvent.h"

#include "AliAnalysisTaskEmcalJetMassResponse.h"

ClassImp(AliAnalysisTaskEmcalJetMassResponse)

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassResponse::AliAnalysisTaskEmcalJetMassResponse() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetMassResponse", kTRUE),
  fContainerBase(0),
  fMinFractionShared(0),
  f1JetMassAvg(0),
  fSingleTrackEmb(kFALSE),
  fSubtractMassless(kFALSE),
  fCreateTree(kFALSE),
  fh2PtJet1DeltaMNoSub(0),
  fh2PtJet2DeltaMNoSub(0),
  fh3PtJet1DeltaPtDeltaMCheat(0),
  fh3PtJet2DeltaPtDeltaMCheat(0),
  fh3PtJet1DeltaPtDeltaM(0),
  fh3PtJet2DeltaPtDeltaM(0),
  fh2PtJet1DeltaE(0),
  fh2PtJet2DeltaE(0),
  fh2PtJet1DeltaP(0),
  fh2PtJet2DeltaP(0),
  fh2PtJet2DeltaM(0),
  fh3PtJet1MJet1MJet2(0),
  fh3PtJet2MJet1MJet2(0),
  fh3PtJet1DeltaPtDeltaMRho(0),
  fh2PtJet1DeltaERho(0),
  fh3PtJet2DeltaPtDeltaMRho(0),
  fh2PtJet2DeltaPxRho(0),
  fh2PtJet2DeltaPyRho(0),
  fh2PtJet2DeltaPzRho(0),
  fh2PtJet2DeltaERho(0),
  fh2PtJet2DeltaMRho(0),
  fh2PtJet2DeltaPtRho(0),
  fh3PtJet2DeltaEDeltaMRho(0),
  fh3PtJet2DeltaPDeltaMRho(0),
  fh2PtJet1DeltaPtVecSub(0),
  fTreeJetBkg(),
  fJet1Vec(new TLorentzVector()),
  fJet2Vec(new TLorentzVector()),
  fBkgVec(new TLorentzVector()),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fJetMassMassless(0)
{
  // Default constructor.

  fh2PtJet1DeltaMNoSub         = new TH2F*[fNcentBins];
  fh2PtJet2DeltaMNoSub         = new TH2F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh2PtJet1DeltaE              = new TH2F*[fNcentBins];
  fh2PtJet2DeltaE              = new TH2F*[fNcentBins];
  fh2PtJet1DeltaP              = new TH2F*[fNcentBins];
  fh2PtJet2DeltaP              = new TH2F*[fNcentBins];
  fh2PtJet2DeltaM              = new TH2F*[fNcentBins];
  fh3PtJet1MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet2MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaMRho    = new TH3F*[fNcentBins];
  fh2PtJet1DeltaERho           = new TH2F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaMRho    = new TH3F*[fNcentBins];
  fh2PtJet2DeltaPxRho          = new TH2F*[fNcentBins];
  fh2PtJet2DeltaPyRho          = new TH2F*[fNcentBins];
  fh2PtJet2DeltaPzRho          = new TH2F*[fNcentBins];
  fh2PtJet2DeltaERho           = new TH2F*[fNcentBins];
  fh2PtJet2DeltaMRho           = new TH2F*[fNcentBins];
  fh2PtJet2DeltaPtRho          = new TH2F*[fNcentBins];
  fh3PtJet2DeltaEDeltaMRho     = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPDeltaMRho     = new TH3F*[fNcentBins];
  fh2PtJet1DeltaPtVecSub       = new TH2F*[fNcentBins];
 
  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1DeltaMNoSub[i]        = 0;
    fh2PtJet2DeltaMNoSub[i]        = 0;
    fh3PtJet1DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet2DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet1DeltaPtDeltaM[i]      = 0; 
    fh3PtJet2DeltaPtDeltaM[i]      = 0;
    fh2PtJet1DeltaE[i]             = 0;
    fh2PtJet2DeltaE[i]             = 0;
    fh2PtJet1DeltaP[i]             = 0;
    fh2PtJet2DeltaP[i]             = 0;
    fh2PtJet2DeltaM[i]             = 0;
    fh3PtJet1MJet1MJet2[i]         = 0;
    fh3PtJet2MJet1MJet2[i]         = 0;
    fh3PtJet1DeltaPtDeltaMRho[i]   = 0; 
    fh2PtJet1DeltaERho[i]          = 0;
    fh3PtJet2DeltaPtDeltaMRho[i]   = 0;
    fh2PtJet2DeltaPxRho[i]         = 0;
    fh2PtJet2DeltaPyRho[i]         = 0;
    fh2PtJet2DeltaPzRho[i]         = 0;
    fh2PtJet2DeltaERho[i]          = 0;
    fh2PtJet2DeltaMRho[i]          = 0;
    fh2PtJet2DeltaPtRho[i]         = 0;
    fh3PtJet2DeltaEDeltaMRho[i]    = 0;
    fh3PtJet2DeltaPDeltaMRho[i]    = 0;
    fh2PtJet1DeltaPtVecSub[i]      = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassResponse::AliAnalysisTaskEmcalJetMassResponse(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fMinFractionShared(0),
  f1JetMassAvg(0),
  fSingleTrackEmb(kFALSE),
  fSubtractMassless(kFALSE),
  fCreateTree(kFALSE),
  fh2PtJet1DeltaMNoSub(0),
  fh2PtJet2DeltaMNoSub(0),
  fh3PtJet1DeltaPtDeltaMCheat(0),
  fh3PtJet2DeltaPtDeltaMCheat(0),
  fh3PtJet1DeltaPtDeltaM(0),
  fh3PtJet2DeltaPtDeltaM(0),
  fh2PtJet1DeltaE(0),
  fh2PtJet2DeltaE(0),
  fh2PtJet1DeltaP(0),
  fh2PtJet2DeltaP(0),
  fh2PtJet2DeltaM(0),
  fh3PtJet1MJet1MJet2(0),
  fh3PtJet2MJet1MJet2(0),
  fh3PtJet1DeltaPtDeltaMRho(0),
  fh2PtJet1DeltaERho(0),
  fh3PtJet2DeltaPtDeltaMRho(0),
  fh2PtJet2DeltaPxRho(0),
  fh2PtJet2DeltaPyRho(0),
  fh2PtJet2DeltaPzRho(0),
  fh2PtJet2DeltaERho(0),
  fh2PtJet2DeltaMRho(0),
  fh2PtJet2DeltaPtRho(0),
  fh3PtJet2DeltaEDeltaMRho(0),
  fh3PtJet2DeltaPDeltaMRho(0),
  fh2PtJet1DeltaPtVecSub(0),
  fTreeJetBkg(0),
  fJet1Vec(new TLorentzVector()),
  fJet2Vec(new TLorentzVector()),
  fBkgVec(new TLorentzVector()),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fJetMassMassless(0)
{
  // Standard constructor.

  fh2PtJet1DeltaMNoSub         = new TH2F*[fNcentBins];
  fh2PtJet2DeltaMNoSub         = new TH2F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaMCheat  = new TH3F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaM       = new TH3F*[fNcentBins];
  fh2PtJet1DeltaE              = new TH2F*[fNcentBins];
  fh2PtJet2DeltaE              = new TH2F*[fNcentBins];
  fh2PtJet1DeltaP              = new TH2F*[fNcentBins];
  fh2PtJet2DeltaP              = new TH2F*[fNcentBins];
  fh2PtJet2DeltaM              = new TH2F*[fNcentBins];
  fh3PtJet1MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet2MJet1MJet2          = new TH3F*[fNcentBins];
  fh3PtJet1DeltaPtDeltaMRho    = new TH3F*[fNcentBins];
  fh2PtJet1DeltaERho           = new TH2F*[fNcentBins];
  fh3PtJet2DeltaPtDeltaMRho    = new TH3F*[fNcentBins];
  fh2PtJet2DeltaPxRho          = new TH2F*[fNcentBins];
  fh2PtJet2DeltaPyRho          = new TH2F*[fNcentBins];
  fh2PtJet2DeltaPzRho          = new TH2F*[fNcentBins];
  fh2PtJet2DeltaERho           = new TH2F*[fNcentBins];
  fh2PtJet2DeltaMRho           = new TH2F*[fNcentBins];
  fh2PtJet2DeltaPtRho          = new TH2F*[fNcentBins];
  fh3PtJet2DeltaEDeltaMRho     = new TH3F*[fNcentBins];
  fh3PtJet2DeltaPDeltaMRho     = new TH3F*[fNcentBins];
  fh2PtJet1DeltaPtVecSub       = new TH2F*[fNcentBins];
 
  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtJet1DeltaMNoSub[i]        = 0;
    fh2PtJet2DeltaMNoSub[i]        = 0;
    fh3PtJet1DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet2DeltaPtDeltaMCheat[i] = 0;
    fh3PtJet1DeltaPtDeltaM[i]      = 0; 
    fh3PtJet2DeltaPtDeltaM[i]      = 0;
    fh2PtJet1DeltaE[i]             = 0;
    fh2PtJet2DeltaE[i]             = 0;
    fh2PtJet1DeltaP[i]             = 0;
    fh2PtJet2DeltaP[i]             = 0;
    fh2PtJet2DeltaM[i]             = 0;
    fh3PtJet1MJet1MJet2[i]         = 0;
    fh3PtJet2MJet1MJet2[i]         = 0;
    fh3PtJet1DeltaPtDeltaMRho[i]   = 0; 
    fh2PtJet1DeltaERho[i]          = 0;
    fh3PtJet2DeltaPtDeltaMRho[i]   = 0;
    fh2PtJet2DeltaPxRho[i]         = 0;
    fh2PtJet2DeltaPyRho[i]         = 0;
    fh2PtJet2DeltaPzRho[i]         = 0;
    fh2PtJet2DeltaERho[i]          = 0;
    fh2PtJet2DeltaMRho[i]          = 0;
    fh2PtJet2DeltaPtRho[i]         = 0;
    fh3PtJet2DeltaEDeltaMRho[i]    = 0;
    fh3PtJet2DeltaPDeltaMRho[i]    = 0;
    fh2PtJet1DeltaPtVecSub[i]      = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(2, TTree::Class());
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

    histName = TString::Format("fh2PtJet1DeltaE_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{E}",histName.Data());
    fh2PtJet1DeltaE[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet1DeltaE[i]);
    
    histName = TString::Format("fh2PtJet2DeltaE_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{E}",histName.Data());
    fh2PtJet2DeltaE[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaE[i]);

    histName = TString::Format("fh2PtJet1DeltaP_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}",histName.Data());
    fh2PtJet1DeltaP[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet1DeltaP[i]);
    
    histName = TString::Format("fh2PtJet2DeltaP_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}",histName.Data());
    fh2PtJet2DeltaP[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaP[i]);

    histName = TString::Format("fh2PtJet2DeltaM_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{M}",histName.Data());
    fh2PtJet2DeltaM[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaM[i]);

    histName = TString::Format("fh3PtJet1MJet1MJet2_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#it{M}_{jet1};#it{M}_{jet2}",histName.Data());
    fh3PtJet1MJet1MJet2[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet1MJet1MJet2[i]);

    histName = TString::Format("fh3PtJet2MJet1MJet2_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet2};#it{M}_{jet1};#it{M}_{jet2}",histName.Data());
    fh3PtJet2MJet1MJet2[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet2MJet1MJet2[i]);

    histName = TString::Format("fh3PtJet1DeltaPtDeltaMRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{T};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet1DeltaPtDeltaMRho[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet1DeltaPtDeltaMRho[i]);

    histName = TString::Format("fh2PtJet1DeltaERho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{E}",histName.Data());
    fh2PtJet1DeltaERho[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet1DeltaERho[i]);

    histName = TString::Format("fh3PtJet2DeltaPtDeltaMRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet2};#delta#it{p}_{T};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet2DeltaPtDeltaMRho[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet2DeltaPtDeltaMRho[i]);

    histName = TString::Format("fh2PtJet2DeltaPxRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{x}",histName.Data());
    fh2PtJet2DeltaPxRho[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaPxRho[i]);

    histName = TString::Format("fh2PtJet2DeltaPyRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{y}",histName.Data());
    fh2PtJet2DeltaPyRho[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaPyRho[i]);

    histName = TString::Format("fh2PtJet2DeltaPzRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{z}",histName.Data());
    fh2PtJet2DeltaPzRho[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaPzRho[i]);

    histName = TString::Format("fh2PtJet2DeltaERho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{E}",histName.Data());
    fh2PtJet2DeltaERho[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaERho[i]);

    histName = TString::Format("fh2PtJet2DeltaMRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{M}",histName.Data());
    fh2PtJet2DeltaMRho[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaMRho[i]);

    histName = TString::Format("fh2PtJet2DeltaPtRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet1};#delta#it{p}_{T}",histName.Data());
    fh2PtJet2DeltaPtRho[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2PtJet2DeltaPtRho[i]);

    histName = TString::Format("fh3PtJet2DeltaEDeltaMRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet2};#delta#it{E};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet2DeltaEDeltaMRho[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet2DeltaEDeltaMRho[i]);

    histName = TString::Format("fh3PtJet2DeltaPDeltaMRho_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,jet2};#delta#it{p};#delta#it{M}_{jet}",histName.Data());
    fh3PtJet2DeltaPDeltaMRho[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh3PtJet2DeltaPDeltaMRho[i]);

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

  // Create a tree.
  if(fCreateTree) {
    fTreeJetBkg = new TTree("fTreeJetBkg", "fTreeJetBkg");
    fTreeJetBkg->Branch("fJet1Vec","TLorentzVector",&fJet1Vec);
    fTreeJetBkg->Branch("fJet2Vec","TLorentzVector",&fJet2Vec);
    fTreeJetBkg->Branch("fBkgVec","TLorentzVector",&fBkgVec);
    fTreeJetBkg->Branch("fArea",&fArea,"fArea/F");
    fTreeJetBkg->Branch("fAreaPhi",&fAreaPhi,"fAreaPhi/F");
    fTreeJetBkg->Branch("fAreaEta",&fAreaEta,"fAreaEta/F");
    fTreeJetBkg->Branch("fRho",&fRho,"fRho/F");
    fTreeJetBkg->Branch("fRhoM",&fRhoM,"fRhoM/F");
    fTreeJetBkg->Branch("fNConst",&fNConst,"fNConst/I");
    fTreeJetBkg->Branch("fJetMassMassless",&fJetMassMassless,"fJetMassMassless/F");
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  if(fCreateTree) PostData(2, fTreeJetBkg);
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

  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      if(fSingleTrackEmb) {
	AliVParticle *vp = GetEmbeddedConstituent(jet1);
	if(!vp) continue;
	fJet2Vec->SetPxPyPzE(vp->Px(),vp->Py(),vp->Pz(),vp->E());
      } else {
	AliEmcalJet *jet2 = jet1->ClosestJet();
	if(!jet2) continue;

	Double_t fraction = jetCont->GetFractionSharedPt(jet1);
	if(fMinFractionShared>0. && fraction<fMinFractionShared) continue;
	fJet2Vec->SetPxPyPzE(jet2->Px(),jet2->Py(),jet2->Pz(),jet2->E());
      }

      Double_t ptJet1 = jet1->Pt() - GetRhoVal(fContainerBase)*jet1->Area();
      Double_t massJet1 = GetJetMass(jet1);//jet1->M();

      Double_t deltaPt = ptJet1 - fJet2Vec->Pt();
      Double_t deltaM  = massJet1 - fJet2Vec->M();
      TLorentzVector vpS = GetSubtractedVector(jet1);
      Double_t deltaE  = vpS.E() - fJet2Vec->E();
      Double_t deltaP  = vpS.P() - fJet2Vec->P();

      fh2PtJet1DeltaMNoSub[fCentBin]->Fill(ptJet1,jet1->M()-fJet2Vec->M());
      fh2PtJet2DeltaMNoSub[fCentBin]->Fill(fJet2Vec->Pt(),jet1->M()-fJet2Vec->M());

      fh3PtJet1DeltaPtDeltaM[fCentBin]->Fill(ptJet1,deltaPt,deltaM);
      fh3PtJet2DeltaPtDeltaM[fCentBin]->Fill(fJet2Vec->Pt(),deltaPt,deltaM);

      fh2PtJet1DeltaE[fCentBin]->Fill(ptJet1,deltaE);
      //  fh2PtJet2DeltaE[fCentBin]->Fill(ptJet2,deltaE);

      fh2PtJet1DeltaP[fCentBin]->Fill(ptJet1,deltaP);
      //      fh2PtJet2DeltaP[fCentBin]->Fill(ptJet2,deltaP);

      fh3PtJet1MJet1MJet2[fCentBin]->Fill(ptJet1,massJet1,fJet2Vec->M());
      fh3PtJet2MJet1MJet2[fCentBin]->Fill(fJet2Vec->Pt(),massJet1,fJet2Vec->M());

      TLorentzVector vpBkgCheat = GetBkgVectorCheat(jet1);
      //      TLorentzVector vpBkgCheat; vpBkgCheat.SetPtEtaPhiM(vpBkgCheatB.Perp(),jet1->Eta(),jet1->Phi(),vpBkgCheatB.M());
      TLorentzVector vpBkg = GetBkgVector(jet1,jetCont);
      fh2PtJet2DeltaE[fCentBin]->Fill(fJet2Vec->Pt(),vpBkg.E()-vpBkgCheat.E());      
      fh2PtJet2DeltaP[fCentBin]->Fill(fJet2Vec->Pt(),vpBkg.P()-vpBkgCheat.P());
      fh2PtJet2DeltaM[fCentBin]->Fill(fJet2Vec->Pt(),vpBkg.M()-vpBkgCheat.M());

      // fh3PtJet1DeltaPtDeltaMRho[fCentBin]->Fill(jet1->Pt()-vpBkg.Perp()*jet1->Area(),jet1->Pt()-vpBkg.Perp()*jet1->Area()-ptJet2,jet1->M()-vpBkg.M()*jet1->Area()-jet2->M());
      // fh3PtJet2DeltaPtDeltaMRho[fCentBin]->Fill(ptJet2,jet1->Pt()-vpBkg.Perp()*jet1->Area()-ptJet2,jet1->M()-vpBkg.M()*jet1->Area()-jet2->M());
      Double_t px = jet1->Px()-vpBkg.Px();
      Double_t py = jet1->Py()-vpBkg.Py();
      Double_t pz = jet1->Pz()-vpBkg.Pz();
      Double_t E  = jet1->E()-vpBkg.E();
      Double_t p2 = px*px + py*py + pz*pz;
      Double_t msub = 0.;
      if((E*E)>p2) msub = TMath::Sqrt(E*E - p2);
      // fh3PtJet1DeltaPtDeltaMRho[fCentBin]->Fill(jet1->Pt()-vpBkg.Perp(),jet1->Pt()-vpBkg.Perp()-ptJet2,jet1->M()-vpBkg.M()-jet2->M());
      // fh3PtJet2DeltaPtDeltaMRho[fCentBin]->Fill(ptJet2,jet1->Pt()-vpBkg.Perp()-ptJet2,jet1->M()-vpBkg.M()-jet2->M());

      fh3PtJet1DeltaPtDeltaMRho[fCentBin]->Fill(jet1->Pt()-vpBkg.Perp(),jet1->Pt()-vpBkg.Perp()-fJet2Vec->Pt(),msub-fJet2Vec->M());
      fh2PtJet1DeltaERho[fCentBin]->Fill(ptJet1,E - fJet2Vec->E());

      fh3PtJet2DeltaPtDeltaMRho[fCentBin]->Fill(fJet2Vec->Pt(),jet1->Pt()-vpBkg.Perp()-fJet2Vec->Pt(),msub-fJet2Vec->M());
      fh2PtJet2DeltaPxRho[fCentBin]->Fill(fJet2Vec->Pt(),px - fJet2Vec->Px());
      fh2PtJet2DeltaPyRho[fCentBin]->Fill(fJet2Vec->Pt(),py - fJet2Vec->Py());
      fh2PtJet2DeltaPzRho[fCentBin]->Fill(fJet2Vec->Pt(),pz - fJet2Vec->Pz());
      fh2PtJet2DeltaERho[fCentBin]->Fill(fJet2Vec->Pt(),E - fJet2Vec->E());
      fh3PtJet2DeltaEDeltaMRho[fCentBin]->Fill(fJet2Vec->Pt(),E-fJet2Vec->E(),msub-fJet2Vec->M());
      fh3PtJet2DeltaPDeltaMRho[fCentBin]->Fill(fJet2Vec->Pt(),TMath::Sqrt(p2)-fJet2Vec->P(),msub-fJet2Vec->M());

      fh2PtJet2DeltaMRho[fCentBin]->Fill(fJet2Vec->Pt(),msub - fJet2Vec->M());
      fh2PtJet2DeltaPtRho[fCentBin]->Fill(fJet2Vec->Pt(),TMath::Sqrt(px*px+py*py) - fJet2Vec->Pt());

      TLorentzVector vpC = GetSubtractedVectorCheat(jet1);
      fh3PtJet1DeltaPtDeltaMCheat[fCentBin]->Fill(vpC.Perp(),vpC.Perp()-fJet2Vec->Pt(),vpC.M()-fJet2Vec->M());
      fh3PtJet2DeltaPtDeltaMCheat[fCentBin]->Fill(fJet2Vec->Pt(),vpC.Perp()-fJet2Vec->Pt(),vpC.M()-fJet2Vec->M());

      if(fJet2Vec->Pt()>20. && fCreateTree) {      
	fBkgVec->SetPxPyPzE(vpBkg.Px(),vpBkg.Py(),vpBkg.Pz(),vpBkg.E());
	fJet1Vec->SetPxPyPzE(px,py,pz,E);                                     //AA jet
	fArea = (Float_t)jet1->Area();
	fAreaPhi = (Float_t)jet1->AreaPhi();
	fAreaEta = (Float_t)jet1->AreaEta();
	fRho  = (Float_t)jetCont->GetRhoVal();
	fRhoM = (Float_t)jetCont->GetRhoMassVal();
	fNConst = (Int_t)jet1->GetNumberOfTracks();
	fJetMassMassless = (Float_t)GetJetMassMasslessConstituents(jet1);
	fTreeJetBkg->Fill();
      }
    }
  }
  return kTRUE;
}

//________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassResponse::GetBkgVector(AliEmcalJet *jet, AliJetContainer *cont) {
  //get background vector

  Double_t rho  = cont->GetRhoVal();
  Double_t rhom = cont->GetRhoMassVal();
  TLorentzVector vpB;
  Double_t aphi = jet->AreaPhi();
  Double_t aeta = jet->AreaEta();
  //  vpB.SetPxPyPzE(rho*TMath::Cos(aphi)*jet->Area(),rho*TMath::Sin(aphi)*jet->Area(),(rho+rhom)*TMath::SinH(aeta)*jet->Area(),(rho+rhom)*TMath::CosH(aeta)*jet->Area());
  vpB.SetPxPyPzE(rho*TMath::Cos(aphi)*jet->Area(),rho*TMath::Sin(aphi)*jet->Area(),(rho+rhom)*TMath::SinH(aeta)*jet->Area(),(rho+rhom)*TMath::CosH(aeta)*jet->Area());
  return vpB;
}

//________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassResponse::GetSubtractedVector(AliEmcalJet *jet) {
  //get subtracted vector
  TLorentzVector vpS;

  // if(f1JetMassAvg) {
  //   Double_t pt = jet->Pt() - GetRhoVal(fContainerBase)*jet->Area();
  //   TLorentzVector vpB; vpB.SetPtEtaPhiE(GetRhoVal(fContainerBase)*jet->Area(),0.,0.,f1JetMassAvg->Eval(pt));
  //   TLorentzVector vpAAboost; vpAAboost.SetPtEtaPhiM(jet->Pt(),0.,0.,jet->M());
  //   TLorentzVector vpSboost = vpAAboost - vpB;
  //   vpS.SetPtEtaPhiM(vpSboost.Perp(),jet->Eta(),jet->Phi(),vpSboost.M());
  // } else {
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  TLorentzVector vpBkg = GetBkgVector(jet,jetCont);
  vpS.SetPxPyPzE(jet->Px()-vpBkg.Px(),jet->Py()-vpBkg.Py(),jet->Pz()-vpBkg.Pz(),jet->E()-vpBkg.E());
  // TLorentzVector vpAAboost; vpAAboost.SetPtEtaPhiE(jet->Pt(),0.,0.,jet->E());
  // TLorentzVector vpBkgboost; vpBkgboost.SetPtEtaPhiE(vpBkg.Perp(),0.,0.,vpBkg.E());
  // TLorentzVector vpSboost = vpAAboost - vpBkgboost;
  //  vpS.SetPtEtaPhiM(vpSboost.Perp(),jet->Eta(),jet->Phi(),vpSboost.M());
  // }
  return vpS;
}

//________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassResponse::GetBkgVectorCheat(AliEmcalJet *jet) {
  //get background vector with cheating
  TLorentzVector vpB; 
  if(fJet2Vec) {
    TLorentzVector vpAAboost; vpAAboost.SetPtEtaPhiM(jet->Pt(),0.,0.,jet->M());
    TLorentzVector vpPPboost; vpPPboost.SetPtEtaPhiM(fJet2Vec->Pt(),0.,0.,fJet2Vec->M());
    Double_t dpt = vpAAboost.Perp()-vpPPboost.Perp();
    Double_t dE = vpAAboost.E()-vpPPboost.E();
    vpB.SetPtEtaPhiE(dpt,0.,0.,dE);
  }
  return vpB;
}

//________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassResponse::GetSubtractedVectorCheat(AliEmcalJet *jet) {
  //get subtracted vector taking pT and energy difference from MC match
  TLorentzVector vpS;
  TLorentzVector vpB = GetBkgVectorCheat(jet);
  TLorentzVector vpAAboost; vpAAboost.SetPtEtaPhiM(jet->Pt(),0.,0.,jet->M());
  TLorentzVector vpSboost = vpAAboost - vpB;
  vpS.SetPtEtaPhiM(vpSboost.Perp(),jet->Eta(),jet->Phi(),vpSboost.M());
  return vpS;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMassResponse::GetJetMass(AliEmcalJet *jet) {

  TLorentzVector vpS = GetSubtractedVector(jet);
  
  AliEmcalJet *jet2 = jet->ClosestJet();
  if(jet2) fh2PtJet1DeltaPtVecSub[fCentBin]->Fill(vpS.Perp(),vpS.Perp()-jet2->Pt());
  
  if(fSubtractMassless) {
    Double_t m = vpS.M()-GetJetMassMasslessConstituents(jet);
    return m;
  } else
    return vpS.M();
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetMassResponse::GetJetMassMasslessConstituents(AliEmcalJet *jet) {
  //Get mass of jet assuming all particles are massless (E==P)
  
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  AliVParticle *vp = 0x0;
  Double_t px = 0.; Double_t py = 0.; Double_t pz = 0.; Double_t E = 0.;
  for(Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
    vp = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
    px+=vp->Px();
    py+=vp->Py();
    pz+=vp->Pz();
    E+=vp->P();
  }

  Double_t m2 = E*E - px*px - py*py - pz*pz;
  if(m2>0.)
    return TMath::Sqrt(E*E - px*px - py*py - pz*pz);
  else
    return 0.;
}

//________________________________________________________________________
AliVParticle* AliAnalysisTaskEmcalJetMassResponse::GetEmbeddedConstituent(AliEmcalJet *jet) {

  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  AliVParticle *vp = 0x0;
  AliVParticle *vpe = 0x0; //embedded particle
  Int_t nc = 0;
  for(Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
    vp = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray())); //check if fTracks is the correct track branch
    if (vp->TestBits(TObject::kBitMask) != (Int_t)(TObject::kBitMask) ) continue;
    if(!vpe) vpe = vp;
    else if(vp->Pt()>vpe->Pt()) vpe = vp;
    nc++;
  }
  AliDebug(11,Form("Found %d embedded particles",nc));
  return vpe;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassResponse::RetrieveEventObjects() {
  //
  // retrieve event objects
  //

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  jetCont->LoadRhoMass(InputEvent());

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetMassResponse::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

