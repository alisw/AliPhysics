//
// Analysis task for angular jet shape G(R) arXiv:1201.2688
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
#include "AliAODEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"

#include "AliAnalysisTaskJetShapeGR.h"

ClassImp(AliAnalysisTaskJetShapeGR)

//________________________________________________________________________
AliAnalysisTaskJetShapeGR::AliAnalysisTaskJetShapeGR() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetShapeGR", kTRUE),
  fContainerBase(0),
  fContainerSub(1),
  fContainerTrue(2),
  fMinFractionShared(0),
  fSingleTrackEmb(kFALSE),
  fCreateTree(kFALSE),
  fTreeJetBkg(),
  fJet1Vec(new TLorentzVector()),
  fJet2Vec(new TLorentzVector()),
  fJetSubVec(new TLorentzVector()),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fMatch(0),
  fDRStep(0.04),
  fMaxR(2.),
  fh2PtTrueDeltaGR(0x0),
  fh2PtTrueDeltaGRRel(0x0),
  fhnGRResponse(0x0),
  fh1PtTrue(0x0),
  fh3DeltaGRNumRPtTrue(0x0),
  fh3DeltaGRDenRPtTrue(0x0),
  fh2DeltaGRNumRPtTrue(0x0),
  fh2DeltaGRDenRPtTrue(0x0),
  fh1PtRaw(0x0),
  fh3DeltaGRNumRPtRaw(0x0),
  fh3DeltaGRDenRPtRaw(0x0),
  fh2DeltaGRNumRPtRaw(0x0),
  fh2DeltaGRDenRPtRaw(0x0),
  fh1PtRawMatch(0x0),
  fh3DeltaGRNumRPtRawMatch(0x0),
  fh3DeltaGRDenRPtRawMatch(0x0),
  fh2DeltaGRNumRPtRawMatch(0x0),
  fh2DeltaGRDenRPtRawMatch(0x0),
  fh1PtMatch(0x0),
  fh3DeltaGRNumRPtMatch(0x0),
  fh3DeltaGRDenRPtMatch(0x0),
  fh2DeltaGRNumRPtMatch(0x0),
  fh2DeltaGRDenRPtMatch(0x0),
  fh2DeltaGRNumRPtTrueMatch(0x0),
  fh2DeltaGRDenRPtTrueMatch(0x0)
{
  // Default constructor.

  fh2PtTrueDeltaGR     = new TH2F*[fNcentBins];
  fh2PtTrueDeltaGRRel  = new TH2F*[fNcentBins];
  fhnGRResponse        = new THnSparse*[fNcentBins];
  fh1PtTrue            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtTrue = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtTrue = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtTrue = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtTrue = new TH2F*[fNcentBins];
  fh1PtRaw            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtRaw = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtRaw = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtRaw = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtRaw = new TH2F*[fNcentBins];
  fh1PtRawMatch            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtRawMatch = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtRawMatch = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtRawMatch = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtRawMatch = new TH2F*[fNcentBins];
  fh1PtMatch            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtMatch = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtMatch = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtMatch = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtMatch = new TH2F*[fNcentBins];
  fh2DeltaGRNumRPtTrueMatch = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtTrueMatch = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtTrueDeltaGR[i]     = 0;
    fh2PtTrueDeltaGRRel[i]  = 0;
    fhnGRResponse[i]        = 0;
    fh1PtTrue[i]            = 0;
    fh3DeltaGRNumRPtTrue[i] = 0;
    fh3DeltaGRDenRPtTrue[i] = 0;
    fh2DeltaGRNumRPtTrue[i] = 0;
    fh2DeltaGRDenRPtTrue[i] = 0;
    fh1PtRaw[i]            = 0;
    fh3DeltaGRNumRPtRaw[i] = 0;
    fh3DeltaGRDenRPtRaw[i] = 0;
    fh2DeltaGRNumRPtRaw[i] = 0;
    fh2DeltaGRDenRPtRaw[i] = 0;
    fh1PtRawMatch[i]            = 0;
    fh3DeltaGRNumRPtRawMatch[i] = 0;
    fh3DeltaGRDenRPtRawMatch[i] = 0;
    fh2DeltaGRNumRPtRawMatch[i] = 0;
    fh2DeltaGRDenRPtRawMatch[i] = 0;
    fh1PtMatch[i]            = 0;
    fh3DeltaGRNumRPtMatch[i] = 0;
    fh3DeltaGRDenRPtMatch[i] = 0;
    fh2DeltaGRNumRPtMatch[i] = 0;
    fh2DeltaGRDenRPtMatch[i] = 0;
    fh2DeltaGRNumRPtTrueMatch[i] = 0;
    fh2DeltaGRDenRPtTrueMatch[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeGR::AliAnalysisTaskJetShapeGR(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fContainerSub(1),
  fContainerTrue(2),
  fMinFractionShared(0),
  fSingleTrackEmb(kFALSE),
  fCreateTree(kFALSE),
  fTreeJetBkg(0),
  fJet1Vec(new TLorentzVector()),
  fJet2Vec(new TLorentzVector()),
  fJetSubVec(new TLorentzVector()),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fMatch(0),
  fDRStep(0.04),
  fMaxR(2.),
  fh2PtTrueDeltaGR(0x0),
  fh2PtTrueDeltaGRRel(0x0),
  fhnGRResponse(0x0),
  fh1PtTrue(0x0),
  fh3DeltaGRNumRPtTrue(0x0),
  fh3DeltaGRDenRPtTrue(0x0),
  fh2DeltaGRNumRPtTrue(0x0),
  fh2DeltaGRDenRPtTrue(0x0),
  fh1PtRaw(0x0),
  fh3DeltaGRNumRPtRaw(0x0),
  fh3DeltaGRDenRPtRaw(0x0),
  fh2DeltaGRNumRPtRaw(0x0),
  fh2DeltaGRDenRPtRaw(0x0),
  fh1PtRawMatch(0x0),
  fh3DeltaGRNumRPtRawMatch(0x0),
  fh3DeltaGRDenRPtRawMatch(0x0),
  fh2DeltaGRNumRPtRawMatch(0x0),
  fh2DeltaGRDenRPtRawMatch(0x0),
  fh1PtMatch(0x0),
  fh3DeltaGRNumRPtMatch(0x0),
  fh3DeltaGRDenRPtMatch(0x0),
  fh2DeltaGRNumRPtMatch(0x0),
  fh2DeltaGRDenRPtMatch(0x0),
  fh2DeltaGRNumRPtTrueMatch(0x0),
  fh2DeltaGRDenRPtTrueMatch(0x0)
{
  // Standard constructor.

  fh2PtTrueDeltaGR     = new TH2F*[fNcentBins];
  fh2PtTrueDeltaGRRel  = new TH2F*[fNcentBins];
  fhnGRResponse        = new THnSparse*[fNcentBins];
  fh1PtTrue            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtTrue = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtTrue = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtTrue = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtTrue = new TH2F*[fNcentBins];
  fh1PtRaw            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtRaw = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtRaw = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtRaw = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtRaw = new TH2F*[fNcentBins];
  fh1PtRawMatch            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtRawMatch = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtRawMatch = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtRawMatch = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtRawMatch = new TH2F*[fNcentBins];
  fh1PtMatch            = new TH1F*[fNcentBins];
  fh3DeltaGRNumRPtMatch = new TH3F*[fNcentBins];
  fh3DeltaGRDenRPtMatch = new TH3F*[fNcentBins];
  fh2DeltaGRNumRPtMatch = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtMatch = new TH2F*[fNcentBins];
  fh2DeltaGRNumRPtTrueMatch = new TH2F*[fNcentBins];
  fh2DeltaGRDenRPtTrueMatch = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtTrueDeltaGR[i]     = 0;
    fh2PtTrueDeltaGRRel[i]  = 0;
    fhnGRResponse[i]        = 0;
    fh1PtTrue[i]            = 0;
    fh3DeltaGRNumRPtTrue[i] = 0;
    fh3DeltaGRDenRPtTrue[i] = 0;
    fh2DeltaGRNumRPtTrue[i] = 0;
    fh2DeltaGRDenRPtTrue[i] = 0;
    fh1PtRaw[i]            = 0;
    fh3DeltaGRNumRPtRaw[i] = 0;
    fh3DeltaGRDenRPtRaw[i] = 0;
    fh2DeltaGRNumRPtRaw[i] = 0;
    fh2DeltaGRDenRPtRaw[i] = 0;
    fh1PtRawMatch[i]            = 0;
    fh3DeltaGRNumRPtRawMatch[i] = 0;
    fh3DeltaGRDenRPtRawMatch[i] = 0;
    fh2DeltaGRNumRPtRawMatch[i] = 0;
    fh2DeltaGRDenRPtRawMatch[i] = 0;
    fh1PtMatch[i]            = 0;
    fh3DeltaGRNumRPtMatch[i] = 0;
    fh3DeltaGRDenRPtMatch[i] = 0;
    fh2DeltaGRNumRPtMatch[i] = 0;
    fh2DeltaGRDenRPtMatch[i] = 0;
    fh2DeltaGRNumRPtTrueMatch[i] = 0;
    fh2DeltaGRDenRPtTrueMatch[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeGR::~AliAnalysisTaskJetShapeGR()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetShapeGR::UserCreateOutputObjects()
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

  const Int_t nBinsR  = 50;
  const Double_t minR = 0.;
  const Double_t maxR = 2.;

  const Int_t nBinsDGR  = 100;
  const Double_t minDGR = 0.;
  const Double_t maxDGR = 10.;

  //Binning for THnSparse
  const Int_t nBinsSparse0 = 4;
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt};

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = Form("fh2PtTrueDeltaGR_%d",i);
    histTitle = Form("fh2PtTrueDeltaGR_%d;#it{p}_{T,true};#it{G(R)}_{sub}-#it{G(R)}_{true}",i);
    fh2PtTrueDeltaGR[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,100,-50.,50.);
    fOutput->Add(fh2PtTrueDeltaGR[i]);

    histName = Form("fh2PtTrueDeltaGRRel_%d",i);
    histTitle = Form("fh2PtTrueDeltaGRRel_%d;#it{p}_{T,true};(#it{G(R)}_{sub}-#it{G(R)}_{true})/#it{G(R)}_{true}",i);
    fh2PtTrueDeltaGRRel[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,200,-1.,1.);
    fOutput->Add(fh2PtTrueDeltaGRRel[i]);

    histName = Form("fhnGRResponse_%d",i);
    histTitle = Form("fhnGRResponse_%d;#it{G(R)}_{sub};#it{G(R)}_{true};#it{p}_{T,sub};#it{p}_{T,true}",i);
    fhnGRResponse[i] = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse0,nBins0,xmin0,xmax0);
    fOutput->Add(fhnGRResponse[i]);

    //Histos for true jets
    histName = Form("fh1PtTrue_%d",i);
    histTitle = Form("%s;#it{p}_{T};#it{N}",histName.Data());
    fh1PtTrue[i] = new TH1F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fh1PtTrue[i]);

    histName = Form("fh3DeltaGRNumRPtTrue_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#delta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRNumRPtTrue[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRNumRPtTrue[i]);

    histName = Form("fh3DeltaGRDenRPtTrue_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#Theta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRDenRPtTrue[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRDenRPtTrue[i]);

    histName = Form("fh2DeltaGRNumRPtTrue_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRNumRPtTrue[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRNumRPtTrue[i]);

    histName = Form("fh2DeltaGRDenRPtTrue_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRDenRPtTrue[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRDenRPtTrue[i]);

    //Histos for raw AA jets
    histName = Form("fh1PtRaw_%d",i);
    histTitle = Form("%s;#it{p}_{T};#it{N}",histName.Data());
    fh1PtRaw[i] = new TH1F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fh1PtRaw[i]);

    histName = Form("fh3DeltaGRNumRPtRaw_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#delta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRNumRPtRaw[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRNumRPtRaw[i]);

    histName = Form("fh3DeltaGRDenRPtRaw_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#Theta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRDenRPtRaw[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRDenRPtRaw[i]);

    histName = Form("fh2DeltaGRNumRPtRaw_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRNumRPtRaw[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRNumRPtRaw[i]);

    histName = Form("fh2DeltaGRDenRPtRaw_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRDenRPtRaw[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRDenRPtRaw[i]);

    //Histos for raw AA jets matched to MC jet
    histName = Form("fh1PtRawMatch_%d",i);
    histTitle = Form("%s;#it{p}_{T};#it{N}",histName.Data());
    fh1PtRawMatch[i] = new TH1F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fh1PtRawMatch[i]);

    histName = Form("fh3DeltaGRNumRPtRawMatch_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#delta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRNumRPtRawMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRNumRPtRawMatch[i]);

    histName = Form("fh3DeltaGRDenRPtRawMatch_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#Theta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRDenRPtRawMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRDenRPtRawMatch[i]);

    histName = Form("fh2DeltaGRNumRPtRawMatch_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRNumRPtRawMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRNumRPtRawMatch[i]);

    histName = Form("fh2DeltaGRDenRPtRawMatch_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRDenRPtRawMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRDenRPtRawMatch[i]);

    //Histos for matched jets
    histName = Form("fh1PtMatch_%d",i);
    histTitle = Form("%s;#it{p}_{T};#it{N}",histName.Data());
    fh1PtMatch[i] = new TH1F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fh1PtMatch[i]);

    histName = Form("fh3DeltaGRNumRPtMatch_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#delta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRNumRPtMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRNumRPtMatch[i]);

    histName = Form("fh3DeltaGRDenRPtMatch_%d",i);
    histTitle = Form("%s;#Sigma#it{p}_{Tk,i}#it{p}_{Tk,j}#Delta#it{R}_{ij}^{2}#Theta_{dR}(#it{r}-#it{R}_{ij});#it{r};#it{p}_{T,true}",histName.Data());
    fh3DeltaGRDenRPtMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsDGR,minDGR,maxDGR,nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh3DeltaGRDenRPtMatch[i]);

    histName = Form("fh2DeltaGRNumRPtMatch_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRNumRPtMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRNumRPtMatch[i]);

    histName = Form("fh2DeltaGRDenRPtMatch_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRDenRPtMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRDenRPtMatch[i]);

    histName = Form("fh2DeltaGRNumRPtTrueMatch_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRNumRPtTrueMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRNumRPtTrueMatch[i]);

    histName = Form("fh2DeltaGRDenRPtTrueMatch_%d",i);
    histTitle = Form("%s;#it{r};#it{p}_{T,true}",histName.Data());
    fh2DeltaGRDenRPtTrueMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsR,minR,maxR,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2DeltaGRDenRPtTrueMatch[i]);

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
    fTreeJetBkg = new TTree("fTreeJetSubGR", "fTreeJetSubGR");
    fTreeJetBkg->Branch("fJet1Vec","TLorentzVector",&fJet1Vec);
    fTreeJetBkg->Branch("fJet2Vec","TLorentzVector",&fJet2Vec);
    fTreeJetBkg->Branch("fJetSubVec","TLorentzVector",&fJetSubVec);
    fTreeJetBkg->Branch("fArea",&fArea,"fArea/F");
    fTreeJetBkg->Branch("fAreaPhi",&fAreaPhi,"fAreaPhi/F");
    fTreeJetBkg->Branch("fAreaEta",&fAreaEta,"fAreaEta/F");
    fTreeJetBkg->Branch("fRho",&fRho,"fRho/F");
    fTreeJetBkg->Branch("fRhoM",&fRhoM,"fRhoM/F");
    fTreeJetBkg->Branch("fMatch",&fMatch,"fMatch/I");
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  if(fCreateTree) PostData(2, fTreeJetBkg);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeGR::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeGR::FillHistograms()
{
  // Fill histograms.
  FillTrueJets();

  AliEmcalJet* jet1 = NULL;
  AliEmcalJet *jet2 = NULL;
  AliEmcalJet *jetS = NULL;
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  AliJetContainer *jetContS = GetJetContainer(fContainerSub);
  AliDebug(11,Form("NJets  Incl: %d  Csub: %d",jetCont->GetNJets(),jetContS->GetNJets()));
  if(jetCont && jetContS) {
    jetCont->ResetCurrentID();
    Double_t rmax = jetCont->GetJetRadius()+0.2;
    Double_t wr = 0.04;
    Int_t nr = TMath::CeilNint(rmax/wr);
    while((jet1 = jetCont->GetNextAcceptJet())) {
      Double_t fraction = 0.;
      fMatch = 0;
      fJet2Vec->SetPtEtaPhiM(0.,0.,0.,0.);
      if(fSingleTrackEmb) {
	AliVParticle *vp = GetEmbeddedConstituent(jet1);
	if(vp) {
	  fJet2Vec->SetPxPyPzE(vp->Px(),vp->Py(),vp->Pz(),vp->E());
	  fMatch = 1;
	}
      } else {
	jet2 = jet1->ClosestJet();
	if(!jet2) continue;

	fraction = jetCont->GetFractionSharedPt(jet1);
	fMatch = 1;
	if(fMinFractionShared>0.) {
	  if(fraction>fMinFractionShared) {
	    fJet2Vec->SetPxPyPzE(jet2->Px(),jet2->Py(),jet2->Pz(),jet2->E());
	    fMatch = 1;
	  } else
	    fMatch = 0;
	}
      }

      //Fill histograms for all AA jets
      Double_t ptcorr = jet1->Pt()-jetCont->GetRhoVal()*jet1->Area();
      fh1PtRaw[fCentBin]->Fill(ptcorr);
      if(fMatch==1) fh1PtRawMatch[fCentBin]->Fill(ptcorr);

      TArrayF numRaw = jet1->GetShapeProperties()->GetGRNumerator();
      TArrayF denRaw = jet1->GetShapeProperties()->GetGRDenominator();
      if(numRaw.GetSize()>0) {
	for(Int_t i = 0; i<nr; i++) {
	  Double_t r = i*wr + 0.5*wr;
	  fh3DeltaGRNumRPtRaw[fCentBin]->Fill(numRaw[i],r,ptcorr);
	  fh3DeltaGRDenRPtRaw[fCentBin]->Fill(denRaw[i],r,ptcorr);
	  fh2DeltaGRNumRPtRaw[fCentBin]->Fill(r,ptcorr,numRaw[i]);
	  fh2DeltaGRDenRPtRaw[fCentBin]->Fill(r,ptcorr,denRaw[i]);
	  if(fMatch==1) {
	    fh3DeltaGRNumRPtRawMatch[fCentBin]->Fill(numRaw[i],r,ptcorr);
	    fh3DeltaGRDenRPtRawMatch[fCentBin]->Fill(denRaw[i],r,ptcorr);
	    fh2DeltaGRNumRPtRawMatch[fCentBin]->Fill(r,ptcorr,numRaw[i]);
	    fh2DeltaGRDenRPtRawMatch[fCentBin]->Fill(r,ptcorr,denRaw[i]);
	  }
	}
      }

      //Fill histograms for matched jets
      if(fMatch==1) {
	fh1PtMatch[fCentBin]->Fill(ptcorr);

	//now get second derivative vs R and do final calculation
	TArrayF num = jet1->GetShapeProperties()->GetGRNumeratorSub();
	TArrayF den = jet1->GetShapeProperties()->GetGRDenominatorSub();
	if(num.GetSize()>0) {
	  for(Int_t i = 0; i<nr; i++) {
	    Double_t r = i*wr + 0.5*wr;
	    fh3DeltaGRNumRPtMatch[fCentBin]->Fill(num[i],r,ptcorr);
	    fh3DeltaGRDenRPtMatch[fCentBin]->Fill(den[i],r,ptcorr);
	    fh2DeltaGRNumRPtMatch[fCentBin]->Fill(r,ptcorr,num[i]);
	    fh2DeltaGRDenRPtMatch[fCentBin]->Fill(r,ptcorr,den[i]);
	    fh2DeltaGRNumRPtTrueMatch[fCentBin]->Fill(r,jet2->Pt(),num[i]);
	    fh2DeltaGRDenRPtTrueMatch[fCentBin]->Fill(r,jet2->Pt(),den[i]);

	    Double_t dGR = 0.;
	    fh2PtTrueDeltaGR[fCentBin]->Fill(jet2->Pt(),dGR);
	  }
	}
      }

      if(fCreateTree) {      
	fJet1Vec->SetPxPyPzE(jet1->Px(),jet1->Py(),jet1->Pz(),jet1->E());
	if(jetS && jetS->Pt()>0.) fJetSubVec->SetPtEtaPhiM(jetS->Pt(),jetS->Eta(),jetS->Phi(),jetS->M());
	else fJetSubVec->SetPtEtaPhiM(0.,0.,0.,0.);
	fArea = (Float_t)jet1->Area();
	fAreaPhi = (Float_t)jet1->AreaPhi();
	fAreaEta = (Float_t)jet1->AreaEta();
	fRho  = (Float_t)jetCont->GetRhoVal();
	fRhoM = (Float_t)jetCont->GetRhoMassVal();
	fNConst = (Int_t)jet1->GetNumberOfTracks();
	fTreeJetBkg->Fill();
      }
    } //jet1 loop
  }//jetCont


  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeGR::FillTrueJets() {

  AliEmcalJet* jet1 = NULL;
  AliJetContainer *jetCont = GetJetContainer(fContainerTrue);
  if(!jetCont)
    return kFALSE;

  AliDebug(11,Form("NJets True: %d",jetCont->GetNJets()));

  //create arrays
  const Int_t nr = TMath::CeilNint(fMaxR/fDRStep);
  TArrayF *fNum = new TArrayF(nr);
  TArrayF *fDen = new TArrayF(nr);

  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      fh1PtTrue[fCentBin]->Fill(jet1->Pt());
      
      //Double_t dev = CalcDeltaGR(jet1,fContainerTrue,fNum,fDen);//num,den);
      for(Int_t i = 0; i<nr; i++) {
	Double_t r = i*fDRStep + 0.5*fDRStep;
	fh3DeltaGRNumRPtTrue[fCentBin]->Fill(fNum->At(i),r,jet1->Pt());
	fh3DeltaGRDenRPtTrue[fCentBin]->Fill(fDen->At(i),r,jet1->Pt());
	fh2DeltaGRNumRPtTrue[fCentBin]->Fill(r,jet1->Pt(),fNum->At(i));
	fh2DeltaGRDenRPtTrue[fCentBin]->Fill(r,jet1->Pt(),fDen->At(i));
      }
    }
  }
  if(fNum)  delete fNum;
  if(fDen)  delete fDen;
  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskJetShapeGR::CalcDeltaGR(AliEmcalJet *jet, Int_t ic, TArrayF *fNum, TArrayF *fDen) { //Double_t *num, Double_t *den) {
  //Calculate G(R)

  //First clear the arrays
  const Int_t nr = TMath::CeilNint(fMaxR/fDRStep);
  for(Int_t i = 0; i<nr; i++) {
    fNum->SetAt(0.,i);
    fDen->SetAt(0.,i);
  }

  AliJetContainer *jetCont = GetJetContainer(ic); 
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  Double_t A = 0.; Double_t B = 0.;
  if(jet->GetNumberOfTracks()<2) return 0.;
  for(Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
    if(!vp1) continue;
    for(Int_t j=i+1; j<jet->GetNumberOfTracks(); j++) {
      vp2 = static_cast<AliVParticle*>(jet->TrackAt(j, jetCont->GetParticleContainer()->GetArray()));
      if(!vp2) continue;
      Double_t dphi = GetDeltaPhi(vp1->Phi(),vp2->Phi());
      Double_t dr2 = (vp1->Eta()-vp2->Eta())*(vp1->Eta()-vp2->Eta()) + dphi*dphi;
      if(dr2>0.) {
	for(Int_t k = 0; k<nr; k++) {
	  Double_t r = k*fDRStep + 0.5*fDRStep;
	  //	Double_t x = jetCont->GetJetRadius()-TMath::Sqrt(dr2);
	  Double_t dr = TMath::Sqrt(dr2);
	  Double_t x = r-dr;
	  //noisy function
	  Double_t noise = TMath::Exp(-x*x/(2*fDRStep*fDRStep))/(TMath::Sqrt(2.*TMath::Pi())*fDRStep);
	  //error function
	  Double_t erf = 0.5*(1.+TMath::Erf(x/(TMath::Sqrt(2.)*fDRStep)));
	  
	  A = vp1->Pt()*vp2->Pt()*dr2*noise;
	  B = vp1->Pt()*vp2->Pt()*dr2*erf;
	  fNum->AddAt(fNum->At(k)+A,k);
	  fDen->AddAt(fDen->At(k)+B,k);
	}
      }
    }
  }

  Double_t deltaGR = 0.;
  if(B>0.) deltaGR = A/B; //useless
  return deltaGR;
}

//________________________________________________________________________
Double_t AliAnalysisTaskJetShapeGR::CalcGR(AliEmcalJet *jet, Int_t ic) {
  //Calculate G(R)
  AliJetContainer *jetCont = GetJetContainer(ic); 
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  Double_t gR = 0.;
  Double_t wr = 0.04;
  const Int_t nr = TMath::CeilNint(jetCont->GetJetRadius()/wr);
  Double_t grArr[999] = {0.};
 
  for(Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
    for(Int_t j=i; j<jet->GetNumberOfTracks(); j++) {
      vp2 = static_cast<AliVParticle*>(jet->TrackAt(j, jetCont->GetParticleContainer()->GetArray()));
      Double_t dphi = GetDeltaPhi(vp1->Phi(),vp2->Phi());
      Double_t dr2 = (vp1->Eta()-vp2->Eta())*(vp1->Eta()-vp2->Eta()) + dphi*dphi;
      Int_t bin = TMath::FloorNint(TMath::Sqrt(dr2)/wr);
      Double_t gr = vp1->Pt()*vp2->Pt()*dr2;
      if(bin<nr) grArr[bin]+=gr;
      
      if(TMath::Sqrt(dr2)<jetCont->GetJetRadius())
	gR += gr;
    }
  }
  return gR;
}

//________________________________________________________________________
AliVParticle* AliAnalysisTaskJetShapeGR::GetEmbeddedConstituent(AliEmcalJet *jet) {

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
Double_t AliAnalysisTaskJetShapeGR::GetDeltaPhi(Double_t phi1,Double_t phi2) {
  //
  // Calculate azimuthal angle difference
  //

  Double_t dPhi = phi1-phi2;

  if(dPhi <-1.*TMath::Pi())  dPhi += TMath::TwoPi();
  if(dPhi > 1.*TMath::Pi())  dPhi -= TMath::TwoPi();

  return dPhi;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeGR::RetrieveEventObjects() {
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
void AliAnalysisTaskJetShapeGR::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

