//
// Do subtraction for jet shapes using derivatives arXiv:1211:2811
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

#include "AliAnalysisTaskJetShapeDeriv.h"

ClassImp(AliAnalysisTaskJetShapeDeriv)

//________________________________________________________________________
AliAnalysisTaskJetShapeDeriv::AliAnalysisTaskJetShapeDeriv() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetShapeDeriv", kTRUE),
  fContainerBase(0),
  fContainerNoEmb(1),
  fMinFractionShared(0),
  fSingleTrackEmb(kFALSE),
  fCreateTree(kFALSE),
  fJetMassVarType(kMass),
  fResponseReference(kDet),
  fUseSumw2(0),
  fTreeJetBkg(),
  fJet1Vec(new TLorentzVector()),
  fJet2Vec(new TLorentzVector()),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fM1st(0),
  fM2nd(0),
  fDeriv1st(0),
  fDeriv2nd(0),
  fMatch(0),
  fh2MSubMatch(0x0),
  fh2MSubPtRawAll(0x0),
  fh3MSubPtRawDRMatch(0x0),
  fh3MSubPtTrueLeadPt(0x0),
  fh3MTruePtTrueLeadPt(0x0),
  fh3PtTrueDeltaMLeadPt(0x0),
  fh3PtTrueDeltaMRelLeadPt(0x0),
  fhnMassResponse(0x0),
  fh2PtTrueSubFacV1(0x0),
  fh2PtRawSubFacV1(0x0),
  fh2PtCorrSubFacV1(0x0),
  fh2NConstSubFacV1(0x0),
  fh2PtTrueSubFacV2(0x0),
  fh2PtRawSubFacV2(0x0),
  fh2PtCorrSubFacV2(0x0),
  fh2NConstSubFacV2(0x0)
{
  // Default constructor.

  fh2MSubMatch             = new TH2F*[fNcentBins];
  fh2MSubPtRawAll          = new TH2F*[fNcentBins];
  fh3MSubPtRawDRMatch      = new TH3F*[fNcentBins];
  fh3MSubPtTrueLeadPt      = new TH3F*[fNcentBins];
  fh3MTruePtTrueLeadPt     = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMLeadPt    = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMRelLeadPt = new TH3F*[fNcentBins];
  fhnMassResponse          = new THnSparse*[fNcentBins];
  fh2PtTrueSubFacV1        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV1         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV1        = new TH2F*[fNcentBins];
  fh2NConstSubFacV1        = new TH2F*[fNcentBins];
  fh2PtTrueSubFacV2        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV2         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV2        = new TH2F*[fNcentBins];
  fh2NConstSubFacV2        = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubMatch[i]             = 0;
    fh2MSubPtRawAll[i]          = 0;
    fh3MSubPtRawDRMatch[i]      = 0;
    fh3MSubPtTrueLeadPt[i]      = 0;
    fh3MTruePtTrueLeadPt[i]     = 0;
    fh3PtTrueDeltaMLeadPt[i]    = 0;
    fh3PtTrueDeltaMRelLeadPt[i] = 0;
    fhnMassResponse[i]          = 0;
    fh2PtTrueSubFacV1[i]        = 0;
    fh2PtRawSubFacV1[i]         = 0;
    fh2PtCorrSubFacV1[i]        = 0;
    fh2NConstSubFacV1[i]        = 0;
    fh2PtTrueSubFacV2[i]        = 0;
    fh2PtRawSubFacV2[i]         = 0;
    fh2PtCorrSubFacV2[i]        = 0;
    fh2NConstSubFacV2[i]        = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  if(fCreateTree) DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeDeriv::AliAnalysisTaskJetShapeDeriv(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fContainerNoEmb(1),
  fMinFractionShared(0),
  fSingleTrackEmb(kFALSE),
  fCreateTree(kFALSE),
  fJetMassVarType(kMass),
  fResponseReference(kDet),
  fUseSumw2(0),
  fTreeJetBkg(0),
  fJet1Vec(new TLorentzVector()),
  fJet2Vec(new TLorentzVector()),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fM1st(0),
  fM2nd(0),
  fDeriv1st(0),
  fDeriv2nd(0),
  fMatch(0),
  fh2MSubMatch(0x0),
  fh2MSubPtRawAll(0x0),
  fh3MSubPtRawDRMatch(0x0),
  fh3MSubPtTrueLeadPt(0x0),
  fh3MTruePtTrueLeadPt(0x0),
  fh3PtTrueDeltaMLeadPt(0x0),
  fh3PtTrueDeltaMRelLeadPt(0x0),
  fhnMassResponse(0x0),
  fh2PtTrueSubFacV1(0x0),
  fh2PtRawSubFacV1(0x0),
  fh2PtCorrSubFacV1(0x0),
  fh2NConstSubFacV1(0x0),
  fh2PtTrueSubFacV2(0x0),
  fh2PtRawSubFacV2(0x0),
  fh2PtCorrSubFacV2(0x0),
  fh2NConstSubFacV2(0x0)
{
  // Standard constructor.

  fh2MSubMatch             = new TH2F*[fNcentBins];
  fh2MSubPtRawAll          = new TH2F*[fNcentBins];
  fh3MSubPtRawDRMatch      = new TH3F*[fNcentBins];
  fh3MSubPtTrueLeadPt      = new TH3F*[fNcentBins];
  fh3MTruePtTrueLeadPt     = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMLeadPt    = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMRelLeadPt = new TH3F*[fNcentBins];
  fhnMassResponse          = new THnSparse*[fNcentBins];
  fh2PtTrueSubFacV1        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV1         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV1        = new TH2F*[fNcentBins];
  fh2NConstSubFacV1        = new TH2F*[fNcentBins];
  fh2PtTrueSubFacV2        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV2         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV2        = new TH2F*[fNcentBins];
  fh2NConstSubFacV2        = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
   fh2MSubMatch[i]             = 0;
    fh2MSubPtRawAll[i]          = 0;
    fh3MSubPtRawDRMatch[i]      = 0;
    fh3MSubPtTrueLeadPt[i]      = 0;
    fh3MTruePtTrueLeadPt[i]     = 0;
    fh3PtTrueDeltaMLeadPt[i]    = 0;
    fh3PtTrueDeltaMRelLeadPt[i] = 0;
    fhnMassResponse[i]          = 0;
    fh2PtTrueSubFacV1[i]        = 0;
    fh2PtRawSubFacV1[i]         = 0;
    fh2PtCorrSubFacV1[i]        = 0;
    fh2NConstSubFacV1[i]        = 0;
    fh2PtTrueSubFacV2[i]        = 0;
    fh2PtRawSubFacV2[i]         = 0;
    fh2PtCorrSubFacV2[i]        = 0;
    fh2NConstSubFacV2[i]        = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  if(fCreateTree) DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeDeriv::~AliAnalysisTaskJetShapeDeriv()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetShapeDeriv::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsPt  = 200;
  const Double_t minPt = -50.;
  const Double_t maxPt = 150.;

  Int_t nBinsM  = 100;
  Double_t minM = -20.;
  Double_t maxM = 80.;
  if(fJetMassVarType==kRatMPt) {
    nBinsM = 100;
    minM   = -0.2;
    maxM   = 0.8;
  }

  Int_t nBinsDM  = 100;
  Double_t minDM = -25.;
  Double_t maxDM = 25.;
  if(fJetMassVarType==kRatMPt) {
    nBinsDM = 100;
    minDM   = -0.5;
    maxDM   = 0.5;
  }

  const Int_t nBinsDRToLJ  = 20; //distance to leading jet in Pb-Pb only event
  const Double_t minDRToLJ = 0.;
  const Double_t maxDRToLJ = 1.;

  const Int_t nBinsPtLead = 20;
  const Double_t minPtLead = 0.;
  const Double_t maxPtLead = 20.;

  const Int_t nBinsV1  = 60;
  const Double_t minV1 = -60.;
  const Double_t maxV1 = 0.;

  const Int_t nBinsV2  = 60;
  const Double_t minV2 = -30.;
  const Double_t maxV2 = 0.;

  //Binning for THnSparse
  const Int_t nBinsSparse0 = 5;
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt,nBinsPtLead};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt, minPtLead};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt, maxPtLead};

  TString histName = "";
  TString histTitle = "";
  TString varName = "#it{M}_{jet}";
  if(fJetMassVarType==kRatMPt) varName = "#it{M}_{jet}/#it{p}_{T,jet}";

  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = Form("fh2MSubMatch_%d",i);
    histTitle = Form("fh2MSubMatch_%d;%s;match",i,varName.Data());
    fh2MSubMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,2,-0.5,1.5);
    fOutput->Add(fh2MSubMatch[i]);

    histName = Form("fh2MSubPtRawAll_%d",i);
    histTitle = Form("fh2MSubPtRawAll_%d;%s;#it{p}_{T}",i,varName.Data());
    fh2MSubPtRawAll[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2MSubPtRawAll[i]);

    histName = Form("fh3MSubPtRawDRMatch_%d",i);
    histTitle = Form("fh3MSubPtRawDRMatch_%d;%s;#it{p}_{T}",i,varName.Data());
    fh3MSubPtRawDRMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt,nBinsDRToLJ,minDRToLJ,maxDRToLJ);
    fOutput->Add(fh3MSubPtRawDRMatch[i]);

    histName = Form("fh3MSubPtTrueLeadPt_%d",i);
    histTitle = Form("fh3MSubPtTrueLeadPt_%d;%s;#it{p}_{T}",i,varName.Data());
    fh3MSubPtTrueLeadPt[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt,nBinsPtLead,minPtLead,maxPtLead);
    fOutput->Add(fh3MSubPtTrueLeadPt[i]);

    histName = Form("fh3MTruePtTrueLeadPt_%d",i);
    histTitle = Form("fh3MTruePtTrueLeadPt_%d;%s;#it{p}_{T}",i,varName.Data());
    fh3MTruePtTrueLeadPt[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt,nBinsPtLead,minPtLead,maxPtLead);
    fOutput->Add(fh3MTruePtTrueLeadPt[i]);

    histName = Form("fh3PtTrueDeltaMLeadPt_%d",i);
    histTitle = Form("fh3PtTrueDeltaMLeadPt_%d;#it{p}_{T,true};#Delta %s",i,varName.Data());
    fh3PtTrueDeltaMLeadPt[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsDM,minDM,maxDM,nBinsPtLead,minPtLead,maxPtLead);
    fOutput->Add(fh3PtTrueDeltaMLeadPt[i]);

    histName = Form("fh3PtTrueDeltaMRelLeadPt_%d",i);
    histTitle = Form("fh3PtTrueDeltaMRelLeadPt_%d;#it{p}_{T,true};Rel #Delta %s",i,varName.Data());
    fh3PtTrueDeltaMRelLeadPt[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,400,-1.,3.,nBinsPtLead,minPtLead,maxPtLead);
    fOutput->Add(fh3PtTrueDeltaMRelLeadPt[i]);

    histName = Form("fhnMassResponse_%d",i);
    histTitle = Form("fhnMassResponse_%d;%s sub;%s true;#it{p}_{T,sub};#it{p}_{T,true};#it{p}_{T,lead trk}",i,varName.Data(),varName.Data());
    fhnMassResponse[i] = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse0,nBins0,xmin0,xmax0);
    fOutput->Add(fhnMassResponse[i]);

    //derivative histograms
    histName = Form("fh2PtTrueSubFacV1_%d",i);
    histTitle = Form("fh2PtTrueSubFacV1_%d;#it{p}_{T,true};-(#rho+#rho_{m})V_{1}",i);
    fh2PtTrueSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV1,minV1,maxV1);

    histName = Form("fh2PtRawSubFacV1_%d",i);
    histTitle = Form("fh2PtRawSubFacV1_%d;#it{p}_{T,raw};-(#rho+#rho_{m})V_{1}",i);
    fh2PtRawSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV1,minV1,maxV1);

    histName = Form("fh2PtCorrSubFacV1_%d",i);
    histTitle = Form("fh2PtCorrSubFacV1_%d;#it{p}_{T,corr};-(#rho+#rho_{m})V_{1}",i);
    fh2PtCorrSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV1,minV1,maxV1);

    histName = Form("fh2NConstSubFacV1_%d",i);
    histTitle = Form("fh2NConstSubFacV1_%d;#it{N}_{const};-(#rho+#rho_{m})V_{1}",i);
    fh2NConstSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,100,0.,200.);

    histName = Form("fh2PtTrueSubFacV2_%d",i);
    histTitle = Form("fh2PtTrueSubFacV2_%d;#it{p}_{T,true};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2PtTrueSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV2,minV2,maxV2);

    histName = Form("fh2PtRawSubFacV2_%d",i);
    histTitle = Form("fh2PtRawSubFacV2_%d;#it{p}_{T,raw};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2PtRawSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV2,minV2,maxV2);

    histName = Form("fh2PtCorrSubFacV2_%d",i);
    histTitle = Form("fh2PtCorrSubFacV2_%d;#it{p}_{T,corr};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2PtCorrSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV2,minV2,maxV2);

    histName = Form("fh2NConstSubFacV2_%d",i);
    histTitle = Form("fh2NConstSubFacV2_%d;#it{N}_{const};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2NConstSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,100,0.,200.);

  }

  if(fUseSumw2) {
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
  }

  TH1::AddDirectory(oldStatus);

  // Create a tree.
  if(fCreateTree) {
    fTreeJetBkg = new TTree("fTreeJetBkg", "fTreeJetBkg");
    fTreeJetBkg->Branch("fJet1Vec","TLorentzVector",&fJet1Vec);
    fTreeJetBkg->Branch("fJet2Vec","TLorentzVector",&fJet2Vec);
    fTreeJetBkg->Branch("fArea",&fArea,"fArea/F");
    fTreeJetBkg->Branch("fAreaPhi",&fAreaPhi,"fAreaPhi/F");
    fTreeJetBkg->Branch("fAreaEta",&fAreaEta,"fAreaEta/F");
    fTreeJetBkg->Branch("fRho",&fRho,"fRho/F");
    fTreeJetBkg->Branch("fRhoM",&fRhoM,"fRhoM/F");
    fTreeJetBkg->Branch("fNConst",&fNConst,"fNConst/I");
    fTreeJetBkg->Branch("fM1st",&fM1st,"fM1st/F");
    fTreeJetBkg->Branch("fM2nd",&fM2nd,"fM2nd/F");
    fTreeJetBkg->Branch("fDeriv1st",&fDeriv1st,"fDeriv1st/F");
    fTreeJetBkg->Branch("fDeriv2nd",&fDeriv2nd,"fDeriv2nd/F");
    fTreeJetBkg->Branch("fMatch",&fMatch,"fMatch/I");
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  if(fCreateTree) PostData(2, fTreeJetBkg);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeDeriv::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeDeriv::FillHistograms()
{
  // Fill histograms.

  AliEmcalJet *jet1  = NULL; //AA jet
  AliEmcalJet *jet2  = NULL; //Embedded Pythia jet
  AliEmcalJet *jetR  = NULL; //true jet for response matrix

  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  fRho  = (Float_t)jetCont->GetRhoVal();
  fRhoM = (Float_t)jetCont->GetRhoMassVal();

  //Get leading jet in Pb-Pb event without embedded objects
  AliJetContainer *jetContNoEmb = GetJetContainer(fContainerNoEmb);
  AliEmcalJet *jetL = NULL;
  if(jetContNoEmb) jetL = jetContNoEmb->GetLeadingJet("rho");

  jetCont->ResetCurrentID();
  while((jet1 = jetCont->GetNextAcceptJet())) {
    jet2 = NULL;

    Double_t mjet1 = jet1->GetSecondOrderSubtracted();
    Double_t ptjet1 = jet1->Pt()-jetCont->GetRhoVal()*jet1->Area();
    Double_t var = mjet1;
    if(fJetMassVarType==kRatMPt) {
      if(ptjet1>0. || ptjet1<0.) var = mjet1/ptjet1;
      else var = -999.;
    }

    //Fill histograms for all AA jets
    fh2MSubPtRawAll[fCentBin]->Fill(var,ptjet1);
    fh2PtRawSubFacV1[fCentBin]->Fill(jet1->Pt(),-1.*(fRho+fRhoM)*jet1->GetFirstDerivative());
    fh2PtCorrSubFacV1[fCentBin]->Fill(jet1->Pt()-fRho*jet1->Area(),-1.*(fRho+fRhoM)*jet1->GetFirstDerivative());
    fh2NConstSubFacV1[fCentBin]->Fill(jet1->GetNumberOfTracks(),-1.*(fRho+fRhoM)*jet1->GetFirstDerivative());
    fh2PtRawSubFacV2[fCentBin]->Fill(jet1->Pt(),0.5*(fRho+fRhoM)*(fRho+fRhoM)*jet1->GetSecondDerivative());
    fh2PtCorrSubFacV2[fCentBin]->Fill(jet1->Pt()-fRho*jet1->Area(),0.5*(fRho+fRhoM)*(fRho+fRhoM)*jet1->GetSecondDerivative());
    fh2NConstSubFacV2[fCentBin]->Fill(jet1->GetNumberOfTracks(),0.5*(fRho+fRhoM)*(fRho+fRhoM)*jet1->GetSecondDerivative());

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
      fraction = jetCont->GetFractionSharedPt(jet1);
      fMatch = 1;
      if(fMinFractionShared>0.) {
        if(fraction>fMinFractionShared) {
          fJet2Vec->SetPxPyPzE(jet2->Px(),jet2->Py(),jet2->Pz(),jet2->E());
          fMatch = 1;

          //choose jet type for true axis of response matrix
          if(fResponseReference==kDet) 
            jetR = jet2;
          else if(fResponseReference==kPart)
            jetR = jet2->GetTaggedJet();
        } else
          fMatch = 0;
      }
    }

    //Fill histograms for matched jets
    fh2MSubMatch[fCentBin]->Fill(var,fMatch);
    if(fMatch==1) {
      Double_t drToLJ = -1.;
      if(jetL) drToLJ = jet1->DeltaR(jetL);
      fh3MSubPtRawDRMatch[fCentBin]->Fill(var,ptjet1,drToLJ);
      if(jetR) {
        Double_t var2 = jetR->M();
        if(fJetMassVarType==kRatMPt) {
          if(jetR->Pt()>0. || jetR->Pt()<0.) var2 = jetR->M()/jetR->Pt();
        }
        fh3MSubPtTrueLeadPt[fCentBin]->Fill(var,jetR->Pt(),jet1->MaxTrackPt());
        fh3MTruePtTrueLeadPt[fCentBin]->Fill(var2,jetR->Pt(),jet1->MaxTrackPt());
        fh3PtTrueDeltaMLeadPt[fCentBin]->Fill(jetR->Pt(),var-var2,jet1->MaxTrackPt());
        if(jetR->M()>0.) fh3PtTrueDeltaMRelLeadPt[fCentBin]->Fill(jetR->Pt(),(var-var2)/var2,jet1->MaxTrackPt());
        Double_t varsp[5] = {var,var2,ptjet1,jetR->Pt(),jet1->MaxTrackPt()};//MRec,MTrue,PtRec,PtTrue,PtLeadRec
        fhnMassResponse[fCentBin]->Fill(varsp);
      }
    }

    if(fCreateTree) {      
      fJet1Vec->SetPxPyPzE(jet1->Px(),jet1->Py(),jet1->Pz(),jet1->E());
      fArea = (Float_t)jet1->Area();
      fAreaPhi = (Float_t)jet1->AreaPhi();
      fAreaEta = (Float_t)jet1->AreaEta();
      fNConst = (Int_t)jet1->GetNumberOfTracks();
      fM1st   = (Float_t)jet1->GetFirstOrderSubtracted();
      fM2nd   = (Float_t)jet1->GetSecondOrderSubtracted();
      fDeriv1st = (Float_t)jet1->GetFirstDerivative();
      fDeriv2nd = (Float_t)jet1->GetSecondDerivative();
      fTreeJetBkg->Fill();
    }
  }

  return kTRUE;
}

//________________________________________________________________________
AliVParticle* AliAnalysisTaskJetShapeDeriv::GetEmbeddedConstituent(AliEmcalJet *jet) {

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
Bool_t AliAnalysisTaskJetShapeDeriv::RetrieveEventObjects() {
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
void AliAnalysisTaskJetShapeDeriv::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

