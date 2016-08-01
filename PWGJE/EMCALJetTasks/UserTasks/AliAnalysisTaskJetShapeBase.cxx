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
#include <TGrid.h>

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

#include "AliAnalysisTaskJetShapeBase.h"

ClassImp(AliAnalysisTaskJetShapeBase)

//________________________________________________________________________
AliAnalysisTaskJetShapeBase::AliAnalysisTaskJetShapeBase() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetShapeBase", kTRUE),
  fContainerBase(0),
  fContainerSub(1),
  fContainerNoEmb(2),
  fContainerOverlap(3),
  fMinFractionShared(0),
  fSingleTrackEmb(kFALSE),
  fCreateTree(kFALSE),
  fJetMassVarType(kMass),
  fResponseReference(kDet),
  fUseSumw2(0),
  fOverlap(0),
  fRadius(0.4),
  fTreeJetBkg(0),
  fJet1Vec(),
  fJet2Vec(),
  fJetSubVec(),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fMatch(0),
  fMinLabelEmb(-kMaxInt),
  fMaxLabelEmb(kMaxInt),
  fSmallSyst(0),
  fh2MSubMatch(0x0),
  fh2MSubPtRawAll(0x0),
  fh3MSubPtRawDRMatch(0x0),
  fh3MSubPtTrueLeadPt(0x0),
  fh3MTruePtTrueLeadPt(0x0),
  fh3PtTrueDeltaMLeadPt(0x0),
  fh3PtTrueDeltaMRelLeadPt(0x0),
  fhnMassResponse(0x0),
  //fhnDeltaMass(0),
  fhRjetTrvspTj(0x0),
  fhNJetsSelEv(0x0),
  fhJetEtaPhiOvl(0x0),
  fhpTTracksJet1(0x0),
  fhpTTracksJetO(0x0),
  fhpTTracksCont(0x0),
  fhptjetSMinusSingleTrack(0x0),
  fhJetSubMatchEtaPhiPt(0x0),
  fTreeEmb(0),
  fFromTree(0),
  fPathTreeinputFile(""),
  fTreeinputName("fTreeJet"),
  fBranchJDetName("fJetDet"),
  fBranchJParName("fJetPar"),
  fThisEntry(0),
  fMaxTreeEntries(0),
  fVecD(0x0),
  fVecP(0x0)
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
  //fhnDeltaMass             = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubMatch[i]             = 0;
    fh2MSubPtRawAll[i]          = 0;
    fh3MSubPtRawDRMatch[i]      = 0;
    fh3MSubPtTrueLeadPt[i]      = 0;
    fh3MTruePtTrueLeadPt[i]     = 0;
    fh3PtTrueDeltaMLeadPt[i]    = 0;
    fh3PtTrueDeltaMRelLeadPt[i] = 0;
    fhnMassResponse[i]          = 0;
    //fhnDeltaMass[i]             = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  if(fCreateTree) DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeBase::AliAnalysisTaskJetShapeBase(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fContainerSub(1),
  fContainerNoEmb(2),
  fContainerOverlap(3),
  fMinFractionShared(0),
  fSingleTrackEmb(kFALSE),
  fCreateTree(kFALSE),
  fJetMassVarType(kMass),
  fResponseReference(kDet),
  fUseSumw2(0),
  fOverlap(0),
  fRadius(0.4),
  fTreeJetBkg(0),
  fJet1Vec(),
  fJet2Vec(),
  fJetSubVec(),
  fArea(0),
  fAreaPhi(0),
  fAreaEta(0),
  fRho(0),
  fRhoM(0),
  fNConst(0),
  fMatch(0),
  fMinLabelEmb(-kMaxInt),
  fMaxLabelEmb(kMaxInt),
  fSmallSyst(0),
  fh2MSubMatch(0x0),
  fh2MSubPtRawAll(0x0),
  fh3MSubPtRawDRMatch(0x0),
  fh3MSubPtTrueLeadPt(0x0),
  fh3MTruePtTrueLeadPt(0x0),
  fh3PtTrueDeltaMLeadPt(0x0),
  fh3PtTrueDeltaMRelLeadPt(0x0),
  fhnMassResponse(0x0),
  //fhnDeltaMass(0),
  fhRjetTrvspTj(0x0),
  fhNJetsSelEv(0x0),
  fhJetEtaPhiOvl(0x0),
  fhpTTracksJet1(0x0),
  fhpTTracksJetO(0x0),
  fhpTTracksCont(0x0),
  fhptjetSMinusSingleTrack(0x0),
  fhJetSubMatchEtaPhiPt(0x0),
  fTreeEmb(0),
  fFromTree(0),
  fPathTreeinputFile(""),
  fTreeinputName("fTreeJet"),
  fBranchJDetName("fJetDet"),
  fBranchJParName("fJetPar"),
  fThisEntry(0),
  fMaxTreeEntries(0),
  fVecD(0x0),
  fVecP(0x0)
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
  //fhnDeltaMass             = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubMatch[i]             = 0;
    fh2MSubPtRawAll[i]          = 0;
    fh3MSubPtRawDRMatch[i]      = 0;
    fh3MSubPtTrueLeadPt[i]      = 0;
    fh3MTruePtTrueLeadPt[i]     = 0;
    fh3PtTrueDeltaMLeadPt[i]    = 0;
    fh3PtTrueDeltaMRelLeadPt[i] = 0;
    fhnMassResponse[i]          = 0;
    //fhnDeltaMass[i]             = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  if(fCreateTree) DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeBase::~AliAnalysisTaskJetShapeBase()
{
  // Destructor.
  delete fTreeEmb;
  delete fVecD;
  delete fVecP;
}

//________________________________________________________________________
void AliAnalysisTaskJetShapeBase::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  //Bool_t oldStatus = TH1::AddDirectoryStatus();
  //TH1::AddDirectory(kFALSE);

  if(!fPathTreeinputFile.IsNull()){
     SetTreeFromFile(fPathTreeinputFile, fTreeinputName);
     if(!fTreeEmb) AliFatal("Something went wrong in setting the tree");
     
     fTreeEmb->SetBranchAddress(fBranchJDetName, &fVecD);
     fTreeEmb->SetBranchAddress(fBranchJParName, &fVecP);
     
     fMaxTreeEntries = fTreeEmb->GetEntries();
     
     fTreeEmb->GetEntry(0);
     fTreeEmb->Show();
     
     
  }
  
  const Int_t nBinsPt  = 200;
  const Double_t minPt = -50.;
  const Double_t maxPt = 150.;

  Int_t nBinsM  = 100;
  Double_t minM = -20.;
  Double_t maxM = 80.;
  if(fSmallSyst) maxM = 40.;
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
  Int_t nBinsDpT  = 100;
  Double_t minDpT = -50.;
  Double_t maxDpT = 50.;

  const Int_t nBinsDRToLJ  = 20; //distance to leading jet in Pb-Pb only event
  const Double_t minDRToLJ = 0.;
  const Double_t maxDRToLJ = 1.;

  const Int_t nBinsPtLead = 20;
  const Double_t minPtLead = 0.;
  const Double_t maxPtLead = 20.;

  //These are good for pPb
  Int_t nBinsRho = 50;
  Double_t minRho = 0.;
  Double_t maxRho = 20.;
  Int_t nBinsRhom = 50;
  Double_t minRhom = 0.;
  Double_t maxRhom = 1.;

  //Binning for THnSparse
  const Int_t nBinsSparse0 = 9;
  //Mass sub;Mass true;#it{p}_{T,sub};#it{p}_{T,true};#it{p}_{T,lead trk}; #rho; #rho_{m}
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt,nBinsPtLead, nBinsRho, nBinsRhom, nBinsM, nBinsPt};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt, minPtLead, minRho, minRhom, minM, minPt};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt, maxPtLead, maxRho, maxRhom, maxM, maxPt};

  const Int_t nBinsSparse0b = 10;
  //Mass sub;Mass true;#it{p}_{T,sub};#it{p}_{T,true};#it{p}_{T,lead trk}; #rho; #rho_{m}; Mass unsub; #it{p}_{T,unsub}
  const Int_t nBins0b[nBinsSparse0b] = {nBinsM, nBinsM, nBinsPt, nBinsPt, nBinsM, nBinsPt, nBinsRho, nBinsRhom, nBinsM, nBinsPt};
  const Double_t xmin0b[nBinsSparse0b]  = { minM, minM, minPt, minPt, minM, minPt, minRho, minRhom, minM, minPt};
  const Double_t xmax0b[nBinsSparse0b]  = { maxM, maxM, maxPt, maxPt, maxM, maxPt, maxRho, maxRhom, maxM, maxPt};
  
  const Int_t nBinsSparse1 = 7;
  // #it{M}_{det,Const} - #it{M}_{part}; #it{p}_{T,det,Const} - #it{p}_{T,part}; #it{M}_{det,Const};  #it{M}_{part}; #it{p}_{T,det,Const}; #it{p}_{T,part}; #it{p}_{T,det,A}
  const Int_t nBins1[nBinsSparse1] = {nBinsDM,nBinsDpT,nBinsM,nBinsM,nBinsPt,nBinsPt,nBinsPt};
  const Double_t xmin1[nBinsSparse1]  = { minDM, minDpT, minM, minM, minPt, minPt, minPt};
  const Double_t xmax1[nBinsSparse1]  = { maxDM, maxDpT, maxM, maxM, maxPt, maxPt, maxPt};

  const Int_t nBinsSparse2 = 10;
  //#it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{unsub} - #it{M}_{part}; #it{p}_{T,unsub} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}
  const Int_t nBins2[nBinsSparse2] = {nBinsDM, nBinsDpT, nBinsDM, nBinsDpT, nBinsM, nBinsM, nBinsPt, nBinsPt, nBinsRho, nBinsRhom};
  const Double_t xmin2[nBinsSparse2]  = {minDM, minDpT, minDM, minDpT, minM, minM, minPt, minPt, minRho, minRhom};
  const Double_t xmax2[nBinsSparse2]  = {maxDM, maxDpT, maxDM, maxDpT, maxM, maxM, maxPt, maxPt, maxRho, maxRhom};

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
    if(fFromTree) {
       histTitle = Form("fhnMassResponse_%d; %s sub; %s true;#it{p}_{T,sub};#it{p}_{T,true};%s (emb, det); #it{p}_{T,emb det}; #rho; #rho_{m}; %s unsub; #it{p}_{T,unsub}", i, varName.Data(),varName.Data(), varName.Data(),varName.Data());
       fhnMassResponse[i] = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse0b, nBins0b, xmin0b, xmax0b);

    } else{
       histTitle = Form("fhnMassResponse_%d;%s sub;%s true;#it{p}_{T,sub};#it{p}_{T,true};#it{p}_{T,lead trk};  #rho; #rho_{m}",i,varName.Data(),varName.Data());
       fhnMassResponse[i] = new THnSparseF(histName.Data(),histTitle.Data(), nBinsSparse0, nBins0, xmin0, xmax0);
    }
    fOutput->Add(fhnMassResponse[i]);
    
    //histName = Form("fhnDeltaMass_%d", i);
    //histTitle = Form("%s; #it{M}_{det,Const} - #it{M}_{part}; #it{p}_{T,det,Const} - #it{p}_{T,part}; #it{M}_{det,Const};  #it{M}_{part}; #it{p}_{T,det,Const}; #it{p}_{T,part}; #it{p}_{T,det,A}",histName.Data());
    //Printf("Nuber of bins %d - write first %d, %f, %f , building %s", nBinsSparse1, nBins1[0], xmin1[0], xmax1[0], histName.Data());
    //fhnDeltaMass[i] = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse1,nBins1,xmin1,xmax1);
    //fOutput->Add(fhnDeltaMass[i]);

  }

  //Chiara's histograms: rho and rhom correlation with pT and mass at reco level with no subtraction
  histName = "fhnDeltaMassAndBkgInfo";
  histTitle = Form("%s; #it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{unsub} - #it{M}_{part}; #it{p}_{T,unsub} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}",histName.Data()); // #it{M}_{unsub} is also deltaM unsub when M_part is zero
  
  fhnDeltaMassAndBkgInfo = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse2,nBins2,xmin2,xmax2);
  fOutput->Add(fhnDeltaMassAndBkgInfo);
  
  if(fOverlap){
     fhRjetTrvspTj = new TH2F("fhRjetTrvspTj", ";R(jet, track);p_{T,jet}", 100, 0., 10., nBinsPt, minPt, maxPt);
     fOutput->Add(fhRjetTrvspTj);
     fhNJetsSelEv = new TH1F("fhNJetsSelEv", "N of jets selected; #it{N}_{jets}/ev;Entries", 20., 0.,19);
     fOutput->Add(fhNJetsSelEv);
     
     fhJetEtaPhiOvl = new TH2F("fhJetEtaPhiOvl", "#eta - #varphi distribution of selected jets; #eta; #varphi", 24., -0.6, 0.6, 50, (-1)*TMath::Pi(), TMath::Pi());
     fOutput->Add(fhJetEtaPhiOvl);
     
     fhpTTracksJetO = new TH1F("hTrackpTO", "Track pT (signal jet); p_{T}", 500,0.,50.);
     fOutput->Add(fhpTTracksJetO);

  }
  
  fhJetSubMatchEtaPhiPt = new TH3F("fhJetSubMatchEtaPhiPt", "Bkg-subtracted matched jet;#eta ; #phi; #it{p}_{T} (GeV/#it{c})", 24., -0.6, 0.6, 50, (-1)*TMath::Pi(), TMath::Pi(), nBinsPt, minPt, maxPt); 
  fOutput->Add(fhJetSubMatchEtaPhiPt);

  fhpTTracksJet1 = new TH1F("hTrackpT1", "Track pT ; p_{T}", 500,0.,50.);
  fOutput->Add(fhpTTracksJet1);
  fhpTTracksCont = new TH1F(Form("fhpTTracksCont"), "Track pT (container) ; p_{T}", 500,0.,50.);
  fOutput->Add(fhpTTracksCont);
  

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

  //TH1::AddDirectory(oldStatus);

  // Create a tree.
  //if(fCreateTree) {
  	  fJet1Vec   = new TLorentzVector();
  	  fJet2Vec   = new TLorentzVector();
  	  fJetSubVec = new TLorentzVector();
  	  Printf("Creating tree");
    fTreeJetBkg = new TTree("fTreeJetSubConst", "fTreeJetSubConst");
    fTreeJetBkg->Branch("fJet1Vec.",fJet1Vec);
    fTreeJetBkg->Branch("fJet2Vec.",fJet2Vec);
    fTreeJetBkg->Branch("fJetSubVec.",fJetSubVec);
    fTreeJetBkg->Branch("fArea",&fArea,"fArea/F");
    fTreeJetBkg->Branch("fAreaPhi",&fAreaPhi,"fAreaPhi/F");
    fTreeJetBkg->Branch("fAreaEta",&fAreaEta,"fAreaEta/F");
    fTreeJetBkg->Branch("fRho",&fRho,"fRho/F");
    fTreeJetBkg->Branch("fRhoM",&fRhoM,"fRhoM/F");
    fTreeJetBkg->Branch("fMatch",&fMatch,"fMatch/I");
    fTreeJetBkg->Branch("fNConst",&fNConst, "fNConst/I");
    
  //}
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
  //if(fCreateTree) 
  PostData(2, fTreeJetBkg);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeBase::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
AliVParticle* AliAnalysisTaskJetShapeBase::GetEmbeddedConstituent(AliEmcalJet *jet) {

  AliParticleContainer *partContEmb = GetParticleContainer(); //the first particle container given is the one with embedded track(s)

  AliVParticle *vpe = 0x0; //embedded particle
  for(Int_t ip = partContEmb->GetNParticles()-1; ip>-1; ip--){
      AliVParticle *vp = partContEmb->GetParticle(ip);
      if(!vp){
      	 AliDebug(2, Form("Particle %d not found", ip));
      	 continue;
      }
      Int_t lab = TMath::Abs(vp->GetLabel());
      if (lab < fMinLabelEmb || lab > fMaxLabelEmb + fMaxTreeEntries)
      	 continue;
      if(!vpe) vpe = vp;
      else if(vp->Pt()>vpe->Pt()) vpe =vp;
      fThisEntry = lab - fMinLabelEmb;
      //Printf("Label of embedded track %d - %d = %d", lab, fMinLabelEmb, fThisEntry);
  }
  
  Double_t deltaR = 99;
  if(vpe && jet) deltaR = jet->DeltaR(vpe);
  if(deltaR < 0.2) return vpe;
  else return 0x0;
}

//__________________________________________________________________________________________________

Int_t AliAnalysisTaskJetShapeBase::MatchEmbeddedConstituentWithParticleLevel() {
   if(!fTreeEmb) return 0;
   
   Int_t entry = fTreeEmb->GetEntry(fThisEntry);
   return entry;
 
}


//________________________________________________________________________
void AliAnalysisTaskJetShapeBase::SetTree(TTree *tree) {
   if(!tree){
      AliError("Null tree");
      return;
   }
   fFromTree = kTRUE;
   fTreeEmb = (TTree*)tree->Clone(Form("%sCpShC", tree->GetName()));
   AliInfo(Form("Input tree set %d (%p -> %p)", fTreeEmb->GetNbranches(), tree, fTreeEmb));
   //fTreeEmb->SetDirectory(0x0);
   
   return;
}

//________________________________________________________________________

void AliAnalysisTaskJetShapeBase::SetTreeFromFile(TString filename, TString treename){
   
   if(filename.Contains("alien")) {
      TGrid::Connect("alien://");
   }
   TFile *f = TFile::Open(filename);
   if(!f->IsOpen()){
      Printf("File %s not found, cannot SetTree", filename.Data());
      return;
   }
   
   TTree *tree = dynamic_cast<TTree*>(f->Get(treename));
   if(!tree){
      Printf("Tree %s not found!!!", treename.Data());
      f->ls();
      return;
   }
   SetTree(tree);
   
   //f->Close();
   //delete f;

   return;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeBase::RetrieveEventObjects() {
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
void AliAnalysisTaskJetShapeBase::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

