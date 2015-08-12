//
// Do subtraction for jet shapes using constituents arXiv:1403.3108
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

#include "AliAnalysisTaskJetShapeConst.h"

ClassImp(AliAnalysisTaskJetShapeConst)

//________________________________________________________________________
AliAnalysisTaskJetShapeConst::AliAnalysisTaskJetShapeConst() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetShapeConst", kTRUE),
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
  fhnDeltaMass(0),
  fhRjetTrvspTj(0x0),
  fhNJetsSelEv(0x0),
  fhJetEtaPhi(0x0),
  fhpTTracksJet1(0x0),
  fhpTTracksJetO(0x0),
  fhpTTracksCont(0x0)
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
  fhnDeltaMass             = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubMatch[i]             = 0;
    fh2MSubPtRawAll[i]          = 0;
    fh3MSubPtRawDRMatch[i]      = 0;
    fh3MSubPtTrueLeadPt[i]      = 0;
    fh3MTruePtTrueLeadPt[i]     = 0;
    fh3PtTrueDeltaMLeadPt[i]    = 0;
    fh3PtTrueDeltaMRelLeadPt[i] = 0;
    fhnMassResponse[i]          = 0;
    fhnDeltaMass[i]             = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  if(fCreateTree) DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeConst::AliAnalysisTaskJetShapeConst(const char *name) : 
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
  fhnDeltaMass(0),
  fhRjetTrvspTj(0x0),
  fhNJetsSelEv(0x0),
  fhJetEtaPhi(0x0),
  fhpTTracksJet1(0x0),
  fhpTTracksJetO(0x0),
  fhpTTracksCont(0x0)
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
  fhnDeltaMass             = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubMatch[i]             = 0;
    fh2MSubPtRawAll[i]          = 0;
    fh3MSubPtRawDRMatch[i]      = 0;
    fh3MSubPtTrueLeadPt[i]      = 0;
    fh3MTruePtTrueLeadPt[i]     = 0;
    fh3PtTrueDeltaMLeadPt[i]    = 0;
    fh3PtTrueDeltaMRelLeadPt[i] = 0;
    fhnMassResponse[i]          = 0;
    fhnDeltaMass[i]             = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  if(fCreateTree) DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetShapeConst::~AliAnalysisTaskJetShapeConst()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetShapeConst::UserCreateOutputObjects()
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
  const Int_t nBinsSparse0 = 5;
  //Mass sub;Mass true;#it{p}_{T,sub};#it{p}_{T,true};#it{p}_{T,lead trk}
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt,nBinsPtLead};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt, minPtLead};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt, maxPtLead};

  const Int_t nBinsSparse1 = 7;
  // #it{M}_{det,Const} - #it{M}_{part}; #it{p}_{T,det,Const} - #it{p}_{T,part}; #it{M}_{det,Const};  #it{M}_{part}; #it{p}_{T,det,Const}; #it{p}_{T,part}; #it{p}_{T,det,A}
  const Int_t nBins1[nBinsSparse1] = {nBinsDM,nBinsDpT,nBinsM,nBinsM,nBinsPt,nBinsPt,nBinsPt};
  const Double_t xmin1[nBinsSparse1]  = { minDM, minDpT, minM, minM, minPt, minPt, minPt};
  const Double_t xmax1[nBinsSparse1]  = { maxDM, maxDpT, maxM, maxM, maxPt, maxPt, maxPt};

  const Int_t nBinsSparse2 = 8;
  //#it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}
  const Int_t nBins2[nBinsSparse2] = {nBinsDM, nBinsDpT, nBinsM, nBinsM, nBinsPt, nBinsPt, nBinsRho, nBinsRhom};
  const Double_t xmin2[nBinsSparse2]  = {minDM, minDpT, minM, minM, minPt, minPt, minRho, minRhom};
  const Double_t xmax2[nBinsSparse2]  = {maxDM, maxDpT, maxM, maxM, maxPt, maxPt, maxRho, maxRhom};

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
    
    histName = Form("fhnDeltaMass_%d", i);
    histTitle = Form("%s; #it{M}_{det,Const} - #it{M}_{part}; #it{p}_{T,det,Const} - #it{p}_{T,part}; #it{M}_{det,Const};  #it{M}_{part}; #it{p}_{T,det,Const}; #it{p}_{T,part}; #it{p}_{T,det,A}",histName.Data());
    Printf("Nuber of bins %d - write first %d, %f, %f , building %s", nBinsSparse1, nBins1[0], xmin1[0], xmax1[0], histName.Data());
    fhnDeltaMass[i] = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse1,nBins1,xmin1,xmax1);
    fOutput->Add(fhnDeltaMass[i]);

  }

  //Chiara's histograms: rho and rhom correlation with pT and mass at reco level with no subtraction
  histName = "fhnDeltaMassAndBkgInfo";
  histTitle = Form("%s; #it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}",histName.Data()); // #it{M}_{unsub} is also deltaM unsub when M_part is zero
  
  fhnDeltaMassAndBkgInfo = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse2,nBins2,xmin2,xmax2);
  fOutput->Add(fhnDeltaMassAndBkgInfo);
  
  if(fOverlap){
     fhRjetTrvspTj = new TH2F("fhRjetTrvspTj", ";R(jet, track);p_{T,jet}", 100, 0., 10., nBinsPt, minPt, maxPt);
     fOutput->Add(fhRjetTrvspTj);
     fhNJetsSelEv = new TH1F("fhNJetsSelEv", "N of jets selected; #it{N}_{jets}/ev;Entries", 20., 0.,19);
     fOutput->Add(fhNJetsSelEv);
     
     fhJetEtaPhi = new TH2F("fhJetEtaPhi", "#eta - #varphi distribution of selected jets; #eta; #varphi", 24., -0.6, 0.6, 50, 0., 2*TMath::Pi());
     fOutput->Add(fhJetEtaPhi);
     
     fhpTTracksJetO = new TH1F("hTrackpTO", "Track pT (signal jet); p_{T}", 500,0.,50.);
     fOutput->Add(fhpTTracksJetO);

  }
  fhpTTracksJet1 = new TH1F("hTrackpT1", "Track pT ; p_{T}", 500,0.,50.);
  fOutput->Add(fhpTTracksJet1);
  fhpTTracksCont = new TH1F(Form("fhpTTracksCont"), "Track pT (container) ; p_{T}", 500,0.,50.);
  fOutput->Add(fhpTTracksCont);

  fhptjetSMinusSingleTrack = new TH1F("fhptjetSMinusSingleTrack", "Subtraction of single track #it{p}_{T}; |#it{p}_{T, jet} - #it{p}_{T, Emb Track}|;Entries", 500,-10.,110.);
  fOutput->Add(fhptjetSMinusSingleTrack);
  
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
    fTreeJetBkg = new TTree("fTreeJetSubConst", "fTreeJetSubConst");
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
Bool_t AliAnalysisTaskJetShapeConst::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeConst::FillHistograms()
{
  // Fill histograms.

  AliEmcalJet* jet1  = NULL; //AA jet
  AliEmcalJet *jet2  = NULL; //Embedded Pythia jet
  AliEmcalJet *jetR  = NULL; //true jet for response matrix
  AliEmcalJet *jetS = NULL;  //subtracted jet
  AliEmcalJet *jetO  = NULL; //hard-ish jet to avoid overlap of single track with
  AliVParticle *vpe  = NULL; //embedded particle
  
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  AliJetContainer *jetContS = GetJetContainer(fContainerSub);
  AliJetContainer *jetContO = GetJetContainer(fContainerOverlap);

  AliParticleContainer *trackCont = GetParticleContainer(0); //track used for the jet finder jet1
  for(Int_t itr=0; itr < trackCont->GetNAcceptedParticles(); itr++){
     fhpTTracksCont->Fill(trackCont->GetParticle(itr)->Pt());
  
  }
  //Get leading jet in Pb-Pb event without embedded objects
  AliJetContainer *jetContNoEmb = GetJetContainer(fContainerNoEmb);
  AliEmcalJet *jetL = NULL;
  if(jetContNoEmb) jetL = jetContNoEmb->GetLeadingJet("rho");
  
  if(fOverlap && !jetContO){
     Printf("Jet Container %d not found, return", fContainerOverlap);
     return kFALSE;
  }

  //rho
  AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(jetContS->GetRhoName()));
  fRho = 0;
  if (!rhoParam) {
     AliError(Form("%s: Could not retrieve rho %s (some histograms will be filled with zero)!", GetName(), jetContS->GetRhoName().Data()));
     
  } else fRho = rhoParam->GetVal();
  //rhom
  AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(jetContS->GetRhoMassName()));
  fRhoM = 0;
  if (!rhomParam) {
     AliError(Form("%s: Could not retrieve rho_m %s (some histograms will be filled with zero)!", GetName(), jetContS->GetRhoMassName().Data()));
     
  } else fRhoM = rhomParam->GetVal();
  
  AliDebug(11,Form("NJets  Incl: %d  Csub: %d",jetCont->GetNJets(),jetContS->GetNJets()));
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      jet2  = NULL;
      jetS  = NULL;
      
      if(jet1->GetTagStatus()<1 || !jet1->GetTaggedJet())
        continue;

     //print constituents of different jet containers
     //jet1
     
     for(Int_t i=0; i<jet1->GetNumberOfTracks(); i++) {
     	AliVParticle *vp = static_cast<AliVParticle*>(jet1->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
     	//    if (vp->TestBits(TObject::kBitMask) != (Int_t)(TObject::kBitMask) ) continue;
     	//Int_t lab = TMath::Abs(vp->GetLabel());
     	if(vp) fhpTTracksJet1 -> Fill(vp->Pt());
     }

      //Get constituent subtacted version of jet
      Int_t ifound = 0;
      Int_t ilab = -1;
      for(Int_t i = 0; i<jetContS->GetNJets(); i++) {
	//if(ifound==1) continue;
	jetS = jetContS->GetJet(i);
        if(!jetS) continue;
	if(jetS->GetLabel()==jet1->GetLabel()) {
	  ifound++;
	  if(ifound==1) ilab = i;
	}
      }
      if(ifound>1) AliDebug(2,Form("Found %d partners",ifound));
      if(ifound==0) jetS = 0x0;
      else jetS = jetContS->GetJet(ilab);
      if(!jetS) continue;

      Double_t mjetS = jetS->M();
      Double_t mUnsubjet1 = jet1->M();
      Double_t ptjet1 = jet1->Pt()-jetCont->GetRhoVal()*jet1->Area();
      Double_t ptjetS = jetS->Pt();
      Double_t ptUnsubjet1 = jet1->Pt();
      Double_t var = mjetS;
      Double_t ptjetSMinusEmbTrpt = ptjetS;
      if(fJetMassVarType==kRatMPt) {
      	 if(ptjetS>0. || ptjetS<0.) var = mjetS/ptjetS;
      	 else var = -999.;
      }
      
      //Fill histograms for all AA jets
      fh2MSubPtRawAll[fCentBin]->Fill(var,ptjetS);
      
      Double_t fraction = 0.;
      fMatch = 0;
      fJet2Vec->SetPtEtaPhiM(0.,0.,0.,0.);
      if(fSingleTrackEmb) {
      	 vpe = GetEmbeddedConstituent(jet1);
      	 if(vpe){
      	    Bool_t reject = kFALSE;
      	    ptjetSMinusEmbTrpt =- vpe->Pt();
      	    if(fOverlap){
      	       Int_t Njets = jetContO->GetNAcceptedJets();
      	       fhNJetsSelEv->Fill(Njets);
      	       jetContO->ResetCurrentID();
      	       while(jetO = jetContO->GetNextAcceptJet()){
      	       	  //print constituents of different jet containers
       	     	  //jetO
       	     	  //Printf("N particle %d",jetO->GetNumberOfTracks());
       	     	  for(Int_t i=0; i<jetO->GetNumberOfTracks(); i++) {
       	     	     AliVParticle* vp = static_cast<AliVParticle*>(jetO->TrackAt(i, jetContO->GetParticleContainer()->GetArray()));
       	     	     //    if (vp->TestBits(TObject::kBitMask) != (Int_t)(TObject::kBitMask) ) continue;
       	     	     //Int_t lab = TMath::Abs(vp->GetLabel());
       	     	     if(vp) fhpTTracksJetO -> Fill(vp->Pt());
       	     	  }
      	       	  Double_t deltaR = jetO->DeltaR(vpe);
      	       	  fhRjetTrvspTj->Fill(deltaR, jetO->Pt());
      	       	  fhJetEtaPhi->Fill(jetO->Eta(), jetO->Phi());
      	       	  if(deltaR < fRadius) {
      	       	     reject = kTRUE;
      	       	     break;
      	       	  }
      	       }
      	    }
      	    if(!reject) {
      	       fJet2Vec->SetPxPyPzE(vpe->Px(),vpe->Py(),vpe->Pz(),vpe->E());
      	       fMatch = 1;
      	    }
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
      
      //      if(fMatch==1 && jet2->GetTagStatus()>0) jet2T = jet2->GetTaggedJet();
      
      //Fill histograms for matched jets
      fh2MSubMatch[fCentBin]->Fill(var,fMatch);
      if(fMatch==1) {
      	 Double_t drToLJ = -1.;
      	 if(jetL) drToLJ = jet1->DeltaR(jetL);
      	 if(fSingleTrackEmb && vpe)
      	    drToLJ = jet1->DeltaR(vpe);
      	 fh3MSubPtRawDRMatch[fCentBin]->Fill(var,ptjet1,drToLJ);
      	 Double_t var2 = 0.;
      	 Double_t mJetR = 0.;
      	 Double_t ptJetR = 0.;
      	 if(jetR) {
      	    mJetR  = jetR->M();
      	    var2 = jetR->M();
      	    ptJetR = jetR->Pt();
      	 }
      	 if(fSingleTrackEmb && vpe) {
      	    mJetR  = vpe->M();
      	    var2   = vpe->M();
      	    ptJetR = vpe->Pt();
      	 }
      	 
      	 if(fJetMassVarType==kRatMPt) {
      	    if(ptJetR>0. || ptJetR<0.) var2 /= ptJetR;
      	 }
      	 fh3MSubPtTrueLeadPt[fCentBin]->Fill(var,ptJetR,jet1->MaxTrackPt());
      	 fh3MTruePtTrueLeadPt[fCentBin]->Fill(var2,ptJetR,jet1->MaxTrackPt());
      	 fh3PtTrueDeltaMLeadPt[fCentBin]->Fill(ptJetR,var-var2,jet1->MaxTrackPt());
      	 if(var2>0.) fh3PtTrueDeltaMRelLeadPt[fCentBin]->Fill(ptJetR,(var-var2)/var2,jet1->MaxTrackPt());
      	 //M sub;M true;#it{p}_{T,sub};#it{p}_{T,true};#it{p}_{T,lead trk}
      	 Double_t varsp[5] = {var,var2,ptjetS,ptJetR,jetS->MaxTrackPt()};
      	 fhnMassResponse[fCentBin]->Fill(varsp);
      	 
      	 Double_t varsp1[7];
      	 //#it{M}_{det,Const} - #it{M}_{part}; #it{p}_{T,det,Const} - #it{p}_{T,part}; #it{M}_{det,Const};  #it{M}_{part}; #it{p}_{T,det,Const}; #it{p}_{T,part}; #it{p}_{T,det,A}
      	 varsp1[0] = var-var2;
      	 varsp1[1] = ptjetS-ptJetR;
      	 varsp1[2] = var;
      	 varsp1[3] = var2;
      	 varsp1[4] = ptjetS;
      	 varsp1[5] = ptJetR;
      	 varsp1[6] = ptjet1;
      	 
      	 fhnDeltaMass[fCentBin]->Fill(varsp1);
      	 
      	 //#it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}
      	 Double_t varsp2[8] = {var-var2, ptjetS-ptJetR, var, mUnsubjet1, ptjetS, ptUnsubjet1, fRho, fRhoM};
      	 fhnDeltaMassAndBkgInfo->Fill(varsp2);
      	 
      	 fhptjetSMinusSingleTrack->Fill(TMath::Abs(ptjetSMinusEmbTrpt));

      }
      
      if(fCreateTree) {      
      	 fJet1Vec->SetPxPyPzE(jet1->Px(),jet1->Py(),jet1->Pz(),jet1->E());
      	 if(jetS->Pt()>0.) fJetSubVec->SetPtEtaPhiM(jetS->Pt(),jetS->Eta(),jetS->Phi(),jetS->M());
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
AliVParticle* AliAnalysisTaskJetShapeConst::GetEmbeddedConstituent(AliEmcalJet *jet) {

  AliParticleContainer *partContEmb = GetParticleContainer(); //the only particle container given is the one with embedded track(s)

  AliVParticle *vpe = 0x0; //embedded particle
  for(Int_t ip = partContEmb->GetNParticles()-1; ip>-1; ip--){
      AliVParticle *vp = partContEmb->GetParticle(ip);
      if(!vp){
      	 AliDebug(2, Form("Particle %d not found", ip));
      	 continue;
      }
      Int_t lab = TMath::Abs(vp->GetLabel());
      if (lab < fMinLabelEmb || lab > fMaxLabelEmb)
      	 continue;
      if(!vpe) vpe = vp;
      else if(vp->Pt()>vpe->Pt()) vpe =vp;
  }
  
  Double_t deltaR = 99;
  if(vpe && jet) deltaR = jet->DeltaR(vpe);
  if(deltaR < 0.2) return vpe;
  else return 0x0;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeConst::RetrieveEventObjects() {
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
void AliAnalysisTaskJetShapeConst::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

