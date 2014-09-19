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
  fh2MSubMatch(0x0),
  fh2MSubPtRawAll(0x0),
  fh3MSubPtRawDRMatch(0x0),
  fh3MSubPtTrueDR(0x0),
  fh3MTruePtTrueDR(0x0),
  fh3PtTrueDeltaMDR(0x0),
  fh3PtTrueDeltaMRelDR(0x0),
  fhnMassResponse(0x0)
{
  // Default constructor.

  fh2MSubMatch         = new TH2F*[fNcentBins];
  fh2MSubPtRawAll      = new TH2F*[fNcentBins];
  fh3MSubPtRawDRMatch  = new TH3F*[fNcentBins];
  fh3MSubPtTrueDR      = new TH3F*[fNcentBins];
  fh3MTruePtTrueDR     = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMDR    = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMRelDR = new TH3F*[fNcentBins];
  fhnMassResponse      = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubMatch[i]         = 0;
    fh2MSubPtRawAll[i]      = 0;
    fh3MSubPtRawDRMatch[i]  = 0;
    fh3MSubPtTrueDR[i]      = 0;
    fh3MTruePtTrueDR[i]     = 0;
    fh3PtTrueDeltaMDR[i]    = 0;
    fh3PtTrueDeltaMRelDR[i] = 0;
    fhnMassResponse[i]      = 0;
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
  fh2MSubMatch(0x0),
  fh2MSubPtRawAll(0x0),
  fh3MSubPtRawDRMatch(0x0),
  fh3MSubPtTrueDR(0x0),
  fh3MTruePtTrueDR(0x0),
  fh3PtTrueDeltaMDR(0x0),
  fh3PtTrueDeltaMRelDR(0x0),
  fhnMassResponse(0x0)
{
  // Standard constructor.

  fh2MSubMatch         = new TH2F*[fNcentBins];
  fh2MSubPtRawAll      = new TH2F*[fNcentBins];
  fh3MSubPtRawDRMatch  = new TH3F*[fNcentBins];
  fh3MSubPtTrueDR      = new TH3F*[fNcentBins];
  fh3MTruePtTrueDR     = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMDR    = new TH3F*[fNcentBins];
  fh3PtTrueDeltaMRelDR = new TH3F*[fNcentBins];
  fhnMassResponse      = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubMatch[i]         = 0;
    fh2MSubPtRawAll[i]      = 0;
    fh3MSubPtRawDRMatch[i]  = 0;
    fh3MSubPtTrueDR[i]      = 0;
    fh3MTruePtTrueDR[i]     = 0;
    fh3PtTrueDeltaMDR[i]    = 0;
    fh3PtTrueDeltaMRelDR[i] = 0;
    fhnMassResponse[i]      = 0;
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

  const Int_t nBinsM  = 100;
  const Double_t minM = -20.;
  const Double_t maxM = 80.;

  // const Int_t nBinsMT  = 50;
  // const Double_t minMT = 0.;
  // const Double_t maxMT = 50.;

  const Int_t nBinsDRToLJ  = 20; //distance to leading jet in Pb-Pb only event
  const Double_t minDRToLJ = 0.;
  const Double_t maxDRToLJ = 1.;

  //Binning for THnSparse
  const Int_t nBinsSparse0 = 5;
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt,nBinsDRToLJ};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt, minDRToLJ};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt, maxDRToLJ};

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = Form("fh2MSubMatch_%d",i);
    histTitle = Form("fh2MSubMatch_%d;#it{M}_{sub};match",i);
    fh2MSubMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,2,-0.5,1.5);
    fOutput->Add(fh2MSubMatch[i]);

    histName = Form("fh2MSubPtRawAll_%d",i);
    histTitle = Form("fh2MSubPtRawAll_%d;#it{M}_{sub};#it{p}_{T}",i);
    fh2MSubPtRawAll[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2MSubPtRawAll[i]);

    histName = Form("fh3MSubPtRawDRMatch_%d",i);
    histTitle = Form("fh3MSubPtRawDRMatch_%d;#it{M}_{sub};#it{p}_{T}",i);
    fh3MSubPtRawDRMatch[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt,nBinsDRToLJ,minDRToLJ,maxDRToLJ);
    fOutput->Add(fh3MSubPtRawDRMatch[i]);

    histName = Form("fh3MSubPtTrueDR_%d",i);
    histTitle = Form("fh3MSubPtTrueDR_%d;#it{M}_{sub};#it{p}_{T}",i);
    fh3MSubPtTrueDR[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt,nBinsDRToLJ,minDRToLJ,maxDRToLJ);
    fOutput->Add(fh3MSubPtTrueDR[i]);

    histName = Form("fh3MTruePtTrueDR_%d",i);
    histTitle = Form("fh3MTruePtTrueDR_%d;#it{M}_{sub};#it{p}_{T}",i);
    fh3MTruePtTrueDR[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt,nBinsDRToLJ,minDRToLJ,maxDRToLJ);
    fOutput->Add(fh3MTruePtTrueDR[i]);

    histName = Form("fh3PtTrueDeltaMDR_%d",i);
    histTitle = Form("fh3PtTrueDeltaMDR_%d;#it{p}_{T,true};#it{M}_{sub}-#it{M}_{true}",i);
    fh3PtTrueDeltaMDR[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,100,-50.,50.,nBinsDRToLJ,minDRToLJ,maxDRToLJ);
    fOutput->Add(fh3PtTrueDeltaMDR[i]);

    histName = Form("fh3PtTrueDeltaMRelDR_%d",i);
    histTitle = Form("fh3PtTrueDeltaMRelDR_%d;#it{p}_{T,true};(#it{M}_{sub}-#it{M}_{true})/#it{M}_{true}",i);
    fh3PtTrueDeltaMRelDR[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,200,-1.,1.,nBinsDRToLJ,minDRToLJ,maxDRToLJ);
    fOutput->Add(fh3PtTrueDeltaMRelDR[i]);

    histName = Form("fhnMassResponse_%d",i);
    histTitle = Form("fhnMassResponse_%d;#it{M}_{sub};#it{M}_{true};#it{p}_{T,sub};#it{p}_{T,true};#it{M}_{sub}^{tagged}",i);
    fhnMassResponse[i] = new THnSparseF(histName.Data(),histTitle.Data(),nBinsSparse0,nBins0,xmin0,xmax0);
    fOutput->Add(fhnMassResponse[i]);
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
  //  AliEmcalJet *jet1T = NULL; //tagged AA jet
  //  AliEmcalJet *jet2T = NULL; //tagged Pythia jet
  AliEmcalJet *jetS = NULL;  //subtracted jet
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  AliJetContainer *jetContS = GetJetContainer(fContainerSub);

  //Get leading jet in Pb-Pb event without embedded objects
  AliJetContainer *jetContNoEmb = GetJetContainer(fContainerNoEmb);
  AliEmcalJet *jetL = NULL;
  if(jetContNoEmb) jetL = jetContNoEmb->GetLeadingJet("rho");

  AliDebug(11,Form("NJets  Incl: %d  Csub: %d",jetCont->GetNJets(),jetContS->GetNJets()));
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      jet2  = NULL;
      //      jet1T = NULL;
      //   jet2T = NULL;
      jetS  = NULL;
      //      if(jet1->GetTagStatus()>0) jet1T = jet1->GetTaggedJet();
 
      //Get constituent subtacted version of jet
      Int_t ifound = 0;
      Int_t ilab = -1;
      for(Int_t i = 0; i<jetContS->GetNJets(); i++) {
	//if(ifound==1) continue;
	jetS = jetContS->GetJet(i);
	if(jetS->GetLabel()==jet1->GetLabel()) { // && jetS->Pt()>0.) {
	  ifound++;
	  if(ifound==1) ilab = i;
	}
      }
      if(ifound>1) AliDebug(2,Form("Found %d partners",ifound));
      if(ifound==0) jetS = 0x0;
      else jetS = jetContS->GetJet(ilab);

      //Fill histograms for all AA jets
      fh2MSubPtRawAll[fCentBin]->Fill(jetS->M(),jet1->Pt()-jetCont->GetRhoVal()*jet1->Area());

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
	  } else
	    fMatch = 0;
	}
      }

      //      if(fMatch==1 && jet2->GetTagStatus()>0) jet2T = jet2->GetTaggedJet();

      //Fill histograms for matched jets
      fh2MSubMatch[fCentBin]->Fill(jetS->M(),fMatch);
      if(fMatch==1) {
	Double_t drToLJ = -1.;
	if(jetL) drToLJ = jet1->DeltaR(jetL);
	fh3MSubPtRawDRMatch[fCentBin]->Fill(jetS->M(),jet1->Pt()-jetCont->GetRhoVal()*jet1->Area(),drToLJ);
	if(jet2) {
	  fh3MSubPtTrueDR[fCentBin]->Fill(jetS->M(),jet2->Pt(),drToLJ);
	  fh3MTruePtTrueDR[fCentBin]->Fill(jet2->M(),jet2->Pt(),drToLJ);
	  fh3PtTrueDeltaMDR[fCentBin]->Fill(jet2->Pt(),jetS->M()-jet2->M(),drToLJ);
	  if(jet2->M()>0.) fh3PtTrueDeltaMRelDR[fCentBin]->Fill(jet2->Pt(),(jetS->M()-jet2->M())/jet2->M(),drToLJ);
	  //	  Double_t mJet1Tagged = -1.;
	  //	  if(jet1T) mJet1Tagged = jet1T->M();
	  Double_t var[5] = {jetS->M(),jet2->M(),jet1->Pt()-jetCont->GetRhoVal()*jet1->Area(),jet2->Pt(),drToLJ};
	  fhnMassResponse[fCentBin]->Fill(var);
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
AliVParticle* AliAnalysisTaskJetShapeConst::GetEmbeddedConstituent(AliEmcalJet *jet) {

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

