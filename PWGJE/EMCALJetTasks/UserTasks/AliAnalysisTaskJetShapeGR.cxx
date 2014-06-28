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
  fh2GRSubMatch(0x0),
  fh2GRSubPtRawAll(0x0),
  fh2GRSubPtRawMatch(0x0),
  fh2GRSubPtTrue(0x0),
  fh2GRTruePtTrue(0x0),
  fh2PtTrueDeltaGR(0x0),
  fh2PtTrueDeltaGRRel(0x0),
  fhnGRResponse(0x0)
{
  // Default constructor.

  fh2GRSubMatch        = new TH2F*[fNcentBins];
  fh2GRSubPtRawAll     = new TH2F*[fNcentBins];
  fh2GRSubPtRawMatch   = new TH2F*[fNcentBins];
  fh2GRSubPtTrue       = new TH2F*[fNcentBins];
  fh2GRTruePtTrue      = new TH2F*[fNcentBins];
  fh2PtTrueDeltaGR     = new TH2F*[fNcentBins];
  fh2PtTrueDeltaGRRel  = new TH2F*[fNcentBins];
  fhnGRResponse     = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2GRSubMatch[i]        = 0;
    fh2GRSubPtRawAll[i]     = 0;
    fh2GRSubPtRawMatch[i]   = 0;
    fh2GRSubPtTrue[i]       = 0;
    fh2GRTruePtTrue[i]      = 0;
    fh2PtTrueDeltaGR[i]     = 0;
    fh2PtTrueDeltaGRRel[i]  = 0;
    fhnGRResponse[i]     = 0;
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
  fh2GRSubMatch(0x0),
  fh2GRSubPtRawAll(0x0),
  fh2GRSubPtRawMatch(0x0),
  fh2GRSubPtTrue(0x0),
  fh2GRTruePtTrue(0x0),
  fh2PtTrueDeltaGR(0x0),
  fh2PtTrueDeltaGRRel(0x0),
  fhnGRResponse(0x0)
{
  // Standard constructor.

  fh2GRSubMatch        = new TH2F*[fNcentBins];
  fh2GRSubPtRawAll     = new TH2F*[fNcentBins];
  fh2GRSubPtRawMatch   = new TH2F*[fNcentBins];
  fh2GRSubPtTrue       = new TH2F*[fNcentBins];
  fh2GRTruePtTrue      = new TH2F*[fNcentBins];
  fh2PtTrueDeltaGR     = new TH2F*[fNcentBins];
  fh2PtTrueDeltaGRRel  = new TH2F*[fNcentBins];
  fhnGRResponse     = new THnSparse*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2GRSubMatch[i]        = 0;
    fh2GRSubPtRawAll[i]     = 0;
    fh2GRSubPtRawMatch[i]   = 0;
    fh2GRSubPtTrue[i]       = 0;
    fh2GRTruePtTrue[i]      = 0;
    fh2PtTrueDeltaGR[i]     = 0;
    fh2PtTrueDeltaGRRel[i]  = 0;
    fhnGRResponse[i]     = 0;
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

  //Binning for THnSparse
  const Int_t nBinsSparse0 = 4;
  const Int_t nBins0[nBinsSparse0] = {nBinsM,nBinsM,nBinsPt,nBinsPt};
  const Double_t xmin0[nBinsSparse0]  = { minM, minM, minPt, minPt};
  const Double_t xmax0[nBinsSparse0]  = { maxM, maxM, maxPt, maxPt};

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = Form("fh2GRSubMatch_%d",i);
    histTitle = Form("fh2GRSubMatch_%d;#it{G(R)}_{sub};match",i);
    fh2GRSubMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,2,-0.5,1.5);
    fOutput->Add(fh2GRSubMatch[i]);

    histName = Form("fh2GRSubPtRawAll_%d",i);
    histTitle = Form("fh2GRSubPtRawAll_%d;#it{G(R)}_{sub};#it{p}_{T}",i);
    fh2GRSubPtRawAll[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2GRSubPtRawAll[i]);

    histName = Form("fh2GRSubPtRawMatch_%d",i);
    histTitle = Form("fh2GRSubPtRawMatch_%d;#it{G(R)}_{sub};#it{p}_{T}",i);
    fh2GRSubPtRawMatch[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2GRSubPtRawMatch[i]);

    histName = Form("fh2GRSubPtTrue_%d",i);
    histTitle = Form("fh2GRSubPtTrue_%d;#it{G(R)}_{sub};#it{p}_{T}",i);
    fh2GRSubPtTrue[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2GRSubPtTrue[i]);

    histName = Form("fh2GRTruePtTrue_%d",i);
    histTitle = Form("fh2GRTruePtTrue_%d;#it{G(R)}_{sub};#it{p}_{T}",i);
    fh2GRTruePtTrue[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2GRTruePtTrue[i]);

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

  AliEmcalJet* jet1 = NULL;
  AliEmcalJet *jet2 = NULL;
  AliEmcalJet *jetS = NULL;
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  AliJetContainer *jetContS = GetJetContainer(fContainerSub);
  AliDebug(11,Form("NJets  Incl: %d  Csub: %d",jetCont->GetNJets(),jetContS->GetNJets()));
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {

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
      fh2GRSubPtRawAll[fCentBin]->Fill(jetS->M(),jet1->Pt()-jetCont->GetRhoVal()*jet1->Area());

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

      //Fill histograms for matched jets
      fh2GRSubMatch[fCentBin]->Fill(jetS->M(),fMatch);
      if(fMatch==1) {
	fh2GRSubPtRawMatch[fCentBin]->Fill(jetS->M(),jet1->Pt()-jetCont->GetRhoVal()*jet1->Area());
	if(jet2) {
	  fh2GRSubPtTrue[fCentBin]->Fill(jetS->M(),jet2->Pt());
	  Double_t gr = CalcGR(jet2,fContainerTrue);
	  Printf("G(R): %f",gr);
	  fh2GRTruePtTrue[fCentBin]->Fill(gr,jet2->Pt());
	  fh2PtTrueDeltaGR[fCentBin]->Fill(jet2->Pt(),jetS->M()-jet2->M());
	  if(jet2->M()>0.) fh2PtTrueDeltaGRRel[fCentBin]->Fill(jet2->Pt(),(jetS->M()-jet2->M())/jet2->M());
	  Double_t var[4] = {0.,0.,jet1->Pt()-jetCont->GetRhoVal()*jet1->Area(),jet2->Pt()};
	  fhnGRResponse[fCentBin]->Fill(var);
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
Double_t AliAnalysisTaskJetShapeGR::CalcGR(AliEmcalJet *jet, Int_t ic) {
  //Calculate G(R)
  AliJetContainer *jetCont = GetJetContainer(ic); 
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  Double_t gR = 0.;
  Double_t wr = 0.04;
  const Int_t nr = TMath::CeilNint(jetCont->GetJetRadius()/wr);
  Double_t grArr[nr];
  for(Int_t i = 0; i<nr; i++) {
    grArr[i] = 0.;
    Printf("bin up edge %d=%f",i,wr+i*wr);
  }
  Printf("jet pt: %f  nconst: %d",jet->Pt(),jet->GetNumberOfTracks());
  for(Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
    for(Int_t j=i; j<jet->GetNumberOfTracks(); j++) {
      vp2 = static_cast<AliVParticle*>(jet->TrackAt(j, jetCont->GetParticleContainer()->GetArray()));
      Double_t dr2 = (vp1->Eta()-vp2->Eta())*(vp1->Eta()-vp2->Eta()) + (vp1->Phi()-vp2->Phi())*(vp1->Phi()-vp2->Phi());
      Int_t bin = TMath::FloorNint(TMath::Sqrt(dr2)/wr);
      Double_t gr = vp1->Pt()*vp2->Pt()*dr2;
      if(bin<nr) grArr[bin]+=gr;
      
      if(TMath::Sqrt(dr2)<jetCont->GetJetRadius())
	gR += gr;
    }
  }
  for(Int_t i = 0; i<nr; i++) {
    Printf("grArr[%d]=%f",i,grArr[i]);
  }
  return gR;
}

//________________________________________________________________________
Double_t AliAnalysisTaskJetShapeGR::CalcDeltaGR(AliEmcalJet *jet, Int_t ic) {
  //Calculate G(R)
  AliJetContainer *jetCont = GetJetContainer(ic); 
  AliVParticle *vp1 = 0x0;
  AliVParticle *vp2 = 0x0;
  Double_t gR = 0.;
  Double_t wr = 0.04;
  const Int_t nr = TMath::CeilNint(jetCont->GetJetRadius()/wr);
  Double_t grArr[nr];
  for(Int_t i = 0; i<nr; i++) {
    grArr[i] = 0.;
    Printf("bin up edge %d=%f",i,wr+i*wr);
  }
  Printf("jet pt: %f  nconst: %d",jet->Pt(),jet->GetNumberOfTracks());
  for(Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
    vp1 = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
    for(Int_t j=i; j<jet->GetNumberOfTracks(); j++) {
      vp2 = static_cast<AliVParticle*>(jet->TrackAt(j, jetCont->GetParticleContainer()->GetArray()));
      Double_t dr2 = (vp1->Eta()-vp2->Eta())*(vp1->Eta()-vp2->Eta()) + (vp1->Phi()-vp2->Phi())*(vp1->Phi()-vp2->Phi());
      Int_t bin = TMath::FloorNint(TMath::Sqrt(dr2)/wr);
      Double_t gr = vp1->Pt()*vp2->Pt()*dr2;
      if(bin<nr) grArr[bin]+=gr;
      
      if(TMath::Sqrt(dr2)<jetCont->GetJetRadius())
	gR += gr;
    }
  }
  for(Int_t i = 0; i<nr; i++) {
    Printf("grArr[%d]=%f",i,grArr[i]);
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

