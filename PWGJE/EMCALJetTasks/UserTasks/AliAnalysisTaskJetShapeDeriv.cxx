#include <TRandom3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TLorentzVector.h>

#include <AliParticleContainer.h>
#include <AliJetContainer.h>
#include <AliRhoParameter.h>
#include <AliEmcalJet.h>

#include "AliAnalysisTaskJetShapeDeriv.h"

ClassImp(AliAnalysisTaskJetShapeDeriv)

//________________________________________________________________________
AliAnalysisTaskJetShapeDeriv::AliAnalysisTaskJetShapeDeriv() : 
  AliAnalysisTaskJetShapeBase("AliAnalysisTaskJetShapeDeriv"),
  fM1st(0),
  fM2nd(0),
  fDeriv1st(0),
  fDeriv2nd(0),
  fPartialExclusion(0),
  fh2MSubPtSubAll(0x0),
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

  fh2MSubPtSubAll          = new TH2F*[fNcentBins];
  fh2PtTrueSubFacV1        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV1         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV1        = new TH2F*[fNcentBins];
  fh2NConstSubFacV1        = new TH2F*[fNcentBins];
  fh2PtTrueSubFacV2        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV2         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV2        = new TH2F*[fNcentBins];
  fh2NConstSubFacV2        = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubPtSubAll[i]          = 0;
    fh2PtTrueSubFacV1[i]        = 0;
    fh2PtRawSubFacV1[i]         = 0;
    fh2PtCorrSubFacV1[i]        = 0;
    fh2NConstSubFacV1[i]        = 0;
    fh2PtTrueSubFacV2[i]        = 0;
    fh2PtRawSubFacV2[i]         = 0;
    fh2PtCorrSubFacV2[i]        = 0;
    fh2NConstSubFacV2[i]        = 0;
  }

}

//________________________________________________________________________
AliAnalysisTaskJetShapeDeriv::AliAnalysisTaskJetShapeDeriv(const char *name) : 
  AliAnalysisTaskJetShapeBase(name),  
  fM1st(0),
  fM2nd(0),
  fDeriv1st(0),
  fDeriv2nd(0),
  fPartialExclusion(0),
  fh2MSubPtSubAll(0x0),
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

  fh2MSubPtSubAll          = new TH2F*[fNcentBins];
  fh2PtTrueSubFacV1        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV1         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV1        = new TH2F*[fNcentBins];
  fh2NConstSubFacV1        = new TH2F*[fNcentBins];
  fh2PtTrueSubFacV2        = new TH2F*[fNcentBins];
  fh2PtRawSubFacV2         = new TH2F*[fNcentBins];
  fh2PtCorrSubFacV2        = new TH2F*[fNcentBins];
  fh2NConstSubFacV2        = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2MSubPtSubAll[i]          = 0;
    fh2PtTrueSubFacV1[i]        = 0;
    fh2PtRawSubFacV1[i]         = 0;
    fh2PtCorrSubFacV1[i]        = 0;
    fh2NConstSubFacV1[i]        = 0;
    fh2PtTrueSubFacV2[i]        = 0;
    fh2PtRawSubFacV2[i]         = 0;
    fh2PtCorrSubFacV2[i]        = 0;
    fh2NConstSubFacV2[i]        = 0;
  }

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

  AliAnalysisTaskJetShapeBase::UserCreateOutputObjects();

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

  const Int_t nBinsV1  = 60;
  const Double_t minV1 = -60.;
  const Double_t maxV1 = 0.;

  const Int_t nBinsV2  = 60;
  const Double_t minV2 = -30.;
  const Double_t maxV2 = 0.;

  const Int_t nBinsdMr  = 200;
  const Double_t mindMr = -2.;
  const Double_t maxdMr = 2.;
  Double_t *binsdMr = new Double_t[nBinsdMr+1];
  for(Int_t i=0; i<=nBinsdMr; i++) binsdMr[i]=(Double_t)mindMr + (maxdMr-mindMr)/nBinsdMr*(Double_t)i ;

  //These are good for pPb
  Int_t nBinsRho = 50;
  Double_t minRho = 0.;
  Double_t maxRho = 20.;
  Int_t nBinsRhom = 50;
  Double_t minRhom = 0.;
  Double_t maxRhom = 1.;
  

  const Int_t nBinsSparse2 = 8;
  //#it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}
  const Int_t nBins2[nBinsSparse2] = {nBinsDM, nBinsDpT, nBinsM, nBinsM, nBinsPt, nBinsPt, nBinsRho, nBinsRhom};
  const Double_t xmin2[nBinsSparse2]  = {minDM, minDpT, minM, minM, minPt, minPt, minRho, minRhom};
  const Double_t xmax2[nBinsSparse2]  = {maxDM, maxDpT, maxM, maxM, maxPt, maxPt, maxRho, maxRhom};

  TString histName = "";
  TString histTitle = "";
  TString varName = "#it{M}_{jet}";
  if(fJetMassVarType==kRatMPt) varName = "#it{M}_{jet}/#it{p}_{T,jet}";
  
  for(Int_t i = 0; i < fNcentBins; i++){
    histName = Form("fh2MSubPtSubAll_%d",i);
    histTitle = Form("fh2MSubPtSubAll_%d;%s;#it{p}_{T, sub}",i,varName.Data());
  
    fh2MSubPtSubAll[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsM,minM,maxM,nBinsPt,minPt,maxPt);
    fOutput->Add(fh2MSubPtSubAll[i]);

    //derivative histograms
    histName = Form("fh2PtTrueSubFacV1_%d",i);
    histTitle = Form("fh2PtTrueSubFacV1_%d;#it{p}_{T,true};-(#rho+#rho_{m})V_{1}",i);
    fh2PtTrueSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV1,minV1,maxV1);
    fOutput->Add(fh2PtTrueSubFacV1[i]);

    histName = Form("fh2PtRawSubFacV1_%d",i);
    histTitle = Form("fh2PtRawSubFacV1_%d;#it{p}_{T,raw};-(#rho+#rho_{m})V_{1}",i);
    fh2PtRawSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV1,minV1,maxV1);
    fOutput->Add(fh2PtRawSubFacV1[i]);

    histName = Form("fh2PtCorrSubFacV1_%d",i);
    histTitle = Form("fh2PtCorrSubFacV1_%d;#it{p}_{T,corr};-(#rho+#rho_{m})V_{1}",i);
    fh2PtCorrSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV1,minV1,maxV1);
    fOutput->Add(fh2PtCorrSubFacV1[i]);

    histName = Form("fh2NConstSubFacV1_%d",i);
    histTitle = Form("fh2NConstSubFacV1_%d;#it{N}_{const};-(#rho+#rho_{m})V_{1}",i);
    fh2NConstSubFacV1[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,100,0.,200.);
    fOutput->Add(fh2NConstSubFacV1[i]);

    histName = Form("fh2PtTrueSubFacV2_%d",i);
    histTitle = Form("fh2PtTrueSubFacV2_%d;#it{p}_{T,true};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2PtTrueSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV2,minV2,maxV2);
    fOutput->Add(fh2PtTrueSubFacV2[i]);

    histName = Form("fh2PtRawSubFacV2_%d",i);
    histTitle = Form("fh2PtRawSubFacV2_%d;#it{p}_{T,raw};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2PtRawSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV2,minV2,maxV2);
    fOutput->Add(fh2PtRawSubFacV2[i]);

    histName = Form("fh2PtCorrSubFacV2_%d",i);
    histTitle = Form("fh2PtCorrSubFacV2_%d;#it{p}_{T,corr};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2PtCorrSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsV2,minV2,maxV2);
    fOutput->Add(fh2PtCorrSubFacV2[i]);

    histName = Form("fh2NConstSubFacV2_%d",i);
    histTitle = Form("fh2NConstSubFacV2_%d;#it{N}_{const};0.5(#rho+#rho_{m})^{2}V_{2}",i);
    fh2NConstSubFacV2[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,100,0.,200.);
    fOutput->Add(fh2NConstSubFacV2[i]);
  }
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
//Bool_t AliAnalysisTaskJetShapeDeriv::Run()
//{
//  // Run analysis code here, if needed. It will be executed before FillHistograms().
//
//  return kTRUE;
//}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeDeriv::FillHistograms()
{
  // Fill histograms.
  AliEmcalJet *jet1  = NULL; //AA jet
  AliEmcalJet *jet2  = NULL; //Embedded Pythia jet
  AliEmcalJet *jetR  = NULL; //true jet for response matrix
  AliEmcalJet *jetO  = NULL; //hard-ish jet to avoid overlap of single track with
  AliVParticle *vpe  = NULL; //embedded particle

  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  if(!jetCont){
     Printf("Jet Container %d not found, return", fContainerBase);
     return kFALSE;
  }
  //Printf("FillHistograms::Jet container %p", jetCont);
  AliJetContainer *jetContO = GetJetContainer(fContainerOverlap);

  if(fOverlap && !jetContO){
     Printf("Jet Container %d not found, return", fContainerOverlap);
     return kFALSE;
  }
  AliParticleContainer *trackCont = GetParticleContainer(0);
  //if(trackCont) Printf("Ci sono");
  
  
  for(Int_t i=0; i<trackCont->GetNParticles(); i++){
     AliVParticle *vp= static_cast<AliVParticle*>(trackCont->GetAcceptParticle(i));
     if(!vp) continue;
     fhpTTracksCont->Fill(vp->Pt());
  }
  
  //rho
  AliRhoParameter* rhoParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(jetCont->GetRhoName()));
  fRho = 0;
  if (!rhoParam) {
     AliError(Form("%s: Could not retrieve rho %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoName().Data()));
      
  } else fRho = rhoParam->GetVal();
  //rhom
  AliRhoParameter* rhomParam = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(jetCont->GetRhoMassName()));
  fRhoM = 0;
  if (!rhomParam) {
     AliError(Form("%s: Could not retrieve rho_m %s (some histograms will be filled with zero)!", GetName(), jetCont->GetRhoMassName().Data()));
      
  } else fRhoM = rhomParam->GetVal();
    
  //Get leading jet in Pb-Pb event without embedded objects
  AliJetContainer *jetContNoEmb = GetJetContainer(fContainerNoEmb);
  AliEmcalJet *jetL = NULL;
  if(jetContNoEmb) jetL = jetContNoEmb->GetLeadingJet("rho");

  jetCont->ResetCurrentID();
  while((jet1 = jetCont->GetNextAcceptJet())) {
    jet2 = NULL;
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
    Double_t mjet1 = jet1->GetShapeProperties()->GetSecondOrderSubtracted();
    Double_t mUnsubjet1 = jet1->M();
    Double_t ptjet1 = jet1->Pt()-fRho*jet1->Area();
    Double_t ptUnsubjet1 = jet1->Pt();
    Double_t var = mjet1;
    
    if(fJetMassVarType==kRatMPt) {
      if(ptjet1>0. || ptjet1<0.) var = mjet1/ptjet1;
      else var = -999.;
    }

    //Fill histograms for all AA jets
    fh2MSubPtRawAll[fCentBin]->Fill(var,ptUnsubjet1);
    fh2MSubPtSubAll[fCentBin]->Fill(var,ptjet1);
    fh2PtRawSubFacV1[fCentBin]->Fill(jet1->Pt(),-1.*(fRho+fRhoM)*jet1->GetShapeProperties()->GetFirstDerivative());
    fh2PtCorrSubFacV1[fCentBin]->Fill(jet1->Pt()-fRho*jet1->Area(),-1.*(fRho+fRhoM)*jet1->GetShapeProperties()->GetFirstDerivative());
    fh2NConstSubFacV1[fCentBin]->Fill(jet1->GetNumberOfTracks(),-1.*(fRho+fRhoM)*jet1->GetShapeProperties()->GetFirstDerivative());
    fh2PtRawSubFacV2[fCentBin]->Fill(jet1->Pt(),0.5*(fRho+fRhoM)*(fRho+fRhoM)*jet1->GetShapeProperties()->GetSecondDerivative());
    fh2PtCorrSubFacV2[fCentBin]->Fill(jet1->Pt()-fRho*jet1->Area(),0.5*(fRho+fRhoM)*(fRho+fRhoM)*jet1->GetShapeProperties()->GetSecondDerivative());
    fh2NConstSubFacV2[fCentBin]->Fill(jet1->GetNumberOfTracks(),0.5*(fRho+fRhoM)*(fRho+fRhoM)*jet1->GetShapeProperties()->GetSecondDerivative());
    
    Double_t fraction = 0.;
    fMatch = 0;
    fJet2Vec->SetPtEtaPhiM(0.,0.,0.,0.);
    if(fSingleTrackEmb) {
       vpe = GetEmbeddedConstituent(jet1);
       if(vpe) {
       	  Bool_t reject = kFALSE;
       	  if(fPartialExclusion) {
       	     
       	     TRandom3 rnd;
       	     rnd.SetSeed(0);
       	     
       	     Double_t ncoll = 6.88; //GetNColl(); //check it out from AliAnalysisTaskDeltaPt and possibly refine
       	     
       	     Double_t prob = 0.;
       	     if(ncoll>0)
       	     	prob = 1./ncoll;
       	     
       	     if(rnd.Rndm()<=prob) reject = kTRUE; //reject cone
       	  }
       	  if(fOverlap){
       	     Int_t Njets = jetContO->GetNAcceptedJets();
       	     fhNJetsSelEv->Fill(Njets);
       	     jetContO->ResetCurrentID();
       	     while((jetO = jetContO->GetNextAcceptJet())){
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
       	     	fhJetEtaPhiOvl->Fill(jetO->Eta(), jetO->Phi());
       	     	if( deltaR < fRadius) {
       	     	   reject = kTRUE;
       	     	   break;
       	     	}
       	     }
       	  }
          if(!reject){
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
    
    //Fill histograms for matched jets
    fh2MSubMatch[fCentBin]->Fill(var,fMatch);
    if(fMatch==1) {
       fhJetSubMatchEtaPhiPt->Fill(jet1->Eta(), jet1->Phi(), ptjet1);
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
        var2   = jetR->M();
        ptJetR = jetR->Pt();
      }
      if(fSingleTrackEmb && vpe) {
      	 if(fFromTree){
      	    Int_t exit = MatchEmbeddedConstituentWithParticleLevel(); //here is GetEntry
      	    
      	    if(exit>0) {
      	       //mJetR  = fVecP->M();
      	       var2   = fVecP->M();
      	       ptJetR = fVecP->Pt();
      	    }
      	 } else{
      	    //mJetR  = vpe->M();
      	    var2   = vpe->M();
      	    ptJetR = vpe->Pt();
      	 }
      }
      
      if(fJetMassVarType==kRatMPt) {
        if(ptJetR>0. || ptJetR<0.) var2 /= ptJetR;
      }
      fh3MSubPtTrueLeadPt[fCentBin]->Fill(var,ptJetR,jet1->MaxTrackPt());
      fh3MTruePtTrueLeadPt[fCentBin]->Fill(var2,ptJetR,jet1->MaxTrackPt());
      fh3PtTrueDeltaMLeadPt[fCentBin]->Fill(ptJetR,var-var2,jet1->MaxTrackPt());
      if(var2>0.) fh3PtTrueDeltaMRelLeadPt[fCentBin]->Fill(ptJetR,(var-var2)/var2,jet1->MaxTrackPt());
      if(fFromTree){
      	 Double_t varsp[10] = {var,var2,ptjet1,ptJetR, fVecD->M(), fVecD->Pt(), fRho, fRhoM, mUnsubjet1, ptUnsubjet1};
      	 fhnMassResponse[fCentBin]->Fill(varsp);
      } else {
      	 Double_t varsp[9] = {var,var2,ptjet1,ptJetR,jet1->MaxTrackPt(), fRho, fRhoM, mUnsubjet1, ptUnsubjet1};//MRec,MTrue,PtRec,PtTrue,PtLeadRec
      	 fhnMassResponse[fCentBin]->Fill(varsp);
      }
      Double_t varsp1[6];
      varsp1[0] = var-var2;
      varsp1[1] = ptjet1-ptJetR;
      varsp1[2] = var;
      varsp1[3] = var2;
      varsp1[4] = ptjet1;
      varsp1[5] = ptJetR;

      fhnDeltaMass[fCentBin]->Fill(varsp1);
      
      //#it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{unsub} - #it{M}_{part}; #it{p}_{T,unsub} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}
      Double_t varsp2[10] = {var-var2, ptjet1-ptJetR, mUnsubjet1 - var2, ptUnsubjet1 - ptJetR, var2, mUnsubjet1, ptjet1, ptUnsubjet1, fRho, fRhoM};
      fhnDeltaMassAndBkgInfo->Fill(varsp2);
    }
    
    if(fCreateTree) {      
      fJet1Vec->SetPxPyPzE(jet1->Px(),jet1->Py(),jet1->Pz(),jet1->E());
      fArea = (Float_t)jet1->Area();
      fAreaPhi = (Float_t)jet1->AreaPhi();
      fAreaEta = (Float_t)jet1->AreaEta();
      fNConst = (Int_t)jet1->GetNumberOfTracks();
      fM1st   = (Float_t)jet1->GetShapeProperties()->GetFirstOrderSubtracted();
      fM2nd   = (Float_t)jet1->GetShapeProperties()->GetSecondOrderSubtracted();
      fDeriv1st = (Float_t)jet1->GetShapeProperties()->GetFirstDerivative();
      fDeriv2nd = (Float_t)jet1->GetShapeProperties()->GetSecondDerivative();
      fTreeJetBkg->Fill();
    }
  }
  return kTRUE;
}

//________________________________________________________________________
//AliVParticle* AliAnalysisTaskJetShapeDeriv::GetEmbeddedConstituent(AliEmcalJet *jet) {
//
//  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
//  //Printf("JEt container %p", jetCont);
//  AliVParticle *vp = 0x0;
//  AliVParticle *vpe = 0x0; //embedded particle
//  Int_t nc = 0;
//  for(Int_t i=0; i<jet->GetNumberOfTracks(); i++) {
//    vp = static_cast<AliVParticle*>(jet->TrackAt(i, jetCont->GetParticleContainer()->GetArray()));
//    //    if (vp->TestBits(TObject::kBitMask) != (Int_t)(TObject::kBitMask) ) continue;
//    //Printf("vp %p", vp);
//    if(!vp) continue;
//    Int_t lab = TMath::Abs(vp->GetLabel());
//    if (lab < fMinLabelEmb || lab > fMaxLabelEmb)
//      continue;
//    if(!vpe) vpe = vp;
//    else if(vp->Pt()>vpe->Pt()) vpe = vp;
//    nc++;
//  }
//
//  AliDebug(11,Form("Found %d embedded particles",nc));
//  return vpe;
//}


