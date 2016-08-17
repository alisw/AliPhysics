#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TTree.h>

#include <AliParticleContainer.h>
#include <AliJetContainer.h>
#include <AliRhoParameter.h>
#include <AliEmcalJet.h>
#include <AliVParticle.h>

#include "AliAnalysisTaskJetShapeConst.h"

ClassImp(AliAnalysisTaskJetShapeConst)

//________________________________________________________________________
AliAnalysisTaskJetShapeConst::AliAnalysisTaskJetShapeConst() : 
  AliAnalysisTaskJetShapeBase("AliAnalysisTaskJetShapeConst"),
  fhptjetSMinusSingleTrack(0x0),
  fhJet1vsJetTag(0x0),
  fhNconstit(0x0),
  fhAreaJet(0x0)
{
  // Default constructor.


}

//________________________________________________________________________
AliAnalysisTaskJetShapeConst::AliAnalysisTaskJetShapeConst(const char *name) : 
  AliAnalysisTaskJetShapeBase(name),  
  fhptjetSMinusSingleTrack(0x0),
  fhJet1vsJetTag(0x0),
  fhNconstit(0x0),
  fhAreaJet(0x0)
{
  // Standard constructor.

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

  AliAnalysisTaskJetShapeBase::UserCreateOutputObjects();
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fhptjetSMinusSingleTrack = new TH1F("fhptjetSMinusSingleTrack", "Subtraction of single track #it{p}_{T}; #it{p}_{T, jet} - #it{p}_{T, Emb Track};Entries", 500,-10.,110.);
  fOutput->Add(fhptjetSMinusSingleTrack);
  
  fhJet1vsJetTag = new TH2F("fhJet1vsJetTag", "Number of jets vs tagged jets; #it{N}_{jet,Tot}; #it{N}_{jet,Tag}", 30, 1., 30., 30, 1., 30.);
  fOutput->Add(fhJet1vsJetTag);
  
  fhNconstit = new TH1F("fhNconstit", "Number of constituents (matched jets); #it{N}_{constituents}", 21, 0., 20.);
  fOutput->Add(fhNconstit);
  
  fhAreaJet = new TH1F("fhAreaJet", "Area (matched jets); Area", 400., 0., 4);
  fOutput->Add(fhAreaJet);
  
  TH1::AddDirectory(oldStatus);
  
  PostData(1, fOutput);
}

//________________________________________________________________________
//Bool_t AliAnalysisTaskJetShapeConst::Run()
//{
//  // Run analysis code here, if needed. It will be executed before FillHistograms().
//
//  return kTRUE;
//}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetShapeConst::FillHistograms()
{
  // Fill histograms.

  AliEmcalJet* jet1  = NULL; //AA jet
  AliEmcalJet *jet2  = NULL; //Embedded Pythia jet
  AliEmcalJet *jetR  = NULL; //true jet for response matrix
  AliEmcalJet *jetS  = NULL; //subtracted jet
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
  
  
  Int_t njet1 = 0, ntagjet2 = 0;
  
  AliDebug(11,Form("NJets  Incl: %d  Csub: %d",jetCont->GetNJets(),jetContS->GetNJets()));
  if(jetCont) {
    jetCont->ResetCurrentID();
    while((jet1 = jetCont->GetNextAcceptJet())) {
      njet1++;
      jet2  = NULL;
      jetS  = NULL;
      if(jet1->GetTagStatus()<1 || !jet1->GetTaggedJet())
        continue;
     
     ntagjet2++;
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
      //Double_t ptjetSMinusEmbTrpt = ptjetS;
      Double_t ptjetSMinusEmbTrpt = ptUnsubjet1;
      
      if(fJetMassVarType==kRatMPt) {
      	 if(ptjetS>0. || ptjetS<0.) var = mjetS/ptjetS;
      	 else var = -999.;
      }
      
      //Fill histograms for all AA jets
      fh2MSubPtRawAll[fCentBin]->Fill(var,ptjetS);
      
      Double_t fraction = 0.;
      fMatch = 0;
      fJet2Vec->SetPtEtaPhiM(0.,0.,0.,0.);
      //Printf("EMB task: pT jet %d = %f", njet1, jet1->Pt());
      if(fSingleTrackEmb) {
      	 vpe = GetEmbeddedConstituent(jet1);
      	 if(vpe){
      	    Bool_t reject = kFALSE;
      	    ptjetSMinusEmbTrpt -= vpe->Pt();
      	    fhptjetSMinusSingleTrack->Fill(ptjetSMinusEmbTrpt);
      	    //Printf("EMB task: pT jet matched = %f -> %f", jet1->Pt(), ptjetSMinusEmbTrpt);
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
      	 fhJetSubMatchEtaPhiPt->Fill(jetS->Eta(), jetS->Phi(), ptjetS);
      	 Double_t drToLJ = -1.;
      	 if(jetL) drToLJ = jet1->DeltaR(jetL);
      	 if(fSingleTrackEmb && vpe)
      	    drToLJ = jet1->DeltaR(vpe);
      	 fh3MSubPtRawDRMatch[fCentBin]->Fill(var,ptjet1,drToLJ);
      	 Double_t var2 = 0.;
      	 //Double_t mJetR = 0.;
      	 Double_t ptJetR = 0.;
      	 if(jetR) {
      	    //mJetR  = jetR->M();
      	    var2 = jetR->M();
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
      	 //M sub;M true;#it{p}_{T,sub};#it{p}_{T,true};#it{p}_{T,lead trk}
      	 
      	 if(!fCreateTree){
      	 	 if(fFromTree){
      	 	 	 // Mass sub; Mass true;#it{p}_{T,sub};#it{p}_{T,true};%s (emb, det); #it{p}_{T,emb det}; #rho; #rho_{m};
      	 	 	 Double_t varsp[10] = {var,var2,ptjetS,ptJetR, fVecD->M(), fVecD->Pt(), fRho, fRhoM, mUnsubjet1, ptUnsubjet1};
      	 	 	 fhnMassResponse[fCentBin]->Fill(varsp);
      	 	 } else {
      	 	 	 Double_t varsp[9] = {var,var2,ptjetS,ptJetR,jetS->MaxTrackPt(), fRho, fRhoM, mUnsubjet1, ptUnsubjet1};
      	 	 	 fhnMassResponse[fCentBin]->Fill(varsp);
      	 	 }
      	 	 
      	 	 Double_t varsp1[8];
      	 	 //#it{M}_{det,Const} - #it{M}_{part}; #it{p}_{T,det,Const} - #it{p}_{T,part}; #it{M}_{det,Const};  #it{M}_{part}; #it{p}_{T,det,Const}; #it{p}_{T,part}; #it{p}_{T,det,A}
      	 	 varsp1[0] = var-var2;
      	 	 varsp1[1] = ptjetS-ptJetR;
      	 	 varsp1[2] = var;
      	 	 varsp1[3] = var2;
      	 	 varsp1[4] = ptjetS;
      	 	 varsp1[5] = ptJetR;
      	 	 varsp1[6] = ptjet1;
      	 	 varsp1[7] = ptjet1 - ptJetR;
      	 	 
      	 	 //fhnDeltaMass[fCentBin]->Fill(varsp1);
      	 	 
      	 	 //#it{M}_{det} - #it{M}_{part}; #it{p}_{T,det} - #it{p}_{T,part}; #it{M}_{unsub} - #it{M}_{part}; #it{p}_{T,unsub} - #it{p}_{T,part}; #it{M}_{det};  #it{M}_{unsub}; #it{p}_{T,det}; #it{p}_{T,unsub}; #rho ; #rho_{m}
      	 	 
      	 	 Double_t varsp2[10] = {var-var2, ptjetS-ptJetR, mUnsubjet1 - var2, ptUnsubjet1 - ptJetR, var2, mUnsubjet1, ptjetS, ptUnsubjet1, fRho, fRhoM};
      	 	 fhnDeltaMassAndBkgInfo->Fill(varsp2);
      	 }
      	 fhNconstit->Fill(jet1->GetNumberOfConstituents());
      	 fhAreaJet ->Fill(jet1->Area());
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
  
  fhJet1vsJetTag->Fill(njet1, ntagjet2);
  
  return kTRUE;
}


