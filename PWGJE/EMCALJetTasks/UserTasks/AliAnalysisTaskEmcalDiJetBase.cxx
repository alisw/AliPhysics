//
// Dijet base analysis task.
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

#include "AliAnalysisTaskEmcalDiJetBase.h"

ClassImp(AliAnalysisTaskEmcalDiJetBase)

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetBase::AliAnalysisTaskEmcalDiJetBase() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalDiJetBase", kTRUE),
  fDebug(kFALSE),
  fJetCorrelationType(kCorrelateAll),
  fJetFullChargedMatchingType(kFraction),
  fTriggerClass(""),
  fContainerCharged(1),
  fContainerFull(0),
  fContainerChargedMC(3),
  fContainerFullMC(2),
  fRhoType(0),
  fRhoChVal(0),
  fRhoFullVal(0),
  fDoChargedCharged(kTRUE),
  fDoFullCharged(kTRUE),
  fDoFullFull(kFALSE),
  fPtMinTriggerJet(10.),
  fDoPtBias(kTRUE),
  fMinFractionShared(0.5),
  fMatchingDone(kFALSE),
  faFullFracIndex(0),
  faFullFracIndexMC(0),
  fhNEvents(0),
  fHistTrialsSelEvents(0)
{
  // Default constructor.

  
  SetMakeGeneralHistograms(kTRUE);

}

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetBase::AliAnalysisTaskEmcalDiJetBase(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fDebug(kFALSE),
  fJetCorrelationType(kCorrelateAll),
  fJetFullChargedMatchingType(kFraction),
  fTriggerClass(""),
  fContainerCharged(1),
  fContainerFull(0),
  fContainerChargedMC(3),
  fContainerFullMC(2),
  fRhoType(0),
  fRhoChVal(0),
  fRhoFullVal(0),
  fDoChargedCharged(kTRUE),
  fDoFullCharged(kTRUE),
  fDoFullFull(kFALSE),
  fPtMinTriggerJet(10.),
  fDoPtBias(kTRUE),
  fMinFractionShared(0.5),
  fMatchingDone(kFALSE),
  faFullFracIndex(0),
  faFullFracIndexMC(0),
  fhNEvents(0),
  fHistTrialsSelEvents(0)
{
  // Standard constructor.

  SetMakeGeneralHistograms(kTRUE);

}

//________________________________________________________________________
AliAnalysisTaskEmcalDiJetBase::~AliAnalysisTaskEmcalDiJetBase()
{
  // Destructor.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetBase::SelectEvent() {
  //
  // Decide if event should be selected for analysis
  //

  fhNEvents->Fill(3.5);

  if(!fTriggerClass.IsNull()) {
    //Check if requested trigger was fired
    TString firedTrigClass = InputEvent()->GetFiredTriggerClasses();

    if(!firedTrigClass.Contains(fTriggerClass))
      return kFALSE;
    else if(fTriggerClass.Contains("J1") && fTriggerClass.Contains("J2")) { //if events with J1&&J2 are requested
      if(!firedTrigClass.Contains("J1") || !firedTrigClass.Contains("J2") ) //check if both are fired
        return kFALSE;
    }
    else if(fTriggerClass.Contains("J1") && firedTrigClass.Contains("J2")) //if J2 is requested also add triggers which have J1&&J2. Reject if J1 is requested and J2 is fired
      return kFALSE;
  }
  
  fhNEvents->Fill(1.5);

  fHistTrialsSelEvents->Fill(fPtHardBin, fNTrials);
  
  return kTRUE;

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetBase::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  fhNEvents = new TH1F("fhNEvents","fhNEvents;selection;N_{evt}",5,0,5);
  fOutput->Add(fhNEvents);

  //Pythia info
  fHistTrialsSelEvents = new TH1F("fHistTrialsSelEvents", "fHistTrialsSelEvents", 12, 0, 12);
  fHistTrialsSelEvents->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistTrialsSelEvents->GetYaxis()->SetTitle("trials");
  fOutput->Add(fHistTrialsSelEvents);

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetBase::IsSameJet(Int_t ijt, Int_t ija, Int_t type, Bool_t isMC) {
   //check if two jets are the same one

   Bool_t bSame = kFALSE;

   if(type<2 && ijt==ija)
     bSame = kTRUE;
   if(type==2) {

     if(fJetFullChargedMatchingType==kFraction) {
       if(isMC && faFullFracIndexMC.At(ijt)==ija)
	 bSame = kTRUE;
       else if(!isMC && faFullFracIndex.At(ijt)==ija)
	 bSame = kTRUE;
     }
     else if(fJetFullChargedMatchingType==kGeo) {
       Int_t contFull = fContainerFull;
       Int_t contChar = fContainerCharged;
       if(isMC) {
	 contFull = fContainerFullMC;
	 contChar = fContainerChargedMC;
       }
       AliEmcalJet *fullJet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ijt, contFull));
       AliEmcalJet *chJet   = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ija, contChar));
       AliEmcalJet *matchJet = fullJet->ClosestJet();
       if(chJet==matchJet)
	 bSame = kTRUE;
     }
     else if(fJetFullChargedMatchingType==kNoMatching) {
       return kFALSE;
     }
     else{
       AliWarning(Form("%s: matching type unknown", GetName()));
     }
   }

   return bSame;
 }


//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDiJetBase::GetJetPt(const AliEmcalJet *jet, Int_t type) {

  if(!jet) return -99;

  if(type==0)
    return jet->Pt() - fRhoFullVal * jet->Area();
  else if(type==1)
    return jet->Pt() - fRhoChVal * jet->Area();
  else
    return jet->Pt();
  
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDiJetBase::GetZ(const AliVParticle *trk, const AliEmcalJet *jet)          const
{  
  // Get Z of constituent trk

  return GetZ(trk->Px(),trk->Py(),trk->Pz(),jet->Px(),jet->Py(),jet->Pz());
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDiJetBase::GetZ(Double_t trkPx, Double_t trkPy, Double_t trkPz, Double_t jetPx, Double_t jetPy, Double_t jetPz) const
{
  // 
  // Get the z of a constituent inside of a jet
  //

  Double_t pJetSq = jetPx*jetPx+jetPy*jetPy+jetPz*jetPz;

  if(pJetSq>0.)
    return (trkPx*jetPx+trkPy*jetPy+trkPz*jetPz)/pJetSq;
  else {
    AliWarning(Form("%s: strange, pjet*pjet seems to be zero pJetSq: %f",GetName(), pJetSq)); 
    return 0;
  }
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDiJetBase::GetDeltaR(const AliEmcalJet* jet1, const AliEmcalJet* jet2) const {
  //
  // Helper function to calculate the distance between two jets
  //

  Double_t dPhi = jet1->Phi() - jet2->Phi();
  if(dPhi>TMath::Pi())dPhi = dPhi - 2.*TMath::Pi();
  if(dPhi<(-1.*TMath::Pi()))dPhi = dPhi + 2.*TMath::Pi();
  Double_t dEta = jet1->Eta() - jet2->Eta();
  Double_t dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);
  return dR;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDiJetBase::GetDeltaPhi(const AliEmcalJet* jet1, const AliEmcalJet* jet2) {
  //
  // Calculate azimuthal angle between the axises of the jets
  //

  return GetDeltaPhi(jet1->Phi(),jet2->Phi());

}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDiJetBase::GetDeltaPhi(Double_t phi1,Double_t phi2) {
  //
  // Calculate azimuthal angle between the axises of the jets
  //

  Double_t dPhi = phi1-phi2;

  if(dPhi <-0.5*TMath::Pi())  dPhi += TMath::TwoPi();
  if(dPhi > 1.5*TMath::Pi())  dPhi -= TMath::TwoPi();

  return dPhi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalDiJetBase::GetFractionSharedPt(const AliEmcalJet *jetFull, const AliEmcalJet *jetCharged) const
{
  //
  // Get fraction of shared pT between matched full and charged jet
  // Uses charged jet pT as baseline: fraction = \Sum_{const,full jet} pT,const,i / pT,jet,ch
  //

  Double_t fraction = 0.;
  Double_t jetPtCh = jetCharged->Pt();
 
  if(jetPtCh>0) {

    AliDebug(11,Form("%s: nConstituents: %d, ch: %d  chne: %d ne: %d",GetName(),jetFull->GetNumberOfConstituents(),jetCharged->GetNumberOfTracks(),jetFull->GetNumberOfTracks(),jetFull->GetNumberOfClusters()));
    
    Double_t sumPt = 0.;
    AliVParticle *vpf = 0x0;
    Int_t iFound = 0;
    for(Int_t icc=0; icc<jetCharged->GetNumberOfTracks(); icc++) {
      Int_t idx = (Int_t)jetCharged->TrackAt(icc);
      iFound = 0;
      for(Int_t icf=0; icf<jetFull->GetNumberOfTracks(); icf++) {
	if(idx == jetFull->TrackAt(icf) && iFound==0 ) {
	  iFound=1;
	  vpf = static_cast<AliVParticle*>(jetFull->TrackAt(icf, fTracks));
	  if(!vpf) continue;
	  if(vpf->Charge()!=0)
	    sumPt += vpf->Pt();
	  continue;
	}
      }
    }
    
    fraction = sumPt/jetPtCh;
  }

  AliDebug(11,Form("%s: charged shared fraction: %.2f",GetName(),fraction));

  return fraction;

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetBase::MatchJetsGeo(Int_t cFull, Int_t cCharged,
						 Int_t iDebug, Float_t maxDist, Int_t type) {

  //
  // Match the full jets to the corresponding charged jets
  // Translation of AliAnalysisHelperJetTasks::GetClosestJets to AliEmcalJet objects
  // type: 0 = use acceptance cuts of container  1 = allow 0.1 one more for cCharged(MC) in eta 2 = allow 0.1 more in eta and phi for cCharged(MC)

  const int kMode = 3;

  const Int_t nChJets   = GetNJets(cCharged);
  const Int_t nFullJets = GetNJets(cFull);

  if(nFullJets==0 || nChJets==0) return;

  TArrayI faChNeMatchIndex;
  faChNeMatchIndex.Set(nChJets+1);
  faChNeMatchIndex.Reset(-1);

  TArrayI faChMatchIndex;
  faChMatchIndex.Set(nFullJets+1);
  faChMatchIndex.Reset(-1);

  static TArrayS iFlag(nChJets*nFullJets);
  if(iFlag.GetSize()<(nChJets*nFullJets)){
    iFlag.Set(nChJets*nFullJets+1);
  }
  iFlag.Reset(0);

  AliJetContainer *contCh = GetJetContainer(cCharged);

  // find the closest distance to the full jet
  for(int ifu = 0;ifu<nFullJets;ifu++){

    AliEmcalJet *fullJet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ifu, cFull));
    if(!fullJet) continue;

    Float_t dist = maxDist;
    
    for(int ich = 0;ich<nChJets;ich++){
      AliEmcalJet *chJet = 0x0;
      if(type==0)
	chJet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ich, cCharged));
      else {
	chJet = static_cast<AliEmcalJet*>(GetJetFromArray(ich, cCharged));
	if(!chJet) continue;
	if(type>0) {
	  if(chJet->Eta()<(contCh->GetJetEtaMin()-0.1) || chJet->Eta()>(contCh->GetJetEtaMax()+0.1))
	    continue;
	  if(type==2) {
	    if(chJet->Phi()<(contCh->GetJetPhiMin()-0.1) || chJet->Phi()>(contCh->GetJetPhiMax()+0.1))
	      continue;
	  }
	}
      }
      if(!chJet)
	continue;

      Double_t frac = GetFractionSharedPt(fullJet,chJet);
      Double_t dR = GetDeltaR(fullJet,chJet);
      if(dR<dist){
	faChMatchIndex[ifu]=ich;
	dist = dR;
      }
      if(iDebug>10) Printf("Distance (%d)--(%d) %3.3f frac:%.2f",ifu,ich,dR,frac);
    }
    if(faChMatchIndex[ifu]>=0) iFlag[ifu*nChJets+faChMatchIndex[ifu]]+=1;//ich closest to ifu
    if(iDebug>10) Printf("Full Distance (%d)--(%d) %3.3f flag[%d] = %d",ifu,faChMatchIndex[ifu],dist,ifu*nChJets+faChMatchIndex[ifu],iFlag[ifu*nChJets+faChMatchIndex[ifu]]);
    
    // reset...
    faChMatchIndex[ifu]=-1;

           
  }//full jet loop


  // other way around
  for(int ich = 0;ich<nChJets;ich++){
    AliEmcalJet *chJet = 0x0;
    if(type==0)
      chJet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ich, cCharged));
    else {
      chJet = static_cast<AliEmcalJet*>(GetJetFromArray(ich, cCharged));
      if(!chJet) continue;;
      if(type>0) {
	if(chJet->Eta()<(contCh->GetJetEtaMin()-0.1) || chJet->Eta()>(contCh->GetJetEtaMax()+0.1))
	  continue;
	if(type==2) {
	  if(chJet->Phi()<(contCh->GetJetPhiMin()-0.1) || chJet->Phi()>(contCh->GetJetPhiMax()+0.1))
	    continue;
	}
      }
    }
    if(!chJet)
      continue;

    Float_t dist = maxDist;
    for(int ifu = 0;ifu<nFullJets;ifu++){
      AliEmcalJet *fullJet = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ifu, cFull));
      if(!fullJet)
	continue;

      Double_t dR = GetDeltaR(fullJet,chJet); 
      if(dR<dist){
	faChNeMatchIndex[ich]=ifu;
        dist = dR;
      }   
    }
    if(faChNeMatchIndex[ich]>=0) iFlag[faChNeMatchIndex[ich]*nChJets+ich]+=2;//ifu closest to ich
    if(iDebug>10) Printf("Other way Distance (%d)--(%d) %3.3f flag[%d] = %d",faChNeMatchIndex[ich],ich,dist,faChNeMatchIndex[ich]*nChJets+ich,iFlag[faChNeMatchIndex[ich]*nChJets+ich]);

    // reset...
    faChNeMatchIndex[ich]=-1;
    
  }
  
  
  // check for "true" correlations
  for(int ifu = 0;ifu<nFullJets;ifu++){
    for(int ich = 0;ich<nChJets;ich++){
      AliDebug(11,Form("%s: Flag[%d][%d] %d ",GetName(),ifu,ich,iFlag[ifu*nChJets+ich]));
      
      if(kMode==3){
        // we have a uniqe correlation
        if(iFlag[ifu*nChJets+ich]==3){

	  AliEmcalJet *chJet = static_cast<AliEmcalJet*>(GetJetFromArray(ich, cCharged));
	  AliEmcalJet *fullJet = static_cast<AliEmcalJet*>(GetJetFromArray(ifu, cFull));
	  Double_t dR = GetDeltaR(fullJet,chJet); 

	  AliDebug(11,Form("closest jets %d  %d  dR =  %f",ich,ifu,dR));

	  chJet->SetClosestJet(fullJet,dR);
	  fullJet->SetClosestJet(chJet,dR);

	}
      }
    }
  }

  fMatchingDone = kTRUE;
  
}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetBase::SetChargedFractionIndex() {

  // take each full jet and set the index of the charged jet with the largest shared charged fraction

  const Int_t nJetsCh   = GetNJets(fContainerCharged);
  const Int_t nJetsFull = GetNJets(fContainerFull);
  faFullFracIndex.Set(nJetsFull+1);
  faFullFracIndex.Reset(-1);

  AliJetContainer  *cont = GetJetContainer(fContainerFull); 
  Float_t radius =  cont->GetJetRadius();

  for(Int_t ifu = 0; ifu<nJetsFull; ifu++) {
    Double_t frac = 0.;
    Double_t dist = 10.;
    AliEmcalJet *jetFull = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ifu, fContainerFull));
    if(!jetFull) {
      faFullFracIndex.SetAt(-1,ifu);
      continue;
    }

    for(Int_t ich = 0; ich<nJetsCh; ich++) {
      AliEmcalJet *jetCh = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ich, fContainerCharged));
      if(!jetCh)
	continue;
      Double_t tmpFrac = GetFractionSharedPt(jetFull,jetCh);
      dist = GetDeltaR(jetFull,jetCh);
      if(tmpFrac>frac && dist<radius) {
	frac = tmpFrac;
	faFullFracIndex.SetAt(ich,ifu);
      }
    }

  }

}

//________________________________________________________________________
void AliAnalysisTaskEmcalDiJetBase::SetChargedFractionIndexMC() {

  // take each full jet and set the index of the charged jet with the largest shared charged fraction

  const Int_t nJetsCh   = GetNJets(fContainerChargedMC);
  const Int_t nJetsFull = GetNJets(fContainerFullMC);
  faFullFracIndexMC.Set(nJetsFull);
  faFullFracIndexMC.Reset(-1);

  AliJetContainer  *cont = GetJetContainer(fContainerFullMC); 
  Float_t radius =  cont->GetJetRadius();

  for(Int_t ifu = 0; ifu<nJetsFull; ifu++) {
    Double_t frac = 0.;
    Double_t dist = 10.;
    AliEmcalJet *jetFull = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ifu, fContainerFullMC));
    if(!jetFull) {
      faFullFracIndexMC.SetAt(-1,ifu);
      continue;
    }

    for(Int_t ich = 0; ich<nJetsCh; ich++) {
      AliEmcalJet *jetCh = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ich, fContainerChargedMC));
      if(!jetCh)
	continue;
      Double_t tmpFrac = GetFractionSharedPt(jetFull,jetCh);
      dist = GetDeltaR(jetFull,jetCh);
      if(tmpFrac>frac && dist<radius) {
	frac = tmpFrac;
	faFullFracIndexMC.SetAt(ich,ifu);
      }
    }

  }

}

//_______________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalDiJetBase::GetLeadingJetOppositeHemisphere(Int_t type, Int_t typea, const AliEmcalJet *jetTrig) {

  // Get leading jet in opposite hemisphere from trigger jet
  // type = correlation type
  // typea = container of associated jets

  Int_t nJetsAssoc = GetNJets(typea);
  Double_t ptLead = -999;
  Int_t    iJetLead = -1;
  for(Int_t ija=0; ija<nJetsAssoc; ija++) {

    AliEmcalJet *jetAssoc = NULL;
    if(type==0) {
      jetAssoc = static_cast<AliEmcalJet*>(GetJetFromArray(ija, typea));
      if(TMath::Abs(jetAssoc->Eta())>0.5)
        jetAssoc = NULL;
    }
    else
      jetAssoc = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ija, typea));

    if(!jetAssoc)
      continue;

    Double_t dPhi = GetDeltaPhi(jetTrig,jetAssoc);
    Double_t phiMin = 0.5*TMath::Pi();
    Double_t phiMax = 1.5*TMath::Pi();
    if(dPhi<phiMin || dPhi>phiMax)
      continue;

    Double_t jetAssocPt = GetJetPt(jetAssoc,typea);

    if(jetAssocPt>ptLead) {
      ptLead = jetAssocPt;
      iJetLead = ija;
    }

  }

  AliEmcalJet *jetAssocLead = NULL;
  if(iJetLead>-1)
    jetAssocLead = static_cast<AliEmcalJet*>(GetJetFromArray(iJetLead, typea));

  return jetAssocLead;

}

//_______________________________________________________________________
AliEmcalJet* AliAnalysisTaskEmcalDiJetBase::GetSecondLeadingJetOppositeHemisphere(Int_t type, Int_t typea, const AliEmcalJet *jetTrig) {

  // Get leading jet in opposite hemisphere from trigger jet
  // type = correlation type
  // typea = container of associated jets

  Int_t nJetsAssoc = GetNJets(typea);
  Double_t ptLead = -999;
  Int_t    iJetLead = -1;
  Int_t    iJetLead2 = -1;
  for(Int_t ija=0; ija<nJetsAssoc; ija++) {

    AliEmcalJet *jetAssoc = NULL;
    if(type==0) {
      jetAssoc = static_cast<AliEmcalJet*>(GetJetFromArray(ija, typea));
      if(TMath::Abs(jetAssoc->Eta())>0.5)
        jetAssoc = NULL;
    }
    else
      jetAssoc = static_cast<AliEmcalJet*>(GetAcceptJetFromArray(ija, typea));

    if(!jetAssoc)
      continue;

    Double_t dPhi = GetDeltaPhi(jetTrig,jetAssoc);
    Double_t phiMin = 0.5*TMath::Pi();
    Double_t phiMax = 1.5*TMath::Pi();
    if(dPhi<phiMin || dPhi>phiMax)
      continue;

    Double_t jetAssocPt = GetJetPt(jetAssoc,typea);

    if(jetAssocPt>ptLead) {
      iJetLead2 = iJetLead;
      ptLead = jetAssocPt;
      iJetLead = ija;
    }

  }

  AliEmcalJet *jetAssocLead2 = NULL;
  if(iJetLead2>-1)
    jetAssocLead2 = static_cast<AliEmcalJet*>(GetJetFromArray(iJetLead2, typea));

  return jetAssocLead2;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalDiJetBase::RetrieveEventObjects() {
  //
  // Retrieve objects from event.
  //

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  if(fRhoType==0) {
    fRhoFullVal = 0.;
    fRhoChVal = 0.;
  }
  if(fRhoType==1) {
    fRhoFullVal = GetRhoVal(fContainerFull);
    fRhoChVal = GetRhoVal(fContainerCharged);
  }

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskEmcalDiJetBase::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}
