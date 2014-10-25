// $Id$

//////////
//Measure Jet-hadron correlations
//Does event Mixing using AliEventPoolManager
/////////

#include "AliAnalysisTaskEmcalJetHMEC.h"

#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include <TClonesArray.h>
#include <TParticle.h>
#include "AliVTrack.h"
#include "TParameter.h"

#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliAODJet.h"
#include "AliEmcalJet.h"
#include "AliESDtrackCuts.h"

#include "TVector3.h"
#include "AliPicoTrack.h"
#include "AliEventPoolManager.h"

ClassImp(AliAnalysisTaskEmcalJetHMEC)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC() : 
  AliAnalysisTaskEmcalJet("HMEC",kFALSE),
  fTracksName(""),
  fJetsName(""),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fTrkBias(5),
  fClusBias(5),
  fTrkEta(0.9),
  fDoEventMixing(0),
  fMixingTracks(50000), fNMIXtracks(5000), fNMIXevents(5),
  fTriggerEventType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDoLessSparseAxes(0), fDoWiderTrackBin(0),
  fCentBinSize(1),
  fESD(0),
  fAOD(0), 
  fPoolMgr(0x0), 
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistJetHEtaPhi(0), 
  fhnMixedEvents(0x0),
  fhnJH(0x0)
{
  // Default Constructor

  for(Int_t ipta=0; ipta<7; ipta++){
    fHistTrackEtaPhi[ipta]=0;
  }

  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistLeadJetPt[icent]=0;
    fHistLeadJetPtBias[icent]=0;
    fHistJetPtTT[icent]=0;
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
  	    fHistJetH[icent][iptjet][ieta]=0;
	    fHistJetHBias[icent][iptjet][ieta]=0;
	    fHistJetHTT[icent][iptjet][ieta]=0;
      }
    }
  }

}
//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC(const char *name) : 
  AliAnalysisTaskEmcalJet(name,kTRUE),
  fTracksName(""),
  fJetsName(""),
  fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),
  fTrkBias(5),
  fClusBias(5),
  fTrkEta(0.9),
  fDoEventMixing(0),
  fMixingTracks(50000), fNMIXtracks(5000), fNMIXevents(5),
  fTriggerEventType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDoLessSparseAxes(0), fDoWiderTrackBin(0),
  fCentBinSize(1),
  fESD(0),
  fAOD(0), 
  fPoolMgr(0x0), 
  fHistTrackPt(0),
  fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistJetHEtaPhi(0),
  fhnMixedEvents(0x0),
  fhnJH(0x0)
{
  // Constructor
  for(Int_t ipta=0; ipta<7; ipta++){
    fHistTrackEtaPhi[ipta]=0;
  }
  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistLeadJetPt[icent]=0;
    fHistLeadJetPtBias[icent]=0;
    fHistJetPtTT[icent]=0;
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	    fHistJetH[icent][iptjet][ieta]=0;
	    fHistJetHBias[icent][iptjet][ieta]=0;
	    fHistJetHTT[icent][iptjet][ieta]=0;
      }
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::UserCreateOutputObjects() {
  // Called once 
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  OpenFile(1);

  // Create histograms
  fHistTrackPt = new TH1F("fHistTrackPt", "P_{T} distribution", 1000, 0.0, 100.0);
  fHistCentrality = new TH1F("fHistCentrality","centrality",100,0,100);
  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet eta-phi",900,-1.8,1.8,720,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron deta-dphi",900,-1.8,1.8,720,-1.6,4.8);

  TString name;

  for(Int_t ipta=0; ipta<7; ++ipta){
    name = Form("fHistTrackEtaPhi_%i", ipta);
    fHistTrackEtaPhi[ipta] = new TH2F(name,name,400,-1,1,720,0.0,2.0*TMath::Pi());
    fOutput->Add(fHistTrackEtaPhi[ipta]);
  }
 
  for(Int_t icent = 0; icent<6; ++icent){
    name = Form("fHistJetPt_%i",icent);   
    fHistJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPt[icent]);

    name = Form("fHistJetPtBias_%i",icent);   
    fHistJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPtBias[icent]);

    name = Form("fHistLeadJetPt_%i",icent);   
    fHistLeadJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistLeadJetPt[icent]);

    name = Form("fHistLeadJetPtBias_%i",icent);   
    fHistLeadJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistLeadJetPtBias[icent]);

    name = Form("fHistJetPtTT_%i",icent);   
    fHistJetPtTT[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPtTT[icent]);

    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
	    name = Form("fHistJetH_%i_%i_%i",icent,iptjet,ieta);   
	    fHistJetH[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	    fOutput->Add(fHistJetH[icent][iptjet][ieta]);

       	name = Form("fHistJetHBias_%i_%i_%i",icent,iptjet,ieta);   
	    fHistJetHBias[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	    fOutput->Add(fHistJetHBias[icent][iptjet][ieta]);

	    name = Form("fHistJetHTT_%i_%i_%i",icent,iptjet,ieta);   
	    fHistJetHTT[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
	    fOutput->Add(fHistJetHTT[icent][iptjet][ieta]);

      }
    }
  }

  UInt_t cifras = 0; // bit coded, see GetDimParams() below 
  if(fDoLessSparseAxes) {
    cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5;
  } else {
    cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; 
  }
  fhnJH = NewTHnSparseF("fhnJH", cifras);
  fhnJH->Sumw2();
  fOutput->Add(fhnJH);

  if(fDoEventMixing){    
    if(fDoLessSparseAxes) { 
      cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5;
    } else {
      cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; 
    }
    fhnMixedEvents = NewTHnSparseF("fhnMixedEvents", cifras);
    fhnMixedEvents->Sumw2();
    fOutput->Add(fhnMixedEvents);
  }
  
  fOutput->Add(fHistTrackPt);
  fOutput->Add(fHistCentrality);
  fOutput->Add(fHistJetEtaPhi);
  fOutput->Add(fHistJetHEtaPhi);

  PostData(1, fOutput);

  //Event Mixing
  Int_t trackDepth = fMixingTracks; 
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
 
  Int_t nZvtxBins  = 5+1+5;
  // bins for second buffer are shifted by 100 cm
  Double_t vertexBins[] = { -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, };
  Double_t* zvtxbin = vertexBins;

//  Int_t nCentralityBins  = 100;
  Int_t nCentralityBins = 100;
  Double_t mult = 1.0;
  if(fCentBinSize==1) { 
    nCentralityBins = 100;
    mult = 1.0;  
  } else if(fCentBinSize==2){
    nCentralityBins = 50;
    mult = 2.0;
  } else if(fCentBinSize==5){
    nCentralityBins = 20;
    mult = 5.0;
  } else if(fCentBinSize==10){
    nCentralityBins = 10;
    mult = 10.0;
  }
  Double_t centralityBins[nCentralityBins];
  for(Int_t ic=0; ic<nCentralityBins; ic++){
   centralityBins[ic]=mult*ic;
  }

  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);

}

//________________________________________________________________________

Double_t AliAnalysisTaskEmcalJetHMEC:: RelativePhi(Double_t mphi,Double_t vphi) {

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  Double_t dphi = vphi-mphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());

  return dphi;//dphi in [-Pi, Pi]                                                                                                    
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetCentBin(Double_t cent) const {
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10) centbin = 0;
  else if (cent>=10 && cent<20) centbin = 1;
  else if (cent>=20 && cent<30) centbin = 2;
  else if (cent>=30 && cent<40) centbin = 3;
  else if (cent>=40 && cent<50) centbin = 4;
  else if (cent>=50 && cent<90) centbin = 5;
  return centbin;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetEtaBin(Double_t eta) const {
  // Get eta bin for histos.

  Int_t etabin = -1;
  if (TMath::Abs(eta)<=0.4) etabin = 0;
  else if (TMath::Abs(eta)>0.4 && TMath::Abs(eta)<0.8) etabin = 1;
  else if (TMath::Abs(eta)>=0.8) etabin = 2;
  return etabin;
}
//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetpTjetBin(Double_t pt) const 
{
  // Get jet pt  bin for histos.

  Int_t ptbin = -1;
  if (pt>=15 && pt<20) ptbin = 0;
  else if (pt>=20 && pt<25) ptbin = 1;
  else if (pt>=25 && pt<30) ptbin = 2;
  else if (pt>=30 && pt<60) ptbin = 3;
  else if (pt>=60) ptbin = 4;

  return ptbin;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::ExecOnce() {
  AliAnalysisTaskEmcalJet::ExecOnce();

}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetHMEC::Run() {
 // Main loop called for each event
  if(!fTracks){
    AliError(Form("No fTracks object!!\n"));
    return kTRUE;
  }
  if(!fJets){
    AliError(Form("No fJets object!!\n"));
    return kTRUE;
  }

  // what kind of event do we have: AOD or ESD?
  Bool_t esdMode = kTRUE; 
  if (dynamic_cast<AliAODEvent*>(InputEvent())) esdMode = kFALSE;

  // if we have ESD event, set up ESD object
  if(esdMode){
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {
      AliError(Form("ERROR: fESD not available\n"));
      return kTRUE;
    }
  }

  // if we have AOD event, set up AOD object
  if(!esdMode){
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) {
      AliError(Form("ERROR: fAOD not available\n"));
      return kTRUE;
    }
  }

  TList *list = InputEvent()->GetList(); 
  if(!list) {
    AliError(Form("ERROR: list not attached\n"));
    return kTRUE;
  }

  // get centrality
  if (fCent<0) {
    AliError(Form("Centrality negative: %f", fCent));
    return kTRUE;
  }

  Int_t centbin = GetCentBin(fCent);
  if(centbin<0) return kTRUE;

  Double_t fvertex[3]={0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);
  Double_t zVtx=fvertex[2];

  if(fabs(zVtx)>10.0) return kTRUE;

  fHistCentrality->Fill(fCent);
    
  TClonesArray *jets = 0;
  TClonesArray *tracks = 0;

  tracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracks));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracks->GetName() ));
    return kTRUE;
  }
  const Int_t Ntracks=tracks->GetEntries();

  jets= dynamic_cast<TClonesArray*>(list->FindObject(fJets));
  if (!jets) {
    AliError(Form("Pointer to tracks %s == 0", fJets->GetName() ));
    return kTRUE;
  }
  const Int_t Njets = jets->GetEntries();
 
  //Leticia's loop to find hardest track
  Int_t iTT=-1;
  Double_t ptmax=-10;

  for (Int_t iTracks = 0; iTracks < Ntracks; iTracks++) {
    AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
      
    if(TMath::Abs(track->Eta())>0.9) continue;
    if(track->Pt()<0.15) continue;
    //iCount++;
    if(track->Pt()>ptmax){
	  ptmax=track->Pt();
	  iTT=iTracks;
    }
  }

  Int_t ijethi=-1;
  Double_t highestjetpt=0.0;
  Int_t passedTTcut=0;

  for (Int_t ijets = 0; ijets < Njets; ijets++){
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijets));
      if (!jet) continue;
      if(!AcceptthisJet(jet)) continue;

      Double_t jetPt = jet->Pt();

      if(highestjetpt<jetPt){
	    ijethi=ijets;
	    highestjetpt=jetPt;
      }
  }

  // see if event is selected
  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  for (Int_t ijet = 0; ijet < Njets; ijet++){    
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));
    if (!jet) continue; 

    // see if event is selected and our jet came from trigger event
    if (!(trig & fTriggerEventType)) continue;
    if (jet->Pt()<0.1) continue;

    if(!AcceptthisJet(jet)) continue;

    Double_t jetphi = jet->Phi();
    Double_t jetPt = jet->Pt();
    Double_t jeteta=jet->Eta();

    Double_t leadjet=0;
    if (ijet==ijethi) leadjet=1;

    fHistJetPt[centbin]->Fill(jet->Pt());
    fHistLeadJetPt[centbin]->Fill(jet->Pt());
    
    if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
      fHistJetPtBias[centbin]->Fill(jet->Pt());
      fHistLeadJetPtBias[centbin]->Fill(jet->Pt());
    }

    fHistJetEtaPhi->Fill(jet->Eta(),jetphi);

    if(iTT>0){
	  AliVTrack* TT = static_cast<AliVTrack*>(tracks->At(iTT));
	  if(TMath::Abs(jetphi-TT->Phi()-TMath::Pi())<0.6) passedTTcut=1;
	  else passedTTcut=0;
    }

    if(passedTTcut)
	  fHistJetPtTT[centbin]->Fill(jet->Pt());

      Int_t iptjet=-1;
      iptjet=GetpTjetBin(jetPt);
      if(iptjet<0) continue;

    if (jetPt > 15) {

      for (Int_t iTracks = 0; iTracks < Ntracks; iTracks++) {
	    AliVTrack* track = static_cast<AliVTrack*>(tracks->At(iTracks));
	    if (!track) {
	      printf("ERROR: Could not receive track %d\n", iTracks);
	      continue;
	    }
	
	    if(TMath::Abs(track->Eta())>fTrkEta) continue;

	    fHistTrackPt->Fill(track->Pt());
	  	  
	    if (track->Pt()<0.15) continue;
	  
	    Double_t trackphi = track->Phi();
	    if (trackphi > TMath::Pi())
	      trackphi = trackphi-2*TMath::Pi();

	    Double_t tracketa=track->Eta();
	    Double_t trackpt=track->Pt();
	    Double_t deta=tracketa-jeteta;
	    Int_t ieta=GetEtaBin(deta);
	    if (ieta<0) {
	      AliError(Form("Eta Bin negative: %f", deta));
	      continue;
	    }
      
	    //Jet pt, track pt, dPhi,deta,fCent
	    Double_t dphijh = RelativePhi(jetphi,trackphi);
	    if (dphijh < -0.5*TMath::Pi())
	      dphijh+= 2*TMath::Pi();
	    if (dphijh > 1.5*TMath::Pi()) 
          dphijh-=2.*TMath::Pi();

  	    fHistJetH[centbin][iptjet][ieta]->Fill(dphijh,track->Pt());
	    fHistJetHEtaPhi->Fill(deta,dphijh);

	    Double_t dR=sqrt(deta*deta+dphijh*dphijh);

	    if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
	      fHistJetHBias[centbin][iptjet][ieta]->Fill(dphijh,trackpt);

          if(fDoLessSparseAxes) { // check if we want all dimensions
            Double_t triggerEntries[6] = {fCent,jetPt,trackpt,deta,dphijh,leadjet};
            fhnJH->Fill(triggerEntries);
          } else { 
	        Double_t triggerEntries[8] = {fCent,jetPt,trackpt,deta,dphijh,leadjet,0.0,dR};                      
            fhnJH->Fill(triggerEntries);
          }
	    }

	    if(passedTTcut)
	      fHistJetHTT[centbin][iptjet][ieta]->Fill(dphijh,trackpt);

	  } //track loop
    }//jet pt cut
  }//jet loop
  
  //Prepare to do event mixing

  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = 0x0;

  if(fDoEventMixing>0){
    
    // event mixing
    
    // 1. First get an event pool corresponding in mult (cent) and
    //    zvertex to the current event. Once initialized, the pool
    //    should contain nMix (reduced) events. This routine does not
    //    pre-scan the chain. The first several events of every chain
    //    will be skipped until the needed pools are filled to the
    //    specified depth. If the pool categories are not too rare, this
    //    should not be a problem. If they are rare, you could lose
    //    statistics.

    // 2. Collect the whole pool's content of tracks into one TObjArray
    //    (bgTracks), which is effectively a single background super-event.

    // 3. The reduced and bgTracks arrays must both be passed into
    //    FillCorrelations(). Also nMix should be passed in, so a weight
    //    of 1./nMix can be applied.

    UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    // if event was not selected (triggered) for any reseason (should never happen) then return
    if (trigger==0)  return kTRUE;

    AliEventPool* pool = fPoolMgr->GetEventPool(fCent, zVtx);
    
    if (!pool){
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCent, zVtx));
      return kTRUE;
    }

if(trigger & fTriggerEventType) {
    //check for a trigger jet
    if (pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {

	for (Int_t ijet = 0; ijet < Njets; ijet++){
	  Double_t leadjet=0;
	  if (ijet==ijethi) leadjet=1;
	  
	  AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));
      if (!jet) continue;

	  if(!AcceptthisJet(jet)) continue;

	  Double_t jetPt = jet->Pt();	
	  Double_t jetphi = jet->Phi();
	  Double_t jeteta=jet->Eta();

      // make sure event contains jet above our threshold (reduce stats of sparse)
      if (jetPt < 15) continue;

	  Int_t nMix = pool->GetCurrentNEvents();
	  
	  //Fill for biased jet triggers only
	  if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){

	    // Fill mixed-event histos here  
	    for (Int_t jMix=0; jMix<nMix; jMix++) {
  		  TObjArray* bgTracks = pool->GetEvent(jMix);
		  const Int_t Nbgtrks = bgTracks->GetEntries();
		  
          for(Int_t ibg=0; ibg<Nbgtrks; ibg++){
		    AliPicoTrack *part = static_cast<AliPicoTrack*>(bgTracks->At(ibg));         
		    if(!part) continue;
            if(TMath::Abs(part->Eta())>0.9) continue;
            if(part->Pt()<0.15) continue;
  
		    Double_t DEta = part->Eta()-jeteta;
		    Double_t DPhi = RelativePhi(jetphi,part->Phi());

		    Double_t DR=TMath::Sqrt(DPhi*DPhi+DEta*DEta);
		    if(DPhi<-0.5*TMath::Pi()) DPhi+=2.*TMath::Pi();
		    if(DPhi>3./2.*TMath::Pi()) DPhi-=2.*TMath::Pi();
            if(fDoLessSparseAxes) {  // check if we want all the axis filled
              Double_t triggerEntries[6] = {fCent,jetPt,part->Pt(),DEta,DPhi,leadjet};
              fhnMixedEvents->Fill(triggerEntries,1./nMix);
            } else {
		      Double_t triggerEntries[8] = {fCent,jetPt,part->Pt(),DEta,DPhi,leadjet,0.0,DR};                      
              fhnMixedEvents->Fill(triggerEntries,1./nMix);
		    }
		    }
	      }
	    }
	  }
    }
}

    if(trigger & fMixingEventType) {
      tracksClone = CloneAndReduceTrackList(tracks);

      //update pool if jet in event or not
      pool->UpdatePool(tracksClone);
    }

  } // end of event mixing

  return kTRUE;
}      

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::Terminate(Option_t *) 
{
  //just terminate
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::AcceptthisJet(AliEmcalJet *jet) 
{
  //applies all jet cuts except pt
  float jetphi = jet->Phi();
  if (jetphi>TMath::Pi())
    jetphi = jetphi-2*TMath::Pi();

  if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
    return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
    return 0;
  if (jet->Area()<fAreacut)
    return 0;
  //prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0)
    return 0;  
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt()>100)
    return 0;

  //passed all above cuts
  return 1;
}

//________________________________________________________________________

THnSparse* AliAnalysisTaskEmcalJetHMEC::NewTHnSparseF(const char* name, UInt_t entries){
   // generate new THnSparseF, axes are defined in GetDimParams()

   Int_t count = 0;
   UInt_t tmp = entries;
   while(tmp!=0){
      count++;
      tmp = tmp &~ -tmp;  // clear lowest bit
   }

   TString hnTitle(name);
   const Int_t dim = count;
   Int_t nbins[dim];
   Double_t xmin[dim];
   Double_t xmax[dim];

   Int_t i=0;
   Int_t c=0;
   while(c<dim && i<32){
      if(entries&(1<<i)){
      
         TString label("");
         GetDimParams(i, label, nbins[c], xmin[c], xmax[c]);
         hnTitle += Form(";%s",label.Data());
         c++;
      }
      
      i++;
   }
   hnTitle += ";";

   return new THnSparseF(name, hnTitle.Data(), dim, nbins, xmin, xmax);
}

void AliAnalysisTaskEmcalJetHMEC::GetDimParams(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
   // stores label and binning of axis for THnSparse

   const Double_t pi = TMath::Pi();
   
   switch(iEntry){
      
   case 0:
      label = "V0 centrality (%)";
         nbins = 10;
         xmin = 0.;
         xmax = 100.;
         break;
      
   case 1:
      label = "corrected jet pt";
      nbins = 20;
      xmin = 0.;
      xmax = 200.;
      break;
      
   case 2:
      if(fDoWiderTrackBin) {
        label = "track pT";
        nbins = 40;
        xmin = 0.;
        xmax = 10.;
      } else {
        label = "track pT";
        nbins = 100;
        xmin = 0.;
        xmax = 10;
      }
      break;
      
   case 3:
      label = "deltaEta";
      nbins = 24;
      xmin = -1.2;
      xmax = 1.2;
      break;

   case 4:
      label = "deltaPhi";
      nbins = 72;
      xmin = -0.5*pi;
      xmax = 1.5*pi;
      break;         

    case 5:
      label = "leading jet";
      nbins = 3;
      xmin = -0.5;
      xmax = 2.5;
      break;
           
    case 6:
      label = "trigger track";
      nbins =10;
      xmin = 0;
      xmax = 50;
      break;

    case 7:
      label = "deltaR";
      nbins = 10;
      xmin = 0.;
      xmax = 5.0;
      break;

    case 8:
      label = "leading track";
      nbins = 13;
      xmin = 0;
      xmax = 50;
      break;
   }
}

//_________________________________________________
// From CF event mixing code PhiCorrelations
TObjArray* AliAnalysisTaskEmcalJetHMEC::CloneAndReduceTrackList(TObjArray* tracks)
{
  // clones a track list by using AliPicoTrack which uses much less memory (used for event mixing)
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<tracks->GetEntriesFast(); i++)
  {
    AliVParticle* particle = (AliVParticle*) tracks->At(i);
    if(TMath::Abs(particle->Eta())>fTrkEta) continue;
    if(particle->Pt()<0.15) continue;

    Double_t trackpt=particle->Pt();

    Int_t hadbin=-1;
    if(trackpt<0.5) hadbin=0;
    else if(trackpt<1) hadbin=1;
    else if(trackpt<2) hadbin=2;
    else if(trackpt<3) hadbin=3;
    else if(trackpt<5) hadbin=4;
    else if(trackpt<8) hadbin=5;
    else if(trackpt<20) hadbin=6;

    if(hadbin>-1) fHistTrackEtaPhi[hadbin]->Fill(particle->Eta(),particle->Phi());

    tracksClone->Add(new AliPicoTrack(particle->Pt(), particle->Eta(), particle->Phi(), particle->Charge(), 0, 0, 0, 0));
  }
  
  return tracksClone;
}
