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
#include "TH3F.h"
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
#include "AliESDVertex.h"
#include "AliCentrality.h"
#include "AliAODJet.h"
#include "AliEmcalJet.h"
#include "AliESDtrackCuts.h"

#include "AliClusterContainer.h"
#include "AliTrackContainer.h"

#include "TVector3.h"
#include "AliPicoTrack.h"
#include "AliEventPoolManager.h"

ClassImp(AliAnalysisTaskEmcalJetHMEC)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC() : 
  AliAnalysisTaskEmcalJet("HMEC",kFALSE),
  /*fTracksName(""),
  fJetsName(""),*/
  /*fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),*/
  fTrkBias(5),
  fClusBias(5),
  fTrkEta(0.9),
  fDoEventMixing(0),
  fMixingTracks(50000), fNMIXtracks(5000), fNMIXevents(5),
  fTriggerEventType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDoEffCorrection(0), fEffFunctionCorrection(0),
  fEmbeddingCorrectionHist(0),
  fDoLessSparseAxes(0), fDoWiderTrackBin(0),
  fCentBinSize(1),
  //fESD(0),
  //fAOD(0), 
  fPoolMgr(0x0), 
  fHistTrackPt(0),
  //fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistClusEtaPhiEn(0), 
  fHistJHPsi(0),
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
  /*fTracksName(""),
  fJetsName(""),*/
  /*fPhimin(-10), 
  fPhimax(10),
  fEtamin(-0.9), 
  fEtamax(0.9),
  fAreacut(0.0),*/
  fTrkBias(5),
  fClusBias(5),
  fTrkEta(0.9),
  fDoEventMixing(0),
  fMixingTracks(50000), fNMIXtracks(5000), fNMIXevents(5),
  fTriggerEventType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDoEffCorrection(0), fEffFunctionCorrection(0),
  fEmbeddingCorrectionHist(0),
  fDoLessSparseAxes(0), fDoWiderTrackBin(0),
  fCentBinSize(1),
  //fESD(0),
  //fAOD(0), 
  fPoolMgr(0x0), 
  fHistTrackPt(0),
  //fHistCentrality(0), 
  fHistJetEtaPhi(0), 
  fHistClusEtaPhiEn(0),  
  fHistJHPsi(0),
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
  //fHistCentrality = new TH1F("fHistCentrality","centrality",100,0,100);
  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet eta-phi",900,-1.8,1.8,720,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron deta-dphi",900,-1.8,1.8,720,-1.6,4.8);

  fHistClusEtaPhiEn = new TH3F("fHistClusEtaPhiEn","Clus eta-phi-energy",900,-1.8,1.8,720,-3.2,3.2,20,0,20);

  fHistJHPsi = new TH3F("fHistJHPsi","Jet-Hadron ntr-trpt-dpsi",20,0,100,200,0,20,120,0,180);

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
  //fOutput->Add(fHistCentrality);
  fOutput->Add(fHistJetEtaPhi);
  fOutput->Add(fHistJetHEtaPhi);
  fOutput->Add(fHistClusEtaPhiEn);
  fOutput->Add(fHistJHPsi);

  PostData(1, fOutput);

  //Event Mixing
  Int_t trackDepth = fMixingTracks; 
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager

  Int_t nZvtxBins  = 10;
  // bins for second buffer are shifted by 100 cm
  Double_t vertexBins[] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};
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

  // for pp we need mult bins for event mixing. Create binning here, to also make a histogram from it
  Int_t nCentralityBins_pp  = 8;
  Double_t centralityBins_pp[9] = {0., 4., 9., 15., 25., 35., 55., 100., 500.};

  if (fForceBeamType != kpp ) {   //all besides pp
    fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);
  }
  else if (fForceBeamType == kpp) { //for pp only
    fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins_pp, centralityBins_pp, nZvtxBins, zvtxbin);
  }

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
/*Int_t AliAnalysisTaskEmcalJetHMEC::GetCentBin(Double_t cent) const {
  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10) centbin = 0;
  else if (cent>=10 && cent<20) centbin = 1;
  else if (cent>=20 && cent<30) centbin = 2;
  else if (cent>=30 && cent<40) centbin = 3;
  else if (cent>=40 && cent<50) centbin = 4;
  else if (cent>=50 && cent<90) centbin = 5;
  return centbin;
}*/

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
  /*if(!fTracks){
    AliError(Form("No fTracks object!!\n"));
    return kTRUE;
  }
  if(!fJets){
    AliError(Form("No fJets object!!\n"));
    return kTRUE;
  }*/

  // what kind of event do we have: AOD or ESD?
  /*Bool_t esdMode = kTRUE; 
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
  }*/

  // TODO: This should be covered by AliAnalysisTaskEmcal
  // get centrality
  if (fCent<0) {
    AliError(Form("Centrality negative: %f", fCent));
    return kTRUE;
  }

  // TODO: Can we just use the fCent = 99, fCentBin = 0 def from AliAnalysisTaskEmcal?
  if (fBeamType == kpp) fCent = 0.0; //put pp centrality to 0.0

  // TODO: Should also be covered here
  //Int_t centbin = GetCentBin(fCent);
  if(fCentBin<0) return kTRUE;

  // Filled by AliAnalysisTaskEmcal
  //Double_t fvertex[3]={0,0,0};
  //InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);
  Double_t zVtx=fVertex[2];

  // Replaced by vz range call
  //if(fabs(zVtx)>10.0) return kTRUE;

  // This is taken care of by hists in AliAnalysisTaskEmcal
  //fHistCentrality->Fill(fCent);

  AliClusterContainer * clusters = GetClusterContainer(0);
  if (!clusters) {
    AliError("Unable to retreive clusters!");
    return kTRUE;
  }

  // TODO: Can we use a QA task for this?
  AliVCluster * cluster = 0;
  TLorentzVector nPart;
  while ((cluster = clusters->GetNextAcceptCluster()))
  {
    cluster->GetMomentum(nPart, fVertex);
    fHistClusEtaPhiEn->Fill( nPart.Eta(), nPart.Phi(), nPart.E() );
  }
  // Reset so that we can iterate again in the future.
  clusters->ResetCurrentID();

  AliTrackContainer * tracks = GetTrackContainer(0);
  if (!tracks) {
    AliError("Unable to retreive tracks!");
    return kTRUE;
  }
  const Int_t Ntracks = tracks->GetNTracks();

  Int_t passedTTcut=0;

  AliJetContainer * jets = GetJetContainer(0);
  if (!jets) {
    AliError("Unable to retreive jets!");
    return kTRUE;
  }

  // TODO: I think this is covered in AliAnalysisTaskEmcal
  // see if event is selected
  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  TVector3 vector_jet, vector_hdr;

  // For comparison below
  AliEmcalJet * leadingJet = jets->GetLeadingJet();

  // Just to be certain that we are interating from the start
  jets->ResetCurrentID();

  AliEmcalJet * jet = 0;
  while ((jet = jets->GetNextAcceptJet())) {
    
    // see if event is selected and our jet came from trigger event
    if (!(trig & fTriggerEventType)) continue;
    // TODO: Can we just move this the jet container?
    //if (jet->Pt()<0.1) continue;

    // TODO: Can this be removed?
    if(!AcceptthisJet(jet)) continue;

    Double_t jetphi = jet->Phi();
    Double_t jetPt = jet->Pt();

    vector_jet.SetXYZ( jet->Px(), jet->Py(), jet->Pz() );

    Double_t leadjet=0;
    if (jet == leadingJet) leadjet=1;

    FillHist(fHistJetPt[fCentBin], jet->Pt());
    FillHist(fHistLeadJetPt[fCentBin], jet->Pt());

    if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
      FillHist(fHistJetPtBias[fCentBin], jet->Pt());
      FillHist(fHistLeadJetPtBias[fCentBin], jet->Pt());
    }

    fHistJetEtaPhi->Fill(jet->Eta(),jetphi);

    // TODO: Would we be better off taking the tracks in the Jet?
    // TODO: Discuss this cut! I don't quite understand it
    AliVTrack * leadingTrack = 0;
    leadingTrack = tracks->GetLeadingTrack();
    if (leadingTrack != 0)
    {
      if(TMath::Abs(jetphi-leadingTrack->Phi()-TMath::Pi())<0.6) passedTTcut=1;
      else passedTTcut=0;
    }

    if(passedTTcut)
      FillHist(fHistJetPtTT[fCentBin], jet->Pt());

    Int_t iptjet=-1;
    iptjet=GetpTjetBin(jetPt);
    if(iptjet<0) continue;

    if (jetPt > 15) {

      AliVTrack * track = 0;
      while ((track = tracks->GetNextAcceptTrack())) {

        // TODO: We should be able to remove this by apply the track eta cut in the run macro
        if(TMath::Abs(track->Eta())>fTrkEta) continue;

        fHistTrackPt->Fill(track->Pt());

        // TODO: Need a second track container - one for jets and one for hadrons
        if (track->Pt()<0.15) continue;

        vector_hdr.SetXYZ( track->Px(), track->Py(), track->Pz() );
        if ( (jetPt>20.) && (jetPt<60.) ) {
          fHistJHPsi->Fill(Ntracks, track->Pt(), vector_hdr.Angle(vector_jet) * TMath::RadToDeg() );
        }

        Double_t trackphi = track->Phi();
        if (trackphi > TMath::Pi())
          trackphi = trackphi-2*TMath::Pi();

        Double_t tracketa=track->Eta();
        Double_t trackpt=track->Pt();
        Double_t deta=tracketa-jet->Eta();
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

        fHistJetH[fCentBin][iptjet][ieta]->Fill(dphijh,track->Pt());
        fHistJetHEtaPhi->Fill(deta,dphijh);

        // calculate single particle tracking efficiency for correlations
        Double_t trefficiency = -999;
        trefficiency = EffCorrection(track->Eta(), track->Pt(), fDoEffCorrection);

        Double_t dR=sqrt(deta*deta+dphijh*dphijh);

        if ((jet->MaxTrackPt()>fTrkBias) || (jet->MaxClusterPt()>fClusBias)){
          fHistJetHBias[fCentBin][iptjet][ieta]->Fill(dphijh,trackpt);

          Double_t eventActivity = 0.0;
          if (fBeamType == kAA || fBeamType == kpA) { //pA and AA
            eventActivity = fCent;
          }
          else if (fBeamType == kpp) {
            eventActivity = static_cast<Double_t>(Ntracks);
          }

          if(fDoLessSparseAxes) { // check if we want all dimensions
            Double_t triggerEntries[6] = {eventActivity,jetPt,trackpt,deta,dphijh,leadjet};
            FillHist(fhnJH, triggerEntries, 1.0/trefficiency);
          } else { 
            Double_t triggerEntries[8] = {eventActivity,jetPt,trackpt,deta,dphijh,leadjet,0.0,dR};
            FillHist(fhnJH, triggerEntries, 1.0/trefficiency);
          }
        }

        if(passedTTcut)
          fHistJetHTT[fCentBin][iptjet][ieta]->Fill(dphijh,trackpt);

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

    AliEventPool *pool = 0;
    if (fBeamType == kAA || fBeamType == kpA) {//everything but pp
      pool = fPoolMgr->GetEventPool(fCent, zVtx);
    }
    else if (fBeamType == kpp) {//pp only
      Double_t Ntrks = (Double_t)Ntracks*1.0;
      pool = fPoolMgr->GetEventPool(Ntrks, zVtx);
    }

    if (!pool){
      if (fBeamType == kAA || fBeamType == kpA) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCent, zVtx));
      else if (fBeamType == kpp) AliFatal(Form("No pool found for ntracks_pp = %d, zVtx = %f", Ntracks, zVtx));
      return kTRUE;
    }

    if(trigger & fTriggerEventType) {
      //check for a trigger jet
      if (pool->IsReady() || pool->NTracksInPool() > fNMIXtracks || pool->GetCurrentNEvents() >= fNMIXevents) {

        jets->ResetCurrentID();
        while ((jet = jets->GetNextAcceptJet())) {
          Double_t leadjet=0;
          if (jet == leadingJet) leadjet=1;

          // TODO: Is this still necessary?
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

                // calculate single particle tracking efficiency of mixed events for correlations
                Double_t mixefficiency = -999;
                mixefficiency = EffCorrection(part->Eta(), part->Pt(), fDoEffCorrection);    

                Double_t DR=TMath::Sqrt(DPhi*DPhi+DEta*DEta);
                if(DPhi<-0.5*TMath::Pi()) DPhi+=2.*TMath::Pi();
                if(DPhi>3./2.*TMath::Pi()) DPhi-=2.*TMath::Pi();

                Double_t eventActivity = 0;
                if (fBeamType == kAA || fBeamType == kpA) { //pA and AA
                  eventActivity = fCent;
                }
                else if (fBeamType == kpp) {
                  eventActivity = static_cast<Double_t>(Ntracks);
                }

                if(fDoLessSparseAxes) {  // check if we want all the axis filled
                  Double_t triggerEntries[6] = {eventActivity,jetPt,part->Pt(),DEta,DPhi,leadjet};
                  FillHist(fhnMixedEvents, triggerEntries, 1./(nMix*mixefficiency), kFALSE);
                } else {
                  Double_t triggerEntries[8] = {eventActivity,jetPt,part->Pt(),DEta,DPhi,leadjet,0.0,DR};
                  FillHist(fhnMixedEvents, triggerEntries, 1./(nMix*mixefficiency), kFALSE);
                }
              }
            }
          }
        }
      }
    }

    if(trigger & fMixingEventType) {
      tracksClone = CloneAndReduceTrackList();

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
  /*float jetphi = jet->Phi();
  if (jetphi>TMath::Pi())
    jetphi = jetphi-2*TMath::Pi();*/

  // Acceptance is taken care of by the jet container
  /*if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax))
    return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax))
    return 0;*/
  /*if (jet->Area()<fAreacut)
    return 0;
  //prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0)
    return 0;  */
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  /*if(jet->MaxTrackPt()>100)
    return 0;*/

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
TObjArray* AliAnalysisTaskEmcalJetHMEC::CloneAndReduceTrackList()
{
  // clones a track list by using AliPicoTrack which uses much less memory (used for event mixing)

  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  AliTrackContainer * tracks = GetTrackContainer(0);
  // Ensure that we start from the beginning
  tracks->ResetCurrentID();
  AliVParticle * particle = 0;
  while ((particle = tracks->GetNextAcceptTrack()))
  {
    // TODO: Probably not necessary since the track container should already apply these cuts. But check it!
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

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetHMEC::EffCorrection(Double_t trackETA, Double_t trackPT, Int_t effSwitch) const {
  // default (current) parameters
  // x-variable = track pt, y-variable = track eta
  Double_t x = trackPT;
  Double_t y = trackETA;
  Double_t TRefficiency = -999;
  Int_t runNUM = fCurrentRunNumber;
  Int_t runSwitchGood = -999;
  //Int_t centbin = -99;

  Double_t etaaxis = 0;
  Double_t ptaxis = 0;

  if(effSwitch < 1) {
    if ((runNUM == 169975 || runNUM == 169981 || runNUM == 170038 || runNUM == 170040 || runNUM == 170083 || runNUM == 170084 || runNUM == 170085 || runNUM == 170088 || runNUM == 170089 || runNUM == 170091 || runNUM == 170152 || runNUM == 170155 || runNUM == 170159 || runNUM == 170163 || runNUM == 170193 || runNUM == 170195 || runNUM == 170203 || runNUM == 170204 || runNUM == 170228 || runNUM == 170230 || runNUM == 170268 || runNUM == 170269 || runNUM == 170270 || runNUM == 170306 || runNUM == 170308 || runNUM == 170309)) runSwitchGood = 0;

    if ((runNUM == 167902 || runNUM == 167903 || runNUM == 167915 || runNUM == 167920 || runNUM == 167987 || runNUM == 167988 || runNUM == 168066 || runNUM == 168068 || runNUM == 168069 || runNUM == 168076 || runNUM == 168104 || runNUM == 168107 || runNUM == 168108 || runNUM == 168115 || runNUM == 168212 || runNUM == 168310 || runNUM == 168311 || runNUM == 168322 || runNUM == 168325 || runNUM == 168341 || runNUM == 168342 || runNUM == 168361 || runNUM == 168362 || runNUM == 168458 || runNUM == 168460 || runNUM == 168461 || runNUM == 168464 || runNUM == 168467 || runNUM == 168511 || runNUM == 168512 || runNUM == 168777 || runNUM == 168826 || runNUM == 168984 || runNUM == 168988 || runNUM == 168992 || runNUM == 169035 || runNUM == 169091 || runNUM == 169094 || runNUM == 169138 || runNUM == 169143 || runNUM == 169144 || runNUM == 169145 || runNUM == 169148 || runNUM == 169156 || runNUM == 169160 || runNUM == 169167 || runNUM == 169238 || runNUM == 169411 || runNUM == 169415 || runNUM == 169417 || runNUM == 169835 || runNUM == 169837 || runNUM == 169838 || runNUM == 169846 || runNUM == 169855 || runNUM == 169858 || runNUM == 169859 || runNUM == 169923 || runNUM == 169956 || runNUM == 170027 || runNUM == 170036 || runNUM == 170081)) runSwitchGood = 1;

    /*if(fCent>=0 && fCent<10) centbin = 0;
    else if (fCent>=10 && fCent<30)	centbin = 1;
    else if (fCent>=30 && fCent<50) centbin = 2;
    else if (fCent>=50 && fCent<90)	centbin = 3;*/

    if(runSwitchGood == 0 && fCentBin == 0) effSwitch = 2;
    if(runSwitchGood == 0 && fCentBin == 1) effSwitch = 3;
    if(runSwitchGood == 0 && fCentBin == 2) effSwitch = 4;
    if(runSwitchGood == 0 && fCentBin == 3) effSwitch = 5;
    if(runSwitchGood == 1 && fCentBin == 0) effSwitch = 6;
    if(runSwitchGood == 1 && fCentBin == 1) effSwitch = 7;
    if(runSwitchGood == 1 && fCentBin == 2) effSwitch = 8;
    if(runSwitchGood == 1 && fCentBin == 3) effSwitch = 9;

  }

  // 0-10% centrality: Semi-Good Runs
  Double_t p0_10SG[17] = {0.906767, 0.0754127, 1.11638, -0.0233078, 0.795454, 0.00935385, -0.000327857, 1.08903, 0.0107272, 0.443252, -0.143411, 0.965822, 0.359156, -0.581221, 1.0739, 0.00632828, 0.706356};
  // 10-30% centrality: Semi-Good Runs
  Double_t p10_30SG[17] = {0.908011, 0.0769254, 1.11912, -0.0249449, 0.741488, 0.0361252, -0.00367954, 1.10424, 0.011472, 0.452059, -0.133282, 0.980633, 0.358222, -0.620256, 1.06871, 0.00564449, 0.753168};
  // 30-50% centrality: Semi-Good Runs
  Double_t p30_50SG[17] = {0.958708, 0.0799197, 1.10817, -0.0357678, 0.75051, 0.0607808, -0.00929713, 0.998801, 0.00692244, 0.615452, -0.0480328, 0.968431, 0.321634, -0.619066, 1.03412, 0.00656201, 0.798666};
  // 50-90% centrality: Semi-Good Runs
  Double_t p50_90SG[17] = {0.944565, 0.0807258, 1.12709, -0.0324746, 0.666452, 0.0842476, -0.00963837, 1.02829, 0.00666852, 0.549625, -0.0603107, 0.981374, 0.309374, -0.619181, 1.05367, 0.005925, 0.744887};

  // 0-10% centrality: Good Runs
  Double_t p0_10G[17] = {0.971679, 0.0767571, 1.13355, -0.0274484, 0.856652, 0.00536795, 3.90795e-05, 1.06889, 0.011007, 0.447046, -0.146626, 0.919777, 0.192601, -0.268515, 1.00243, 0.00620849, 0.709477};
  // 10-30% centrality: Good Runs
  Double_t p10_30G[17] = {0.97929, 0.0776039, 1.12213, -0.0300645, 0.844722, 0.0134788, -0.0012333, 1.07955, 0.0116835, 0.456608, -0.132743, 0.930964, 0.174175, -0.267154, 0.993118, 0.00574892, 0.765256};
  // 30-50% centrality: Good Runs
  Double_t p30_50G[17] = {0.997696, 0.0816769, 1.14341, -0.0353734, 0.752151, 0.0744259, -0.0102926, 1.01561, 0.00713274, 0.57203, -0.0640248, 0.947747, 0.102007, -0.194698, 0.999164, 0.00568476, 0.7237};
  // 50-90% centrality: Good Runs
  Double_t p50_90G[17] = {0.97041, 0.0813559, 1.12151, -0.0368797, 0.709327, 0.0701501, -0.00784043, 1.06276, 0.00676173, 0.53607, -0.0703117, 0.982534, 0.0947881, -0.18073, 1.03229, 0.00580109, 0.737801};

  // set up a switch for different parameter values...
  switch(effSwitch) {
    case 1 :
      // first switch value - TRefficiency not used so = 1
      TRefficiency = 1.0;
      break;

    case 2 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (0-10%):
      ptaxis = (x<2.9)*(p0_10SG[0]*exp(-pow(p0_10SG[1]/x,p0_10SG[2])) + p0_10SG[3]*x) + (x>=2.9)*(p0_10SG[4] + p0_10SG[5]*x + p0_10SG[6]*x*x);
      etaaxis = (y<-0.07)*(p0_10SG[7]*exp(-pow(p0_10SG[8]/TMath::Abs(y+0.91),p0_10SG[9])) + p0_10SG[10]*y) + (y>=-0.07 && y<=0.4)*(p0_10SG[11] + p0_10SG[12]*y + p0_10SG[13]*y*y) + (y>0.4)*(p0_10SG[14]*exp(-pow(p0_10SG[15]/TMath::Abs(-y+0.91),p0_10SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 3 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (10-30%):
      ptaxis = (x<2.9)*(p10_30SG[0]*exp(-pow(p10_30SG[1]/x,p10_30SG[2])) + p10_30SG[3]*x) + (x>=2.9)*(p10_30SG[4] + p10_30SG[5]*x + p10_30SG[6]*x*x);
      etaaxis = (y<-0.07)*(p10_30SG[7]*exp(-pow(p10_30SG[8]/TMath::Abs(y+0.91),p10_30SG[9])) + p10_30SG[10]*y) + (y>=-0.07 && y<=0.4)*(p10_30SG[11] + p10_30SG[12]*y + p10_30SG[13]*y*y) + (y>0.4)*(p10_30SG[14]*exp(-pow(p10_30SG[15]/TMath::Abs(-y+0.91),p10_30SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 4 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (30-50%):
      ptaxis = (x<2.9)*(p30_50SG[0]*exp(-pow(p30_50SG[1]/x,p30_50SG[2])) + p30_50SG[3]*x) + (x>=2.9)*(p30_50SG[4] + p30_50SG[5]*x + p30_50SG[6]*x*x);
      etaaxis = (y<-0.07)*(p30_50SG[7]*exp(-pow(p30_50SG[8]/TMath::Abs(y+0.91),p30_50SG[9])) + p30_50SG[10]*y) + (y>=-0.07 && y<=0.4)*(p30_50SG[11] + p30_50SG[12]*y + p30_50SG[13]*y*y) + (y>0.4)*(p30_50SG[14]*exp(-pow(p30_50SG[15]/TMath::Abs(-y+0.91),p30_50SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 5 :
      // Parameter values for Semi-GOOD TPC (LHC11h) runs (50-90%):
      ptaxis = (x<2.9)*(p50_90SG[0]*exp(-pow(p50_90SG[1]/x,p50_90SG[2])) + p50_90SG[3]*x) + (x>=2.9)*(p50_90SG[4] + p50_90SG[5]*x + p50_90SG[6]*x*x);
      etaaxis = (y<-0.07)*(p50_90SG[7]*exp(-pow(p50_90SG[8]/TMath::Abs(y+0.91),p50_90SG[9])) + p50_90SG[10]*y) + (y>=-0.07 && y<=0.4)*(p50_90SG[11] + p50_90SG[12]*y + p50_90SG[13]*y*y) + (y>0.4)*(p50_90SG[14]*exp(-pow(p50_90SG[15]/TMath::Abs(-y+0.91),p50_90SG[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 6 :
      // Parameter values for GOOD TPC (LHC11h) runs (0-10%):
      ptaxis = (x<2.9)*(p0_10G[0]*exp(-pow(p0_10G[1]/x,p0_10G[2])) + p0_10G[3]*x) + (x>=2.9)*(p0_10G[4] + p0_10G[5]*x + p0_10G[6]*x*x);
      etaaxis = (y<0.0)*(p0_10G[7]*exp(-pow(p0_10G[8]/TMath::Abs(y+0.91),p0_10G[9])) + p0_10G[10]*y) + (y>=0.0 && y<=0.4)*(p0_10G[11] + p0_10G[12]*y + p0_10G[13]*y*y) + (y>0.4)*(p0_10G[14]*exp(-pow(p0_10G[15]/TMath::Abs(-y+0.91),p0_10G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 7 :
      // Parameter values for GOOD TPC (LHC11h) runs (10-30%):
      ptaxis = (x<2.9)*(p10_30G[0]*exp(-pow(p10_30G[1]/x,p10_30G[2])) + p10_30G[3]*x) + (x>=2.9)*(p10_30G[4] + p10_30G[5]*x + p10_30G[6]*x*x);
      etaaxis = (y<0.0)*(p10_30G[7]*exp(-pow(p10_30G[8]/TMath::Abs(y+0.91),p10_30G[9])) + p10_30G[10]*y) + (y>=0.0 && y<=0.4)*(p10_30G[11] + p10_30G[12]*y + p10_30G[13]*y*y) + (y>0.4)*(p10_30G[14]*exp(-pow(p10_30G[15]/TMath::Abs(-y+0.91),p10_30G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 8 :
      // Parameter values for GOOD TPC (LHC11h) runs (30-50%):
      ptaxis = (x<2.9)*(p30_50G[0]*exp(-pow(p30_50G[1]/x,p30_50G[2])) + p30_50G[3]*x) + (x>=2.9)*(p30_50G[4] + p30_50G[5]*x + p30_50G[6]*x*x);
      etaaxis = (y<0.0)*(p30_50G[7]*exp(-pow(p30_50G[8]/TMath::Abs(y+0.91),p30_50G[9])) + p30_50G[10]*y) + (y>=0.0 && y<=0.4)*(p30_50G[11] + p30_50G[12]*y + p30_50G[13]*y*y) + (y>0.4)*(p30_50G[14]*exp(-pow(p30_50G[15]/TMath::Abs(-y+0.91),p30_50G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    case 9 :
      // Parameter values for GOOD TPC (LHC11h) runs (50-90%):
      ptaxis = (x<2.9)*(p50_90G[0]*exp(-pow(p50_90G[1]/x,p50_90G[2])) + p50_90G[3]*x) + (x>=2.9)*(p50_90G[4] + p50_90G[5]*x + p50_90G[6]*x*x);
      etaaxis = (y<0.0)*(p50_90G[7]*exp(-pow(p50_90G[8]/TMath::Abs(y+0.91),p50_90G[9])) + p50_90G[10]*y) + (y>=0.0 && y<=0.4)*(p50_90G[11] + p50_90G[12]*y + p50_90G[13]*y*y) + (y>0.4)*(p50_90G[14]*exp(-pow(p50_90G[15]/TMath::Abs(-y+0.91),p50_90G[16])));
      TRefficiency = ptaxis*etaaxis;
      break;

    default :
      // no Efficiency Switch option selected.. therefore don't correct, and set eff = 1
      TRefficiency = 1.0;

  }

  return TRefficiency;
}

void AliAnalysisTaskEmcalJetHMEC::FillHist(TH1 * hist, Double_t fillValue, Double_t weight, Bool_t noCorrection)
{
  if (fEmbeddingCorrectionHist == 0 || noCorrection == kTRUE)
  {
    hist->Fill(fillValue, weight);
  }
  else
  {
    // Determine where to get the values in the correction hist
    Int_t xBin = fEmbeddingCorrectionHist->GetXaxis()->FindBin(fillValue);

    std::vector <Double_t> yBinsContent;
    accessSetOfYBinValues(fEmbeddingCorrectionHist, xBin, yBinsContent);

    // Loop over all possible bins to contribute.
    // If content is 0 then calling Fill won't make a difference
    for (Int_t index = 1; index <= fEmbeddingCorrectionHist->GetYaxis()->GetNbins(); index++)
    {
      // Determine the value to fill based on the center of the bins.
      // This in principle allows the binning between the correction and hist to be different
      Double_t fillLocation = fEmbeddingCorrectionHist->GetYaxis()->GetBinCenter(index); 
      Printf("fillLocation: %f, weight: %f", fillLocation, yBinsContent.at(index-1));
      // minus 1 since loop starts at 1
      hist->Fill(fillLocation, weight*yBinsContent.at(index-1));
    }

    //TEMP
    //hist->Draw();
    //END TEMP
  }
}

void AliAnalysisTaskEmcalJetHMEC::FillHist(THnSparse * hist, Double_t *fillValue, Double_t weight, Bool_t noCorrection)
{
  if (fEmbeddingCorrectionHist == 0 || noCorrection == kTRUE)
  {
    hist->Fill(fillValue, weight);
  }
  else
  {
    // Jet pt is always located in the second position
    Double_t jetPt = fillValue[1];

    // Determine where to get the values in the correction hist
    Int_t xBin = fEmbeddingCorrectionHist->GetXaxis()->FindBin(jetPt);

    std::vector <Double_t> yBinsContent;
    accessSetOfYBinValues(fEmbeddingCorrectionHist, xBin, yBinsContent);

    // Loop over all possible bins to contribute.
    // If content is 0 then calling Fill won't make a difference
    for (Int_t index = 1; index <= fEmbeddingCorrectionHist->GetYaxis()->GetNbins(); index++)
    {
      // Determine the value to fill based on the center of the bins.
      // This in principle allows the binning between the correction and hist to be different
      fillValue[1]  = fEmbeddingCorrectionHist->GetYaxis()->GetBinCenter(index); 
      Printf("fillValue[1]: %f, weight: %f", fillValue[1], yBinsContent.at(index-1));
      // minus 1 since loop starts at 1
      hist->Fill(fillValue, weight*yBinsContent.at(index-1));
    }
  }
}

void AliAnalysisTaskEmcalJetHMEC::accessSetOfYBinValues(TH2F * hist, Int_t xBin, std::vector <Double_t> & yBinsContent, Double_t scaleFactor)
{
  for (Int_t index = 1; index <= hist->GetYaxis()->GetNbins(); index++)
  {
    //yBinsContent[index-1] = hist->GetBinContent(hist->GetBin(xBin,index));
    yBinsContent.push_back(hist->GetBinContent(hist->GetBin(xBin,index)));

    if (scaleFactor >= 0)
    {
      // -1 since index starts at 1
      hist->SetBinContent(hist->GetBin(xBin,index), yBinsContent.at(index-1)/scaleFactor);
    }
  }
}

