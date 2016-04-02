// $Id$

//////////
//Measure Jet-hadron correlations
//Does event Mixing using AliEventPoolManager
/////////

#include "AliAnalysisTaskEmcalJetHMEC.h"

// TODO: Clean up includes
//#include "TChain.h"
//#include "TTree.h"
//#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
//#include "TCanvas.h"
//#include <TClonesArray.h>
//#include <TParticle.h>
#include "AliVTrack.h"
//#include "TParameter.h"

#include "AliAODEvent.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDVertex.h"
//#include "AliCentrality.h"
//#include "AliAODJet.h"
#include "AliEmcalJet.h"
//#include "AliESDtrackCuts.h"

#include "AliClusterContainer.h"
#include "AliTrackContainer.h"

#include "TVector3.h"
#include "AliBasicParticle.h"
#include "AliEventPoolManager.h"

ClassImp(AliAnalysisTaskEmcalJetHMEC)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetHMEC",kFALSE),
  fTrkBias(5),
  fClusBias(5),
  fDoEventMixing(0),
  fMixingTracks(50000), fNMIXtracks(5000), fNMIXevents(5),
  fTriggerEventType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDoEffCorrection(0), fEffFunctionCorrection(0),
  fEmbeddingCorrectionHist(0),
  fDoLessSparseAxes(0), fDoWiderTrackBin(0),
  fCentBinSize(1),
  fPoolMgr(0x0), 
  fHistTrackPt(0),
  fHistJetEtaPhi(0), 
  fHistClusEtaPhiEn(0), 
  fHistJHPsi(0),
  fHistJetHEtaPhi(0), 
  fhnMixedEvents(0x0),
  fhnJH(0x0)
{
  // Default Constructor
  InitializeArraysToZero();
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetHMEC::AliAnalysisTaskEmcalJetHMEC(const char *name) : 
  AliAnalysisTaskEmcalJet(name,kTRUE),
  fClusBias(5),
  fDoEventMixing(0),
  fMixingTracks(50000), fNMIXtracks(5000), fNMIXevents(5),
  fTriggerEventType(AliVEvent::kEMCEJE), fMixingEventType(AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral),
  fDoEffCorrection(0), fEffFunctionCorrection(0),
  fEmbeddingCorrectionHist(0),
  fDoLessSparseAxes(0), fDoWiderTrackBin(0),
  fCentBinSize(1),
  fPoolMgr(0x0), 
  fHistTrackPt(0),
  fHistJetEtaPhi(0), 
  fHistClusEtaPhiEn(0),  
  fHistJHPsi(0),
  fHistJetHEtaPhi(0),
  fhnMixedEvents(0x0),
  fhnJH(0x0)
{
  // Constructor
  InitializeArraysToZero();
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::InitializeArraysToZero()
{
  // TODO: Replace loop limits with enum!
  for(Int_t ipta=0; ipta<7; ipta++){
    fHistTrackEtaPhi[ipta]=0;
  }
  for(Int_t icent = 0; icent<6; ++icent){
    fHistJetPt[icent]=0;
    fHistJetPtBias[icent]=0;
    fHistLeadJetPtBias[icent]=0;
    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
        fHistJetH[icent][iptjet][ieta]=0;
        fHistJetHBias[icent][iptjet][ieta]=0;
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
  fHistJetEtaPhi = new TH2F("fHistJetEtaPhi","Jet eta-phi",900,-1.8,1.8,720,-3.2,3.2);
  fHistJetHEtaPhi = new TH2F("fHistJetHEtaPhi","Jet-Hadron deta-dphi",900,-1.8,1.8,720,-1.6,4.8);

  fHistClusEtaPhiEn = new TH3F("fHistClusEtaPhiEn","Clus eta-phi-energy",900,-1.8,1.8,720,-3.2,3.2,20,0,20);

  fHistJHPsi = new TH3F("fHistJHPsi","Jet-Hadron ntr-trpt-dpsi",20,0,100,200,0,20,120,0,180);

  TString name;

  // TODO: Replace loop limits with enum!
  for(Int_t ipta=0; ipta<7; ++ipta){
    name = Form("fHistTrackEtaPhi_%i", ipta);
    fHistTrackEtaPhi[ipta] = new TH2F(name,name,400,-1,1,720,0.0,2.0*TMath::Pi());
    fOutput->Add(fHistTrackEtaPhi[ipta]);
  }

  // TODO: Replace loop limits with enum!
  for(Int_t icent = 0; icent<6; ++icent){
    name = Form("fHistJetPt_%i",icent);   
    fHistJetPt[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPt[icent]);

    name = Form("fHistJetPtBias_%i",icent);   
    fHistJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistJetPtBias[icent]);

    name = Form("fHistLeadJetPtBias_%i",icent);   
    fHistLeadJetPtBias[icent] = new TH1F(name,name,200,0,200);
    fOutput->Add(fHistLeadJetPtBias[icent]);

    for(Int_t iptjet = 0; iptjet<5; ++iptjet){
      for(Int_t ieta = 0; ieta<3; ++ieta){	
        name = Form("fHistJetH_%i_%i_%i",icent,iptjet,ieta);   
        fHistJetH[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
        fOutput->Add(fHistJetH[icent][iptjet][ieta]);

        name = Form("fHistJetHBias_%i_%i_%i",icent,iptjet,ieta);   
        fHistJetHBias[icent][iptjet][ieta]=new TH2F(name,name,72,-0.5*TMath::Pi(),1.5*TMath::Pi(),300,0,30);
        fOutput->Add(fHistJetHBias[icent][iptjet][ieta]);
      }
    }
  }

  // TODO: Remove variable that is always filled as 0!
  UInt_t cifras = 0; // bit coded, see GetDimParams() below 
  if(fDoLessSparseAxes) {
    cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5;
  } else {
    cifras = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7; 
  }
  fhnJH = NewTHnSparseF("fhnJH", cifras);
  fhnJH->Sumw2();
  fOutput->Add(fhnJH);

  // TODO: Remove variable that is always filled as 0!
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
  fOutput->Add(fHistJetEtaPhi);
  fOutput->Add(fHistJetHEtaPhi);
  fOutput->Add(fHistClusEtaPhiEn);
  fOutput->Add(fHistJHPsi);

  PostData(1, fOutput);

  // TODO: Refactor this section. Cleanup variable names
  //Event Mixing
  Int_t trackDepth = fMixingTracks; 
  Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager

  Int_t nZvtxBins  = 10;
  // bins for second buffer are shifted by 100 cm
  Double_t vertexBins[] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};
  Double_t* zvtxbin = vertexBins;

  // TODO: Refactor here in particular
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
Int_t AliAnalysisTaskEmcalJetHMEC::GetEtaBin(Double_t eta) const {
  // Get eta bin for histos.

  Int_t etabin = -1;
  eta = TMath::Abs(eta);
  if      (eta <= 0.4)              etabin = 0;
  else if (eta >  0.4 && eta < 0.8) etabin = 1;
  else if (eta >= 0.8)              etabin = 2;
  return etabin;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetTrackPtBin(Double_t pt) const
{
  Int_t ptBin = -1;
  if      (pt <  0.5) ptBin = 0;
  else if (pt <  1  ) ptBin = 1;
  else if (pt <  2  ) ptBin = 2;
  else if (pt <  3  ) ptBin = 3;
  else if (pt <  5  ) ptBin = 4;
  else if (pt <  8  ) ptBin = 5;
  else if (pt < 20  ) ptBin = 6;

  return ptBin;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHMEC::GetJetPtBin(Double_t pt) const 
{
  // Get jet pt  bin for histos.

  Int_t ptBin = -1;
  if      (pt >= 15 && pt < 20) ptBin = 0;
  else if (pt >= 20 && pt < 25) ptBin = 1;
  else if (pt >= 25 && pt < 30) ptBin = 2;
  else if (pt >= 30 && pt < 60) ptBin = 3;
  else if (pt >= 60)            ptBin = 4;

  return ptBin;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::ExecOnce() {
  AliAnalysisTaskEmcalJet::ExecOnce();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetHMEC::Run() {
  // Main loop called for each event

  // Retrieve clusters
  AliClusterContainer * clusters = GetClusterContainer(0);
  if (!clusters) {
    AliError(Form("%s: Unable to retrieve clusters!", GetName()));
    return kTRUE;
  }

  // Retrieve tracks
  AliTrackContainer * tracks = GetTrackContainer("tracksForCorrelations");
  if (!tracks) {
    AliError(Form("%s: Unable to retrieve tracks!", GetName()));
    return kTRUE;
  }

  // Retrieve jets
  AliJetContainer * jets = GetJetContainer(0);
  if (!jets) {
    AliError(Form("%s: Unable to retrieve jets!", GetName()));
    return kTRUE;
  }


  // Used to calculate the angle betwene the jet and the hadron
  TVector3 jetVector, trackVector;
  // Get z vertex
  Double_t zVertex=fVertex[2];
  // Flags
  Bool_t biasedJet = kFALSE;
  Bool_t leadJet = kFALSE;
  // Relative angles and distances
  Double_t deltaPhi = 0;
  Double_t deltaEta = 0;
  Double_t deltaR = 0;
  // Event activity (centrality or multipilicity)
  Double_t eventActivity = 0;
  // Efficiency correction
  Double_t efficiency = -999;
  // Determining bins for histogram indices
  Int_t jetPtBin = -1;
  Int_t etaBin = -1;
  // For the jets below
  AliEmcalJet * jet = 0;
  // For comparison to the current jet
  AliEmcalJet * leadingJet = jets->GetLeadingJet();

  // Determine the trigger for the current event
  UInt_t trigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  // Just to be certain that we are interating from the start
  jets->ResetCurrentID();
  while ((jet = jets->GetNextAcceptJet())) {
    
    // Selects only events that we are interested in (ie triggered)
    if (!(trigger & fTriggerEventType)) continue;

    // Jet properties
    // Determine if we have the lead jet
    leadJet = kFALSE;
    if (jet == leadingJet) leadJet = kTRUE;

    // Determine if the jet is biased
    biasedJet = BiasedJet(jet);

    // Calculate vector
    jetVector.SetXYZ(jet->Px(), jet->Py(), jet->Pz());

    // Fill jet properties
    FillHist(fHistJetPt[fCentBin], jet->Pt());
    if (biasedJet == kTRUE) {
      FillHist(fHistJetPtBias[fCentBin], jet->Pt());
    }

    fHistJetEtaPhi->Fill(jet->Eta(), jet->Phi());

    if (jet->Pt() > 15) {

      AliVTrack * track = 0;
      while ((track = tracks->GetNextAcceptTrack())) {

        // Determine relative angles and distances and set the respective variables
        GetDeltaEtaDeltaPhiDeltaR(track, jet, deltaEta, deltaPhi, deltaR);

        // Determine bins for filling histograms
        // jet Pt
        jetPtBin = GetJetPtBin(jet->Pt());
        if (jetPtBin < 0)
        {
          AliError(Form("Jet Pt Bin negative: %f", jet->Pt()));
          continue;
        }
        // eta
        etaBin = GetEtaBin(deltaEta);
        if (etaBin < 0) {
          AliError(Form("Eta Bin negative: %f", deltaEta));
          continue;
        }

        // Fill track properties
        fHistTrackPt->Fill(track->Pt());

        trackVector.SetXYZ( track->Px(), track->Py(), track->Pz() );
        if ( (jet->Pt() > 20.) && (jet->Pt() < 60.) ) {
          fHistJHPsi->Fill(tracks->GetNTracks(), track->Pt(), trackVector.Angle(jetVector) * TMath::RadToDeg() );
        }

        fHistJetH[fCentBin][jetPtBin][etaBin]->Fill(deltaPhi ,track->Pt());
        fHistJetHEtaPhi->Fill(deltaEta, deltaPhi);

        // Calculate single particle tracking efficiency for correlations
        efficiency = EffCorrection(track->Eta(), track->Pt(), fDoEffCorrection);

        if (biasedJet == kTRUE) {
          fHistJetHBias[fCentBin][jetPtBin][etaBin]->Fill(deltaPhi ,track->Pt());

          if (fBeamType == kAA || fBeamType == kpA) { //pA and AA
            eventActivity = fCent;
          }
          else if (fBeamType == kpp) {
            eventActivity = static_cast<Double_t>(tracks->GetNTracks());
          }

          if(fDoLessSparseAxes) { // check if we want all dimensions
            Double_t triggerEntries[6] = {eventActivity, jet->Pt(), track->Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet)};
            FillHist(fhnJH, triggerEntries, 1.0/efficiency);
          } else { 
            Double_t triggerEntries[8] = {eventActivity, jet->Pt(), track->Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet), 0.0, deltaR};
            FillHist(fhnJH, triggerEntries, 1.0/efficiency);
          }
        }

      } //track loop
    }//jet pt cut
  }//jet loop

  //Prepare to do event mixing

  // create a list of reduced objects. This speeds up processing and reduces memory consumption for the event pool
  TObjArray* tracksClone = 0;

  if(fDoEventMixing > 0){

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

    AliEventPool *pool = 0;
    if (fBeamType == kAA || fBeamType == kpA) {//everything but pp
      pool = fPoolMgr->GetEventPool(fCent, zVertex);
    }
    else if (fBeamType == kpp) {//pp only
      pool = fPoolMgr->GetEventPool(static_cast<Double_t>(tracks->GetNTracks()), zVertex);
    }

    if (!pool){
      if (fBeamType == kAA || fBeamType == kpA) AliFatal(Form("No pool found for centrality = %f, zVertex = %f", fCent, zVertex));
      else if (fBeamType == kpp) AliFatal(Form("No pool found for ntracks_pp = %d, zVertex = %f", tracks->GetNTracks(), zVertex));
      return kTRUE;
    }

    // The number of events in the pool
    Int_t nMix = pool->GetCurrentNEvents();

    if(trigger & fTriggerEventType) {
      // check for a trigger jet
      if (pool->IsReady() || pool->NTracksInPool() >= fNMIXtracks || nMix >= fNMIXevents) {

        jets->ResetCurrentID();
        while ((jet = jets->GetNextAcceptJet())) {
          // Jet properties
          // Determine if we have the lead jet
          leadJet = kFALSE;
          if (jet == leadingJet) { leadJet = kTRUE; }

          // Determine if the jet is biased
          biasedJet = BiasedJet(jet);

          // Make sure event contains jet above our threshold (reduce stats of sparse)
          if (jet->Pt() < 15) continue;

          // Fill for biased jet triggers only
          if (biasedJet == kTRUE) {

            // Fill mixed-event histos here  
            for (Int_t jMix=0; jMix < nMix; jMix++) {
              TObjArray* bgTracks = pool->GetEvent(jMix);

              for(Int_t ibg=0; ibg < bgTracks->GetEntries(); ibg++){

                AliBasicParticle *bgTrack = static_cast<AliBasicParticle*>(bgTracks->At(ibg));
                if(!bgTrack)
                {
                  AliError(Form("%s:Failed to retrieve tracks from mixed events", GetName()));
                }

                // Calculate single particle tracking efficiency of mixed events for correlations
                efficiency = EffCorrection(bgTrack->Eta(), bgTrack->Pt(), fDoEffCorrection);

                // Phi is [-0.5*TMath::Pi(), 3*TMath::Pi()/2.]
                GetDeltaEtaDeltaPhiDeltaR(bgTrack, jet, deltaEta, deltaPhi, deltaR);

                if (fBeamType == kAA || fBeamType == kpA) { //pA and AA
                  eventActivity = fCent;
                }
                else if (fBeamType == kpp) {
                  eventActivity = static_cast<Double_t>(tracks->GetNTracks());
                }

                if(fDoLessSparseAxes) {  // check if we want all the axis filled
                  Double_t triggerEntries[6] = {eventActivity, jet->Pt(), bgTrack->Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet)};
                  FillHist(fhnMixedEvents, triggerEntries, 1./(nMix*efficiency), kTRUE);
                } else {
                  Double_t triggerEntries[8] = {eventActivity, jet->Pt(), bgTrack->Pt(), deltaEta, deltaPhi, static_cast<Double_t>(leadJet), 0.0, deltaR};
                  FillHist(fhnMixedEvents, triggerEntries, 1./(nMix*efficiency), kTRUE);
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
Bool_t AliAnalysisTaskEmcalJetHMEC::BiasedJet(AliEmcalJet * jet)
{
  if ((jet->MaxTrackPt() > fTrkBias) || (jet->MaxClusterPt() > fClusBias))
  {
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHMEC::GetDeltaEtaDeltaPhiDeltaR(AliVParticle * particleOne, AliVParticle * particleTwo, Double_t & deltaPhi, Double_t & deltaEta, Double_t & deltaR)
{
  // TODO: Understand order of arguments to DeltaPhi vs DeltaEta
  // Returns deltaPhi in symmetric range so that we can calculate DeltaR.
  deltaPhi = DeltaPhi(particleTwo->Phi(), particleOne->Phi(), -1.0*TMath::Pi(), TMath::Pi());
  deltaEta = particleOne->Eta() - particleTwo->Eta();
  deltaR = TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);

  // Adjust to the normal range after the DeltaR caluclation
  deltaPhi = DeltaPhi(particleTwo->Phi(), particleOne->Phi(), -0.5*TMath::Pi(), 3*TMath::Pi()/2.);
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

//________________________________________________________________________
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
  // clones a track list by using AliBasicTrack which uses much less memory (used for event mixing)
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);

  // Loop over all tracks
  AliVParticle * particle = 0;
  AliBasicParticle * clone = 0;
  AliTrackContainer * tracks = GetTrackContainer("tracksForCorrelations");
  // Ensure that we start from the beginning
  tracks->ResetCurrentID();

  while ((particle = tracks->GetNextAcceptTrack()))
  {
    // Fill some QA information about the tracks
    Int_t trackPtBin = GetTrackPtBin(particle->Pt());
    if(trackPtBin > -1) fHistTrackEtaPhi[trackPtBin]->Fill(particle->Eta(),particle->Phi());

    // Create new particle
    clone = new AliBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
    // Set so that we can do comparisons using the IsEqual() function.
    clone ->SetUniqueID(particle->GetUniqueID());

    tracksClone->Add(clone);
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

    /*if(runSwitchGood == 0 && fCentBin == 0) effSwitch = 2;
    if(runSwitchGood == 0 && fCentBin == 1) effSwitch = 3;
    if(runSwitchGood == 0 && fCentBin == 2) effSwitch = 4;
    if(runSwitchGood == 0 && fCentBin == 3) effSwitch = 5;
    if(runSwitchGood == 1 && fCentBin == 0) effSwitch = 6;
    if(runSwitchGood == 1 && fCentBin == 1) effSwitch = 7;
    if(runSwitchGood == 1 && fCentBin == 2) effSwitch = 8;
    if(runSwitchGood == 1 && fCentBin == 3) effSwitch = 9;*/

    // Determine which efficiency to use.
    // This is just a way to map all possible values of the cent bin and runSwitchGood to a unique flag.
    // 4 is the number of cent bins, and we want to index the effSwitch starting at 2.
    effSwitch = 2 + runSwitchGood*4 + fCentBin;
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

//________________________________________________________________________
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

//________________________________________________________________________
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
      fillValue[1] = fEmbeddingCorrectionHist->GetYaxis()->GetBinCenter(index); 
      Printf("fillValue[1]: %f, weight: %f", fillValue[1], yBinsContent.at(index-1));
      // minus 1 since loop starts at 1
      hist->Fill(fillValue, weight*yBinsContent.at(index-1));
    }
  }
}

//________________________________________________________________________
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

