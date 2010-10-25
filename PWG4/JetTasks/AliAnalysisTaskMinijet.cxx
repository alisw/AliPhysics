#include <TChain.h>
#include <TList.h>

#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include "TRandom.h"

#include "AliVEvent.h"
#include "AliVParticle.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"

#include "AliAnalysisManager.h"

#include "AliAnalysisTaskMinijet.h"

// Analysis Task for mini jet activity analysis
// Authors: Eva Sicking

ClassImp(AliAnalysisTaskMinijet)

//________________________________________________________________________
AliAnalysisTaskMinijet::AliAnalysisTaskMinijet(const char *name)
  : AliAnalysisTaskSE(name),
    fUseMC(kFALSE),
    fCuts(0),
    fRadiusCut(0.7),
    fTriggerPtCut(0.7),
    fAssociatePtCut(0.4),
    fLeadingOrRandom(1),
    fMode(0),
    fVertexZCut(10.),
    fESDEvent(0),
    fAODEvent(0),
    fHists(0),
    fHistPt(0),
    fHistPtMC(0)

{
  for(Int_t i = 0;i< 4;i++){
    fVertexZ[i]               = 0;
    fPt[i]                    = 0;
    fEta[i]                   = 0;
    fPhi[i]                   = 0;
    
    fDPhiDEtaEventAxis[i]     = 0;
    fTriggerNch[i]            = 0;
    fTriggerTracklet[i]       = 0;
    
    fNch07Nch[i]              = 0;
    pNch07Nch[i]              = 0;
    fNch07Tracklet[i]         = 0;
    fNchTracklet[i]         = 0;
    pNch07Tracklet[i]         = 0;
    for(Int_t j=0;j<150;j++){
      fDPhiEventAxisNchBin[i][j]   = 0;
      fDPhiEventAxisNchBinTrig[i][j]   = 0;

      fDPhiEventAxisTrackletBin[i][j]   = 0;
      fDPhiEventAxisTrackletBinTrig[i][j]   = 0;
    }
  }
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMinijet::~AliAnalysisTaskMinijet()
{
  // Destructor

  if (fHists && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHists;
}

//________________________________________________________________________
void AliAnalysisTaskMinijet::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

   
  fHistPt = new TH1F("fHistPt", "P_{T} distribution REC", 15, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);

  
  if (fUseMC) {
    fHistPtMC = new TH1F("fHistPtMC", "P_{T} distribution MC", 15, 0.1, 3.1);
    fHistPtMC->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPtMC->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPtMC->SetMarkerStyle(kFullCircle);
  }

  TString labels[4]={"ESDrec", "ESDmc", "AODrec", "AODmc"};

  for(Int_t i=2*fMode;i<1+2*fMode+fUseMC;i++){
  
    fVertexZ[i]                  = new TH1F(Form("fVertexZ%s",labels[i].Data()),
						Form("fVertexZ%s",labels[i].Data()) ,  
					    200, -10., 10);
    fPt[i]                       = new TH1F(Form("fPt%s",labels[i].Data()),
					    Form("fPt%s",labels[i].Data()) ,  
					    500, 0., 50);
    fEta[i]                      = new TH1F (Form("fEta%s",labels[i].Data()),
						 Form("fEta%s",labels[i].Data()) ,  
						 100, -1., 1);
    fPhi[i]                      = new TH1F(Form("fPhi%s",labels[i].Data()),
						Form("fPhi%s",labels[i].Data()) ,  
						360, 0.,2*TMath::Pi());
    fDPhiDEtaEventAxis[i]        = new TH2F(Form("fDPhiDEtaEventAxis%s",labels[i].Data()),
						Form("fDPhiDEtaEventAxis%s",labels[i].Data()) ,  
						180, -1., 2*TMath::Pi()-1, 200, -2.,2.);
    fTriggerNch[i]               = new TH1F(Form("fTriggerNch%s",labels[i].Data()),
						Form("fTriggerNch%s",labels[i].Data()) ,  
						250, -0.5, 249.5);
    fTriggerTracklet[i]          = new TH1F(Form("fTriggerTracklet%s",labels[i].Data()),
						Form("fTriggerTracklet%s",labels[i].Data()) ,  
						250, -0.5, 249.5);
    fNch07Nch[i]                 = new TH2F(Form("fNch07Nch%s",labels[i].Data()),
						Form("fNch07Nch%s",labels[i].Data()) ,  
						250, -2.5, 247.5,250, -2.5, 247.5);
    pNch07Nch[i]                 = new TProfile(Form("pNch07Nch%s",labels[i].Data()),
						    Form("pNch07Nch%s",labels[i].Data()) ,  
						    250, -2.5, 247.5);
    fNch07Tracklet[i]            = new TH2F(Form("fNch07Tracklet%s",labels[i].Data()),
						Form("fNch07Tracklet%s",labels[i].Data()) ,  
						250, -2.5, 247.5,250, -2.5, 247.5);
    fNchTracklet[i]              = new TH2F(Form("fNchTracklet%s",labels[i].Data()),
						Form("fNchTracklet%s",labels[i].Data()) ,  
						250, -2.5, 247.5,250, -2.5, 247.5);
    pNch07Tracklet[i]            = new TProfile(Form("pNch07Tracklet%s",labels[i].Data()),
						    Form("pNch07Tracklet%s",labels[i].Data()) ,  
						    250, -2.5, 247.5);
    for(Int_t j=0;j<150;j++){
      fDPhiEventAxisNchBin[i][j]          = new TH1F(Form("fDPhiEventAxisNchBin%s%02d",
							      labels[i].Data(),j),
							 Form("fDPhiEventAxisNchBin%s%02d",
							      labels[i].Data(),j) ,  
							 180, 0., TMath::Pi());
      fDPhiEventAxisNchBinTrig[i][j]          = new TH1F(Form("fDPhiEventAxisNchBinTrig%s%02d",
							      labels[i].Data(),j),
							 Form("fDPhiEventAxisNchBinTrig%s%02d",
							      labels[i].Data(),j) ,  
							 180, 0., TMath::Pi());
      
      fDPhiEventAxisTrackletBin[i][j]     = new TH1F(Form("fDPhiEventAxisTrackletBin%s%02d",
							      labels[i].Data(),j),
							 Form("fDPhiEventAxisTrackletBin%s%02d",
							      labels[i].Data(),j) ,  
							 180, 0., TMath::Pi());
      fDPhiEventAxisTrackletBinTrig[i][j]     = new TH1F(Form("fDPhiEventAxisTrackletBinTrig%s%02d",
							      labels[i].Data(),j),
							 Form("fDPhiEventAxisTrackletBinTrig%s%02d",
							      labels[i].Data(),j) ,  
							 180, 0., TMath::Pi());
    }
  }

  fHists = new TList();
  fHists->SetOwner();

  fHists->Add(fHistPt);
  if(fUseMC)fHists->Add(fHistPtMC); 

  for(Int_t i=2*fMode;i<1+2*fMode+fUseMC;i++){
    fHists->Add(fVertexZ[i]);
    fHists->Add(fPt[i]);
    fHists->Add(fEta[i]);
    fHists->Add(fPhi[i]);
    fHists->Add(fDPhiDEtaEventAxis[i]);
    fHists->Add(fTriggerNch[i]);
    fHists->Add(fTriggerTracklet[i]);
    fHists->Add(fNch07Nch[i]);
    fHists->Add(pNch07Nch[i]);
    fHists->Add(fNch07Tracklet[i]);
    fHists->Add(fNchTracklet[i]);
    fHists->Add(pNch07Tracklet[i]);
    for(Int_t j=0;j<150;j++){
      fHists->Add(fDPhiEventAxisNchBin[i][j]);
      fHists->Add(fDPhiEventAxisNchBinTrig[i][j]);

      fHists->Add(fDPhiEventAxisTrackletBin[i][j]);
      fHists->Add(fDPhiEventAxisTrackletBinTrig[i][j]);
    }
  }
  PostData(1, fHists);
 
}

//________________________________________________________________________
void AliAnalysisTaskMinijet::UserExec(Option_t *)
{
  // Main loop
  // Called for each event

  AliVEvent *event = InputEvent();
  if (!event) {
     Error("UserExec", "Could not retrieve event");
     return;
  }
  
  //get events, either ESD or AOD
  fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  

  //arrays for track properties
  Float_t *pt  = 0;
  Float_t *eta = 0;
  Float_t *phi = 0;
  //number of accepted tracks and tracklets
  Int_t ntracks = 0;  //return value for reading functions for ESD and AOD
  Int_t *numbers = 0; // numbers[0]=nAcceptedTracks, numbers[1]=ntracklets, 


  //read data and analyse. Implemented for ESD, ESDmc, AOD, AODmc
  //-------------------------------------------------------------
  if (fESDEvent && fMode==0) {//ESDs
    //ESD reading and analysis
    //------
    ntracks = LoopESD(&pt, &eta, &phi, &numbers); //read tracks
    if(ntracks>0){
      if(fDebug){
	Printf("ntracks=%d", numbers[0]);
	Printf("ntracklets=%d", numbers[1]);
      }
      Analyse(pt, eta, phi, ntracks, numbers[1], 0); //analyse tracks
    }
    CleanArrays(pt, eta, phi, numbers); // clean up array memory
    //------

    //ESD MC reading and analysis
    //------
    if (fUseMC){
      ntracks = LoopESDMC(&pt, &eta, &phi); //read mc particles
      if(ntracks>0){
	Analyse(pt, eta, phi, ntracks, 0, 1);//analyse
      }
      CleanArrays(pt, eta, phi);// clean up array memory
    }
    //------
  }
  
  if (fAODEvent && fMode==1) {//AODs

    //AOD reading and analysis
    //------
    ntracks = LoopAOD(&pt, &eta, &phi, &numbers);//read tracks
    if(ntracks>0){
      if(fDebug){
	Printf("ntracks=%d", numbers[0]);
	Printf("ntracklets=%d", numbers[1]);
      }
      Analyse(pt, eta, phi, ntracks, numbers[1], 2);//analyse
    }
    CleanArrays(pt, eta, phi, numbers);// clean up array memory
    //------

    //AOD MCreading and analysis
    //------
    if (fUseMC){
      ntracks = LoopAODMC(&pt, &eta, &phi);//read tracks
      if(ntracks>0){
	Analyse(pt, eta, phi, ntracks, 0, 3);//analyse
      }
      CleanArrays(pt, eta, phi);// clean up array memory
    }
    //------
  }
  //-------------------------------------------------------------
  
  
  // Post output data.
  //  PostData(1, fHists);

}      


//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopESD(Float_t **ptArray, Float_t ** etaArray, 
				      Float_t ** phiArray, Int_t **numbers )
{
  // gives back the number of esd tracks and pointer to arrays with track
  // properties (pt, eta, phi)
  
  // Retreive the number of all tracks for this event.
  Int_t ntracks = fESDEvent->GetNumberOfTracks();
  if(fDebug)Printf("all ESD tracks: %d", ntracks);


  //first loop to check how many tracks are accepted
  //------------------
  Int_t nAcceptedTracks=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliESDtrack *track = (AliESDtrack *)fESDEvent->GetTrack(iTracks);
    if (!track) {
      Error("UserExec", "Could not receive track %d", iTracks);
      continue;
    }
    if (fCuts->AcceptTrack(track)) ++nAcceptedTracks;
  }
  //------------------
  
  //generate arrays
  *ptArray = new Float_t[nAcceptedTracks]; 
  *etaArray = new Float_t[nAcceptedTracks]; 
  *phiArray = new Float_t[nAcceptedTracks]; 
  *numbers = new Int_t[2]; //ntracksAccepted, ntracklets

  //check if event is pile up or no tracks are accepted, return to user exec
  if(fESDEvent->IsPileupFromSPD(3,0,8)) return 0;  
  if(nAcceptedTracks==0) return 0;

  //accept event only, if vertex is good and is within fVertexZcut region
  const AliESDVertex*	vertexESD   = fESDEvent->GetPrimaryVertex();
  if(vertexESD->GetNContributors()==0)return 0;
  Float_t fVz= vertexESD->GetZ();
  if(TMath::Abs(fVz)>fVertexZCut) return 0;
  fVertexZ[0]->Fill(fVz);
  
  
  // Track loop
  Int_t iAcceptedTrack=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliESDtrack *track = (AliESDtrack *)fESDEvent->GetTrack(iTracks);
    if (!track) {
      Error("UserExec", "Could not receive track %d", iTracks);
      continue;
    }
    if (fCuts->AcceptTrack(track)){
      (*ptArray)[iAcceptedTrack]  = track->Pt();
      (*etaArray)[iAcceptedTrack] = track->Eta();
      (*phiArray)[iAcceptedTrack++] = track->Phi();
      fHistPt->Fill(track->Pt());
    }
  }


  //tracklet loop
  Int_t ntrackletsAccept=0;
  AliMultiplicity * mult = (AliMultiplicity*)(fESDEvent->GetMultiplicity());
  Int_t ntracklets = mult->GetNumberOfTracklets();
  for(Int_t i=0;i< ntracklets;i++){
    if(mult->GetDeltaPhi(i)<0.05){
      ntrackletsAccept++;
    }
  }

  (*numbers)[0] = nAcceptedTracks;
  (*numbers)[1] = ntrackletsAccept;

  return nAcceptedTracks;

}   

//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopESDMC(Float_t **ptArray, Float_t ** etaArray, 
					Float_t ** phiArray)
{
  // gives back the number of charged prim MC particle and pointer to arrays 
  // with particle properties (pt, eta, phi)

  AliMCEvent *mcEvent = (AliMCEvent*) MCEvent();
  if (!mcEvent) {
    Error("UserExec", "Could not retrieve MC event");
    return 0;
  }

  AliStack* stack = 0x0;
  if(fUseMC) stack = MCEvent()->Stack();
  
  Int_t ntracks = mcEvent->GetNumberOfTracks();
  if(fDebug)Printf("MC particles: %d", ntracks);


  //----------------------------------
  //first track loop to check how many chared primary tracks are available
  Int_t nChargedPrimaries=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliMCParticle *track = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iTracks));
    if (!track) {
      Error("UserExec", "Could not receive track %d", iTracks);
      continue;
    }
    if (!(stack->IsPhysicalPrimary(track->Label()))) continue;
    if (track->Particle()->GetPDG()->Charge() == 0) continue;
    //same cuts as on ESDtracks
    if(TMath::Abs(track->Eta())>0.9 || track->Pt()<0.2 || track->Pt()>200) continue;
    ++nChargedPrimaries;
  }
  //----------------------------------


  *ptArray = new Float_t[nChargedPrimaries]; 
  *etaArray = new Float_t[nChargedPrimaries]; 
  *phiArray = new Float_t[nChargedPrimaries]; 

  if(nChargedPrimaries==0) return 0;  

  AliGenEventHeader*  header = MCEvent()->GenEventHeader();
  TArrayF mcV;
  header->PrimaryVertex(mcV);
  Float_t vzMC = mcV[2];
  if(TMath::Abs(vzMC)>fVertexZCut) return 0;
  fVertexZ[1]->Fill(vzMC);
  


  //track loop
  Int_t iChargedPrimaries=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliMCParticle *track = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iTracks));
    if (!track) {
      Error("UserExec", "Could not receive track %d", iTracks);
      continue;
    }
    if (!(stack->IsPhysicalPrimary(track->Label()))) continue;
    if (track->Particle()->GetPDG()->Charge() == 0) continue;
    //same cuts as on ESDtracks
    if(TMath::Abs(track->Eta())>0.9 || track->Pt()<0.2 || track->Pt()>200) continue;
   
    
    fHistPtMC->Fill(track->Pt());
    //fills arrays with track properties
    (*ptArray)[iChargedPrimaries]  = track->Pt();
    (*etaArray)[iChargedPrimaries] = track->Eta();
    (*phiArray)[iChargedPrimaries++] = track->Phi();

  } //track loop

  return nChargedPrimaries;
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopAOD(Float_t **ptArray, Float_t ** etaArray, 
				      Float_t ** phiArray, Int_t ** numbers)
{
  // gives back the number of AOD tracks and pointer to arrays with track 
  // properties (pt, eta, phi)

  
  // Retreive the number of tracks for this event.
  Int_t ntracks = fAODEvent->GetNumberOfTracks();
  if(fDebug) Printf("AOD tracks: %d", ntracks);
  

  Int_t nAcceptedTracks=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliAODTrack *track = (AliAODTrack *)fAODEvent->GetTrack(iTracks);
    if (!track) {
      Error("UserExec", "Could not receive track %d", iTracks);
      continue;
    }
    if(track->TestFilterBit(16))nAcceptedTracks++;
  }
  
  *ptArray = new Float_t[nAcceptedTracks];
  *etaArray = new Float_t[nAcceptedTracks]; 
  *phiArray = new Float_t[nAcceptedTracks]; 
  *numbers = new Int_t[2]; 

 
  if(nAcceptedTracks==0) return 0;
  AliAODVertex*	vertex= fAODEvent->GetPrimaryVertex();
  if(!vertex) return 0;
  Double_t vzAOD=vertex->GetZ();
  if(vertex->GetNContributors()==0) return 0;
  if(TMath::Abs(vzAOD)>fVertexZCut) return 0;
  fVertexZ[2]->Fill(vzAOD);

  // Track loop to fill a pT spectrum
  Int_t iAcceptedTracks=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliAODTrack *track = (AliAODTrack *)fAODEvent->GetTrack(iTracks);
    if (!track) {
     Error("UserExec", "Could not receive track %d", iTracks);
     continue;
    }
    if(!track->TestFilterBit(16))continue;
    fHistPt->Fill(track->Pt());

    //fills arrays with track properties
    (*ptArray)[iAcceptedTracks]  = track->Pt();
    (*etaArray)[iAcceptedTracks] = track->Eta();
    (*phiArray)[iAcceptedTracks++] = track->Phi();

  } //track loop 

  //tracklet loop
  Int_t ntrackletsAccept=0;
  AliAODTracklets * mult= (AliAODTracklets*)fAODEvent->GetTracklets();
  for(Int_t i=0;i<mult->GetNumberOfTracklets();++i){
    if(TMath::Abs(mult->GetDeltaPhi(i))<0.05){
      ++ntrackletsAccept;
    }
  }

  (*numbers)[0] = nAcceptedTracks;
  (*numbers)[1] = ntrackletsAccept;
  
  return nAcceptedTracks;

}   

//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopAODMC(Float_t **ptArray, Float_t ** etaArray, 
					Float_t ** phiArray)
{
  // gives back the number of AOD MC particles and pointer to arrays with particle 
  // properties (pt, eta, phi)
  
  //retreive MC particles from event
  TClonesArray *mcArray = (TClonesArray*)fAODEvent->
    FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    Printf("%s:%d No MC particle branch found",(char*)__FILE__,__LINE__);
    return kFALSE;
  }
  
  Int_t ntracks = mcArray->GetEntriesFast();
  if(fDebug)Printf("MC particles: %d", ntracks);


  // Track loop: chek how many particles will be accepted
  Float_t vzMC=0.;
  Int_t nAcceptedTracks=0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
    if (!track) {
      Error("UserExec", "Could not receive track %d", it);
      continue;
    }
    if(!track->IsPhysicalPrimary())continue;
    if (track->Charge() == 0) continue;
    if(TMath::Abs(track->Eta())>1. || track->Pt()<0.2 || track->Pt()>200) continue; //same cuts as in ESD filter
    if(track->TestBit(16))nAcceptedTracks++;
    if(nAcceptedTracks==1) vzMC= track->Zv(); // check only one time. (only one vertex per event allowed)
  }

  //generate array with size of number of accepted tracks
  *ptArray = new Float_t[nAcceptedTracks]; 
  *etaArray = new Float_t[nAcceptedTracks]; 
  *phiArray = new Float_t[nAcceptedTracks]; 
  
  if(nAcceptedTracks==0) return 0;

  //check vertex
  if(TMath::Abs(vzMC)>fVertexZCut) return 0;
  fVertexZ[3]->Fill(vzMC);
  

  // Track loop: fill arrays for accepted tracks 
  Int_t iAcceptedTracks=0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
    if (!track) {
      Error("UserExec", "Could not receive track %d", it);
      continue;
    }
    if(!track->IsPhysicalPrimary())continue;
    if (track->Charge() == 0) continue;
    if(TMath::Abs(track->Eta())>0.9 || track->Pt()<0.2 || track->Pt()>200) continue;

    if(track->TestBit(16)){
      fHistPtMC->Fill(track->Pt());
      (*ptArray)[iAcceptedTracks]  = track->Pt();
      (*etaArray)[iAcceptedTracks] = track->Eta();
      (*phiArray)[iAcceptedTracks++] = track->Phi();
    }
  }
  
  return nAcceptedTracks;

} 

//________________________________________________________________________
void AliAnalysisTaskMinijet::Analyse(Float_t *pt, Float_t *eta, Float_t *phi, Int_t ntracks, 
				     Int_t ntracklets, Int_t mode)
{

  // analyse track properties (comming from either ESDs or AODs) in order to compute 
  // mini jet activity (chared tracks) as function of charged multiplicity
  
  // ntracks and ntracklets are already the number of accepted tracks and tracklets

  if(fDebug) Printf("In Analysis\n");
  
  Float_t ptEventAxis=0;  // pt leading
  Float_t etaEventAxis=0; // eta leading
  Float_t phiEventAxis=0; // phi leading
  
  Float_t ptOthers  = 0; // pt others // for all other tracks around event axis -> see loop
  Float_t etaOthers = 0; // eta others
  Float_t phiOthers = 0; // phi others
  
  Int_t *pindex  = new Int_t[ntracks];//index needed for sorting of pt array

  for (Int_t i = 0; i < ntracks; i++) {
    //filling of simple check plots
    fPt[mode]  ->Fill( pt[i]);
    fEta[mode] ->Fill(eta[i]);
    fPhi[mode] ->Fill(phi[i]);
    pindex[i]=0; //set all values to zero
  }

  //sort all tracks according to their pt
  if(ntracks) TMath::Sort(ntracks, pt, pindex, kTRUE);  
  
  
  // define event axis: leading or random track (with pt>fTriggerPtCut) 
  // ---------------------------------------
  Int_t highPtTracks=0;
  for(Int_t i=0;i<ntracks;i++){//show just the filled entries, skip empty ones.
    if(fDebug) Printf("%d:  pt = %f, number %i \n",mode, pt[pindex[i]],i);
    if(pt[pindex[i]]>fTriggerPtCut) highPtTracks++;
  }

  //  plot of multiplicity distributions
  fNch07Nch[mode]->Fill(ntracks, highPtTracks);     
  pNch07Nch[mode]->Fill(ntracks, highPtTracks);     
  if(ntracklets){
    fNch07Tracklet[mode]->Fill(ntracklets, highPtTracks);     
    fNchTracklet[mode]->Fill(ntracklets, ntracks);     
    pNch07Tracklet[mode]->Fill(ntracklets, highPtTracks);     
  }
 
  //analysis can only be performed with event axis, defined by high pt track
  if(highPtTracks>0){

    //check setter of event axis 
    //default option: random=1, 
    //second option:leading=0
    Int_t axis=-1;
    if(fLeadingOrRandom==0)axis=0;
    else if (fLeadingOrRandom==1)axis= Int_t((highPtTracks)*gRandom->Rndm());
    else Printf("Wrong settings for event axis.");
    if(fDebug)Printf("Axis tracks has pT=%f", pt[pindex[axis]]);

    //---------------------------------------


    if(ntracks>1){ // require at least two tracks (leading and prob. accosicates)
      
      //EventAxisRandom track properties
      ptEventAxis  = pt [pindex[axis]];
      etaEventAxis = eta[pindex[axis]];
      phiEventAxis = phi[pindex[axis]];

           
      //track loop for event propoerties around event axis with pt>triggerPtCut
      //loop only over already accepted tracks except event axis 
      if(ptEventAxis>fTriggerPtCut){

	for (Int_t iTrack = 0; iTrack < ntracks; iTrack++) {
	  
	  if(axis==iTrack)continue; // no double counting
	  
	  ptOthers   = pt [pindex[iTrack]];
	  etaOthers  = eta[pindex[iTrack]];
	  phiOthers  = phi[pindex[iTrack]];


	  if(ptOthers>fAssociatePtCut){ // only tracks which fullfill associate pt cut

	    Float_t dPhi=TMath::Abs(phiOthers-phiEventAxis);
	    if(dPhi>TMath::Pi())      dPhi=2*TMath::Pi()-dPhi;
	    Float_t dEta=etaOthers-etaEventAxis;
	    
	    Float_t dphiplot = phiOthers-phiEventAxis;
	    if(dphiplot>2*TMath::Pi()-1) dphiplot = dphiplot-2*TMath::Pi();
	    else if(dphiplot<-1)dphiplot=dphiplot+2*TMath::Pi();
	    fDPhiDEtaEventAxis[mode]->Fill(dphiplot, dEta);
	    
	    if(ntracks<150){
	      fDPhiEventAxisNchBin[mode][ntracks]->Fill(dPhi);
	      if(ptOthers>fTriggerPtCut)
		fDPhiEventAxisNchBinTrig[mode][ntracks]->Fill(dPhi);
	    }

	    if(ntracklets<150){
	      fDPhiEventAxisTrackletBin[mode][ntracklets]->Fill(dPhi);
	      if(ptOthers>fTriggerPtCut)
		fDPhiEventAxisTrackletBinTrig[mode][ntracklets]->Fill(dPhi);
	    }

	  }//tracks fulfill assoc track cut
	  
	}// end track loop


	// fill histogram with number of tracks (pt>fAssociatePtCut) around event axis
	// how often is there a trigger particle at a certain Nch bin
	fTriggerNch[mode]->Fill(ntracks); 
	fTriggerTracklet[mode]->Fill(ntracklets); 

      }//if track pt is at least trigger pt

    } //if there are more than 1 track

  }//if there is at least one high pt track

  if(pindex){// clean up array memory used for TMath::Sort
    delete[] pindex; 
    pindex=0;
  }

  PostData(1, fHists);

}

//________________________________________________________________________
void AliAnalysisTaskMinijet::Terminate(Option_t*)
{
  //terminate function is called at the end
  //can be used to draw histograms etc.

}

//________________________________________________________________________
void AliAnalysisTaskMinijet::CleanArrays(Float_t* pt, Float_t* eta, Float_t* phi,Int_t* numbers)
{
  //clean up of memory used for arrays of track properties

  if(pt){
    delete[] pt; 
    pt=0; 
  }
  if(eta){
    delete[] eta; 
    eta=0; 
  }
  if(phi){
    delete[] phi; 
    phi=0; 
  }
  if(numbers){
    delete[] numbers; 
    numbers=0; 
  }

}
