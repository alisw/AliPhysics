#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
// #include "THnF.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliESDEvent.h"
#include "AlidNdPtEventCuts.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliPhysicsSelection.h"

#include "AliVMultiplicity.h"

#include "AliMCEvent.h"
#include "AliStack.h"

#include "AlidNdPtUnifiedAnalysisTask.h"

/// \cond CLASSIMP
ClassImp(AlidNdPtUnifiedAnalysisTask);
/// \endcond

//________________________________________________________________________
AlidNdPtUnifiedAnalysisTask::AlidNdPtUnifiedAnalysisTask(const char *name) 
  : AliAnalysisTaskSE(name),
    fEvent(0),
    fMCEvent(0),
    fMCStack(0),
    fOutputList(0),
    fHistTrack(0),
    fHistEvent(0),
    fESDtrackCuts(0),
    fESDtrackCutsInit(kFALSE),
    fEventCuts(0),
    fEventCutsInit(kFALSE),
    fHistMCGenPrimTrack(0),
    fHistMCRecTrack(0),
    fHistMCRecPrimTrack(0),
    fHistMCRecSecTrack(0),    
    fHistMCRecEvent(0),
    fHistMCTrigEvent(0),
    fHistMCGenEvent(0),
    fHistMCGenINEL0Event(0),
    fHistMCTrigINEL0Event(0),
    fHistMCRecINEL0Event(0),
//     fHistMCGenTrackINEL0(0),
    fBinsPt(0),
    fBinsEta(0),
    fBinsMultCent(0),
    fBinsZv(0),
    fIsMC(kFALSE),
    fIsESD(kTRUE),
    fIsEventINEL0(kFALSE),
    fUseMultiplicity(kTRUE),
    fMinEta(-10),
    fMaxEta(10),
    fMinPt(0),
    fMaxPt(999),
    fSigmaMeanXYZv(),
    fMeanXYZv(),
    fEventTriggerRequired(kTRUE),
    fZvtx(10),
    fTPCRefit(kFALSE),
    fITSRefit(kFALSE),
    fAcceptKinks(kTRUE),
    fminNCrossedRowsTPC(0),
    fminRatioCrossedRowsOverFindableClustersTPC(0),
    fmaxFractionSharedTPCCluster(0),
    fmaxchi2perTPCcl(0),
    fmaxchi2perITScl(0),
    fdcatoVertex2D(kFALSE),
    fsigmatovertex(kFALSE),
    fmaxdcazITSTPC(0),
    fmaxchi2TPCconstrained(0),
    fminActiveLength(0),
    fUseGeomCut(kFALSE),
    fParticleLabel(kPrimary)
{
  
  // Set default binning
  Double_t binsMultCentDefault[2] = {0,10000};
  Double_t binsPtDefault[69] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
  Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
  Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};

  SetBinsPt(68,binsPtDefault);
  SetBinsEta(30,binsEtaDefault);
  SetBinsMultCent(1,binsMultCentDefault);
  SetBinsZv(12,binsZvDefault);
  
  SetMeanXYZv(0.0,0.0,0.0);
  SetSigmaMeanXYZv(1.0,1.0,10.0);
  SetZvtx(10.);
  SetEventTriggerRequired(kTRUE);
  // Constructor
  
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::UserCreateOutputObjects()
{
  // Create histograms
 
  // Called once
  OpenFile(1,"recreate");
  fOutputList = new TList();
  fOutputList -> SetOwner();
  
  Int_t nBinsTrack[4]={fBinsPt->GetSize()-1,fBinsEta->GetSize()-1,fBinsZv->GetSize()-1,fBinsMultCent->GetSize()-1};
  Double_t minTrack[4]={fBinsPt->GetAt(0),fBinsEta->GetAt(0),fBinsZv->GetAt(0),fBinsMultCent->GetAt(0)};
  Double_t maxTrack[4]={fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsEta->GetSize()-1),fBinsZv->GetAt(fBinsZv->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};
  
  Int_t nBinsEvent[2]={fBinsZv->GetSize()-1,fBinsMultCent->GetSize()-1};
  Double_t minEvent[2]={fBinsZv->GetAt(0),fBinsMultCent->GetAt(0)};
  Double_t maxEvent[2]={fBinsZv->GetAt(fBinsZv->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};  

  fHistTrack = new THnF("fHistTrack", "Histogram for Tracks",4,nBinsTrack,minTrack,maxTrack);
  fHistTrack -> SetBinEdges(0,fBinsPt->GetArray());
  fHistTrack -> SetBinEdges(1,fBinsEta->GetArray());
  fHistTrack -> SetBinEdges(2,fBinsZv->GetArray());
  fHistTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
  fHistTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
  fHistTrack->GetAxis(1)->SetTitle("#eta");
  fHistTrack->GetAxis(2)->SetTitle("Zv (cm)");
  fHistTrack->GetAxis(3)->SetTitle("multiplicity (multCuts)");
  fHistTrack -> Sumw2();
  
  fHistEvent = new THnF("fHistEvent", "Histogram for Events",2,nBinsEvent,minEvent,maxEvent);
  fHistEvent -> SetBinEdges(0,fBinsZv->GetArray());
  fHistEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
  fHistEvent->GetAxis(0)->SetTitle("Zv (cm)");
  fHistEvent->GetAxis(1)->SetTitle("multiplicity (multCuts)");
  fHistEvent -> Sumw2();
  
  if(fIsMC){
    
    fHistMCGenPrimTrack = new THnF("fHistMCGenPrimTrack", "Histogram for generated MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCGenPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCGenPrimTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCGenPrimTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCGenPrimTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCGenPrimTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCGenPrimTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCGenPrimTrack -> Sumw2();
    
    fHistMCRecTrack = new THnF("fHistMCRecTrack", "Histogram for reconstructed MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCRecTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCRecTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCRecTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCRecTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCRecTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCRecTrack -> Sumw2();
    
    fHistMCRecPrimTrack = new THnF("fHistMCRecPrimTrack", "Histogram for reconstructed primary MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCRecPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCRecPrimTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCRecPrimTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecPrimTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCRecPrimTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCRecPrimTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCRecPrimTrack -> Sumw2();
    
    fHistMCRecSecTrack = new THnF("fHistMCRecSecTrack", "Histogram for reconstructed secondary MC Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistMCRecSecTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistMCRecSecTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistMCRecSecTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistMCRecSecTrack->GetAxis(1)->SetTitle("#eta");
    fHistMCRecSecTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistMCRecSecTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
    fHistMCRecSecTrack -> Sumw2();

    fHistMCGenEvent = new THnF("fHistMCGenEvent", "Histogram for generated MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCGenEvent -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCGenEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCGenEvent->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCGenEvent->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCGenEvent -> Sumw2();
    
    fHistMCTrigEvent = new THnF("fHistMCTrigEvent", "Histogram for triggered MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCTrigEvent -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCTrigEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCTrigEvent->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCTrigEvent->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCTrigEvent -> Sumw2();
        
    fHistMCRecEvent = new THnF("fHistMCRecEvent", "Histogram for reconstructed MC  Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCRecEvent -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCRecEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCRecEvent->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCRecEvent->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCRecEvent -> Sumw2();
    
    fHistMCGenINEL0Event = new THnF("fHistMCGenINEL0Event", "Histogram for generated INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCGenINEL0Event -> SetBinEdges(0,fBinsZv->GetArray());
    fHistMCGenINEL0Event -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCGenINEL0Event->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCGenINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCGenINEL0Event -> Sumw2();
    
    fHistMCTrigINEL0Event = new THnF("fHistMCTrigINEL0Event","Histogram for triggered INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
    fHistMCTrigINEL0Event->SetBinEdges(0,fBinsZv->GetArray());
    fHistMCTrigINEL0Event->SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCTrigINEL0Event->GetAxis(0)->SetTitle("Zv (cm)");
    fHistMCTrigINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCTrigINEL0Event->Sumw2();
    
    fHistMCRecINEL0Event = new THnF("fHistMCRecINEL0Event","Histogram for reconstructed INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent); 
    fHistMCRecINEL0Event->SetBinEdges(0,fBinsZv->GetArray());
    fHistMCRecINEL0Event->SetBinEdges(1,fBinsMultCent->GetArray());
    fHistMCRecINEL0Event->GetAxis(0)->SetTitle("mcZv (cm)");
    fHistMCRecINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
    fHistMCRecINEL0Event->Sumw2();
    
//     fHistMCGenTrackINEL0 = new THnF("fHistMCGenTrackINEL0","Histogram for generated tracks for INEL>0 MC Events",4,nBinsTrack,minTrack,maxTrack);
//     fHistMCGenTrackINEL0->SetBinEdges(0,fBinsPt->GetArray());
//     fHistMCGenTrackINEL0->SetBinEdges(1,fBinsEta->GetArray());
//     fHistMCGenTrackINEL0->SetBinEdges(2,fBinsZv->GetArray());
//     fHistMCGenTrackINEL0->SetBinEdges(3,fBinsMultCent->GetArray());
//     fHistMCGenTrackINEL0->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
//     fHistMCGenTrackINEL0->GetAxis(1)->SetTitle("#eta");
//     fHistMCGenTrackINEL0->GetAxis(2)->SetTitle("Zv (cm)");
//     fHistMCGenTrackINEL0->GetAxis(3)->SetTitle("true multiplicity (MC)");
//     fHistMCGenTrackINEL0 -> Sumw2();
  }
    
  fOutputList->Add(fHistTrack);
  fOutputList->Add(fHistEvent);
  
  if(fIsMC){
    fOutputList->Add(fHistMCGenPrimTrack); 
    fOutputList->Add(fHistMCRecTrack);
    fOutputList->Add(fHistMCRecPrimTrack);
    fOutputList->Add(fHistMCRecSecTrack);
    fOutputList->Add(fHistMCGenEvent);
    fOutputList->Add(fHistMCTrigEvent);
    fOutputList->Add(fHistMCRecEvent);
    fOutputList->Add(fHistMCGenINEL0Event);
    fOutputList->Add(fHistMCTrigINEL0Event);
    fOutputList->Add(fHistMCRecINEL0Event);
//     fOutputList->Add(fHistMCGenTrackINEL0);
  }
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::UserExec(Option_t *) 
{
  /// Main loop (called for each event)
  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fEvent) {printf("ERROR: fEvent not available\n"); return;} 

  if(fIsMC)
   {
    fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
    if (!fMCEvent) {printf("ERROR: fMCEvent not available\n"); return;}
   }

  Double_t eventValues[2] = {fEvent->GetPrimaryVertex()->GetZ(),GetEventMultCent(fEvent)};

  // Generated Events in MC
  if(fIsMC) fHistMCGenEvent->Fill(eventValues);
   
  // check if mc event is in inel0 class
  // 1 charged particle in abs(eta)<1.0, pt>0
  Bool_t fIsEventINEL0 = SelectMCEventINEL0(fMCEvent,0,1.0);
  if(fIsEventINEL0) fHistMCGenINEL0Event->Fill(eventValues);
  
  // trigger selection
  Bool_t isEventTriggered = kTRUE;
  AliPhysicsSelection *physicsSelection = NULL;
  AliTriggerAnalysis* triggerAnalysis = NULL;
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (!inputHandler){ Printf("ERROR: Could not receive input handler"); return; }
  
  
  // always MB
  isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;

  physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
  if(!physicsSelection) {Printf("ERROR: Could not receive physicsSelection"); return;}
  if(!isEventTriggered) {Printf("isEventTriggered:   kFALSE"); return;}
  
  // Triggered Events in MC
  if(fIsMC)
  {
   fHistMCTrigEvent->Fill(eventValues);
   if(fIsEventINEL0) fHistMCTrigINEL0Event->Fill(eventValues);
  }  
    /// Apply event cuts
  if(!IsEventAcceptedGeometrics(fEvent)) return;
  if(!IsEventAcceptedQuality(fEvent)) return;
  
  // Reconstructed Events in DATA and MC, events selected for analysis 
  fHistEvent->Fill(eventValues);
  if(fIsMC)
  {
   fHistMCRecEvent->Fill(eventValues);
   if(fIsEventINEL0) fHistMCRecINEL0Event->Fill(eventValues);
  }
  
  /// Generated paricle loop for MC
  fMCStack=0;
  if (fIsMC){
//     fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
//     if (!fMCEvent) {printf("ERROR: fMCEvent not available\n"); return;}
    
    fMCStack = fMCEvent->Stack();
    if (!fMCStack) {printf("ERROR: fMCStack not available\n"); return;}
       
    for (Int_t iParticle = 0; iParticle < fMCStack->GetNtrack(); iParticle++){
       
      TParticle *mcGenParticle = fMCStack->Particle(iParticle);
  
      
      /// \li Acceptance cuts for generated particles
      
      if(!mcGenParticle) {printf("ERROR: mcGenParticle  not available\n"); continue;}
      if(!IsTrackAcceptedKinematics(mcGenParticle)) continue;     
      
      if(IsSelectedParticle(iParticle)){						/// \li TODO MC label???
        Double_t mcGenPrimTrackValue[4] = {mcGenParticle->Pt(), mcGenParticle->Eta(),fEvent->GetPrimaryVertex()->GetZ(),GetEventMultCent(fEvent)};
        fHistMCGenPrimTrack->Fill(mcGenPrimTrackValue);
// 	if(fIsEventINEL0) fHistMCGenTrackINEL0->Fill(mcGenPrimTrackValue);
      }
  
    }    
  }

  
  /// Track loop for Data and MC
  for (Int_t iTracks = 0; iTracks < fEvent ->GetNumberOfTracks(); iTracks++) {
   AliVTrack *track = fEvent->GetVTrack(iTracks);				/// \li TODO GetTrack? GetVTrack?    
   
   if (!track){ printf("ERROR: Could not receive track %d\n", iTracks); continue; }
   
    /// \li Trackcuts
    if(!IsTrackAcceptedKinematics(track)) continue;
    if(!IsTrackAcceptedQuality(track)) continue;
    
    if( fUseGeomCut == kTRUE && !(IsTrackAcceptedGeometricalCut(track, fEvent->GetMagneticField())) ) continue;

    
    /// \li Fill Track Histograms
    Double_t trackValues[4] = {track->Pt(), track->Eta(), fEvent->GetPrimaryVertex()->GetZ(), GetEventMultCent(fEvent)}; 
    fHistTrack->Fill(trackValues);
    
    /// Monte-Carlo reconstructed loop
    if(fIsMC){
      Int_t mcLabel = TMath::Abs(track->GetLabel());				///  \li TODO Protection do negative values have a meening?
      
      if (!fMCStack) {return;}
      TParticle *mcParticle = fMCStack->Particle(mcLabel);
      if(!mcParticle) {printf("ERROR: mcParticle not available\n"); continue;}
      
      if(!IsTrackAcceptedKinematics(mcParticle)) continue;			/// \li TODO In principal, this should not be needed.
      Double_t mcRecTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta() ,fEvent->GetPrimaryVertex()->GetZ(), GetEventMultCent(fEvent)};
      fHistMCRecTrack->Fill(mcRecTrackValue);
      
      if(IsSelectedParticle(mcLabel)){
	Double_t mcPrimTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), fEvent->GetPrimaryVertex()->GetZ() ,GetEventMultCent(fEvent)};
        fHistMCRecPrimTrack->Fill(mcPrimTrackValue);
      }
//       SetMCParticleType(AlidNdPtUnifiedAnalysisTask::kPrimary);
      if(IsSecondary(mcLabel))
      {
	Double_t mcSecTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), fEvent->GetPrimaryVertex()->GetZ() ,GetEventMultCent(fEvent)};
        fHistMCRecSecTrack->Fill(mcSecTrackValue);
      }
    }
    
  } //track loop 
  PostData(1, fOutputList);
} 

//________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::Terminate(Option_t *) 
{

}

/// Track Acceptance cuts, for tracks.
///
/// \param AliVTrack Input track
///
/// \return Is track accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedKinematics(AliVTrack *track) 
{

  if(!track) return kFALSE;

  Float_t eta = track->Eta();
  Float_t pt = track->Pt();

  if(eta < fMinEta) return kFALSE;
  if(eta > fMaxEta) return kFALSE;
  if(pt < fMinPt) return kFALSE;
  if(pt > fMaxPt) return kFALSE;

  return kTRUE;
}
/// Track Acceptance cuts, for MC particles.
///
/// \param TParticle Input particle
///
/// \return Is particle accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedKinematics(TParticle *mcTrack) 
{
  if(!mcTrack) return kFALSE;

  Float_t eta = mcTrack->Eta();
  Float_t pt = mcTrack->Pt();

  if(eta < fMinEta) return kFALSE;
  if(eta > fMaxEta) return kFALSE;
  if(pt < fMinPt) return kFALSE;
  if(pt > fMaxPt) return kFALSE;
  return kTRUE;
}
/// Track Quality cuts.
///
/// \param AliVTrack Input track
///
/// \return Is track accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedQuality(AliVTrack *track) 
{
  if(fIsESD){
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*> (track);
    if(!fESDtrackCutsInit) InitESDTrackCuts();			/// Check if track cuts were already initialized
    if(!fESDtrackCuts->AcceptTrack(esdTrack)) return kFALSE;
  }
  return kTRUE;
}


/// Event Acceptance cuts.
///
/// \param AliVEvent Input event
///
/// \return Is event accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAcceptedGeometrics(AliVEvent *event) 
{
  if(!fEventCutsInit) InitdNdPtEventCuts();
  if(TMath::Abs(event->GetPrimaryVertex()->GetZ())>fZvtx) return kFALSE;
  return kTRUE;
}
/// Event Quality cuts.
///
/// \param AliVEvent Input event
///
/// \return Is event accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAcceptedQuality(AliVEvent *event) 
{
//   AliVVertex *vertex = event->GetPrimaryVertexTracks();
  if(!event) return kFALSE;
  if(!IsVertexOK(event)) return kFALSE;
  return kTRUE;
}
/// Select primary charged particles
///
/// \param Int_t MC lable of particle
///
/// \return Is particle accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsSelectedParticle(Int_t mcLabel)
{
  if(!fMCStack->IsPhysicalPrimary(mcLabel)) return kFALSE;				/// Reject secondary particles
  if ( TMath::Abs(fMCStack->Particle(mcLabel)->GetPDG()->Charge()/3.) < 0.01 )return kFALSE;	/// Reject neutral particles
  Int_t ipdg = TMath::Abs(fMCStack->Particle(mcLabel)->GetPdgCode());
  switch(fParticleLabel)
  {
    case AlidNdPtUnifiedAnalysisTask::kPrimary: /// Accept all primary particles
     break;
    case AlidNdPtUnifiedAnalysisTask::kPion: /// accept pions
     if(ipdg!=211) return kFALSE; break;				/// Reject other particles then pions
    case AlidNdPtUnifiedAnalysisTask::kKaon: /// accept kaons
      if(ipdg!=321) return kFALSE; break;				/// Reject other particles then kaons
    case AlidNdPtUnifiedAnalysisTask::kProtons: /// accept protons
      if(ipdg!=2212) return kFALSE; break;				/// Reject other particles then protons
    case AlidNdPtUnifiedAnalysisTask::kRest:
    if(ipdg==211 || ipdg==321 || ipdg==2212) return kFALSE;  break;			/// Select all other particles
  }
//   printf("PDG Code: %d\n", ipdg);
  return kTRUE;
}

/// Select secondary charged particles
///
/// \param Int_t MC lable of particle
///
/// \return Is particle accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsSecondary(Int_t mcLabel)
{
  if(fMCStack->IsPhysicalPrimary(mcLabel)) return kFALSE;				/// Reject primary particles
  if(TMath::Abs(fMCStack->Particle(mcLabel)->GetPDG()->Charge()/3.) < 0.01 )return kFALSE;	/// Reject neutral particles
  return kTRUE;
}

/// Function to either fill Multiplicity or centrality
///
/// \param AliVEvent event to be analised
///
/// \return Double_t with centrality or multiplicity
Double_t AlidNdPtUnifiedAnalysisTask::GetEventMultCent(AliVEvent *event)
{
  if(fUseMultiplicity){
    AliVMultiplicity* multiplicity = event->GetMultiplicity();
    if(!multiplicity) {printf("ERROR: multiplicity not available\n"); return 999;}
    Int_t mult = multiplicity->GetNumberOfTracklets();
    return mult;
  }else{
    return 999;
    
  }
}
/// Function to initialize the ESD track cuts, if needed
void AlidNdPtUnifiedAnalysisTask::InitESDTrackCuts(){

  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
  if(!fESDtrackCuts) {printf("ERROR: fESDtrackCuts not available\n"); return;}

  fESDtrackCuts->SetRequireTPCRefit(fTPCRefit);
  fESDtrackCuts->SetRequireITSRefit(fITSRefit);
  fESDtrackCuts->SetAcceptKinkDaughters(fAcceptKinks);
  if(fminNCrossedRowsTPC >0) fESDtrackCuts->SetMinNCrossedRowsTPC(fminNCrossedRowsTPC);
  if(fminRatioCrossedRowsOverFindableClustersTPC >0) fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(fminRatioCrossedRowsOverFindableClustersTPC);
  if(fmaxFractionSharedTPCCluster >0) fESDtrackCuts->SetMaxFractionSharedTPCClusters(fmaxFractionSharedTPCCluster);
  if(fmaxchi2perTPCcl >0) fESDtrackCuts->SetMaxChi2PerClusterTPC(fmaxchi2perTPCcl);
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);///TODO
  if(fmaxchi2perITScl >0) fESDtrackCuts->SetMaxChi2PerClusterITS(fmaxchi2perITScl);

  if(fdcatoVertex2D >0) fESDtrackCuts->SetDCAToVertex2D(fdcatoVertex2D);
  if(fsigmatovertex >0) fESDtrackCuts->SetRequireSigmaToVertex(fsigmatovertex);
  if(fmaxdcazITSTPC >0) fESDtrackCuts->SetMaxDCAToVertexZ(fmaxdcazITSTPC);
  fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01"); /// TODO
  if(fmaxchi2TPCconstrained >0) fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(fmaxchi2TPCconstrained);
  if(fminActiveLength>0) fESDtrackCuts->SetMinLengthActiveVolumeTPC(fminActiveLength);
  
  fESDtrackCutsInit = kTRUE;
}
/// Function to initialize the dNdPtEventCuts
void AlidNdPtUnifiedAnalysisTask::InitdNdPtEventCuts(){
  
  fEventCuts = new AlidNdPtEventCuts("AlidNdPtEventCuts","Event cuts");
  if(!fEventCuts) {printf("ERROR: fEventCuts not available\n"); return;}
  
  fEventCuts->SetMeanXYZv(fMeanXYZv[0],fMeanXYZv[1],fMeanXYZv[2]);
  fEventCuts->SetSigmaMeanXYZv(fSigmaMeanXYZv[0],fSigmaMeanXYZv[1],fSigmaMeanXYZv[2]);
  fEventCuts->SetTriggerRequired(fEventTriggerRequired);
  
  fEventCutsInit = kTRUE;
   
}

///Function to check vertex quality
Bool_t AlidNdPtUnifiedAnalysisTask::IsVertexOK(AliVEvent *event){
  
  if(fIsESD){
    Float_t requiredZResolution = 1000;
    AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(event);
    const AliESDVertex *esdVertex = ESDevent->GetPrimaryVertexTracks();
    if(!esdVertex){printf("ERROR: vertex not available\n"); return kFALSE;}
    if(esdVertex->GetNContributors()<1) {
    // SPD vertex
      esdVertex = ESDevent->GetPrimaryVertexSPD();
    }
//     AliESDVertex *esdVertex = dynamic_cast<AliESDVertex*> (vertex);
    if(!esdVertex->GetStatus()){return kFALSE;}
    Double_t zRes = esdVertex->GetZRes();
    if (zRes > requiredZResolution) return kFALSE;

    const AliESDVertex *vertexSPD = ESDevent->GetPrimaryVertexSPD();
    // always check for SPD vertex
    if(!vertexSPD) return kFALSE;
    if(!vertexSPD->GetStatus()) return kFALSE;
    if (vertexSPD->IsFromVertexerZ() && vertexSPD->GetDispersion() > 0.02) return kFALSE;
  
  }else{
    //AOD code goes here
  }
  return kTRUE;
}

Bool_t AlidNdPtUnifiedAnalysisTask::SelectMCEventINEL0(AliMCEvent* mcEvent, Double_t ptmin, Double_t etarange){
//
// select INEL>0 events with at least
// one prompt (MC primary) particle in acceptance
// pT>0, |eta|<1.0 for normalization
//

if(!mcEvent) return kFALSE;
// AliMCEvent* fmcEvent = static_cast<AliMCEvent*>(mcEvent);
AliStack* stack = mcEvent->Stack(); 
if(!stack) return kFALSE;

  Int_t count = 0;
  for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc) 
  {
    TParticle* particle = stack->Particle(iMc);
    if (!particle) continue;

    // only charged particles
    if(!particle->GetPDG()) continue;
    Double_t charge = particle->GetPDG()->Charge()/3.;
    if(charge == 0) continue;

    // physical primary
    Bool_t prim = stack->IsPhysicalPrimary(iMc);
    if(!prim) continue;

    if(particle->Pt() < ptmin) continue;
    if(TMath::Abs(particle->Eta()) > etarange) continue;

    count++;
  }

  if(count > 0) return kTRUE;
  else return kFALSE;

}

/// Function to apply the TPC geometrical cut
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedGeometricalCut(AliVTrack *track, Double_t bMagZ){

  
  if(!track) return kFALSE;
  
  if(track->Charge()==0) { return kFALSE; }
  
  //
  // as done in AliAnalysisTaskFragmentationFunction
  //
  
  Short_t sign = track->Charge();
  Double_t xyz[50];
  Double_t pxpypz[50];
  Double_t cv[21];
  
  for(Int_t i = 0; i < 21; i++) cv[i] = 0;
  for(Int_t i = 0; i < 50; i++) xyz[i] = 0;
  for(Int_t i = 0; i < 50; i++) pxpypz[i] = 0;
  
  track->GetXYZ(xyz);
  track->GetPxPyPz(pxpypz);
  track->GetCovarianceXYZPxPyPz(cv);
  
  AliExternalTrackParam par(xyz, pxpypz, cv, sign);
  static AliESDtrack dummy;
 
  Double_t dLengthInTPC = 0;
  dLengthInTPC = dummy.GetLengthInActiveZone(&par,3,236, bMagZ ,0,0);

  
  // cut on length
  if( dLengthInTPC < (130-5*TMath::Abs(1./track->Pt())) ) { return kFALSE; }
  
  return kTRUE;
}