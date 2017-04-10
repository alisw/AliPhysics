#include "TCanvas.h"
#include "TParticle.h"
#include "TH1I.h"
#include <TDatabasePDG.h>

#include "AliAnalysisDataSlot.h"
#include "AliMultSelection.h"
#include "AliAnalysisDataContainer.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAODTracklets.h"
#include "AliMultiplicity.h"
#include "AliCFManager.h"
#include "AliCFCutBase.h"
#include "AliCFContainer.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliLog.h"
#include "AliGenEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliCentrality.h"

#include "AliSingleTrackEffCuts.h"
#include "AliCFSingleTrackEfficiencyTask.h"


ClassImp(AliCFSingleTrackEfficiencyTask)

//__________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask() :
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fTrackCuts(0x0),
  fMCCuts(0x0),
  fTriggerMask(AliVEvent::kAny),
  fSetFilterBit(kFALSE),
  fbit(0),
  fRemoveNegativeLabelTracks(kTRUE),
  fMatchToKinematicTrack(kTRUE),
  fUseGeneratedKine(kFALSE),
  fEvalCentrality(kFALSE),
  fCentralityEstimator("V0M"),
  fConfiguration(kFast), // default  use the minimal configuration
  fHistEventsProcessed(0x0)
{
  //
  //Default constructor
  //
}

//___________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask(const Char_t* name,AliESDtrackCuts *trackcuts, AliSingleTrackEffCuts * mccuts) :
  AliAnalysisTaskSE(name),
  fReadAODData(0),
  fCFManager(0x0),
  fQAHistList(0x0),
  fTrackCuts(trackcuts),
  fMCCuts(mccuts),
  fTriggerMask(AliVEvent::kAny),
  fSetFilterBit(kFALSE),
  fbit(0),
  fRemoveNegativeLabelTracks(kTRUE),
  fMatchToKinematicTrack(kTRUE),
  fUseGeneratedKine(kFALSE),
  fEvalCentrality(kFALSE),
  fCentralityEstimator("V0M"),
  fConfiguration(kFast), // default  use the minimal configuration
  fHistEventsProcessed(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliCFSingleTrackEfficiencyTask","Calling Constructor");

  // Output slot #1 writes into a TList container (nevents histogran)
  DefineOutput(1,TH1I::Class());
  // Output slot #2 writes into a TList container (distributions)
  DefineOutput(2,AliCFContainer::Class());
  // Output slot #3 writes the QA list
  DefineOutput(3,TList::Class());
  // Output slot #3 writes the ESD track cuts used
  DefineOutput(4,AliESDtrackCuts::Class());
  // Output slot #4 writes the particle and event selection object
  DefineOutput(5,AliSingleTrackEffCuts::Class());
}

//_________________________________________________________________________________________________________________
AliCFSingleTrackEfficiencyTask& AliCFSingleTrackEfficiencyTask::operator=(const AliCFSingleTrackEfficiencyTask& c) 
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;

    fReadAODData = c.fReadAODData;
    fCFManager  = c.fCFManager;
    fQAHistList = c.fQAHistList;

    if(c.fTrackCuts) { delete fTrackCuts; fTrackCuts = new AliESDtrackCuts(*(c.fTrackCuts)); }
    if(c.fMCCuts) { delete fMCCuts; fMCCuts = new AliSingleTrackEffCuts(*(c.fMCCuts)); }
    fTriggerMask = c.fTriggerMask;

    fSetFilterBit  = c.fSetFilterBit;
    fbit = c.fbit;
    fRemoveNegativeLabelTracks = c.fRemoveNegativeLabelTracks;
    fMatchToKinematicTrack = c.fMatchToKinematicTrack;
    fUseGeneratedKine = c.fUseGeneratedKine;

    fEvalCentrality = c.fEvalCentrality;
    fCentralityEstimator = c.fCentralityEstimator;
      
    fConfiguration = c.fConfiguration;

    fHistEventsProcessed = c.fHistEventsProcessed;
  }
  return *this;
}

//________________________________________________________________________________________________________
AliCFSingleTrackEfficiencyTask::AliCFSingleTrackEfficiencyTask(const AliCFSingleTrackEfficiencyTask& c) :
  AliAnalysisTaskSE(c),
  fReadAODData(c.fReadAODData),
  fCFManager(c.fCFManager),
  fQAHistList(c.fQAHistList),
  fTrackCuts(c.fTrackCuts),
  fMCCuts(c.fMCCuts),
  fTriggerMask(c.fTriggerMask),
  fSetFilterBit(c.fSetFilterBit),
  fbit(c.fbit),
  fRemoveNegativeLabelTracks(c.fRemoveNegativeLabelTracks),
  fMatchToKinematicTrack(c.fMatchToKinematicTrack),
  fUseGeneratedKine(c.fUseGeneratedKine),
  fEvalCentrality(c.fEvalCentrality),
  fCentralityEstimator(c.fCentralityEstimator),
  fConfiguration(c.fConfiguration),
  fHistEventsProcessed(c.fHistEventsProcessed)
{
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliCFSingleTrackEfficiencyTask::~AliCFSingleTrackEfficiencyTask()
{
  //
  // Destructor
  //
  Info("~AliCFSingleTrackEfficiencyTask","Calling Destructor");

  if (fCFManager)           delete fCFManager;
  if (fHistEventsProcessed) delete fHistEventsProcessed;
  if (fQAHistList) { fQAHistList->Clear(); delete fQAHistList; }
  if (fTrackCuts)           delete fTrackCuts;
  if (fMCCuts)              delete fMCCuts;
}

//_________________________________________________
void AliCFSingleTrackEfficiencyTask::Init() 
{
  //
  // Initialization, checks + copy cuts
  //
  if(!fMCCuts) {
    AliFatal(" MC Cuts not defined");
    return;
  }
  if(!fTrackCuts) {
    AliFatal(" Track Cuts not defined");
    return;
  }

  AliESDtrackCuts* copyfTrackCuts = new AliESDtrackCuts(*fTrackCuts);
  const char* nameoutput = GetOutputSlot(4)->GetContainer()->GetName();
  copyfTrackCuts->SetName(nameoutput);
  // Post the data
  PostData(4,copyfTrackCuts);

  AliSingleTrackEffCuts* copyfMCCuts = new AliSingleTrackEffCuts(*fMCCuts);
  nameoutput = GetOutputSlot(5)->GetContainer()->GetName();
  copyfMCCuts->SetName(nameoutput);
  // Post the data
  PostData(5,copyfMCCuts);


  if(fEvalCentrality) {
    Bool_t isCentEstimatorOk = kFALSE;
    TString validEstimators[9] = { "V0M", "V0A", "V0C", "TRK", "TKL", "CL1", "ZNA", "ZNC" "ZPA" };
    for(Int_t i=0; i<9; i++ ) { 
      if(fCentralityEstimator==validEstimators[i]) isCentEstimatorOk = kTRUE; 
    }
    if(!isCentEstimatorOk) {
      AliFatal(Form("Chosen centrality estimator %s is not valid\n",fCentralityEstimator.Data()));
      return;
    }
  }

  return;
}

//_________________________________________________________________
void AliCFSingleTrackEfficiencyTask::UserExec(Option_t *)
{
  //
  // User Exec
  //

  Info("UserExec","Start of method") ;

  AliVEvent* event = fInputEvent;

  if(!fInputEvent) {
    AliFatal("NO EVENT FOUND!");
    return;
  }
  if(!fMCEvent) {
    AliFatal("NO MC INFO FOUND");
    return;
  }

  fHistEventsProcessed->Fill(0.5); // # of Event proceed        
  Bool_t IsEventMCSelected = kFALSE;
  Bool_t isAOD = fInputEvent->IsA()->InheritsFrom("AliAODEvent");

  //
  // Step 0: MC Gen Event Selection
  //
  if(isAOD) {
    if(!event && AODEvent() && IsStandardAOD()) {
      // In case there is an AOD handler writing a standard AOD, use the AOD 
      // event in memory rather than the input (ESD) event.    
      event = dynamic_cast<AliAODEvent*> (AODEvent());
    }
    IsEventMCSelected = fMCCuts->IsMCEventSelected(event);//AODs
  } else {
    IsEventMCSelected = fMCCuts->IsMCEventSelected(fMCEvent);//ESDs
  }

  // pass to the manager the event info to the cuts that need it 
  if(!isAOD) fCFManager->SetMCEventInfo(fMCEvent);
  else fCFManager->SetMCEventInfo(event);
  fCFManager->SetRecEventInfo(event);

  if(!IsEventMCSelected) {
    AliDebug(3,"MC event not passing the quality criteria \n");
    PostData(1,fHistEventsProcessed);
    PostData(2,fCFManager->GetParticleContainer());
    PostData(3,fQAHistList);
    return;
  }
  fHistEventsProcessed->Fill(1.5); // # of Event after passing MC cuts
       

  //
  // Step 1-3: Check the MC generated particles
  // 
  if(!isAOD) CheckESDMCParticles();
  else CheckAODMCParticles();

  //
  // Step 4-7: Reconstructed event and track selection
  //
  Bool_t isRecoEventOk = fMCCuts->IsRecoEventSelected(event);

  if(isRecoEventOk) {
    fHistEventsProcessed->Fill(2.5); // # of Event after passing all cuts
    CheckReconstructedParticles();
  }
  else AliDebug(3,"Event not passing quality criteria\n");


  PostData(1,fHistEventsProcessed);
  PostData(2,fCFManager->GetParticleContainer());
  PostData(3,fQAHistList);

  return;
}


//___________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::Terminate(Option_t*)
{
  //
  // Terminate
  //

  Info("Terminate","Start and end of Method");
  AliAnalysisTaskSE::Terminate();

  /*
  //draw some example plots....
  AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
       
  TH1D* h00 =   cont->ShowProjection(5,0) ;
  TH1D* h01 =   cont->ShowProjection(5,1) ;
  TH1D* h02 =   cont->ShowProjection(5,2) ;
  TH1D* h03 =   cont->ShowProjection(5,3) ;
  TH1D* h04 =   cont->ShowProjection(5,4) ;
  TH1D* h05 =   cont->ShowProjection(5,5) ;
  TH1D* h06 =   cont->ShowProjection(5,6) ;
  TH1D* h07 =   cont->ShowProjection(5,7) ;
       
  h00->SetMarkerStyle(23) ;
  h01->SetMarkerStyle(24) ;
  h02->SetMarkerStyle(25) ;
  h03->SetMarkerStyle(26) ;
  h04->SetMarkerStyle(27) ;
  h05->SetMarkerStyle(28) ;
  h06->SetMarkerStyle(28) ;
       
       
       
  TCanvas * c =new TCanvas("c","",1400,800);
  c->Divide(4,2);
       
  c->cd(1);
  h00->Draw("p");
  c->cd(2);
  h01->Draw("p");
  c->cd(3);
  h02->Draw("p");
  c->cd(4);
  h03->Draw("p");
  c->cd(5);
  h04->Draw("p");
  c->cd(6);
  h05->Draw("p");
  c->cd(7);
  h06->Draw("p");
  c->cd(8);
  h07->Draw("p");
       
  c->SaveAs("plots.eps");
  */
}


//___________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::UserCreateOutputObjects() 
{
  //
  // UserCreateOutputObjects
  //
  
  Info("CreateOutputObjects","CreateOutputObjects of task %s", GetName());

  //slot #1
  OpenFile(1);
       
  const char* nameoutput=GetOutputSlot(1)->GetContainer()->GetName();
  //       fHistEventsProcessed = new TH1I(nameoutput"fHistEventsProcessed","fHistEventsProcessed",3,0,3) ;
  fHistEventsProcessed = new TH1I(nameoutput,"fHistEventsProcessed",3,0,3) ;
  fHistEventsProcessed->GetXaxis()->SetBinLabel(1,"All events");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(2,"Good MC events");
  fHistEventsProcessed->GetXaxis()->SetBinLabel(3,"Good Reconstructed events");
    
  PostData(1,fHistEventsProcessed) ;
  PostData(2,fCFManager->GetParticleContainer()) ;
  PostData(3,fQAHistList) ;
       
  return;
}


//_________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::CheckESDMCParticles()
{
  //
  // Check ESD generated particles
  //
  if (!fMCEvent) {
    AliFatal("NO MC INFO FOUND");
    return;
  }

  Double_t containerInput[7] = { 0., 0., 0., 0., 0., 0., 0. };

  TArrayF vtxPos(3);
  AliGenEventHeader *genHeader = NULL;
  genHeader = fMCEvent->GenEventHeader();
  genHeader->PrimaryVertex(vtxPos);
  containerInput[4] = vtxPos[2]; // z-vtx position

  Double_t multiplicity = (Double_t)GetNumberOfTrackletsInEtaRange(-1.0,1.0);
  containerInput[5] = multiplicity; //reconstructed number of tracklets

  Double_t centrality = -1.;
  if(fEvalCentrality) centrality = GetCentrality();
  containerInput[6] = centrality;

  // loop on the MC Generated Particles
  for (Int_t ipart=0; ipart<fMCEvent->GetNumberOfTracks(); ipart++) {

    AliMCParticle *mcPart = (AliMCParticle*)fMCEvent->GetTrack(ipart);
    containerInput[0] = (Float_t)mcPart->Pt();
    containerInput[1] = mcPart->Eta();
    containerInput[2] = mcPart->Phi();
    containerInput[3] = mcPart->Theta();

    // Step 1. Particle passing through Generation criteria and filling
    if( !fMCCuts->IsMCParticleGenerated(mcPart) ) {
      AliDebug(3,"MC Particle not passing through genetations criteria\n");
      continue;
    }
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCGenCut);

    // Step 2. Particle passing through Kinematic criteria and filling
    if( !fMCCuts->IsMCParticleInKineAcceptance(mcPart) ) {
      AliDebug(3,"MC Particle not in the kine acceptance\n");
      continue;
    }
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCKineCut);

    // Step 3. Particle passing through Track ref criteria and filling
    // did leave signal (enough clusters) on the detector
    if( !fMCCuts->IsMCParticleInReconstructable(mcPart) ) {
      AliDebug(3,"MC Particle not in the reconstructible\n");
      continue;
    }
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCAccpCut);

  }// end of particle loop

  return;
}


//________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::CheckAODMCParticles()
{
  //
  // Check AOD generated particles
  //
  if (!fInputEvent) {
    AliFatal("NO EVENT FOUND!");
    return;
  }

  AliAODEvent* event = dynamic_cast<AliAODEvent*>(fInputEvent);
  if(!event && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    event = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  TClonesArray* mcArray = dynamic_cast<TClonesArray*>(event->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcArray) {
    AliError("Could not find Monte-Carlo in AOD");
    return;
  }
  AliAODMCHeader *mcHeader=NULL;
  mcHeader = dynamic_cast<AliAODMCHeader*>(event->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  if (!mcHeader) {
    AliError("Could not find MC Header in AOD");
    return;
  }

  Double_t containerInput[7] = { 0., 0., 0., 0., 0., 0., 0. };
  // Set the z-vertex position
  containerInput[4]  = mcHeader->GetVtxZ();
  // Multiplicity of the event defined as Ntracklets |eta|<1.0
  Double_t multiplicity = (Double_t)GetNumberOfTrackletsInEtaRange(-1.0,1.0);
  containerInput[5] = multiplicity;
  // Determine the event centrality
  Double_t centrality = 0.;
  if(fEvalCentrality) centrality = GetCentrality();
  containerInput[6] = centrality;

  for (Int_t ipart=0; ipart<mcArray->GetEntriesFast(); ipart++) {

    AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(ipart));
    if (!mcPart){
      AliError("Failed casting particle from MC array!, Skipping particle");
      continue;
    }
    containerInput[0] = (Float_t)mcPart->Pt();
    containerInput[1] = mcPart->Eta();
    containerInput[2] = mcPart->Phi();
    containerInput[3] = mcPart->Theta();

    // Step 1. Particle passing through Generation criteria and filling
    if( !fMCCuts->IsMCParticleGenerated(mcPart) ) {
      AliDebug(3,"MC Particle not passing quality criteria\n");
      continue;
    }
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCGenCut);

    // Step 2. Particle passing through Kinematic criteria and filling
    if( !fMCCuts->IsMCParticleInKineAcceptance(mcPart) ) {
      AliDebug(3,"MC Particle not in the acceptance\n");
      continue;
    }
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCKineCut);

    // Step 3. Particle passing through Track Ref criteria and filling
    // but no info available for Track ref in AOD fillng same as above
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInput,kStepMCAccpCut);
      
  }

  return;
}


//_______________________________________________________________________________
AliESDtrack * AliCFSingleTrackEfficiencyTask::ConvertTrack(AliAODTrack *track)
{
  //
  // Convert an AOD track to an  ESD track to apply ESDtrackCuts
  //

  AliAODEvent *aodEvent = (AliAODEvent*)fInputEvent;
  const AliAODVertex *primary = aodEvent->GetPrimaryVertex();
  Double_t pos[3],cov[6];
  primary->GetXYZ(pos);
  primary->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  AliESDtrack *esdTrack =  new AliESDtrack(track);
  // set the TPC cluster info
  esdTrack->SetTPCClusterMap(track->GetTPCClusterMap());
  esdTrack->SetTPCSharedMap(track->GetTPCSharedMap());
  esdTrack->SetTPCPointsF(track->GetTPCNclsF());
  // needed to calculate the impact parameters
  esdTrack->RelateToVertex(&vESD,0.,3.);
  //  std::cout << " primary vtx "<< primary << std::endl;
  //  std::cout << " esdvtx "<< vESD.GetName() << std::endl;
  //  std::cout<< " esdtrack pt "<< esdTrack.Pt() << " and status " << esdTrack.GetStatus() <<endl;
  //  std::cout << " aod track "<< track<< " and status " << track->GetStatus() << std::endl;
  //  std::cout << " esd track "<< esdTrack<< " and status " << esdTrack->GetStatus() << std::endl;
  return esdTrack;
}


//___________________________________________________________________________
void AliCFSingleTrackEfficiencyTask::CheckReconstructedParticles()
{
  //
  // Check reconstructed particles
  //

  AliVEvent* event = fInputEvent;
  Bool_t isAOD = fInputEvent->IsA()->InheritsFrom("AliAODEvent");
  if(!event && AODEvent() && IsStandardAOD()) {
   event  = dynamic_cast<AliAODEvent*> (AODEvent());
  }
  if (!event) return;

  Double_t containerInput[7] = { 0, 0, 0, 0, 0, 0, 0};   // contains reconstructed quantities
  Double_t containerInputMC[7] = { 0, 0, 0, 0, 0, 0, 0}; // contains generated quantities
  
  const AliVVertex *vertex = event->GetPrimaryVertex();
  // set the z-vertex position
  containerInput[4] = vertex->GetZ();
  containerInputMC[4] = containerInput[4];
  // set the event multiplicity as Ntracklets in |eta|<1.0
  Double_t multiplicity = (Double_t)GetNumberOfTrackletsInEtaRange(-1.0,1.0);
  containerInput[5] = multiplicity;
  containerInputMC[5] = multiplicity;
  // Determine the event centrality
  Double_t centrality = 0.;
  if(fEvalCentrality) centrality = GetCentrality();
  containerInput[6] = centrality;
  containerInputMC[6] = centrality;

  // Reco tracks track loop
  AliVParticle* track = NULL;
  for (Int_t iTrack = 0; iTrack<event->GetNumberOfTracks(); iTrack++) {

    track = event->GetTrack(iTrack);
    if(!track) {
      AliDebug(3,Form("Track %d not found",iTrack));
      continue;
    }

    Double_t mom[3];
    track->PxPyPz(mom);
    Double_t pt = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
    containerInput[0] = pt;
    containerInput[1] = track->Eta();
    containerInput[2] = track->Phi();
    containerInput[3] = track->Theta();

    //
    // Step 4. Track that are recostructed, filling
    //
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed);

    //
    // Step 5. Track that are recostructed and pass acceptance cuts, filling
    //
    if (!fMCCuts->IsRecoParticleKineAcceptance(track)) continue;
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoKineCuts);

    // is track associated to particle ? if yes + implimenting the physical primary..
    Int_t label = TMath::Abs(track->GetLabel());
    if(fRemoveNegativeLabelTracks) {
        if (label<=0) {
            AliDebug(3,"Particle not matching MC label \n");
            continue;
        }
    }

    //
    // Step 6-7. Track that are recostructed + Kine + Quality criteria, filling
    //

    // check particle selections at MC level
    if(fMatchToKinematicTrack){
        AliVParticle *mcPart  = (AliVParticle*)fMCEvent->GetTrack(TMath::Abs(label));
        if(!mcPart) continue;
        containerInputMC[0] = (Float_t)mcPart->Pt();
        containerInputMC[1] = mcPart->Eta();
        containerInputMC[2] = mcPart->Phi();
        containerInputMC[3] = mcPart->Theta();

        if (!fMCCuts->IsMCParticleGenerated(mcPart)) continue;
        //    cout<< "MC matching did work"<<endl;
	if(fUseGeneratedKine){
	  for(Int_t ivar=0; ivar<=3; ivar++) containerInput[ivar] = containerInputMC[ivar];
	}
    }

    // for filter bit selection
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(track);
    if(isAOD && !aodTrack) continue;
    if(isAOD && fSetFilterBit) if (!aodTrack->TestFilterMask(BIT(fbit))) continue;
    //    cout<<" Filter bit check passed"<<endl;

    Bool_t isESDtrack = track->IsA()->InheritsFrom("AliESDtrack");
    AliESDtrack *tmptrack = NULL;
    if(isESDtrack) {
      tmptrack = dynamic_cast<AliESDtrack*>(track); // ESDs
    } else {
      if (aodTrack) tmptrack = ConvertTrack(aodTrack); // AODs
    }
    if (!tmptrack) continue;

    // exclude global constrained and TPC only tracks (filter bits 128 and 512)
    Int_t id = tmptrack->GetID();
    if(isAOD && id<0) {
      AliDebug(3,"Track removed bc corresponds to either filter bit 128 or 512 (TPC only tracks)\n");
      delete tmptrack; tmptrack=NULL;
      continue;
    }

    // Apply ESD track cuts
    if( !fTrackCuts->IsSelected(tmptrack) ){
      AliDebug(3,"Reconstructed track not passing quality criteria\n");
      if(isAOD) { delete tmptrack; tmptrack=NULL; }
      continue;
    }
    //    cout<<" analysis cuts passed"<<endl;
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepReconstructedMC);
    fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoQualityCuts);
    //    cout << " checking particle pid"<<endl;

    //
    // Step8, PID check
    //
    if( !fMCCuts->IsRecoParticlePID(track) ){
      AliDebug(3,"Reconstructed track not passing PID criteria\n");
      if(isAOD) { delete tmptrack; tmptrack=NULL; }
      continue;
    }
    if(fConfiguration!=kFast) fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoPID);
    //    cout << " all steps filled"<<endl;

    if(isAOD) { delete tmptrack; tmptrack=NULL; }
  }
  return;
}

//______________________________________________________________________
Int_t AliCFSingleTrackEfficiencyTask::GetNumberOfTrackletsInEtaRange(Double_t mineta, Double_t maxeta)
{
  //
  // counts tracklets in given eta range
  //

  AliAODEvent* event = dynamic_cast<AliAODEvent*>(fInputEvent);
  Bool_t isAOD = fInputEvent->IsA()->InheritsFrom("AliAODEvent");
  if(!event && AODEvent() && IsStandardAOD()) {
   event  = dynamic_cast<AliAODEvent*> (AODEvent());
   if (!event) return -1;
  }
  Int_t count=0;

  if(isAOD) {
    AliAODTracklets* tracklets = event->GetTracklets();
    Int_t nTr=tracklets->GetNumberOfTracklets();
    for(Int_t iTr=0; iTr<nTr; iTr++){
      Double_t theta=tracklets->GetTheta(iTr);
      Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
      if(eta>mineta && eta<maxeta) count++;
    }
  } else {  
    AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
    if (!esdEvent) return -1;
    const AliMultiplicity *mult = esdEvent->GetMultiplicity();
    if (mult) {
      if (mult->GetNumberOfTracklets()>0) {	
	for (Int_t n=0; n<mult->GetNumberOfTracklets(); n++) {
	  Double_t eta = -TMath::Log( TMath::Tan( mult->GetTheta(n) / 2.0 ) );
	  if(eta>mineta && eta<maxeta) count++;
	}
      }
    }
  }

  return count;
}

//______________________________________________________________________
Double_t AliCFSingleTrackEfficiencyTask::GetCentrality()
{
  //
  // Get centrality
  //
  if(fInputEvent->GetRunNumber()<244824) return GetCentralityOldFramework();
  Double_t cent=-1;
  AliMultSelection *multSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  if(!multSelection){
    AliWarning("AliMultSelection could not be found in the aod event list of objects");
    return cent;
  }
  cent=multSelection->GetMultiplicityPercentile(fCentralityEstimator.Data());
  Int_t qual = multSelection->GetEvSelCode();
  if(qual == 199 ) cent=-1;
  return cent;

}

//______________________________________________________________________
Double_t AliCFSingleTrackEfficiencyTask::GetCentralityOldFramework()
{
  //
  // Get centrality, Run1 framework
  //
  Bool_t isAOD = fInputEvent->IsA()->InheritsFrom("AliAODEvent");
  Double_t cent = -1;

  if(isAOD) {
    AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
    if(!aodEvent) return cent;
    AliAODHeader* header = dynamic_cast<AliAODHeader*>(aodEvent->GetHeader());
    if(!header) AliFatal("Not a standard AOD");
    if(!header) return cent;
    AliCentrality *centrality = header->GetCentralityP();
    if(!centrality) return cent;
    //    cout<<" about to get cent perc AOD"<<endl;
    cent = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
  } else {
    AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
    if(!esdEvent) return cent;
    AliCentrality *centrality = esdEvent->GetCentrality();
    if(!centrality) return cent;
    cent = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
  }

  return cent;
}
