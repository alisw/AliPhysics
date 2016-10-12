#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "TParticle.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"


#include "AliAnalysisTaskSimpleTreeMaker.h"

ClassImp(AliAnalysisTaskSimpleTreeMaker)

Int_t numEvents = 0;

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker():
  AliAnalysisTaskSE(),
  fTree(0),
  fStream(0),
  fQAhist(0),
  esdEvent(0), 
  mcEvent(0), 
  fESDtrackCuts(0),
  fPIDResponse(0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100), 
  fPtMin(0.4),
  fPtMax(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fESigITSMin(-3.),
  fESigITSMax(3.),
  fESigTPCMin(-5.),
  fESigTPCMax(5.),
  fESigTOFMin(-3.),
  fESigTOFMax(3.),
  fPSigTPCMin(-99.),
  fPSigTPCMax(-3.),
  fPIDcutITS(kFALSE),
  fPIDcutTOF(kFALSE),
  fPionPIDcutTPC(kFALSE),
  isIonColl(kFALSE),
  fIsMC(kTRUE)
{

}

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker(const char *name) :
  AliAnalysisTaskSE(name),
  fTree(0),
  fStream(0),
  fQAhist(0),
  esdEvent(0), 
  mcEvent(0), 
  fESDtrackCuts(0),
  fPIDResponse(0),
  fCentralityPercentileMin(0),
  fCentralityPercentileMax(100), 
  fPtMin(0.4),
  fPtMax(10),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fESigITSMin(-3.),
  fESigITSMax(3.),
  fESigTPCMin(-5.),
  fESigTPCMax(5.),
  fESigTOFMin(-3.),
  fESigTOFMax(3.),
  fPSigTPCMin(-99.),
  fPSigTPCMax(-3.),
  fPIDcutITS(kFALSE),
  fPIDcutTOF(kFALSE),
  fPionPIDcutTPC(kFALSE),
  isIonColl(kFALSE),
  fIsMC(kTRUE)
{
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
  fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
  Printf("sup");
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class()); //will be connected to fTree
  DefineOutput(2, TH1F::Class());
}

//________________________________________________________________________

//~ AliAnalysisTaskSimpleTreeMaker::~AliAnalysisTaskSimpleTreeMaker() {

  //~ // Destructor

  //~ // ... not implemented

//~ }


//________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::UserCreateOutputObjects() {

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  inputHandler->SetNeedField();
     
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse){
	  return;} 

  AliInfo("Init");

  fStream = new TTreeStream("tracks");
  fTree   = (TTree*)fStream->GetTree();
   
  fQAhist = new TH1F("h1", "Event and track QA", 6, 0, 1);
  PostData(1, fTree);
  PostData(2, fQAhist);
  
  AliInfo("Finished setting up the output");
}

//________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  esdEvent = (AliESDEvent*)InputEvent();
  if(!esdEvent) {
    AliError("No AliESDEvent");
    return;
  } 
  fQAhist->Fill("Events_ESDcheck",1);
  numEvents += 1;
  
  // Process also MC truth  
  mcEvent = MCEvent();
  if (!mcEvent) {
    AliError("Could not retrieve MC event");
    return;
  }
  fQAhist->Fill("Events_MCcheck",1);
  

  // check event cuts
  if( IsEventAccepted(esdEvent) == 0){ 
    return;
  }
  fQAhist->Fill("Events_accepted",1);
  
  // PID Response task active?
  fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();

  if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  
  AliESDVertex* vertex = (AliESDVertex*)esdEvent->GetPrimaryVertex();

  Double_t primaryVertex[3];
  primaryVertex[0] = vertex->GetX();
  primaryVertex[1] = vertex->GetY();
  primaryVertex[2] = vertex->GetZ();

  Double_t impactParameter = -999;
  if(fIsMC){
    AliGenHijingEventHeader* hHijing = 0;
    AliGenEventHeader*  mcGenH  = mcEvent->GenEventHeader();

    if (mcGenH->InheritsFrom(AliGenHijingEventHeader::Class())){
    //Option 1: Just Hijing
      hHijing = (AliGenHijingEventHeader*)mcGenH;
    } else if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
      //Option 2: cocktail involving Hijing
      TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
      hHijing = dynamic_cast<AliGenHijingEventHeader*>(headers->FindObject("hijing_0"));
    }   
    if(hHijing){
      impactParameter = hHijing->ImpactParameter();
    }
  }

  Int_t eventTracks = esdEvent->GetNumberOfTracks();
  Int_t runNumber = esdEvent->GetRunNumber();
  Int_t numTracks = 0;
  
  //Loop over tracks for event
  for (Int_t iTrack = 0; iTrack < eventTracks; iTrack++){ 
    AliESDtrack* track = esdEvent->GetTrack(iTrack);
	  if (!track) {
		  AliError(Form("Could not receive track %d", iTrack));
	    continue;
	  }

    fQAhist->Fill("Tracks_all",1);
    //Apply global track filter
    if(!fESDtrackCuts->AcceptTrack(track)){ continue;}
    //Apply momentum and eta cuts
    Double_t pt   = track->Pt();
    if( pt < fPtMin || pt > fPtMax ) { continue; }
    Double_t eta  = track->Eta();
    if( eta < fEtaMin || eta > fEtaMax ) { continue; } 
    fQAhist->Fill("Tracks_KineCuts", 1);

    //Get electron nSigma in TPC for cut (inclusive cut)
    Double_t EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron);
    if( EnSigmaTPC > fESigTPCMax || EnSigmaTPC < fESigTPCMin) { continue; }
    
    fQAhist->Fill("Tracks_PIDcuts",1); 
    numTracks += 1;
    
    
    //Get rest of electron nSigma values and apply cuts if requested (inclusive cuts)
    Double_t EnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron);
    if(fPIDcutITS){
      if(EnSigmaITS < fESigITSMin || EnSigmaITS > fESigITSMax){ continue; }
    }
    Double_t EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron);
    if(fPIDcutTOF){
      if(EnSigmaTOF < fESigTOFMin || EnSigmaTOF > fESigTOFMax){ continue; }
    }

    //Get pion nSigma for TPC and apply cut if requested (exclusive cut)
    Double_t PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion);
    if(fPionPIDcutTPC){
      if(fPionPIDcutTPC > fPSigTPCMin && fPionPIDcutTPC < fPSigTPCMax){ continue; }
    }
    
    //Get rest of nSigma values for pion and kaon
    Double_t PnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kPion);
    Double_t PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kPion);
    
    Double_t KnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kKaon);
    Double_t KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon);
    Double_t KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kKaon);
    
    Double_t phi  = track->Phi();

    //Get ITS and TPC signals
    Double_t ITSsignal = track->GetITSsignal();
    Double_t TPCsignal = track->GetTPCsignal();
    Double_t TOFsignal = track->GetTOFsignal();
    
        
    //DCA values
    Float_t DCAxy = 0.;
    Float_t DCAz = 0.;
    track->GetImpactParameters( &DCAxy, &DCAz);
    
    //ITS clusters and shared clusters
    Int_t nITS = track->GetNumberOfITSClusters();
    Double_t nITS_shared = 0.;
    for(Int_t d = 0; d < 6; d++){
            nITS_shared += (Double_t) track->HasSharedPointOnITSLayer(d);
      }
    nITS_shared /= nITS;
   
    //Get chi2 values 
    Double_t chi2ITS = track->GetITSchi2();
    Double_t chi2TPC = track->GetTPCchi2();

    Int_t fCutMaxChi2TPCConstrainedVsGlobalVertexType = fESDtrackCuts->kVertexTracks | fESDtrackCuts->kVertexSPD;
    const AliESDVertex* vertex = 0;
    if (fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTracks){
      vertex = track->GetESDEvent()->GetPrimaryVertexTracks();}

    if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexSPD){
      vertex = track->GetESDEvent()->GetPrimaryVertexSPD();}

    if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTPC){
      vertex = track->GetESDEvent()->GetPrimaryVertexTPC();}

    Double_t goldenChi2 = track->GetChi2TPCConstrainedVsGlobal(vertex);

    
    Int_t charge = -10;
    charge = track->Charge();
    Int_t label = -999;
    label = track->GetLabel();

    //Declare MC variables
    Double_t mcEta   = -99;
    Double_t mcPhi   = -99;
    Double_t mcPt    = -99;
    Int_t iPdg       = 0;
    Int_t iPdgMother = 0;
    if(fIsMC){
      AliMCParticle* mcTrack = (AliMCParticle*) mcEvent->GetTrack(TMath::Abs(label));
      TParticle* mcParticle = (TParticle*)mcTrack->Particle();
      mcEta = mcTrack->Eta();
      mcPhi = mcTrack->Phi();
      mcPt = mcTrack->Pt();

      iPdg = mcTrack->PdgCode();

      Int_t gMotherIndex = mcTrack->GetMother();
      AliMCParticle* motherTrack = (AliMCParticle*)(mcEvent->GetTrack(gMotherIndex));
      iPdgMother = motherTrack->PdgCode();
    }

    //Stream values into tree
    (*fStream)    << "tracks" <<
    "pt="         << pt << 
    "eta="        << eta << 
    "phi="        << phi << 
    "EsigITS="    << EnSigmaITS <<
    "EsigTPC="    << EnSigmaTPC <<
    "EsigTOF="    << EnSigmaTOF <<
    "PsigITS="    << PnSigmaITS <<
    "PsigTPC="    << PnSigmaTPC <<
    "PsigTOF="    << PnSigmaTOF <<
    "KsigITS="    << KnSigmaITS <<
    "KsigTPC="    << KnSigmaTPC <<
    "KsigTOF="    << KnSigmaTOF <<
    "charge="     << charge <<
    "ITSsignal="  << ITSsignal <<
    "TPCsignal="  << TPCsignal << 
    "TOFsignal="  << TOFsignal <<
    "vertexX="    << primaryVertex[0] <<
    "vertexY="    << primaryVertex[1] <<
    "vertexZ="    << primaryVertex[2] <<
    "nITS="       << nITS <<
    "nITSshared=" << nITS_shared << 
    "DCAxy="      << DCAxy <<
    "DCAz="       << DCAz <<
    "chi2ITS="    << chi2ITS <<
    "chi2TPC="    << chi2TPC <<
    "goldenChi2=" << goldenChi2 <<
    "mcEta="      << mcEta <<
    "mcPhi="      << mcPhi <<
    "mcPt="       << mcPt <<
    "pdg="        << iPdg <<
    "pdgMother="  << iPdgMother <<
    "runNumber="  << runNumber << 
    "numEvents="  << numEvents <<
    "\n";
  }

  Printf("Accepted- Event: %i, Tracks: %i", numEvents, numTracks);
}

//~ //________________________________________________________________________

void  AliAnalysisTaskSimpleTreeMaker::FinishTaskOutput(){
  // Finish task output

  // not implemented ...

}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::Terminate(Option_t *) {
  // Draw result to the screen

  // Called once at the end of the query

  // not implemented ...

}
//~ 


//________________________________________________________________________

Double_t AliAnalysisTaskSimpleTreeMaker::IsEventAccepted(AliESDEvent *event){

  
    if (TMath::Abs(event->GetVertex()->GetZ()) < 10  &&  event->GetPrimaryVertexSPD() ){
      if (event->GetPrimaryVertexSPD()->GetNContributors() >0) return 1.;
      else return 0.;}
    return 0.;
}



