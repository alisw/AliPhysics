#include "AliTOFAnalysisTaskCalibTree.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "TTree.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliPhysicsSelection.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliTOFT0v1.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"

ClassImp(AliTOFAnalysisTaskCalibTree)
  
//_______________________________________________________
  
AliTOFAnalysisTaskCalibTree::AliTOFAnalysisTaskCalibTree(const Char_t* name) :
AliAnalysisTaskSE(name),
  fInitFlag(kFALSE),                   
  fEventSelectionFlag(kFALSE),        
  fVertexSelectionFlag(kFALSE),       
  fVertexCut(1000.),              
  fDiscardPileupEventFlag(kFALSE),   
  fCalibrateTOFsignal(kFALSE),       
  fComputeT0TOF(kFALSE),             
  fUseT0TOF(kFALSE),                 
  fUseLHCClockPhase(kFALSE),         
  fPrimaryDCASelectionFlag(kFALSE),  
  fRunNumber(0),                   
  fESDEvent(0),             
  fEventCuts(new AliPhysicsSelection()),    
  fTrackCuts(new AliESDtrackCuts()),  
  fESDpid(new AliESDpid()),
  fStartTime(0),
  fEndTime(0),
  fElapsedTime(0),
  fIsCollisionCandidate(kFALSE),
  fHasVertex(kFALSE),
  fVertex(NULL),
  fGRPManager(new AliGRPManager()),         
  fGRPObject(NULL),     
  fSpecificStorageParOffline(), 
  fSpecificStorageRunParams(),  
  fTimeResolution(80.),
  fTOFcalib(new AliTOFcalib()),             
  fTOFT0maker(new AliTOFT0maker(fESDpid, fTOFcalib)),         
  fTOFT0v1(new AliTOFT0v1(fESDpid))             

{
  /* 
   * default constructor 
   */

}

//_______________________________________________________

AliTOFAnalysisTaskCalibTree::~AliTOFAnalysisTaskCalibTree()
{
  /*
   * default destructor
   */
  delete fEventCuts;
  delete fTrackCuts;
  delete fESDpid;
  delete fTOFcalib;
  delete fTOFT0maker;
  delete fTOFT0v1;

}

//_______________________________________________________

void
AliTOFAnalysisTaskCalibTree::UserCreateOutputObjects()
{
  /*
   * user create output objects
   */

  /* setup output tree */
  OutputTree()->Branch("run", &fRunNumber, "run/I");
  OutputTree()->Branch("timestamp", &ftimestamp, "timestamp/i");
  OutputTree()->Branch("timezero", &ftimezero, "timezero/F");
  OutputTree()->Branch("vertex", &fVertexZ, "vertex/F");
  OutputTree()->Branch("nhits", &fnhits, "nhits/I");
  OutputTree()->Branch("momentum", &fmomentum, "momentum[nhits]/F");
  OutputTree()->Branch("length", &flength, "length[nhits]/F");
  OutputTree()->Branch("index", &findex, "index[nhits]/I");
  OutputTree()->Branch("time", &ftime, "time[nhits]/F");
  OutputTree()->Branch("tot", &ftot, "tot[nhits]/F");
  OutputTree()->Branch("texp", &ftexp, "texp[nhits]/F");

}

//_______________________________________________________

void AliTOFAnalysisTaskCalibTree::UserExec(Option_t *option) {
  //
  // user exec
  //

  // unset fill AOD
  ((AliAODHandler*)(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()))->SetFillAOD(kFALSE);

  //init run
  if (!InitRun()) return;

  //init event 
  if (!InitEvent()) return;

  /*** ACCEPTED EVENT ***/
  
  // set vertex
  fVertexZ = fVertex->GetZ();

  // compute T0-TOF using all availeble tracks 
  fTOFT0v1->DefineT0("all", 0.0, 0.5);
  Float_t timeZeroTOF = -1000. * fTOFT0v1->GetResult(0);
  Float_t timeZeroTOF_sigma = 1000. * fTOFT0v1->GetResult(1);
  Int_t timeZeroTOF_tracks = fTOFT0v1->GetResult(3);

  // check T0-TOF sigma 
  if (timeZeroTOF_sigma >= 250.) ftimezero = 999999.; 
  else ftimezero = timeZeroTOF;

  // reset
  fnhits = 0;
  
  // loop over ESD tracks
  Int_t nTracks = fESDEvent->GetNumberOfTracks();
  AliESDtrack *track;
  Int_t index;
  Double_t momentum, length, time, tot, timei[AliPID::kSPECIES];
  for (Int_t itrk = 0; itrk < nTracks; itrk++) {
    // get track
    track = fESDEvent->GetTrack(itrk);
    if (!track) continue;
    // check accept track
    if (!fTrackCuts->AcceptTrack(track)) continue;
    // check primary DCA
    if (fPrimaryDCASelectionFlag && !HasPrimaryDCA(track)) continue;
    // check TOF measurement
    if (!HasTOFMeasurement(track)) continue;

    /*** ACCEPTED TRACK WITH TOF MEASUREMENT ***/

    // get track info 
    momentum = track->P();
    length = track->GetIntegratedLength();
    index = track->GetTOFCalChannel();
    time = track->GetTOFsignal();
    tot = track->GetTOFsignalToT();
    track->GetIntegratedTimes(timei);

    // add hit to array (if there is room)
    if (fnhits > MAXHITS) continue;
    fmomentum[fnhits] = momentum;
    flength[fnhits] = length;
    findex[fnhits] = index;
    ftime[fnhits] = time;
    ftot[fnhits] = tot;
    ftexp[fnhits] = timei[AliPID::kPion];
    fnhits++;

  } /* end of loop over ESD tracks */

  // check number of hits and set fill output tree
  if (fnhits > 0)
    ((AliAODHandler*)(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler()))->SetFillAOD(kTRUE);

}

//_______________________________________________________

Bool_t AliTOFAnalysisTaskCalibTree::InitRun() {

  //
  // init run
  //
  
  // get ESD event 
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESDEvent) {
    AliError("cannot get ESD event");
    return kFALSE;
  }

  // get run number 
  Int_t runNb = fESDEvent->GetRunNumber();

  // check run already initialized 
  if (fInitFlag && fRunNumber == runNb) return kTRUE;
  fInitFlag = kFALSE;

  // init cdb 
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  if (!fSpecificStorageParOffline.IsNull())
    cdb->SetSpecificStorage("TOF/Calib/ParOffline", fSpecificStorageParOffline.Data());
  if (!fSpecificStorageRunParams.IsNull())
    cdb->SetSpecificStorage("TOF/Calib/RunParams", fSpecificStorageRunParams.Data());
  cdb->SetRun(runNb);
  
  // init TOF calib
  if (!fTOFcalib->Init()) {
    AliError("cannot init TOF calib");
    return kFALSE;
  }

  // get GRP data
  if (!fGRPManager->ReadGRPEntry()) {
    AliError("error while reading \"GRPEntry\" from OCDB");
    return kFALSE;
  }
  fGRPObject = fGRPManager->GetGRPData();
  if (!fGRPObject) {
    AliError("cannot get \"GRPData\" from GRP manager");
    return kFALSE;
  }
  fStartTime = fGRPObject->GetTimeStart();
  fEndTime = fGRPObject->GetTimeEnd();
  AliInfo(Form("got \"GRPData\": startTime=%d, endTime=%d", fStartTime, fEndTime));
  
  AliInfo(Form("initialized for run %d", runNb));
  fInitFlag = kTRUE;
  fRunNumber = runNb;
  return kTRUE;
}

//_______________________________________________________

Bool_t AliTOFAnalysisTaskCalibTree::InitEvent() {

  //
  // init event
  //

  // get ESD event 
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESDEvent) return kFALSE;

  // set event time and elapsed time
  ftimestamp = fESDEvent->GetTimeStamp();
  fElapsedTime = fESDEvent->GetTimeStamp() - fStartTime;

  // event selection
  fIsCollisionCandidate = fEventCuts->IsCollisionCandidate(fESDEvent);
  if (fEventSelectionFlag && !fIsCollisionCandidate) return kFALSE;

  // vertex selection
  fVertex = fESDEvent->GetPrimaryVertexTracks();
  if (fVertex->GetNContributors() < 1) {
    fVertex = fESDEvent->GetPrimaryVertexSPD();
    if (fVertex->GetNContributors() < 1) fHasVertex = kFALSE;
    else fHasVertex = kTRUE;
  }
  else fHasVertex = kTRUE;
  if (fVertexSelectionFlag && (!fHasVertex || TMath::Abs(fVertex->GetZ()) > fVertexCut)) return kFALSE;

  // discard pile-up event is requested
  if (fDiscardPileupEventFlag) {
    if (fESDEvent->IsPileupFromSPD()) {
      printf("PILE-UP event, will be discarded\n");
      return kFALSE;
    }
  }

  // init TOF response 
  fESDpid->GetTOFResponse().SetTimeResolution(fTimeResolution);

  // init TOF-T0 maker 
  fTOFT0maker->SetTimeResolution(fTimeResolution);
  
  // calibrate ESD if requested
  if (fCalibrateTOFsignal)
    fTOFcalib->CalibrateESD(fESDEvent);

  // compute T0-TOF and apply it if requested
  if (fComputeT0TOF) 
    fTOFT0maker->ComputeT0TOF(fESDEvent);
  if (fUseT0TOF) {
    fTOFT0maker->ApplyT0TOF(fESDEvent);
    //fESDpid->MakePID(fESDEvent, kFALSE, 0.);
  }

  // init T0-TOF
  fTOFT0v1->Init(fESDEvent);

  return kTRUE;

}

//_______________________________________________________

Bool_t AliTOFAnalysisTaskCalibTree::HasTOFMeasurement(AliESDtrack *track) {

  //
  // has TOF measurement
  //

  //check TOF status flags 
  if (!(track->GetStatus() & AliESDtrack::kTOFout) ||
      !(track->GetStatus() & AliESDtrack::kTIME)) return kFALSE;

  // check integrated length
  if (track->GetIntegratedLength() < 350.) return kFALSE;

  // TOF measurement ok
  return kTRUE;
}

//_______________________________________________________

Bool_t AliTOFAnalysisTaskCalibTree::HasPrimaryDCA(AliESDtrack *track) {

  //
  // has primary DCA
  //

  // cut on transverse impact parameter
  Float_t d0z0[2],covd0z0[3];
  track->GetImpactParameters(d0z0, covd0z0);
  Float_t sigma= 0.0050 + 0.0060 / TMath::Power(track->Pt(), 0.9);
  Float_t d0max = 7. * sigma;
  //
  Float_t sigmaZ = 0.0146 + 0.0070 / TMath::Power(track->Pt(), 1.114758);
  if (track->Pt() > 1.) sigmaZ = 0.0216;
  Float_t d0maxZ = 5. * sigmaZ;
  //
  if(TMath::Abs(d0z0[0]) > d0max || TMath::Abs(d0z0[1]) > d0maxZ) 
    return kFALSE;
  
  /* primary DCA ok */
  return kTRUE;
}
