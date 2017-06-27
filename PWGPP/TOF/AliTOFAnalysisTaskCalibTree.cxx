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
  fSpecificStorageFineSlewing(),
  fTimeResolution(80.),
  fTOFcalib(new AliTOFcalib()),             
  fTOFT0maker(new AliTOFT0maker(fESDpid, fTOFcalib)),         
  fTOFT0v1(new AliTOFT0v1(fESDpid)),
  fMaxHits(MAXHITS),
  ftimestamp(0),
  fVertexZ(99999),
  ftimezero(9999999),
  fnhits(-1),
  fmomentum(0),
  flength(0),
  findex(0),
  ftime(0),
  ftot(0),
  ftexp(0),
  fDeltax(0),
  fDeltaz(0),
  fDeltat(0),
  fDeltaraw(0),
  fSaveCoordinates(kFALSE),
  fOutputTree(0x0)             

{
  /* 
   * default constructor 
   */

  fmomentum = new Float_t[fMaxHits];
  flength = new Float_t[fMaxHits];
  findex = new Int_t[fMaxHits];
  ftime = new Float_t[fMaxHits];
  ftot = new Float_t[fMaxHits];
  ftexp = new Float_t[fMaxHits];
  fDeltax = new Float_t[fMaxHits];
  fDeltaz = new Float_t[fMaxHits];
  fDeltat = new Float_t[fMaxHits];
  fDeltaraw = new Float_t [fMaxHits]; 
  for (Int_t i = 0; i < fMaxHits; i++){
      fmomentum[i] = 999999;
      flength[i] = 999999;
      findex[i] = -1;
      ftime[i] = 999999;
      ftot[i] = 999999;
      ftexp[i] = 999999;
      fDeltax[i] = 999999;
      fDeltaz[i] = 999999;
      fDeltat[i] = 999999;
      fDeltaraw[i] = 999999;

  }

  DefineOutput(1, TTree::Class());  

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
  delete fGRPManager;
  delete fGRPObject;
  delete fESDEvent;
  delete fVertex;

  delete fmomentum;
  delete flength;
  delete findex;
  delete ftime;
  delete ftot;
  delete ftexp;
  delete fDeltax;
  delete fDeltaz;
  delete fDeltat;
  delete fDeltaraw;
  if (fOutputTree) {
    delete fOutputTree;
    fOutputTree = 0x0;
  }

}

//_______________________________________________________

void
AliTOFAnalysisTaskCalibTree::UserCreateOutputObjects()
{
  /*
   * user create output objects
   */

  fOutputTree = new TTree("aodTree", "Tree with TOF calib output"); 

  /* setup output tree */
  
if (fSaveCoordinates) { //save reduced tree for coordinate study
   
  fOutputTree->Branch("nhits", &fnhits, "nhits/I");
  fOutputTree->Branch("index", findex, "index[nhits]/I");
  fOutputTree->Branch("tot", ftot, "tot[nhits]/F");
  fOutputTree->Branch("deltax", fDeltax, "deltax[nhits]/F");
  fOutputTree->Branch("deltaz", fDeltaz, "deltaz[nhits]/F");
  fOutputTree->Branch("deltat", fDeltat, "deltat[nhits]/F");
  fOutputTree->Branch("deltaraw", fDeltaraw, "deltaraw[nhits]/F");
} 
   else { //default tree
  fOutputTree->Branch("run", &fRunNumber, "run/I");
  fOutputTree->Branch("timestamp", &ftimestamp, "timestamp/i");
  fOutputTree->Branch("timezero", &ftimezero, "timezero/F");
  fOutputTree->Branch("vertex", &fVertexZ, "vertex/F");
  fOutputTree->Branch("nhits", &fnhits, "nhits/I");
  fOutputTree->Branch("momentum", fmomentum, "momentum[nhits]/F");
  fOutputTree->Branch("length", flength, "length[nhits]/F");
  fOutputTree->Branch("index", findex, "index[nhits]/I");
  fOutputTree->Branch("time", ftime, "time[nhits]/F");
  fOutputTree->Branch("tot", ftot, "tot[nhits]/F");
  fOutputTree->Branch("texp", ftexp, "texp[nhits]/F");

   }

  PostData(1, fOutputTree);

}

//_______________________________________________________

void AliTOFAnalysisTaskCalibTree::UserExec(Option_t *) {
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
  //Int_t timeZeroTOF_tracks = fTOFT0v1->GetResult(3);  // not used for the time being

  // check T0-TOF sigma 
  if (timeZeroTOF_sigma >= 250.) { ftimezero = 999999.; 
   timeZeroTOF = 0.;}
  else ftimezero = timeZeroTOF;
  
  // reset
  fnhits = 0;
  
  // loop over ESD tracks
  Int_t nTracks = fESDEvent->GetNumberOfTracks();
  AliESDtrack *track;
  Int_t index;
  Double_t momentum, length, time, tot, timei[AliPID::kSPECIES], deltax, deltaz, deltat, deltaraw;
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
    deltax = track->GetTOFsignalDx();
    deltaz = track->GetTOFsignalDz();
    deltat = (time - timei[AliPID::kPion] - timeZeroTOF);
    deltaraw = (track->GetTOFsignalRaw() - timei[AliPID::kPion] - timeZeroTOF);
    // add hit to array (if there is room)
    if (fnhits > fMaxHits) continue;
    fmomentum[fnhits] = momentum;
    flength[fnhits] = length;
    ftexp[fnhits] = timei[AliPID::kPion];
    findex[fnhits] = index;
    ftime[fnhits] = time;
    ftot[fnhits] = tot;
    fDeltax[fnhits] = deltax;
    fDeltaz[fnhits] = deltaz;
    fDeltat[fnhits] = deltat;
    fDeltaraw[fnhits] = deltaraw;
    fnhits++;

  } /* end of loop over ESD tracks */

  // check number of hits and set fill output tree
  if (fnhits > 0)
    fOutputTree->Fill();

  PostData(1, fOutputTree);

}

//_______________________________________________________

Bool_t AliTOFAnalysisTaskCalibTree::InitRun() {

  //¯
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
  if (!fSpecificStorageFineSlewing.IsNull())
    cdb->SetSpecificStorage("TOF/Calib/FineSlewing", fSpecificStorageFineSlewing.Data());
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
  if (fEventSelectionFlag){
    fIsCollisionCandidate = fEventCuts->IsCollisionCandidate(fESDEvent);
    if (!fIsCollisionCandidate) return kFALSE;
  }

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
