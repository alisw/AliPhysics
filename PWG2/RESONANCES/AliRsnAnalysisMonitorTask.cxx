//
// Implementation file for implementation of data analysis aft 900 GeV
//
// Author: A. Pulvirenti
//

#include "Riostream.h"
#include <iomanip>

#include "TH1.h"
#include "TTree.h"
#include "TParticle.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "AliLog.h"
#include "AliESDpid.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliITSPIDResponse.h"
#include "AliRsnMonitorTrack.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"

#include "AliRsnAnalysisMonitorTask.h"

ClassImp(AliRsnAnalysisMonitorTask)

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorTask::AliRsnAnalysisMonitorTask(const char *name) :
  AliAnalysisTaskSE(name),
  fOut(0x0),
  fTrack(0x0),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0),
  fEventCuts("eventCuts", AliRsnCut::kEvent),
  fTrackCuts("trackCuts", AliRsnCut::kDaughter)
  
{
//
// Constructor
//

  DefineOutput(1, TTree::Class());
}

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorTask::AliRsnAnalysisMonitorTask(const AliRsnAnalysisMonitorTask& copy) :
  AliAnalysisTaskSE(copy),
  fOut(0x0),
  fTrack(0x0),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0),
  fEventCuts(copy.fEventCuts),
  fTrackCuts(copy.fTrackCuts)
{
//
// Copy constructor
//
}

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorTask& AliRsnAnalysisMonitorTask::operator=(const AliRsnAnalysisMonitorTask& copy)
{
//
// Assignment operator
//

  fTOFcalibrateESD = copy.fTOFcalibrateESD;
  fTOFcorrectTExp = copy.fTOFcorrectTExp;
  fTOFuseT0 = copy.fTOFuseT0;
  fTOFtuneMC = copy.fTOFtuneMC;
  fTOFresolution = copy.fTOFresolution;

  return (*this);
}

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorTask::~AliRsnAnalysisMonitorTask()
{
//
// Destructor
//

  if (fOut)    delete fOut;
  if (fESDpid) delete fESDpid;
  if (fTrack)  delete fTrack;
}

//__________________________________________________________________________________________________
void AliRsnAnalysisMonitorTask::UserCreateOutputObjects()
{
//
// Create the output data container
//

  // setup TPC response
  fESDpid = new AliESDpid;
  fESDpid->GetTPCResponse().SetBetheBlochParameters(fTPCpar[0], fTPCpar[1], fTPCpar[2], fTPCpar[3], fTPCpar[4]);

  // setup TOF maker & calibration
  fTOFcalib = new AliTOFcalib;
  fTOFmaker = new AliTOFT0maker(fESDpid, fTOFcalib);
  fTOFmaker->SetTimeResolution(fTOFresolution);
  
  // create output branch object
  fTrack = new AliRsnMonitorTrack;
    
  // create output tree
  OpenFile(1);
  fOut = new TTree("rsnMonitor", "Informations on single tracks for cut checking");
  fOut->Branch("tracks", "AliRsnMonitorTrack", &fTrack);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisMonitorTask::UserExec(Option_t *)
{
//
// Main execution function.
// Fills the fHEvents data member with the following legenda:
// 0 -- event OK, prim vertex with tracks
// 1 -- event OK, prim vertex with SPD
// 2 -- event OK but vz large
// 3 -- event bad
//

  // retrieve ESD event and related stack (if available)
  AliESDEvent *esd   = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliStack    *stack = 0x0;
  
  // skip NULL events
  if (!esd) return;
  if (fMCEvent) stack = fMCEvent->Stack();
  
  // create interface objects to AliRsnEvent to check event cuts
  AliRsnEvent event;
  event.SetRef(esd);
  event.SetRefMC(fMCEvent);
  if (!fEventCuts.IsSelected(&event)) return;
  
  // check the event
  Int_t type = EventEval(esd);
    
  // if processable, then process it
  if      (type == 0) ProcessESD(esd, esd->GetPrimaryVertexTracks(), stack);
  else if (type == 1) ProcessESD(esd, esd->GetPrimaryVertexSPD()   , stack);
  else return;
  
  // update histogram container
  PostData(1, fOut);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisMonitorTask::Terminate(Option_t *)
{
//
// Terminate
//
}

//__________________________________________________________________________________________________
Int_t AliRsnAnalysisMonitorTask::EventEval(AliESDEvent *esd)
{
//
// Checks if the event is good for analysis.
// Returns:
// ---> 0 if a good primary vertex with tracks was found,
// ---> 1 if a good SPD primary vertex was found
// ---> 2 otherwise (event to be rejected)
// In any case, adds an entry to the TTree, to keep trace of all events.
//

  // get the best primary vertex:
  // first try that with tracks, then the SPD one
  const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
  const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
  if(vTrk->GetNContributors() > 0)
  {
    return 0;
  }
  else if (vSPD->GetNContributors() > 0)
  {
    return 1;
  }
  else
  {
    return 2;
  }
}

//__________________________________________________________________________________________________
void AliRsnAnalysisMonitorTask::ProcessESD
(AliESDEvent *esd, const AliESDVertex *v, AliStack *stack)
{
//
// Process the ESD container, to read all tracks and copy their useful values.
// All info are stored into an AliRsnMonitorTrack object and saved into the
// TClonesArray which is one of the branches of the output TTree.
//

  // create interfacr objects
  AliRsnEvent    event;
  AliRsnDaughter daughter;
  event.SetRef(esd);
  event.SetRefMC(fMCEvent);

  // ITS stuff #1 
  // create the response function and initialize it to MC or not
  // depending if the AliStack object is there or not
  Bool_t isMC = (stack != 0x0);
  AliITSPIDResponse itsrsp(isMC);

  // TOF stuff #1: init OCDB
  Int_t run = esd->GetRunNumber();
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(run);
  // TOF stuff #2: init calibration
  fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
  fTOFcalib->Init();
  // TOF stuff #3: calibrate
  if (fTOFcalibrateESD) fTOFcalib->CalibrateESD(esd);
  if (fTOFtuneMC) fTOFmaker->TuneForMC(esd);
  if (fTOFuseT0) 
  {
    fTOFmaker->ComputeT0TOF(esd);
    fTOFmaker->ApplyT0TOF(esd);
    fESDpid->MakePID(esd, kFALSE, 0.);
  }
  
  // loop on all tracks
  Int_t               i, k, nITS, ntracks = esd->GetNumberOfTracks();;
  Bool_t              isTPC, isITSSA, isTOF;
  Float_t             b[2], bCov[3];
  Double_t            time[10];
  
  for (i = 0; i < ntracks; i++)
  {
    AliESDtrack *track = esd->GetTrack(i);
    event.SetDaughter(daughter, i, AliRsnDaughter::kTrack);
    
    // reset the output object
    // 'usable' flag will need to be set to 'ok'
    fTrack->Reset();
    
    // check cuts
    fTrack->CutsPassed() = fTrackCuts.IsSelected(&daughter);
        
    // skip NULL pointers, kink daughters and tracks which
    // cannot be propagated to primary vertex
    if (!track) continue;
    if ((Int_t)track->GetKinkIndex(0) > 0) continue;
    if (!track->RelateToVertex(v, esd->GetMagneticField(), kVeryBig)) continue;
    
    // get MC info if possible
    if (stack) fTrack->AdoptMC(TMath::Abs(track->GetLabel()), stack);
    
    // copy general info
    fTrack->Status() = (UInt_t)track->GetStatus();
    fTrack->Length() = (Double_t)track->GetIntegratedLength();
    fTrack->Charge() = (Int_t)track->Charge();
    fTrack->PrecX()  = (Double_t)track->Px();
    fTrack->PrecY()  = (Double_t)track->Py();
    fTrack->PrecZ()  = (Double_t)track->Pz();
    
    // evaluate some flags from the status to decide what to do next in some points
    isTPC   = ((fTrack->Status() & AliESDtrack::kTPCin)  != 0);
    isITSSA = ((fTrack->Status() & AliESDtrack::kTPCin)  == 0 && (fTrack->Status() & AliESDtrack::kITSrefit) != 0 && (fTrack->Status() & AliESDtrack::kITSpureSA) == 0 && (fTrack->Status() & AliESDtrack::kITSpid) != 0);
    isTOF   = ((fTrack->Status() & AliESDtrack::kTOFout) != 0 && (fTrack->Status() & AliESDtrack::kTIME) != 0);
    
    // accept only tracks which are TPC+ITS or ITS standalone
    if (isITSSA)
    {
      fTrack->ITSsa() = kTRUE;
      fTrack->TOFok() = kFALSE;
    }
    else if (isTPC)
    {
      fTrack->ITSsa() = kFALSE;
      fTrack->TOFok() = isTOF;
    }
    else
      continue;

    // get DCA to primary vertex
    track->GetImpactParameters(b, bCov);
    fTrack->DCAr() = (Double_t)b[0];
    fTrack->DCAz() = (Double_t)b[1];
    
    // get ITS info
    for (k = 0; k < 6; k++)
    {
      fTrack->ITSmap(k) = track->HasPointOnITSLayer(k);
    }
    if (isITSSA)
    {
      fTrack->ITSchi2() = track->GetITSchi2();
      fTrack->ITSsignal() = track->GetITSsignal();
      nITS = fTrack->SSDcount() + fTrack->SDDcount();
      for (k = 0; k < AliPID::kSPECIES; k++)
      {
        fTrack->ITSnsigma(k) = itsrsp.GetNumberOfSigmas(fTrack->Prec(), fTrack->ITSsignal(), (AliPID::EParticleType)k, nITS, kTRUE);
      }
    }
    
    // get TPC info
    if (isTPC)
    {
      fTrack->TPCcount()  = (Int_t)track->GetTPCclusters(0);
      fTrack->TPCchi2()   = (Double_t)track->GetTPCchi2();
      fTrack->TPCsignal() = (Double_t)track->GetTPCsignal();
      fTrack->PtpcX()     = fTrack->PtpcY() = fTrack->PtpcZ() = 1E10;
      if (track->GetInnerParam())
      {
        fTrack->PtpcX() = track->GetInnerParam()->Px();
        fTrack->PtpcY() = track->GetInnerParam()->Py();
        fTrack->PtpcZ() = track->GetInnerParam()->Pz();
        for (k = 0; k < AliPID::kSPECIES; k++) 
        {
          fTrack->TPCnsigma(k) = fESDpid->NumberOfSigmasTPC(track, (AliPID::EParticleType)k);
        }
      }
    }
    
    // get TOF info
    if (isTOF)
    {
      track->GetIntegratedTimes(time);
      fTrack->TOFsignal() = (Double_t)track->GetTOFsignal();
      for (k = 0; k < AliPID::kSPECIES; k++)
      {
        fTrack->TOFref(k)   = time[k];
        fTrack->TOFsigma(k) = (Double_t)fTOFmaker->GetExpectedSigma(fTrack->Prec(), time[k], AliPID::ParticleMass(k));
      }
    }
    
    // add entry to TTree
    fOut->Fill();
  }
}
