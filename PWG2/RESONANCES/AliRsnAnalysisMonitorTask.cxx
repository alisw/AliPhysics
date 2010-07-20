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

#include "AliRsnAnalysisMonitorTask.h"

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorTask::AliRsnAnalysisMonitorTask(const char *name) :
  AliAnalysisTaskSE(name),
  fEventType(2),
  fNTracks(0),
  fOut(0x0),
  fTracks(0x0),
  fMaxITSband(1E6),
  fTPCpLimit(0.35),
  fLargeTPCband(-1E6),
  fSmallTPCband( 1E6),
  fESDtrackCutsTPC(),
  fESDtrackCutsITS(),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0)
  
{
//
// Constructor
//

  DefineOutput(1, TTree::Class());
}

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorTask::AliRsnAnalysisMonitorTask(const AliRsnAnalysisMonitorTask& copy) :
  AliAnalysisTaskSE(copy),
  fEventType(2),
  fNTracks(0),
  fOut(0x0),
  fTracks(0x0),
  fMaxITSband(copy.fMaxITSband),
  fTPCpLimit(copy.fTPCpLimit),
  fLargeTPCband(copy.fLargeTPCband),
  fSmallTPCband(copy.fSmallTPCband),
  fESDtrackCutsTPC(copy.fESDtrackCutsTPC),
  fESDtrackCutsITS(copy.fESDtrackCutsITS),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0)
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

  fMaxITSband = copy.fMaxITSband;
  
  fTPCpLimit    = copy.fTPCpLimit;
  fLargeTPCband = copy.fSmallTPCband;
  fSmallTPCband = copy.fLargeTPCband;
  
  fESDtrackCutsTPC = copy.fESDtrackCutsTPC;
  fESDtrackCutsITS = copy.fESDtrackCutsITS;
  
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
  
  // initialize all the branch arrays
  fTracks = new TClonesArray("AliRsnMonitorTrack", 0);
  
  // create output tree
  OpenFile(1);
  fOut = new TTree("rsnMonitor", "Informations on single tracks for cut checking");
  
  // link branches
  fOut->Branch("ntracks", &fNTracks     , "ntracks/I"   );
  fOut->Branch("evtype" , &fEventType   , "evtype/I"    );
  fOut->Branch("vertex" , &fVertex      , "vertex[3]/F");
  fOut->Branch("tracks" , "TClonesArray", &fTracks      );
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

  static Int_t evNum = 0;
  evNum++;

  // retrieve ESD event and related stack (if available)
  AliESDEvent *esd   = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliStack    *stack = (fMCEvent ? fMCEvent->Stack() : 0x0);
  
  // check the event
  EventEval(esd);
    
  // if processable, then process it
  if      (fEventType == 0) ProcessESD(esd, esd->GetPrimaryVertexTracks(), stack);
  else if (fEventType == 1) ProcessESD(esd, esd->GetPrimaryVertexSPD()   , stack);
  else
  {
    fTracks->Delete();
    fTracks->Clear();
    fNTracks = 0;
  }
  
  // add a new entry in the TTree
  fOut->Fill();
  
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
void AliRsnAnalysisMonitorTask::EventEval(AliESDEvent *esd)
{
//
// Checks if the event is good for analysis.
// Sets the 'fEventType' flag to:
// ---> 0 if a good primary vertex with tracks was found,
// ---> 1 if a good SPD primary vertex was found
// ---> 2 otherwise (event to be rejected)
// In any case, adds an entry to the TTree, to keep trace of all events.
//
  
  // get number of tracks
  fNTracks = esd->GetNumberOfTracks();

  // get the best primary vertex:
  // first try that with tracks, then the SPD one
  const AliESDVertex *vTrk  = esd->GetPrimaryVertexTracks();
  const AliESDVertex *vSPD  = esd->GetPrimaryVertexSPD();
  if(vTrk->GetNContributors() > 0)
  {
    fVertex[0] = vTrk->GetXv();
    fVertex[1] = vTrk->GetYv();
    fVertex[2] = vTrk->GetZv();
    fEventType = 0;
  }
  else if (vSPD->GetNContributors() > 0)
  {
    fVertex[0] = vSPD->GetXv();
    fVertex[1] = vSPD->GetYv();
    fVertex[2] = vSPD->GetZv();
    fEventType = 1;
  }
  else
  {
    fNTracks = 0;
    fEventType = 2;
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

  // clear array
  fTracks->Delete();
  fTracks->Clear();
  
  // reject empty events
  if (!fNTracks) return;

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
  // TOF stuff #4: define fixed functions for compatibility range
  Double_t a1 = 0.01, a2 = -0.03;
  Double_t b1 = 0.25, b2 =  0.25;
  Double_t c1 = 0.05, c2 = -0.03;
  Double_t ymax, ymin;
  
  // loop on all tracks
  Int_t               i, k, size, nITS;
  Double_t            tpcMaxNSigma, itsdedx[4], tofRef, tofRel;
  Bool_t              okTOF, okTrack, isTPC, isITSSA, matchedTOF;
  UChar_t             itsCluMap;
  Float_t             b[2], bCov[3];
  AliRsnMonitorTrack  mon;
  
  for (i = 0; i < fNTracks; i++)
  {
    AliESDtrack *track = esd->GetTrack(i);
    
    // skip NULL pointers, kink daughters and tracks which
    // cannot be propagated to primary vertex
    if (!track) continue;
    if ((Int_t)track->GetKinkIndex(0) > 0) continue;
    if (!track->RelateToVertex(v, esd->GetMagneticField(), kVeryBig)) continue;
    
    // reset the output object
    // 'usable' flag will need to be set to 'ok'
    mon.Reset();
    
    // copy general info
    mon.Status() = (UInt_t)track->GetStatus();
    mon.Length() = (Double_t)track->GetIntegratedLength();
    mon.Charge() = (Int_t)track->Charge();
    mon.PrecX()  = (Double_t)track->Px();
    mon.PrecY()  = (Double_t)track->Py();
    mon.PrecZ()  = (Double_t)track->Pz();
    
    // evaluate some flags from the status to decide what to do next in some points
    isTPC      = ((mon.Status() & AliESDtrack::kTPCin)  != 0);
    isITSSA    = ((mon.Status() & AliESDtrack::kTPCin)  == 0 && (mon.Status() & AliESDtrack::kITSrefit) != 0 && (mon.Status() & AliESDtrack::kITSpureSA) == 0 && (mon.Status() & AliESDtrack::kITSpid) != 0);
    matchedTOF = ((mon.Status() & AliESDtrack::kTOFout) != 0 && (mon.Status() & AliESDtrack::kTIME) != 0);
    
    // accept only tracks which are TPC+ITS or ITS standalone
    if (!isTPC && !isITSSA) continue;
    
    // set the track type in the output object
    mon.ITSsa() = isITSSA;

    // get DCA to primary vertex
    track->GetImpactParameters(b, bCov);
    mon.DCAr() = (Double_t)b[0];
    mon.DCAz() = (Double_t)b[1];
    
    // get ITS info
    track->GetITSdEdxSamples(itsdedx);
    mon.ITSchi2() = track->GetITSchi2();
    mon.ITSsignal() = track->GetITSsignal();
    for (k = 0; k < 6; k++)
    {
      mon.ITSmap(k) = track->HasPointOnITSLayer(k);
      if (k < 4) mon.ITSdedx(k) = itsdedx[k];
    }
    
    // get TPC info
    if (isTPC)
    {
      mon.TPCcount()  = (Int_t)track->GetTPCclusters(0);
      mon.TPCdedx()   = (Double_t)track->GetTPCsignal();
      mon.TPCchi2()   = (Double_t)track->GetTPCchi2();
      mon.TPCnsigma() = fESDpid->NumberOfSigmasTPC(track, AliPID::kKaon);
      mon.PtpcX()     = mon.PtpcY() = mon.PtpcZ() = 1E10;
      if (track->GetInnerParam())
      {
        mon.PtpcX() = track->GetInnerParam()->Px();
        mon.PtpcY() = track->GetInnerParam()->Py();
        mon.PtpcZ() = track->GetInnerParam()->Pz();
        for (k = 0; k < AliPID::kSPECIES; k++) mon.TPCref(k) = fESDpid->GetTPCResponse().GetExpectedSignal(mon.Ptpc(), (AliPID::EParticleType)k);
      }
    }
    
    // get TOF info
    Double_t time[10];
    track->GetIntegratedTimes(time);
    mon.TOFsignal() = (Double_t)track->GetTOFsignal();
    for (k = 0; k < AliPID::kSPECIES; k++)
    {
      mon.TOFref(k)   = time[k];
      mon.TOFsigma(k) = (Double_t)fTOFmaker->GetExpectedSigma(mon.Prec(), time[k], AliPID::ParticleMass(k));
    }
    
    // if we are here, the track is usable
    mon.SetUsable();
    
    // now check the track against its cuts
    // and update the flag related to it
    // first, assume that cuts were passed
    // and if they aren't, just update the flag accordingly
    mon.CutsPassed() = kTRUE;
    if (TMath::Abs(mon.DCAz()) > 3.0) mon.CutsPassed() = kFALSE;
    
    if (isTPC)
    {
      // check standard ESD cuts
      if (!fESDtrackCutsTPC.IsSelected(track)) mon.CutsPassed() = kFALSE;
      
      // check TPC dE/dx
      if (mon.Ptpc() > fTPCpLimit) tpcMaxNSigma = fSmallTPCband; else tpcMaxNSigma = fLargeTPCband;
      if (TMath::Abs(mon.TPCnsigma()) > tpcMaxNSigma) mon.CutsPassed() = kFALSE;
      
      // check TOF (only if momentum is large than function asymptote and flags are OK)
      okTOF = kTRUE;
      if (matchedTOF && mon.Prec() > TMath::Max(b1, b2))
      {
        tofRef = mon.TOFref(AliPID::kKaon);
        if (tofRef > 0.0)
        {
          tofRel   = (mon.TOFsignal() - tofRef) / tofRef;
          ymax     = a1 / (mon.Prec() - b1) + c1;
          ymin     = a2 / (mon.Prec() - b2) + c2;
          okTOF    = (tofRel >= ymin && tofRel <= ymax);
        }
      }
      if (!okTOF) mon.CutsPassed() = kFALSE;
    }
    else
    {
      // check standard ESD cuts
      if (!fESDtrackCutsITS.IsSelected(track)) mon.CutsPassed() = kFALSE;
      
      // check dE/dx
      itsCluMap = track->GetITSClusterMap();
      nITS      = 0;
      for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITS;
      if (nITS < 3) okTrack = kFALSE; // track not good for PID
      mon.ITSnsigma() = itsrsp.GetNumberOfSigmas(mon.Prec(), mon.ITSsignal(), AliPID::kKaon, nITS, kTRUE);
      if (TMath::Abs(mon.ITSnsigma()) > fMaxITSband) mon.CutsPassed() = kFALSE;
    }
    
    // collect only tracks which are declared usable
    if (mon.IsUsable())
    {
      size = (Int_t)fTracks->GetEntriesFast();
      new ((*fTracks)[size]) AliRsnMonitorTrack(mon);
    }
  }
}
