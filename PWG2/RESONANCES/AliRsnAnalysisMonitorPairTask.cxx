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

#include "AliRsnAnalysisMonitorPairTask.h"

ClassImp(AliRsnAnalysisMonitorPairTask)

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorPairTask::AliRsnAnalysisMonitorPairTask(const char *name) :
  AliAnalysisTaskSE(name),
  fOut(0x0),
  fInvMass(0.0),
  fRangeMin(0.0),
  fRangeMax(100.0),
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
AliRsnAnalysisMonitorPairTask::AliRsnAnalysisMonitorPairTask(const AliRsnAnalysisMonitorPairTask& copy) :
  AliAnalysisTaskSE(copy),
  fOut(0x0),
  fInvMass(0.0),
  fRangeMin(copy.fRangeMin),
  fRangeMax(copy.fRangeMax),
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
AliRsnAnalysisMonitorPairTask& AliRsnAnalysisMonitorPairTask::operator=(const AliRsnAnalysisMonitorPairTask& copy)
{
//
// Assignment operator
//

  fMass[0] = copy.fMass[0];
  fMass[1] = copy.fMass[1];
  
  fRangeMin = copy.fRangeMin;
  fRangeMax = copy.fRangeMax;

  fTOFcalibrateESD = copy.fTOFcalibrateESD;
  fTOFcorrectTExp = copy.fTOFcorrectTExp;
  fTOFuseT0 = copy.fTOFuseT0;
  fTOFtuneMC = copy.fTOFtuneMC;
  fTOFresolution = copy.fTOFresolution;

  return (*this);
}

//__________________________________________________________________________________________________
AliRsnAnalysisMonitorPairTask::~AliRsnAnalysisMonitorPairTask()
{
//
// Destructor
//

  if (fOut)    delete fOut;
  if (fESDpid) delete fESDpid;
  if (fTrack)  delete fTrack;
}

//__________________________________________________________________________________________________
void AliRsnAnalysisMonitorPairTask::UserCreateOutputObjects()
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
  fTrack[0] = new AliRsnMonitorTrack;
  fTrack[1] = new AliRsnMonitorTrack;
    
  // create output tree
  OpenFile(1);
  fOut = new TTree("rsnPairMonitor", "Informations on pairs for cut checking");
  fOut->Branch("track1", "AliRsnMonitorTrack", &fTrack[0]);
  fOut->Branch("track2", "AliRsnMonitorTrack", &fTrack[1]);
  fOut->Branch("minv"  , &fInvMass, "minv/F");
}

//__________________________________________________________________________________________________
void AliRsnAnalysisMonitorPairTask::UserExec(Option_t *)
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
  event.SetRef(esd, fMCEvent);
  if (!fEventCuts.IsSelected(&event)) return;
  
  // check the event
  Int_t type = EventEval(esd);
    
  // if processable, then process it
  if      (type == 0) ProcessESD(esd, esd->GetPrimaryVertexTracks());
  else if (type == 1) ProcessESD(esd, esd->GetPrimaryVertexSPD());
  else return;
  
  // update histogram container
  PostData(1, fOut);
}

//__________________________________________________________________________________________________
void AliRsnAnalysisMonitorPairTask::Terminate(Option_t *)
{
//
// Terminate
//
}

//__________________________________________________________________________________________________
Int_t AliRsnAnalysisMonitorPairTask::EventEval(AliESDEvent *esd)
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
void AliRsnAnalysisMonitorPairTask::ProcessESD(AliESDEvent *esd, const AliESDVertex *v)
{
//
// Process the ESD container, to read all tracks and copy their useful values.
// All info are stored into an AliRsnMonitorTrack object and saved into the
// TClonesArray which is one of the branches of the output TTree.
//

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
  Int_t          i0 , i1, ch0, ch1, ntracks = esd->GetNumberOfTracks();
  TLorentzVector v0, v1, sum;
  
  // two nested loops on tracks
  for (i0 = 0; i0 < ntracks; i0++)
  {
    if (!ProcessTrack(0, i0, esd, v)) continue;
    ch0 = fTrack[0]->Charge();
    v0.SetXYZM(fTrack[0]->PrecX(), fTrack[0]->PrecY(), fTrack[0]->PrecZ(), fMass[0]);
    
    for (i1 = i0 + 1; i1 < ntracks; i1++)
    {
      if (!ProcessTrack(1, i1, esd, v)) continue;
      ch1 = fTrack[1]->Charge();
      
      // skip like-sign pairs
      if (ch1 == ch0) continue;
      
      // check invmass range
      v1.SetXYZM(fTrack[1]->PrecX(), fTrack[1]->PrecY(), fTrack[1]->PrecZ(), fMass[1]);  
      sum = v0 + v1;
      fInvMass = sum.M();
      if (fInvMass < fRangeMin || fInvMass > fRangeMax) continue;
      
      // if here, add an entry to the tree
      fOut->Fill();
    }
  }
}

//__________________________________________________________________________________________________
Bool_t AliRsnAnalysisMonitorPairTask::ProcessTrack(Int_t myIndex, Int_t esdIndex, AliESDEvent *esd, const AliESDVertex *v)
{
//
// Process a single track, and stores all its info into one of the two branch objects
// defined by the first index, while the second chooses the track in the owner ESD
//

  if (myIndex < 0 || myIndex > 1) return kFALSE;
  
  // ITS stuff #1 
  // create the response function and initialize it to MC or not
  // depending if the AliStack object is there or not
  AliStack *stack = 0x0;
  if (fMCEvent) stack = fMCEvent->Stack();
  Bool_t isMC = (stack != 0x0);
  AliITSPIDResponse itsrsp(isMC);

  // create interfacr objects
  AliRsnEvent    event;
  event.SetRef(esd, fMCEvent);
  AliRsnDaughter daughter;
  AliESDtrack *track = esd->GetTrack(esdIndex);
  event.SetDaughter(daughter, esdIndex, AliRsnDaughter::kTrack);
  
  // skip NULL pointers, kink daughters and tracks which
  // cannot be propagated to primary vertex
  if (!track) return kFALSE;
  if ((Int_t)track->GetKinkIndex(0) > 0) return kFALSE;
  if (!track->RelateToVertex(v, esd->GetMagneticField(), kVeryBig)) return kFALSE;
  
  // useful variables
  Int_t     k, nITS;
  Bool_t    isTPC, isITSSA, isTOF;
  Float_t   b[2], bCov[3];
  Double_t  time[10];
  
  // reset the output object
  // 'usable' flag will need to be set to 'ok'
  fTrack[myIndex]->Reset();
  
  // check cuts
  fTrack[myIndex]->CutsPassed() = fTrackCuts.IsSelected(&daughter);
  
  // get MC info if possible
  if (stack) fTrack[myIndex]->AdoptMC(TMath::Abs(track->GetLabel()), stack);
  
  // copy general info
  fTrack[myIndex]->Status() = (UInt_t)track->GetStatus();
  fTrack[myIndex]->Length() = (Double_t)track->GetIntegratedLength();
  fTrack[myIndex]->Charge() = (Int_t)track->Charge();
  fTrack[myIndex]->PrecX()  = (Double_t)track->Px();
  fTrack[myIndex]->PrecY()  = (Double_t)track->Py();
  fTrack[myIndex]->PrecZ()  = (Double_t)track->Pz();
  
  // evaluate some flags from the status to decide what to do next in some points
  isTPC   = ((fTrack[myIndex]->Status() & AliESDtrack::kTPCin)  != 0);
  isITSSA = ((fTrack[myIndex]->Status() & AliESDtrack::kTPCin)  == 0 && (fTrack[myIndex]->Status() & AliESDtrack::kITSrefit) != 0 && (fTrack[myIndex]->Status() & AliESDtrack::kITSpureSA) == 0 && (fTrack[myIndex]->Status() & AliESDtrack::kITSpid) != 0);
  isTOF   = ((fTrack[myIndex]->Status() & AliESDtrack::kTOFout) != 0 && (fTrack[myIndex]->Status() & AliESDtrack::kTIME) != 0);
  
  // accept only tracks which are TPC+ITS or ITS standalone
  if (isITSSA)
  {
    fTrack[myIndex]->ITSsa() = kTRUE;
    fTrack[myIndex]->TOFok() = kFALSE;
  }
  else if (isTPC)
  {
    fTrack[myIndex]->ITSsa() = kFALSE;
    fTrack[myIndex]->TOFok() = isTOF;
  }
  else
    return kFALSE;

  // get DCA to primary vertex
  track->GetImpactParameters(b, bCov);
  fTrack[myIndex]->DCAr() = (Double_t)b[0];
  fTrack[myIndex]->DCAz() = (Double_t)b[1];
  
  // get ITS info
  for (k = 0; k < 6; k++)
  {
    fTrack[myIndex]->ITSmap(k) = track->HasPointOnITSLayer(k);
  }
  if (isITSSA)
  {
    fTrack[myIndex]->ITSchi2() = track->GetITSchi2();
    fTrack[myIndex]->ITSsignal() = track->GetITSsignal();
    nITS = fTrack[myIndex]->SSDcount() + fTrack[myIndex]->SDDcount();
    for (k = 0; k < AliPID::kSPECIES; k++)
    {
      fTrack[myIndex]->ITSnsigma(k) = itsrsp.GetNumberOfSigmas(fTrack[myIndex]->Prec(), fTrack[myIndex]->ITSsignal(), (AliPID::EParticleType)k, nITS, kTRUE);
    }
  }
  
  // get TPC info
  if (isTPC)
  {
    fTrack[myIndex]->TPCcount()  = (Int_t)track->GetTPCclusters(0);
    fTrack[myIndex]->TPCchi2()   = (Double_t)track->GetTPCchi2();
    fTrack[myIndex]->TPCsignal() = (Double_t)track->GetTPCsignal();
    fTrack[myIndex]->PtpcX()     = fTrack[myIndex]->PtpcY() = fTrack[myIndex]->PtpcZ() = 1E10;
    if (track->GetInnerParam())
    {
      fTrack[myIndex]->PtpcX() = track->GetInnerParam()->Px();
      fTrack[myIndex]->PtpcY() = track->GetInnerParam()->Py();
      fTrack[myIndex]->PtpcZ() = track->GetInnerParam()->Pz();
      for (k = 0; k < AliPID::kSPECIES; k++) 
      {
        fTrack[myIndex]->TPCnsigma(k) = fESDpid->NumberOfSigmasTPC(track, (AliPID::EParticleType)k);
      }
    }
  }
  
  // get TOF info
  if (isTOF)
  {
    track->GetIntegratedTimes(time);
    fTrack[myIndex]->TOFsignal() = (Double_t)track->GetTOFsignal();
    for (k = 0; k < AliPID::kSPECIES; k++)
    {
      fTrack[myIndex]->TOFref(k)   = time[k];
      fTrack[myIndex]->TOFsigma(k) = (Double_t)fTOFmaker->GetExpectedSigma(fTrack[myIndex]->Prec(), time[k], AliPID::ParticleMass(k));
    }
  }
  
  return kTRUE;
}
