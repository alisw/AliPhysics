//
// Class AliRsnCutPIDTOF
//
// General implementation of a single cut strategy, which can be:
// - a value contained in a given interval  [--> IsBetween()   ]
// - a value equal to a given reference     [--> MatchesValue()]
//
// In all cases, the reference value(s) is (are) given as data members
// and each kind of cut requires a given value type (Int, UInt, Double),
// but the cut check procedure is then automatized and chosen thanks to
// an enumeration of the implemented cut types.
// At the end, the user (or any other point which uses this object) has
// to use the method IsSelected() to check if this cut has been passed.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <Riostream.h>

#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDpid.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliVTrack.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutPIDTOF.h"

ClassImp(AliRsnCutPIDTOF)

//Bool_t         AliRsnCutPIDTOF::fgTOFcalibrateESD = kTRUE;
Bool_t         AliRsnCutPIDTOF::fgTOFcorrectTExp  = kTRUE;
Bool_t         AliRsnCutPIDTOF::fgTOFuseT0        = kTRUE;
Bool_t         AliRsnCutPIDTOF::fgTOFtuneMC       = kFALSE;
Double_t       AliRsnCutPIDTOF::fgTOFresolution   = 100.0;
AliTOFT0maker* AliRsnCutPIDTOF::fgTOFmaker        = 0x0;
AliTOFcalib*   AliRsnCutPIDTOF::fgTOFcalib        = 0x0;
Int_t          AliRsnCutPIDTOF::fgLastRun         = -1;
Int_t          AliRsnCutPIDTOF::fgLastEventID     = -1;
AliESDEvent*   AliRsnCutPIDTOF::fgLastEvent       = 0x0;


//_________________________________________________________________________________________________
AliRsnCutPIDTOF::AliRsnCutPIDTOF(const char *name, AliPID::EParticleType pid, Bool_t isMC, Double_t min, Double_t max, Bool_t forceMatching) :
  AliRsnCut(name, AliRsnCut::kDaughter, min, max),
  fIsMC(isMC),
  fForceMatching(forceMatching),
  fPIDtype(pid),
  fESDpid(),
  fAODpid()
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDTOF::AliRsnCutPIDTOF(const AliRsnCutPIDTOF& copy) :
  AliRsnCut(copy),
  fIsMC(copy.fIsMC),
  fForceMatching(copy.fForceMatching),
  fPIDtype(copy.fPIDtype),
  fESDpid(copy.fESDpid),
  fAODpid(copy.fAODpid)
{
//
// Copy constructor
//
}

//_________________________________________________________________________________________________
AliRsnCutPIDTOF& AliRsnCutPIDTOF::operator=(const AliRsnCutPIDTOF& copy)
{
//
// Assignment operator
//

  fIsMC = copy.fIsMC;
  fForceMatching = copy.fForceMatching;
  fPIDtype = copy.fPIDtype;
  fESDpid = copy.fESDpid;
  fAODpid = copy.fAODpid;
  
  return (*this);
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDTOF::IsSelected(TObject *object)
{
//
// Cut checker.
//

  // coherence check
  if (!TargetOK(object)) return kFALSE;
  
  // reject always non-track objects
  AliVTrack *vtrack = dynamic_cast<AliVTrack*>(fDaughter->GetRef());
  if (!vtrack)
  {
    AliDebug(AliLog::kDebug + 2, Form("This object is not either an ESD nor AOD track, it is an %s", fDaughter->GetRef()->ClassName()));
    return kFALSE;
  }
  
  // for non TOF-matched tracks, the TOF PID check cannot be done:
  // -- if 'fForceMatching' is kTRUE
  //    all non matched tracks are rejected as if they didn't pass the cut
  // -- if 'fForceMatching' is kFALSE
  //    all non matched tracks are ignored, as if they did pass the cut
  ULong_t status = (ULong_t)vtrack->GetStatus();
  if ((status & AliESDtrack::kTOFout) == 0 || (status & AliESDtrack::kTIME) == 0)
  {
    AliDebug(AliLog::kDebug + 2, "Track is not matched with TOF");
    if (fForceMatching)
      return kFALSE;
    else
      return kTRUE;
  }
  
  // retrieve real object type
  AliESDtrack *esdTrack = fDaughter->GetRefESDtrack();
  AliAODTrack *aodTrack = fDaughter->GetRefAODtrack();
  if (esdTrack) 
  {
    AliDebug(AliLog::kDebug + 2, "Checking an ESD track");
    ProcessCurrentEvent();
    return CheckESD(esdTrack);
  }
  else if (aodTrack)
  {
    AliDebug(AliLog::kDebug + 2, "Checking an AOD track");
    return CheckAOD(aodTrack);
  }
  else
  {
    AliDebug(AliLog::kDebug + 2, Form("This object is not either an ESD nor AOD track, it is an %s", fDaughter->GetRef()->ClassName()));
    return kFALSE;
  }
}

//_________________________________________________________________________________________________
void AliRsnCutPIDTOF::ProcessCurrentEvent()
{
//
// Repeats the PID for current event.
// In order to avoid to repeat the operation for each track
// this function uses the data-member pointers to check
// if current event was already processed or not.
//

  // get current event and skip if it is the same as before
  AliESDEvent *esd = fgCurrentEvent->GetRefESD();
  if (esd == fgLastEvent) return;
  
  // if event has changed, get run number and 
  // reinitialize the calib and time maker
  // if this has changed
  Int_t run = esd->GetRunNumber();
  if (run != fgLastRun)
  {
    AliInfo("============================================================================================");
    AliInfo(Form("*** CHANGING RUN NUMBER: PREVIOUS = %d --> CURRENT = %d ***", fgLastRun, run));
    AliInfo("============================================================================================");
    fgLastRun = run;
  
    AliCDBManager::Instance()->SetDefaultStorage("raw://");
    AliCDBManager::Instance()->SetRun(fgLastRun);
    
    if (fgTOFmaker) delete fgTOFmaker;
    if (fgTOFcalib) delete fgTOFcalib;
    
    fgTOFcalib = new AliTOFcalib();
    if (fIsMC)
    {
      fgTOFcalib->SetRemoveMeanT0(kFALSE);
      fgTOFcalib->SetCalibrateTOFsignal(kFALSE);
    }
    else
    {
      fgTOFcalib->SetRemoveMeanT0(kTRUE);
      fgTOFcalib->SetCalibrateTOFsignal(kTRUE);
    }
    if (fgTOFcorrectTExp) fgTOFcalib->SetCorrectTExp(kTRUE);
    fgTOFcalib->Init();
    
    fgTOFmaker = new AliTOFT0maker(&fESDpid, fgTOFcalib);
    fgTOFmaker->SetTimeResolution(fgTOFresolution);
  }

  // repeat the calibration and PID computations
  /*if (fgTOFcalibrateESD)*/ fgTOFcalib->CalibrateESD(esd);
  if (fgTOFtuneMC) fgTOFmaker->TuneForMC(esd);
  if (fgTOFuseT0)
  {
    fgTOFmaker->ComputeT0TOF(esd);
    fgTOFmaker->ApplyT0TOF(esd);
    fESDpid.MakePID(esd, kFALSE, 0.);
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDTOF::CheckESD(AliESDtrack *track)
{
//
// Check TOF particle identification for ESD tracks.
// Uses the AlifESDpid odata member.
//

  // require a minimum length to have meaningful match
  if (track->GetIntegratedLength() < 350.) return kFALSE;
  
  // setup TOF PID response
  AliTOFPIDResponse &tofrsp = fESDpid.GetTOFResponse();
  
  // get info for computation
  Double_t momentum = track->P();
  Double_t time     = track->GetTOFsignal();
  Double_t timeint[AliPID::kSPECIES];
  tofrsp.GetStartTime(momentum);
  track->GetIntegratedTimes(timeint);
  Double_t timeDiff = time - timeint[(Int_t)fPIDtype];
  Double_t sigmaRef = tofrsp.GetExpectedSigma(momentum, timeint[(Int_t)fPIDtype], AliPID::ParticleMass(fPIDtype));
  
  // port values to standard AliRsnCut checker
  fCutValueD = timeDiff / sigmaRef;
  return OkRangeD();
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutPIDTOF::CheckAOD(AliAODTrack *track) 
{
//
// check TOF particle identification for AOD tracks.
// Uses the AliAODpidUtil data member.
//

  fCutValueD = (Double_t)fAODpid.NumberOfSigmasTOF(track, fPIDtype);
  return OkRangeD();
}

//_________________________________________________________________________________________________
void AliRsnCutPIDTOF::Print(const Option_t *) const
{
//
// Print information on this cut
//

  AliInfo(Form("Cut name, type            : %s %s", GetName(), ClassName()));
  AliInfo(Form("TOF PID cut range (sigmas): %.3f %.3f", fMinD, fMaxD));
  AliInfo(Form("Unmatched tracks are      : %s", (fForceMatching ? "rejected" : "accepted")));
}
