//
// Class AliRsnCutESD2010
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

#include "AliESDpid.h"
#include "AliTOFT0maker.h"
#include "AliTOFcalib.h"
#include "AliCDBManager.h"
#include "AliITSPIDResponse.h"

#include "AliRsnEvent.h"
#include "AliRsnDaughter.h"
#include "AliRsnCutESD2010.h"

ClassImp(AliRsnCutESD2010)

//_________________________________________________________________________________________________
AliRsnCutESD2010::AliRsnCutESD2010() :
  AliRsnCut(AliRsnCut::kDaughter),
  fIsMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fMaxITSband(1E6),
  fTPCpLimit(0.35),
  fMinTPCband(-1E6),
  fMaxTPCband( 1E6),
  fESDtrackCutsTPC(),
  fESDtrackCutsITS(),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0),
  fLastRun(-1)
{
//
// Default constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutESD2010::AliRsnCutESD2010
(const char *name) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),
  fIsMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fMaxITSband(1E6),
  fTPCpLimit(0.35),
  fMinTPCband(-1E6),
  fMaxTPCband( 1E6),
  fESDtrackCutsTPC(),
  fESDtrackCutsITS(),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(kFALSE),
  fTOFcorrectTExp(kFALSE),
  fTOFuseT0(kFALSE),
  fTOFtuneMC(kFALSE),
  fTOFresolution(0.0),
  fLastRun(-1)
{
//
// Main constructor.
//
}

//_________________________________________________________________________________________________
AliRsnCutESD2010::AliRsnCutESD2010
(const AliRsnCutESD2010& copy) :
  AliRsnCut(copy),
  fIsMC(kFALSE),
  fCheckITS(copy.fCheckITS),
  fCheckTPC(copy.fCheckTPC),
  fCheckTOF(copy.fCheckTOF),
  fMaxITSband(copy.fMaxITSband),
  fTPCpLimit(copy.fTPCpLimit),
  fMinTPCband(copy.fMinTPCband),
  fMaxTPCband(copy.fMaxTPCband),
  fESDtrackCutsTPC(copy.fESDtrackCutsTPC),
  fESDtrackCutsITS(copy.fESDtrackCutsITS),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(copy.fTOFcalibrateESD),
  fTOFcorrectTExp(copy.fTOFcorrectTExp),
  fTOFuseT0(copy.fTOFuseT0),
  fTOFtuneMC(copy.fTOFtuneMC),
  fTOFresolution(copy.fTOFresolution),
  fLastRun(-1)
{
//
// Main constructor.
//

  Initialize();
}

//_________________________________________________________________________________________________
void AliRsnCutESD2010::Initialize()
{
//
// Main constructor.
//

  // setup TPC response
  fESDpid = new AliESDpid;
  fESDpid->GetTPCResponse().SetBetheBlochParameters(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);
  
  // setup TOF maker & calibration
  fTOFcalib = new AliTOFcalib;
  fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
  fTOFmaker = new AliTOFT0maker(fESDpid, fTOFcalib);
  fTOFmaker->SetTimeResolution(fTOFresolution);
  fLastRun = -1;
}

//_________________________________________________________________________________________________
void AliRsnCutESD2010::SetEvent(AliRsnEvent *event)
{
  // don't do anything if the event is the same as before
  if (fEvent != 0x0 && fEvent == event) return;
  
  // retrieve the ESD event
  AliESDEvent *esd = event->GetRefESD();
  if (!esd)
  {
    fEvent = 0x0;
    return;
  }
  else
  {
    fEvent = event;
  }
  
  // initialize DB to current run
  Int_t run = esd->GetRunNumber();
  if (run != fLastRun)
  {
    cout << "Run = " << run << " -- LAST = " << fLastRun << endl;
    fLastRun = run;
    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->SetDefaultStorage("raw://");
    cdb->SetRun(run);
    fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
    fTOFcalib->Init();
  }
  
  // if required, calibrate the TOF t0 maker with current event
  if (fTOFcalibrateESD) fTOFcalib->CalibrateESD(esd);
  if (fTOFtuneMC) fTOFmaker->TuneForMC(esd);
  if (fTOFuseT0) 
  {
    fTOFmaker->ComputeT0TOF(esd);
    fTOFmaker->ApplyT0TOF(esd);
    fESDpid->MakePID(esd, kFALSE, 0.);
  }
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutESD2010::IsSelected(TObject *obj1, TObject* /*obj2*/)
{
//
// Cut checker.
//

  // coherence check: require an ESD track
  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(obj1);
  if (!daughter) return kFALSE;
  AliESDtrack *track = daughter->GetRefESDtrack();
  if (!track) return kFALSE;
  
  // if no reference event, skip
  if (!fEvent) return kFALSE;
  
  // ITS stuff #1 create the response function
  AliITSPIDResponse itsrsp(fIsMC);
  
  // TOF: define fixed function for compatibility range
  Double_t a1 = 0.01, a2 = -0.03;
  Double_t b1 = 0.25, b2 =  0.25;
  Double_t c1 = 0.05, c2 = -0.03;
  Double_t ymax, ymin;
  
  ULong_t  status;
  Int_t    k, nITS;
  Double_t times[10], tpcNSigma, tpcMaxNSigma, itsSignal, itsNSigma, mom, tofTime, tofSigma, tofRef, tofRel;
  Bool_t   okTOF, isTPC, isITSSA;
  UChar_t  itsCluMap;
  
  // get commonly used variables
  status  = (ULong_t)track->GetStatus();
  mom     = track->P();
  isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
  isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
  
  // accept only tracks which are TPC+ITS or ITS standalone
  if (!isTPC && !isITSSA) return kFALSE;
  
  if (isTPC)
  {
    // check standard ESD cuts
    if (!fESDtrackCutsTPC.IsSelected(track)) return kFALSE;
    
    // check TPC dE/dx
    if (fCheckTPC)
    {
      tpcNSigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kKaon));
      if (track->GetInnerParam()->P() > fTPCpLimit) tpcMaxNSigma = fMinTPCband; else tpcMaxNSigma = fMaxTPCband;
      if (tpcNSigma > tpcMaxNSigma) return kFALSE;
    }
    
    // check TOF (only if momentum is large than function asymptote and flags are OK)
    if (fCheckTOF)
    {
      okTOF = kTRUE;
      if (((status & AliESDtrack::kTOFout) != 0) && ((status & AliESDtrack::kTIME) != 0) && mom > TMath::Max(b1, b2))
      {
        track->GetIntegratedTimes(times);
        tofTime  = (Double_t)track->GetTOFsignal();
        tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kKaon], AliPID::ParticleMass(AliPID::kKaon));
        tofRef   = times[AliPID::kKaon];
        if (tofRef > 0.0)
        {
          tofRel   = (tofTime - tofRef) / tofRef;
          ymax     = a1 / (mom - b1) + c1;
          ymin     = a2 / (mom - b2) + c2;
          okTOF    = (tofRel >= ymin && tofRel <= ymax);
        }
      }
      if (!okTOF) return kFALSE;
    }
  }
  else
  {
    // check standard ESD cuts
    if (!fESDtrackCutsITS.IsSelected(track)) return kFALSE;
    
    // check dE/dx
    if (fCheckITS)
    {
      if ((status & AliESDtrack::kITSpid)  != 0) return kFALSE;
      itsSignal = track->GetITSsignal();
      itsCluMap = track->GetITSClusterMap();
      nITS      = 0;
      for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITS;
      if (nITS < 3) return kFALSE; // track not good for PID
      itsNSigma = itsrsp.GetNumberOfSigmas(mom, itsSignal, AliPID::kKaon, nITS, kTRUE);
      if (TMath::Abs(itsNSigma) > fMaxITSband) return kFALSE;
    }
  }
  
  return kTRUE;
}
