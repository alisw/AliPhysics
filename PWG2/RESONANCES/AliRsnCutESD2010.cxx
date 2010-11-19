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

Int_t AliRsnCutESD2010::fgLastRun = -1;

//_________________________________________________________________________________________________
AliRsnCutESD2010::AliRsnCutESD2010
(const char *name, Bool_t isMC) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),
  fIsMC(isMC),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fUseITSTPC(kTRUE),
  fUseITSSA(kTRUE),
  fPID(AliPID::kKaon),
  fMaxITSPIDmom(0.0),
  fMaxITSband(3.0),
  fTPCpLimit(0.35),
  fMinTPCband(3.0),
  fMaxTPCband(5.0),
  fESDpid(0x0),
  fESDtrackCutsTPC(),
  fESDtrackCutsITS(),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(!isMC),
  fTOFcorrectTExp(kTRUE),
  fTOFuseT0(kTRUE),
  fTOFtuneMC(isMC),
  fTOFresolution(100.0),
  fMinTOF(-2.5),
  fMaxTOF( 3.5)
{
//
// Main constructor.
//

  SetMC(isMC);
}

//_________________________________________________________________________________________________
AliRsnCutESD2010::AliRsnCutESD2010
(const AliRsnCutESD2010& copy) :
  AliRsnCut(copy),
  fIsMC(copy.fIsMC),
  fCheckITS(copy.fCheckITS),
  fCheckTPC(copy.fCheckTPC),
  fCheckTOF(copy.fCheckTOF),
  fUseITSTPC(copy.fUseITSTPC),
  fUseITSSA(copy.fUseITSSA),
  fPID(copy.fPID),
  fMaxITSPIDmom(copy.fMaxITSPIDmom),
  fMaxITSband(copy.fMaxITSband),
  fTPCpLimit(copy.fTPCpLimit),
  fMinTPCband(copy.fMinTPCband),
  fMaxTPCband(copy.fMaxTPCband),
  fESDpid(0x0),
  fESDtrackCutsTPC(copy.fESDtrackCutsTPC),
  fESDtrackCutsITS(copy.fESDtrackCutsITS),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(copy.fTOFcalibrateESD),
  fTOFcorrectTExp(copy.fTOFcorrectTExp),
  fTOFuseT0(copy.fTOFuseT0),
  fTOFtuneMC(copy.fTOFtuneMC),
  fTOFresolution(copy.fTOFresolution),
  fMinTOF(copy.fMinTOF),
  fMaxTOF(copy.fMaxTOF)
{
//
// Copy constructor.
//

  Int_t i = 0;
  for (i = 0; i < 5; i++) fTPCpar[i] = copy.fTPCpar[i];
}

//_________________________________________________________________________________________________
AliRsnCutESD2010& AliRsnCutESD2010::operator=(const AliRsnCutESD2010& copy)
{
//
// Assignment operator
//

  AliRsnCut::operator=(copy);

  fIsMC = copy.fIsMC;
  fCheckITS = copy.fCheckITS;
  fCheckTPC = copy.fCheckTPC;
  fCheckTOF = copy.fCheckTOF;
  fUseITSTPC = copy.fUseITSTPC;
  fUseITSSA = copy.fUseITSSA;
  fPID = copy.fPID;
  fMaxITSPIDmom = copy.fMaxITSPIDmom;
  fMaxITSband = copy.fMaxITSband;
  fTPCpLimit = copy.fTPCpLimit;
  fMinTPCband = copy.fMinTPCband;
  fMaxTPCband = copy.fMaxTPCband;
  fTOFcalibrateESD = copy.fTOFcalibrateESD;
  fTOFcorrectTExp = copy.fTOFcorrectTExp;
  fTOFuseT0 = copy.fTOFuseT0;
  fTOFtuneMC = copy.fTOFtuneMC;
  fTOFresolution = copy.fTOFresolution;
  fMinTOF = copy.fMinTOF;
  fMaxTOF = copy.fMaxTOF;
  
  Int_t i = 0;
  for (i = 0; i < 5; i++) fTPCpar[i] = copy.fTPCpar[i];
  
  delete fESDpid;
  fESDpid = 0x0;
  
  copy.GetCutsTPC()->Copy(fESDtrackCutsTPC);
  copy.GetCutsITS()->Copy(fESDtrackCutsITS);
  
  return (*this);
}

//_________________________________________________________________________________________________
void AliRsnCutESD2010::SetMC(Bool_t isMC)
{
//
// Sets some aspects of cuts depending on the fact that runs on MC or not
//

  fIsMC = isMC;
  
  if (isMC)
  {
    SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
    SetTOFcalibrateESD(kFALSE);
    SetTOFtuneMC(kTRUE);
  }
  else
  {
    SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
    SetTOFcalibrateESD(kTRUE);
    SetTOFtuneMC(kFALSE);
  }
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
  
  // if the run number has changed,
  // reset the CDB in order to point to the correct one
  Int_t run = esd->GetRunNumber();
  if (run != fgLastRun)
  {
    AliInfo("============================================================================================");
    AliInfo(Form("*** CHANGING RUN NUMBER: PREVIOUS = %d --> CURRENT = %d ***", fgLastRun, run));
    AliInfo("============================================================================================");
    fgLastRun = run;

    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->ClearCache(); // suggestion by Annalisa
    cdb->Clear();      // suggestion by Annalisa
    cdb->SetDefaultStorage("raw://");
    cdb->SetRun(run);
  }
  
  // if absent, initialize ESD pid response
  if (!fESDpid)
  {
    AliITSPIDResponse itsresponse(fIsMC);
    
    fESDpid = new AliESDpid;
    fESDpid->GetTPCResponse().SetBetheBlochParameters(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);
    fESDpid->GetITSResponse() = itsresponse;
  }
  
  // if absent, initialize TOF objects for calibration
  if (!fTOFcalib) fTOFcalib = new AliTOFcalib;
  fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
  fTOFcalib->Init();
  if (!fTOFmaker) fTOFmaker = new AliTOFT0maker(fESDpid, fTOFcalib);
  fTOFmaker->SetTimeResolution(fTOFresolution);

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
Bool_t AliRsnCutESD2010::IsSelected(TObject *object)
{
//
// Cut checker.
//

  // coherence check: require an ESD track
  AliRsnDaughter *daughter = dynamic_cast<AliRsnDaughter*>(object);
  if (!daughter) return kFALSE;
  AliESDtrack *track = daughter->GetRefESDtrack();
  if (!track) return kFALSE;

  // if no reference event, skip
  if (!fEvent) return kFALSE;

  // check quality and track type and reject tracks not passing this step
  if (!OkQuality(track))
  {
    AliDebug(AliLog::kDebug+2, "Failed quality cut");
    return kFALSE;
  }
  
  // ITS PID can be checked always
  // if PID is not required, the flag is sed as
  // if the cut was alsways passed 
  Bool_t okITSpid = OkITSPID(track);
  if (!fCheckITS) okITSpid = kTRUE;
  
  // TPC PID can be checked only for TPC+ITS tracks
  // if PID is not required, the flag is sed as
  // if the cut was alsways passed
  Bool_t okTPCpid = kFALSE;
  if (IsITSTPC(track)) okTPCpid = OkTPCPID(track);
  if (!fCheckTPC) okTPCpid = kTRUE;
  
  // TOF PID can be checked only if TOF is matched
  // if PID is not required, the flag is sed as
  // if the cut was alsways passed
  Bool_t okTOFpid = kFALSE;
  if (IsITSTPC(track) && MatchTOF(track)) okTOFpid = OkTOFPID(track);
  if (!fCheckTOF) okTOFpid = kTRUE;
  
  // now combine all outcomes according to the different possibilities:
  // -- ITS standalone:
  //    --> only ITS PID, always
  // -- ITS + TPC:
  //    --> ITS PID, only for momenta lower than 'fMaxITSPIDmom' and when the ITSpid flag is active
  //    --> TPC PID, always --> MASTER (first to be checked, if fails, track is rejected)
  //    --> TOF PID, only if matched
  if (IsITSSA(track))
  {
    if (!okITSpid)
    {
      AliDebug(AliLog::kDebug+2, "ITS standalone track --> ITS PID failed");
      return kFALSE;
    }
  }
  else // checking IsITSTPC() is redundant due to OkQuality() cut check
  {
    if (!okTPCpid)
    {
      AliDebug(AliLog::kDebug+2, "ITS+TPC track --> TPC PID failed");
      return kFALSE;
    }
    else if (MatchTOF(track) && !okTOFpid)
    {
      AliDebug(AliLog::kDebug+2, "ITS+TPC track --> TOF matched but TOF PID failed");
      return kFALSE;
    }
    else if (track->IsOn(AliESDtrack::kITSpid) && track->P() <= fMaxITSPIDmom && !okITSpid)
    {
      AliDebug(AliLog::kDebug+2, Form("ITS+TPC track --> Momentum lower than limit (%.2f) and ITS PID failed", fMaxITSPIDmom));
      return kFALSE;
    }
  }
  
  // arriving here, the track has survived all checks
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliRsnCutESD2010::OkQuality(AliESDtrack *track)
{
//
// Check track quality parameters.
// Rejects all tracks which are not either TPC+ITS nor ITS standalone.
// If tracks of any type are not flagged to be used, they are rejected anyway.
//

  if (IsITSTPC(track)) return (fUseITSTPC && fESDtrackCutsTPC.IsSelected(track));
  if (IsITSSA (track)) return (fUseITSSA  && fESDtrackCutsITS.IsSelected(track));

  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliRsnCutESD2010::OkITSPID (AliESDtrack *track)
{
//
// Check ITS particle identification with 3sigma cut
//

  // count PID layers and reject if they are too few
  Int_t   k, nITSpidLayers = 0;
  UChar_t itsCluMap = track->GetITSClusterMap();
  for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITSpidLayers;
  if (nITSpidLayers < 3)
  {
    AliDebug(AliLog::kDebug+2, "Rejecting track with too few ITS pid layers");
    return kFALSE;
  }
  
  // check the track type (ITS+TPC or ITS standalone)
  // and reject it if it is of none of the allowed types
  Bool_t isSA = kFALSE;
  if (IsITSTPC(track)) isSA = kFALSE;
  else if (IsITSSA(track)) isSA = kTRUE;
  else
  {
    AliWarning("Track is neither ITS+TPC nor ITS standalone");
    return kFALSE;
  }
  
  // create the PID response object and compute nsigma
  AliITSPIDResponse &itsrsp = fESDpid->GetITSResponse();
  Double_t mom    = track->P();
  Double_t nSigma = itsrsp.GetNumberOfSigmas(mom, track->GetITSsignal(), fPID, nITSpidLayers, isSA);
  
  // evaluate the cut
  Bool_t ok = (TMath::Abs(nSigma) <= fMaxITSband);
  
  // debug message
  AliDebug(AliLog::kDebug + 2, Form("ITS nsigma = %f -- max = %f -- cut %s", nSigma, fMaxITSband, (ok ? "passed" : "failed")));
  
  // outcome
  return ok;
}

//______________________________________________________________________________
Bool_t AliRsnCutESD2010::OkTPCPID (AliESDtrack *track)
{
//
// Check TPC particle identification with {3|5}sigmacut,
// depending on the track total momentum.
//

  // setup TPC PID response
  AliTPCPIDResponse &tpcrsp = fESDpid->GetTPCResponse();
  tpcrsp.SetBetheBlochParameters(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);
  
  // get momentum and number of sigmas and choose the reference band
  Double_t mom       = track->GetInnerParam()->P();
  Double_t nSigma    = tpcrsp.GetNumberOfSigmas(mom, track->GetTPCsignal(), track->GetTPCsignalN(), fPID);
  Double_t maxNSigma = fMaxTPCband;
  if (mom < fTPCpLimit) maxNSigma = fMinTPCband;
  
  // evaluate the cut
  Bool_t ok = (TMath::Abs(nSigma) <= maxNSigma);
  
  // debug message
  AliDebug(AliLog::kDebug + 2, Form("TPC nsigma = %f -- max = %f -- cut %s", nSigma, maxNSigma, (ok ? "passed" : "failed")));
  
  // outcome
  return ok;
}

//______________________________________________________________________________
Bool_t AliRsnCutESD2010::OkTOFPID (AliESDtrack *track)
{
//
// Check TOF particle identification if matched there.
//

  // check integrated length
  if (track->GetIntegratedLength() < 350.) return kFALSE;
  
  // setup TOF PID response
  AliTOFPIDResponse &tofrsp = fESDpid->GetTOFResponse();
  
  // get info for computation
  Double_t momentum = track->P();
  Double_t time     = track->GetTOFsignal();
  Double_t timeint[AliPID::kSPECIES];
  tofrsp.GetStartTime(momentum);
  track->GetIntegratedTimes(timeint);

  // check the cut
  Double_t timeDiff = time - timeint[(Int_t)fPID];
  Double_t sigmaRef = tofrsp.GetExpectedSigma(momentum, timeint[(Int_t)fPID], AliPID::ParticleMass(fPID));
  Double_t nSigma   = timeDiff / sigmaRef;
  
  // evaluate the cut
  Bool_t ok = (nSigma >= fMinTOF && nSigma <= fMaxTOF);
  
  // debug message
  AliDebug(AliLog::kDebug + 2, Form("TOF nsigma = %f -- range = %f - %f -- cut %s", nSigma, fMinTOF, fMaxTOF, (ok ? "passed" : "failed")));
  
  // outcome
  return ok;
}
