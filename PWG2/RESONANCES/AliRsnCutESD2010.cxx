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
AliRsnCutESD2010::AliRsnCutESD2010
(const char *name, Bool_t isMC) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),
  fIsMC(isMC),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fUseGlobal(kTRUE),
  fUseITSSA(kTRUE),
  fMaxEta(1E6),
  fMaxITSband(3.0),
  fTPCpLimit(0.35),
  fMinTPCband(3.0),
  fMaxTPCband(5.0),
  fESDtrackCutsTPC(),
  fESDtrackCutsITS(),
  fESDpid(0x0),
  fTOFmaker(0x0),
  fTOFcalib(0x0),
  fTOFcalibrateESD(!isMC),
  fTOFcorrectTExp(kTRUE),
  fTOFuseT0(kTRUE),
  fTOFtuneMC(isMC),
  fTOFresolution(100.0),
  fMinTOF(-2.5),
  fMaxTOF( 3.5),
  fLastRun(-1)
{
//
// Main constructor.
//

  SetMC(isMC);
  
  // set default quality cuts for TPC+ITS tracks
  // TPC  
  fESDtrackCutsTPC.SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
  fESDtrackCutsTPC.SetMinNClustersTPC(70);
  fESDtrackCutsTPC.SetMaxChi2PerClusterTPC(4);
  fESDtrackCutsTPC.SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsTPC.SetRequireTPCRefit(kTRUE);
  // ITS
  fESDtrackCutsTPC.SetRequireITSRefit(kTRUE);
  fESDtrackCutsTPC.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fESDtrackCutsTPC.SetMaxDCAToVertexXYPtDep("0.0350+0.0490/pt^1.0");
  fESDtrackCutsTPC.SetMaxDCAToVertexZ(1.e6);
  fESDtrackCutsTPC.SetDCAToVertex2D(kFALSE);
  fESDtrackCutsTPC.SetRequireSigmaToVertex(kFALSE);
  fESDtrackCutsTPC.SetEtaRange(-fMaxEta, fMaxEta);
  
  // set default quality cuts for ITS standalone tracks
  fESDtrackCutsITS.SetRequireITSStandAlone(kTRUE);
  fESDtrackCutsITS.SetRequireITSPureStandAlone(kFALSE);
  fESDtrackCutsITS.SetRequireITSRefit(kTRUE); 
  fESDtrackCutsITS.SetMinNClustersITS(4);
  fESDtrackCutsITS.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fESDtrackCutsITS.SetMaxChi2PerClusterITS(2.5);
  fESDtrackCutsITS.SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55");
  fESDtrackCutsITS.SetEtaRange(-fMaxEta, fMaxEta);
}

//_________________________________________________________________________________________________
AliRsnCutESD2010::AliRsnCutESD2010
(const AliRsnCutESD2010& copy) :
  AliRsnCut(copy),
  fIsMC(copy.fIsMC),
  fCheckITS(copy.fCheckITS),
  fCheckTPC(copy.fCheckTPC),
  fCheckTOF(copy.fCheckTOF),
  fUseGlobal(copy.fUseGlobal),
  fUseITSSA(copy.fUseITSSA),
  fMaxEta(copy.fMaxEta),
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
  fMinTOF(copy.fMinTOF),
  fMaxTOF(copy.fMaxTOF),
  fLastRun(-1)
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
  fUseGlobal = copy.fUseGlobal;
  fUseITSSA = copy.fUseITSSA;
  fMaxEta = copy.fMaxEta;
  fMaxITSband = copy.fMaxITSband;
  fTPCpLimit = copy.fTPCpLimit;
  fMinTPCband = copy.fMinTPCband;
  fMaxTPCband = copy.fMaxTPCband;
  fESDtrackCutsTPC = copy.fESDtrackCutsTPC;
  fESDtrackCutsITS = copy.fESDtrackCutsITS;
  fTOFcalibrateESD = copy.fTOFcalibrateESD;
  fTOFcorrectTExp = copy.fTOFcorrectTExp;
  fTOFuseT0 = copy.fTOFuseT0;
  fTOFtuneMC = copy.fTOFtuneMC;
  fTOFresolution = copy.fTOFresolution;
  fMinTOF = copy.fMinTOF;
  fMaxTOF = copy.fMaxTOF;
  fLastRun = copy.fLastRun;
  
  Int_t i = 0;
  for (i = 0; i < 5; i++) fTPCpar[i] = copy.fTPCpar[i];
  
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

  // if absent, initialize ESD pid responst
  if (!fESDpid)
  {
    fESDpid = new AliESDpid;
    fESDpid->GetTPCResponse().SetBetheBlochParameters(fTPCpar[0],fTPCpar[1],fTPCpar[2],fTPCpar[3],fTPCpar[4]);
  }

  // initialize DB to current run
  Int_t run = esd->GetRunNumber();
  if (run != fLastRun)
  {
    cout << "Run = " << run << " -- LAST = " << fLastRun << endl;
    fLastRun = run;

    // setup TOF maker & calibration
    if (!fTOFcalib) fTOFcalib = new AliTOFcalib;
    fTOFcalib->SetCorrectTExp(fTOFcorrectTExp);
    if (!fTOFmaker) fTOFmaker = new AliTOFT0maker(fESDpid, fTOFcalib);
    fTOFmaker->SetTimeResolution(fTOFresolution);

    AliCDBManager *cdb = AliCDBManager::Instance();
    cdb->ClearCache(); // suggestion by Annalisa
    cdb->Clear();      // suggestion by Annalisa
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

  // ITS: create the response function
  AliITSPIDResponse itsrsp(fIsMC);

  // TOF: define fixed function for compatibility range
  //Double_t a1 = 0.01, a2 = -0.03;
  //Double_t b1 = 0.25, b2 =  0.25;
  //Double_t c1 = 0.05, c2 = -0.03;
  //Double_t ymax, ymin;

  ULong_t  status;
  Int_t    k, nITS;
  Double_t times[10], tpcNSigma, tpcMaxNSigma, itsSignal, itsNSigma, mom, tofTime, tofSigma, tofRef, tofRel;
  Bool_t   okQuality, okTOF, okTPC, okITS, isTPC, isITSSA, isTOF;
  UChar_t  itsCluMap;

  // get commonly used variables
  status  = (ULong_t)track->GetStatus();
  mom     = track->P();
  isTPC   = ((status & AliESDtrack::kTPCin) != 0);
  isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
  isTOF   = (((status & AliESDtrack::kTOFout) != 0) && ((status & AliESDtrack::kTIME) != 0) /* && mom > TMath::Max(b1, b2)*/);
  
  // check if the track type matches what is required
  if (!isTPC && !isITSSA) 
  {
    AliDebug(AliLog::kDebug + 2, "Track is not either a TPC track or a ITS standalone. Rejected");
    return kFALSE;
  }
  else if (isTPC && !fUseGlobal)
  {
    AliDebug(AliLog::kDebug + 2, "Global tracks not used. Rejected");
    return kFALSE;
  }
  else if (isITSSA && !fUseITSSA)
  {
    AliDebug(AliLog::kDebug + 2, "ITS standalone not used. Rejected");
    return kFALSE;
  }
  
  // does a preliminary check on TOF values, if necessary
  // then, if the reference time or TOF signal are meaningless
  // even if the 'isTOF' flag is true, switch it to false
  if (isTOF)
  {
    track->GetIntegratedTimes(times);
    tofTime  = (Double_t)track->GetTOFsignal();
    tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kKaon], AliPID::ParticleMass(AliPID::kKaon));
    tofRef   = times[AliPID::kKaon];
    if (tofRef <= 0.0 && tofSigma <= 0.0) isTOF = kFALSE;
  }
  
  // check quality (eta range must be adapted)
  AliESDtrackCuts *cuts = 0x0;
  if (isTPC)   cuts = &fESDtrackCutsTPC;
  if (isITSSA) cuts = &fESDtrackCutsITS;
  if (!cuts)   return kFALSE;
  cuts->SetEtaRange(-fMaxEta, fMaxEta);
  okQuality = cuts->IsSelected(track);
  AliDebug(AliLog::kDebug + 2, Form("Global quality cut = %s", (okQuality ? "GOOD" : "BAD")));
  if (!okQuality) return kFALSE;
  
  if (isTPC) // this branch is entered by all global tracks
  {
    // check TPC dE/dx:
    if (fCheckTPC)
    {
      tpcNSigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kKaon));
      if (track->GetInnerParam()->P() > fTPCpLimit) tpcMaxNSigma = fMinTPCband; else tpcMaxNSigma = fMaxTPCband;
      okTPC = (tpcNSigma <= tpcMaxNSigma);
      AliDebug(AliLog::kDebug + 2, Form("TPC nsigma = %f -- max = %f -- cut %s", tpcNSigma, tpcMaxNSigma, (okTPC ? "passed" : "failed")));
    }
    else
    {
      // if TPC is not checked, it is as if all tracks do pass the cut
      okTPC = kTRUE;
      AliDebug(AliLog::kDebug + 2, "TPC not checked, track accepted");
    }

    // check TOF (only if flags are OK)
    if (fCheckTOF)
    {
      if (isTOF)
      {
        // TOF can be checked only when track is matched there
        track->GetIntegratedTimes(times);
        tofTime  = (Double_t)track->GetTOFsignal();
        tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kKaon], AliPID::ParticleMass(AliPID::kKaon));
        tofRef   = times[AliPID::kKaon];
        /*
        tofRel   = (tofTime - tofRef) / tofRef;
        ymax     = a1 / (mom - b1) + c1;
        ymin     = a2 / (mom - b2) + c2;
        okTOF    = (tofRel >= ymin && tofRel <= ymax);
        */
        tofRel   = (tofTime - tofRef) / tofSigma;
        okTOF    = (tofRel >= fMinTOF && tofRel <= fMaxTOF);
        AliDebug(AliLog::kDebug + 2, Form("TOF nsigma = %f -- range = %f %f -- cut %s", tofRel, fMinTOF, fMaxTOF, (okTOF ? "passed" : "failed")));
      }
      else
      {
        // if TOF is not matched, the answer depends on TPC:
        // - if TPC is required, track is checked only there and TOF simply ignored
        // - if TPC is not required, track is rejected when TOF does not match it, if TOF check is required
        if (fCheckTPC) okTOF = kTRUE; else okTOF = kFALSE;
      }
    }
    else
    {
      okTOF = kTRUE;
    }
    
    // properly combine the outcome of TPC and TOF cuts
    return okTPC && okTOF;
  }
  else if (isITSSA) // this branch is entered by all ITS standalone tracks
  {
    // check dE/dx only if this is required, otherwise ITS standalone are just added but not checked for PID
    if (fCheckITS)
    {
      itsSignal = track->GetITSsignal();
      itsCluMap = track->GetITSClusterMap();
      nITS      = 0;
      for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITS;
      if (nITS < 3) return kFALSE;
      itsNSigma = itsrsp.GetNumberOfSigmas(mom, itsSignal, AliPID::kKaon, nITS, kTRUE);
      okITS = (TMath::Abs(itsNSigma) <= fMaxITSband);
      AliDebug(AliLog::kDebug + 2, Form("ITS nsigma = %f -- max = %f -- cut %s", itsNSigma, fMaxITSband, (okITS ? "passed" : "failed")));
    }
    else
    {
      okITS = kTRUE;
    }

    return okITS;
  }
  else
  {
    // if we are here, the track is surely bad
    return kFALSE;
  }
}
