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
  fUseGlobal(kTRUE),
  fUseITSSA(kTRUE),
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
  fMaxTOFband(1E6),
  fLastRun(-1)
{
//
// Default constructor.
//

  Int_t i = 0;
  for (i = 0; i < 5; i++) fTPCpar[i] = 0.0;
}

//_________________________________________________________________________________________________
AliRsnCutESD2010::AliRsnCutESD2010
(const char *name) :
  AliRsnCut(name, AliRsnCut::kDaughter, 0.0, 0.0),
  fIsMC(kFALSE),
  fCheckITS(kTRUE),
  fCheckTPC(kTRUE),
  fCheckTOF(kTRUE),
  fUseGlobal(kTRUE),
  fUseITSSA(kTRUE),
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
  fMaxTOFband(1E6),
  fLastRun(-1)
{
//
// Main constructor.
//

  Int_t i = 0;
  for (i = 0; i < 5; i++) fTPCpar[i] = 0.0;
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
  fMaxTOFband(copy.fMaxTOFband),
  fLastRun(-1)
{
//
// Main constructor.
//

  Int_t i = 0;
  for (i = 0; i < 5; i++) fTPCpar[i] = copy.fTPCpar[i];
}

//_________________________________________________________________________________________________
void AliRsnCutESD2010::InitializeToDefaults(Bool_t isSim)
{
//
// Main constructor.
//

  // ----> set TPC range for PID and calibration
  SetTPCrange(5.0, 3.0);
  SetTPCpLimit(0.35);
  
  // ----> set ITS range for PID
  SetITSband(3.0);
  
  // ----> set TOF range for PID
  SetTOFband(3.0);
  
  // ----> set TPC calibration
  if (isSim) SetTPCpar(2.15898 / 50.0, 1.75295E1, 3.40030E-9, 1.96178, 3.91720);
  else       SetTPCpar(1.41543 / 50.0, 2.63394E1, 5.0411E-11, 2.12543, 4.88663);
  
  // ----> set standard quality cuts for TPC global tracks
  fESDtrackCutsTPC.SetRequireTPCStandAlone(kTRUE); // require to have the projection at inner TPC wall
  fESDtrackCutsTPC.SetMinNClustersTPC(80);
  fESDtrackCutsTPC.SetMaxChi2PerClusterTPC(4);
  fESDtrackCutsTPC.SetAcceptKinkDaughters(kFALSE);
  fESDtrackCutsTPC.SetRequireTPCRefit(kTRUE);
  fESDtrackCutsTPC.SetRequireITSRefit(kTRUE);
  fESDtrackCutsTPC.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  //fESDtrackCutsTPC.SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9"); // DCA pt dependent: 7*(0.0050+0.0060/pt0.9)
  fESDtrackCutsTPC.SetMaxDCAToVertexXYPtDep("0.035+0.03/pt^0.9"); // DCA pt dependent: 5*(0.0050+0.0060/pt0.9)
  fESDtrackCutsTPC.SetMaxDCAToVertexZ(1e6); // disabled
  fESDtrackCutsTPC.SetDCAToVertex2D(kFALSE); // each DCA is checked separately
  fESDtrackCutsTPC.SetRequireSigmaToVertex(kFALSE);
  
  // ----> set standard quality cuts for ITS standalone tracks
  fESDtrackCutsITS.SetRequireITSStandAlone(kTRUE);
  fESDtrackCutsITS.SetRequireITSRefit(kTRUE);
  fESDtrackCutsITS.SetMinNClustersITS(4);
  fESDtrackCutsITS.SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  fESDtrackCutsITS.SetMaxChi2PerClusterITS(2.);
  // fESDtrackCutsITS.SetMaxDCAToVertexXYPtDep("0.0595+0.0182/pt^1.55"); // DCA pt dependent (7 sigma)
  fESDtrackCutsITS.SetMaxDCAToVertexXYPtDep("0.0425+0.013/pt^1.55"); // DCA pt dependent (5 sigma)
  fESDtrackCutsITS.SetMaxDCAToVertexZ(1e6); // disabled
  fESDtrackCutsITS.SetDCAToVertex2D(kFALSE); // each DCA is checked separately
  
  // ----> set the TOF calibration depending on type of input (sim/data)
  SetTOFcorrectTExp(kTRUE);
  SetTOFuseT0(kTRUE);
  SetTOFresolution(100.0);
  if (isSim)
  {
    SetTOFcalibrateESD(kFALSE);
    SetTOFtuneMC(kTRUE);
  }
  else
  {
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
  Bool_t   okQuality, okTOF, okTPC, okITS, isTPC, isITSSA, isTOF;
  UChar_t  itsCluMap;
  
  // get commonly used variables
  status  = (ULong_t)track->GetStatus();
  mom     = track->P();
  isTPC   = ((status & AliESDtrack::kTPCin)  != 0);
  isITSSA = ((status & AliESDtrack::kTPCin)  == 0 && (status & AliESDtrack::kITSrefit) != 0 && (status & AliESDtrack::kITSpureSA) == 0 && (status & AliESDtrack::kITSpid) != 0);
  isTOF   = (((status & AliESDtrack::kTOFout) != 0) && ((status & AliESDtrack::kTIME) != 0) /* && mom > TMath::Max(b1, b2)*/);
  
  // accept only tracks which are TPC+ITS or ITS standalone if their flag is switched on
  if (!isTPC   && !isITSSA) return kFALSE;
  
  // this branch is entered by all global tracks
  // is their usage is required in the initialization
  if (isTPC && fUseGlobal)
  {
    // check standard ESD cuts
    okQuality = fESDtrackCutsTPC.IsSelected(track);
    AliDebug(AliLog::kDebug + 2, Form("Global quality cut = %s", (okQuality ? "GOOD" : "BAD")));
    if (!okQuality) return kFALSE;
  
    // check TPC dE/dx
    if (fCheckTPC)
    {
      tpcNSigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track, AliPID::kKaon));
      if (track->GetInnerParam()->P() > fTPCpLimit) tpcMaxNSigma = fMinTPCband; else tpcMaxNSigma = fMaxTPCband;
      okTPC = (tpcNSigma <= tpcMaxNSigma);
      AliDebug(AliLog::kDebug + 2, Form("TPC nsigma = %f -- max = %f -- cut %s", tpcNSigma, tpcMaxNSigma, (okTPC ? "passed" : "failed")));
      if (!okTPC) return kFALSE;
    }
    
    // check TOF (only if momentum is large than function asymptote and flags are OK)
    if (fCheckTOF && isTOF)
    {
      track->GetIntegratedTimes(times);
      tofTime  = (Double_t)track->GetTOFsignal();
      tofSigma = fTOFmaker->GetExpectedSigma(mom, times[AliPID::kKaon], AliPID::ParticleMass(AliPID::kKaon));
      tofRef   = times[AliPID::kKaon];
      if (tofRef > 0.0)
      {
        /*
        tofRel   = (tofTime - tofRef) / tofRef;
        ymax     = a1 / (mom - b1) + c1;
        ymin     = a2 / (mom - b2) + c2;
        okTOF    = (tofRel >= ymin && tofRel <= ymax);
        */
        tofRel   = TMath::Abs(tofTime - tofRef) / tofSigma;
        okTOF    = (tofRel <= fMaxTOFband);
        AliDebug(AliLog::kDebug + 2, Form("TOF nsigma = %f -- max = %f -- cut %s", tofRel, fMaxTOFband, (okTOF ? "passed" : "failed")));
        if (!okTOF) return kFALSE;
      }
    }
    else
    {
      //
      // the opposite of previous condition (fCheckTOF && isTOF) is
      // NOT fCheckTOF OR NOT isTOF
      
      // if TOF must not be checked, TPC is assumed to be checked
      // otherwise, just quality cuts have been checked, but anyway
      // all other cuts have already been checked and no further check 
      // is needed here;
      // Just print a debug message
      if (!isTOF) AliDebug(AliLog::kDebug + 2, "TOF not matched");
      
      // instead, if track has not a match in TOF, we must reject it
      // in case it was not already checked in the TPC, since we cannot
      // say anything about not matched track if TOF is the only detector
      if (!fCheckTPC && !isTOF) return kFALSE;
    }
      
    // if we arrive here, the cut is passed,
    // since in all points where something goew wrong,
    // and exit point is implemented which returns kFALSE
    return kTRUE;
  }
  
  // this branch is entered by all ITS standalone tracks
  // is their usage is required in the initialization
  if (isITSSA && fUseITSSA)
  {
    // check standard ESD cuts
    okQuality = fESDtrackCutsITS.IsSelected(track);
    AliDebug(AliLog::kDebug + 2, Form("Global quality cut = %s", (okQuality ? "GOOD" : "BAD")));
    if (!okQuality) return kFALSE;
    
    // check dE/dx
    itsSignal = track->GetITSsignal();
    itsCluMap = track->GetITSClusterMap();
    nITS      = 0;
    for(k = 2; k < 6; k++) if(itsCluMap & (1 << k)) ++nITS;
    if (nITS < 3) return kFALSE;
    itsNSigma = itsrsp.GetNumberOfSigmas(mom, itsSignal, AliPID::kKaon, nITS, kTRUE);
    okITS = (TMath::Abs(itsNSigma) <= fMaxITSband);
    AliDebug(AliLog::kDebug + 2, Form("ITS nsigma = %f -- max = %f -- cut %s", itsNSigma, fMaxITSband, (okITS ? "passed" : "failed")));
    if (!okITS) return kFALSE;
    
    // if we arrive here, the cut is passed,
    // since in all points where something goew wrong,
    // and exit point is implemented which returns kFALSE
    return kTRUE;
  }
  
  // if we are here, the track is surely bad
  return kFALSE;
}
