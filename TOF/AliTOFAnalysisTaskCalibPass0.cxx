/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
***************************************************************************/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides TOF pass0/passX calibration tools   //
//                                                           //
///////////////////////////////////////////////////////////////

#include "AliTOFAnalysisTaskCalibPass0.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliPhysicsSelection.h"
#include "AliESDtrackCuts.h"
#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliGRPObject.h"
#include "AliTOFcalib.h"
#include "TList.h"
#include "TH2F.h"
#include "TFile.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TString.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TF1.h"
#include "TTime.h"
#include "TGrid.h"
#include "TTimeStamp.h"
#include "AliTOFRunParams.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "TSystem.h"
#include "AliLog.h"

ClassImp(AliTOFAnalysisTaskCalibPass0)
  
//_______________________________________________________

const Char_t *AliTOFAnalysisTaskCalibPass0::fgkStatusCodeName[AliTOFAnalysisTaskCalibPass0::kNStatusCodes] = {
  "ok",
  "input error",
  "data error",
  "not active",
  "low statistics",
  "no measurement",
  "store error"
};

const Int_t AliTOFAnalysisTaskCalibPass0::fgkMaxNumberOfPoints = 10000; // max number of points
Double_t AliTOFAnalysisTaskCalibPass0::fgMinVertexIntegral = 100.;
Double_t AliTOFAnalysisTaskCalibPass0::fgMinDeltatIntegral = 2000.;
Double_t AliTOFAnalysisTaskCalibPass0::fgMinVertexIntegralSample = 1000.;
Double_t AliTOFAnalysisTaskCalibPass0::fgMinDeltatIntegralSample = 20000.;

//_______________________________________________________
  
AliTOFAnalysisTaskCalibPass0::AliTOFAnalysisTaskCalibPass0() :
  AliAnalysisTaskSE("TOFCalib-Pass0"),
  fStatus(kOk),
  fInitFlag(kFALSE),
  fEventSelectionFlag(kFALSE),
  fVertexSelectionFlag(kFALSE),
  fVertexCut(10.),
  fRunNumber(0),
  fESDEvent(NULL),
  fEventCuts(new AliPhysicsSelection()),
  fTrackCuts(new AliESDtrackCuts()),
  fStartTime(0),
  fEndTime(0),
  fEventTime(0),
  fElapsedTime(0),
  fkVertex(NULL),
  fGRPManager(new AliGRPManager()),
  fkGRPObject(NULL),
  fTOFcalib(new AliTOFcalib()),
  fHistoList(new TList()),
  fHistoVertexTimestamp(NULL),
  fHistoDeltatTimestamp(NULL),
  fHistoDeltazEta(NULL),
  fHistoDeltazCosTheta(NULL),
  fHistoAcceptedTracksEtaPt(NULL),
  fHistoMatchedTracksEtaPt(NULL)
{
  /* 
   * default constructor 
   */

  /* define output */
  fHistoList->SetOwner(kTRUE);
  DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliTOFAnalysisTaskCalibPass0::~AliTOFAnalysisTaskCalibPass0()
{
  /*
   * default destructor
   */

  if (fHistoList) delete fHistoList;
}

//_______________________________________________________

Bool_t
AliTOFAnalysisTaskCalibPass0::InitRun()
{
  /*
   * init run
   */

  /* get ESD event */
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESDEvent) {
    AliError("cannot get ESD event");
    return kFALSE;
  }
  /* get run number */
  Int_t runNb = fESDEvent->GetRunNumber();
  /* check run already initialized */
  if (fInitFlag && fRunNumber == runNb) return kTRUE;
  if (fInitFlag && fRunNumber != runNb) {
    AliFatal(Form("run number has changed within same job and this is not allowed: %d -> %d", fRunNumber, runNb));
    return kTRUE;
  }
  fInitFlag = kFALSE;
  /* init cdb if not yet done */
  AliCDBManager *cdb = AliCDBManager::Instance();
  if (!cdb->IsDefaultStorageSet())
    cdb->SetDefaultStorage("raw://");
  if (cdb->GetRun() < 0)
    cdb->SetRun(runNb);
    /* init TOF calib */
  if (!fTOFcalib->Init()) {
    AliError("cannot init TOF calib");
    return kFALSE;
  }
  /* get GRP data */
  if (!fGRPManager->ReadGRPEntry()) {
    AliError("error while reading \"GRPEntry\" from OCDB");
    return kFALSE;
  }
  fkGRPObject = fGRPManager->GetGRPData();
  if (!fkGRPObject) {
    AliError("cannot get \"GRPData\" from GRP manager");
    return kFALSE;
  }
  fStartTime = fkGRPObject->GetTimeStart();
  fEndTime = fkGRPObject->GetTimeEnd();
  AliInfo(Form("got \"GRPData\": startTime=%d, endTime=%d", fStartTime, fEndTime));
  
  AliInfo(Form("initialized for run %d", runNb));
  fInitFlag = kTRUE;
  fRunNumber = runNb;
  
  /* set histo title with run-number and start-time */
  fHistoVertexTimestamp->SetTitle(Form("run: %d, startTimestamp: %u", fRunNumber, fStartTime));
  fHistoDeltatTimestamp->SetTitle(Form("run: %d, startTimestamp: %u", fRunNumber, fStartTime));
  fHistoDeltazEta->SetTitle(Form("run: %d, startTimestamp: %u", fRunNumber, fStartTime));
  fHistoDeltazCosTheta->SetTitle(Form("run: %d, startTimestamp: %u", fRunNumber, fStartTime));
  fHistoAcceptedTracksEtaPt->SetTitle(Form("run: %d, startTimestamp: %u", fRunNumber, fStartTime));
  fHistoMatchedTracksEtaPt->SetTitle(Form("run: %d, startTimestamp: %u", fRunNumber, fStartTime));
  
  return kTRUE;
}

//_______________________________________________________

Bool_t
AliTOFAnalysisTaskCalibPass0::InitEvent()
{
  /*
   * init event
   */

  /* get ESD event */
  fESDEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  if (!fESDEvent) return kFALSE;
  /* set event time and elapsed time */
  fEventTime = fESDEvent->GetTimeStamp();
  fElapsedTime = fESDEvent->GetTimeStamp() - fStartTime;
  /* event selection */
  if (fEventSelectionFlag && !fEventCuts->IsCollisionCandidate(fESDEvent)) return kFALSE;
  /* vertex selection */
  fkVertex = fESDEvent->GetPrimaryVertexTracks();
  if (fVertexSelectionFlag && fkVertex->GetNContributors() < 1) {
    fkVertex = fESDEvent->GetPrimaryVertexSPD();
    if (fkVertex->GetNContributors() < 1) return kFALSE;
  }
  if (fVertexSelectionFlag && TMath::Abs(fkVertex->GetZ()) > fVertexCut) return kFALSE;
  /* calibrate ESD if requested */
  fTOFcalib->CalibrateESD(fESDEvent);

  return kTRUE;
}

//_______________________________________________________

Bool_t
AliTOFAnalysisTaskCalibPass0::HasTOFMeasurement(const AliESDtrack *track) const
{
  /*
   * has TOF measurement
   */

  /* check TOF status flags */
  if (!(track->GetStatus() & AliESDtrack::kTOFout) ||
      !(track->GetStatus() & AliESDtrack::kTIME)) return kFALSE;
  /* check integrated length */
  if (track->GetIntegratedLength() < 350.) return kFALSE;

  /* TOF measurement ok */
  return kTRUE;
}

//_______________________________________________________

void
AliTOFAnalysisTaskCalibPass0::UserCreateOutputObjects()
{
  /*
   * user create output objects
   */

  /* time binning */
  Int_t timeBins = 288;
  Float_t timeMin = 0.;
  Float_t timeMax = 24. * 3600.;
  /* vertex binning */
  Int_t vertexBins = 200;
  Float_t vertexMin = -50.;
  Float_t vertexMax = 50.;
  /* deltat binning */
  Int_t deltatBins = 2000;
  Float_t deltatMin = -24400.;
  Float_t deltatMax = 24400.;
  /* deltaz binning */
  Int_t deltazBins = 200;
  Float_t deltazMin = -10.;
  Float_t deltazMax = 10.;
  /* eta binning */
  Int_t etaBins = 200;
  Float_t etaMin = -1.;
  Float_t etaMax = 1.;
  /* p binning */
  Int_t pBins = 500;
  Float_t pMin = 0.;
  Float_t pMax = 5.;
  
  fHistoVertexTimestamp = new TH2F("hHistoVertexTimestamp", "Vertex position;elapsed time (s);z (cm);", timeBins, timeMin, timeMax, vertexBins, vertexMin, vertexMax);
  fHistoList->Add(fHistoVertexTimestamp);

  fHistoDeltatTimestamp = new TH2F("hHistoDeltatTimestamp", "Global time shift (T0-fill);elapsed time (s);t - t_{exp}^{(#pi)} (ps);", timeBins, timeMin, timeMax, deltatBins, deltatMin, deltatMax);
  fHistoList->Add(fHistoDeltatTimestamp);

  fHistoDeltazEta = new TH2F("hHistoDeltazEta", "Matching residuals (longitudinal);#eta;#Deltaz (cm);", etaBins, etaMin, etaMax, deltazBins, deltazMin, deltazMax);
  fHistoList->Add(fHistoDeltazEta);

  fHistoDeltazCosTheta = new TH2F("hHistoDeltazCosTheta", "Matching residuals (longitudinal);cos #theta;#Deltaz (cm);", etaBins, etaMin, etaMax, deltazBins, deltazMin, deltazMax);
  fHistoList->Add(fHistoDeltazCosTheta);

  fHistoAcceptedTracksEtaPt = new TH2F("hHistoAcceptedTracksEtaPt", "Accepted tracks;#eta;#p_{T} (GeV/c);", etaBins, etaMin, etaMax, pBins, pMin, pMax);
  fHistoList->Add(fHistoAcceptedTracksEtaPt);

  fHistoMatchedTracksEtaPt = new TH2F("hHistoMatchedTracksEtaPt", "Matched tracks;#eta;p_{T} (GeV/c);", etaBins, etaMin, etaMax, pBins, pMin, pMax);
  fHistoList->Add(fHistoMatchedTracksEtaPt);

  /* post data */
  PostData(1, fHistoList);
}

//_______________________________________________________

void
AliTOFAnalysisTaskCalibPass0::UserExec(Option_t *)
{
  /*
   * user exec
   */

  /* init run */
  if (!InitRun()) return;
  /* init event */
  if (!InitEvent()) return;

  /*** ACCEPTED EVENT ***/
  
  /* fill vertex histo */
  fHistoVertexTimestamp->Fill(fElapsedTime, fkVertex->GetZ());

  /* loop over ESD tracks */
  Int_t nTracks = fESDEvent->GetNumberOfTracks();
  AliESDtrack *track;
  Double_t eta, costheta, pt, time, timei[AliPID::kSPECIES], deltat, deltaz;
  for (Int_t itrk = 0; itrk < nTracks; itrk++) {
    /* get track */
    track = fESDEvent->GetTrack(itrk);
    if (!track) continue;
    /* check accept track */
    if (!fTrackCuts->AcceptTrack(track)) continue;
    /* get track info */
    eta = track->Eta();
    costheta = TMath::Cos(track->Theta());
    pt = track->Pt();
    /* fill accepted tracks histo */
    fHistoAcceptedTracksEtaPt->Fill(eta, pt);
    /* check TOF measurement */
    if (!HasTOFMeasurement(track)) continue;

    /*** ACCEPTED TRACK WITH TOF MEASUREMENT ***/

    /* fill matched tracks histo */
    fHistoMatchedTracksEtaPt->Fill(eta, pt);
    /* get TOF info */
    time = track->GetTOFsignal();
    track->GetIntegratedTimes(timei);
    deltat = time - timei[AliPID::kPion];
    deltaz = track->GetTOFsignalDz();
    
    /* fill histos */
    fHistoDeltatTimestamp->Fill(fElapsedTime, deltat);
    fHistoDeltazEta->Fill(eta, deltaz);
    fHistoDeltazCosTheta->Fill(costheta, deltaz);
    
  } /* end of loop over ESD tracks */

  /* post data */
  PostData(1, fHistoList);
  
}

//_______________________________________________________

Int_t
AliTOFAnalysisTaskCalibPass0::GetStatus()
{
  /*
   * get status
   */

  switch (fStatus) {

    /* OK, return zero */
  case kOk:
    return 0;
    break;
    
    /* non-fatal error, return negative status */
  case kNotActive:
  case kLowStatistics:
  case kNoMeasurement:
    return -fStatus;
    break;
    
    /* fatal error, return positive status */
  case kInputError: 
  case kDataError:
  case kStoreError:
    return fStatus;
    break;

    /* anything else, return negative large number */
  default:
    return -999;
    break;
  }
  
  /* should never arrive here, anyway return negative large number */
  return -999;
}

//_______________________________________________________

Bool_t
AliTOFAnalysisTaskCalibPass0::ProcessOutput(const Char_t *filename, const Char_t *dbString)
{
  /*
   * process output
   */

  Int_t ret = DoProcessOutput(filename, dbString);
  Int_t status = GetStatus();
  if (status == 0) {
    AliInfo(Form("TOF calibration successful: %s (status=%d)", fgkStatusCodeName[fStatus], status));
  }
  else if (status > 0) {
    AliInfo(Form("TOF calibration failed: %s (status=%d)", fgkStatusCodeName[fStatus], status));
  }
  else if (status < 0) {
    AliInfo(Form("TOF calibration failed (expected): %s (status=%d)", fgkStatusCodeName[fStatus], status));
  }
  
  return ret;
}

//_______________________________________________________

Bool_t
AliTOFAnalysisTaskCalibPass0::DoProcessOutput(const Char_t *filename, const Char_t *dbString)
{
  /*
   * do process output
   */

  /* reset status to OK */
  fStatus = kOk;

  /* open file */
  TFile *file = TFile::Open(filename);
  if (!file || !file->IsOpen()) {
    AliError(Form("cannot open output file %s", filename));
    fStatus = kInputError;
    return kFALSE;
  }
  /* get histograms */
  TList *list = (TList *)file->Get("TOFHistos");
  TH2F *histoVertexTimestamp = NULL;
  TH2F *histoDeltatTimestamp = NULL;
  TH2F *histoDeltazEta = NULL;
  TH2F *histoDeltazCosTheta = NULL;
  TH2F *histoAcceptedTracksEtaPt = NULL;
  TH2F *histoMatchedTracksEtaPt = NULL;
  if (list) {
    AliInfo(Form("getting histograms from \"Histos\" list from file %s", filename));
    histoVertexTimestamp = (TH2F *)list->FindObject("hHistoVertexTimestamp");
    histoDeltatTimestamp = (TH2F *)list->FindObject("hHistoDeltatTimestamp");
    histoDeltazEta = (TH2F *)list->FindObject("hHistoDeltazEta");
    histoDeltazCosTheta = (TH2F *)list->FindObject("hHistoDeltazCosTheta");
    histoAcceptedTracksEtaPt = (TH2F *)list->FindObject("hHistoAcceptedTracksEtaPt");
    histoMatchedTracksEtaPt = (TH2F *)list->FindObject("hHistoMatchedTracksEtaPt");
  }
  else {
    AliInfo(Form("getting histograms directly from file %s", filename));
    histoVertexTimestamp = (TH2F *)file->Get("hHistoVertexTimestamp");
    histoDeltatTimestamp = (TH2F *)file->Get("hHistoDeltatTimestamp");
    histoDeltazEta = (TH2F *)file->Get("hHistoDeltazEta");
    histoDeltazCosTheta = (TH2F *)file->Get("hHistoDeltazCosTheta");
    histoAcceptedTracksEtaPt = (TH2F *)file->Get("hHistoAcceptedTracksEtaPt");
    histoMatchedTracksEtaPt = (TH2F *)file->Get("hHistoMatchedTracksEtaPt");
  }
  /* check histos */ 
  if (!histoVertexTimestamp) {
    AliError(Form("cannot get \"hHistoVertexTimestamp\" object from file %s", filename));
    fStatus = kInputError;
    return kFALSE;
  }
  if (!histoDeltatTimestamp) {
    AliError(Form("cannot get \"hHistoDeltatTimestamp\" object from file %s", filename));
    fStatus = kInputError;
    return kFALSE;
  }
  if (!histoDeltazEta) {
    AliError(Form("cannot get \"hHistoDeltazEta\" object from file %s", filename));
    fStatus = kInputError;
    return kFALSE;
  }
  if (!histoDeltazCosTheta) {
    AliError(Form("cannot get \"hHistoDeltazCosTheta\" object from file %s", filename));
    fStatus = kInputError;
    return kFALSE;
  }
  if (!histoAcceptedTracksEtaPt) {
    AliError(Form("cannot get \"hHistoAccptedTracksEtaPt\" object from file %s", filename));
    fStatus = kInputError;
    return kFALSE;
  }
  if (!histoMatchedTracksEtaPt) {
    AliError(Form("cannot get \"hHistoMatchedTracksEtaPt\" object from file %s", filename));
    fStatus = kInputError;
    return kFALSE;
  }

  /* check matching performance */
  if (!CheckMatchingPerformance(histoDeltazEta, histoAcceptedTracksEtaPt, histoMatchedTracksEtaPt)) {
    AliError("error while checking matching efficiency");
    return kFALSE;
  }
  /* calibrate and store */
  if (!CalibrateAndStore(histoVertexTimestamp, histoDeltatTimestamp, dbString)) {
    AliError("error while calibrating and storing");
    return kFALSE;
  }

  /* success */
  return kTRUE;
}

//_______________________________________________________

Bool_t
AliTOFAnalysisTaskCalibPass0::CheckMatchingPerformance(const TH2F *histoDeltazCosTheta, const TH2F *histoAcceptedTracksEtaPt, const TH2F *histoMatchedTracksEtaPt) const
{
  /*
   * check matching performance
   */

  /* check pointers */
  if (!histoDeltazCosTheta || !histoAcceptedTracksEtaPt || !histoMatchedTracksEtaPt)
    return kFALSE;
  /* dummy for the time being */
  return kTRUE;
}

//_______________________________________________________

Bool_t
AliTOFAnalysisTaskCalibPass0::CalibrateAndStore(TH2F *histoVertexTimestamp, TH2F *histoDeltatTimestamp, const Char_t *dbString)
{
  /*
   * calibrate and store
   */

  /* check pointers */
  if (!histoVertexTimestamp || !histoDeltatTimestamp)
    return kFALSE;

  /*** GET RUN-NUMBER AND START-TIMESTAMP ***/

  TString str;
  TObjArray *strarr = NULL;
  TObjString *ostr = NULL;
  str = histoVertexTimestamp->GetTitle();
  strarr = str.Tokenize(",");
  if (!strarr) {
    AliError("problems whith tokenize histogram title");
    fStatus = kDataError;
    return kFALSE;
  }
  
  /* run number */
  ostr = (TObjString *)strarr->At(0);
  if (!ostr) {
    AliError("problems while getting run number from histogram title");
    fStatus = kDataError;
    return kFALSE;
  }
  str = ostr->GetString();
  if (!str.BeginsWith("run:")) {
    AliError("problems while getting run number from histogram title");
    fStatus = kDataError;
    return kFALSE;
  }
  str.Remove(0, 5);
  Int_t runNb = atoi(str.Data());
  if (runNb <= 0) {
    AliError(Form("bad run number: %d", runNb));
    fStatus = kDataError;
    return kFALSE;
  }
  AliInfo(Form("got run number: %d", runNb));

  /* start timestamp */
  ostr = (TObjString *)strarr->At(1);
  if (!ostr) {
    AliError("problems while getting start timestamp from histogram title");
    fStatus = kDataError;
    return kFALSE;
  }
  str = ostr->GetString();
  str.Remove(0, 1); /* remove empty space at the beginning */
  if (!str.BeginsWith("startTimestamp:")) {
    AliError("problems while getting start timestamp from histogram title");
    fStatus = kDataError;
    return kFALSE;
  }
  str.Remove(0, 16);
  UInt_t startTimestamp = atoi(str.Data());
  if (startTimestamp <= 0) {
    AliError(Form("bad start timestamp: %d", startTimestamp));
    fStatus = kDataError;
    return kFALSE;
  }
  TTimeStamp ts = startTimestamp;
  AliInfo(Form("got start timestamp: %d (%s)", startTimestamp, ts.AsString()));

  /*** CALIBRATION STAGE ***/

  /* get fit function */
  TF1 *fitFunc = (TF1 *)gROOT->GetFunction("gaus");

  /* projection-x */
  TH1D *histoVertexTimestamppx = histoVertexTimestamp->ProjectionX("histoVertexTimestamppx");
  TH1D *histoDeltatTimestamppx = histoDeltatTimestamp->ProjectionX("histoDeltatTimestamppx");

  /* check statistics */
  if (histoVertexTimestamppx->Integral() < fgMinVertexIntegral ||
      histoDeltatTimestamppx->Integral() < fgMinDeltatIntegral) {
    fStatus = kLowStatistics;
    return kFALSE;
  }

  /* define mix and max time bin */
  Int_t minBin = histoVertexTimestamppx->FindFirstBinAbove(0);
  Int_t maxBin = histoVertexTimestamppx->FindLastBinAbove(0);
  Float_t minTime = histoVertexTimestamppx->GetBinLowEdge(minBin);
  Float_t maxTime = histoVertexTimestamppx->GetBinLowEdge(maxBin + 1);
  AliInfo(Form("min/max time defined: %d < t < %d s [%d, %d]", (Int_t)minTime, (Int_t)maxTime, minBin, maxBin));

  /* loop over time bins */
  Int_t nPoints = 0;
  Float_t time[fgkMaxNumberOfPoints], timeerr[fgkMaxNumberOfPoints];
  Float_t vertexMean[fgkMaxNumberOfPoints], vertexMeanerr[fgkMaxNumberOfPoints];
  Float_t vertexSigma[fgkMaxNumberOfPoints], vertexSigmaerr[fgkMaxNumberOfPoints];
  Float_t timeZeroMean[fgkMaxNumberOfPoints], timeZeroMeanerr[fgkMaxNumberOfPoints];
  Float_t timeZeroSigma[fgkMaxNumberOfPoints], timeZeroSigmaerr[fgkMaxNumberOfPoints];
  Float_t averageTimeZero = 0.;
  Float_t averageTimeSigma = 0.;
  Float_t averageVertexSpread = 0.;
  for (Int_t ibin = minBin; ibin <= maxBin; ibin++) {

    /* define time window */
    Int_t startBin = ibin;
    Int_t endBin = ibin;
    while(histoVertexTimestamppx->Integral(startBin, endBin) < fgMinVertexIntegralSample ||
	  histoDeltatTimestamppx->Integral(startBin, endBin) < fgMinDeltatIntegralSample) {
      if (endBin < maxBin) endBin++;
      else if (startBin > minBin) startBin--;
      else break;
    }
    if (histoVertexTimestamppx->Integral(startBin, endBin) < fgMinVertexIntegral ||
        histoDeltatTimestamppx->Integral(startBin, endBin) < fgMinDeltatIntegral) continue;
    Float_t startTime = histoVertexTimestamppx->GetBinLowEdge(startBin);
    Float_t endTime = histoVertexTimestamppx->GetBinLowEdge(endBin + 1);
    Float_t vertexIntegral = histoVertexTimestamppx->Integral(startBin, endBin);
    Float_t deltatIntegral = histoDeltatTimestamppx->Integral(startBin, endBin);
    AliInfo(Form("time window defined: %d < t < %d s [%d, %d]: %d vertices, %d tracks", (Int_t)startTime, (Int_t)endTime, startBin, endBin, (Int_t)vertexIntegral, (Int_t)deltatIntegral));

    /* projection-y */
    TH1D *histoVertexTimestamppy = histoVertexTimestamp->ProjectionY("histoVertexTimestamppy", startBin, endBin);
    TH1D *histoDeltatTimestamppy = histoDeltatTimestamp->ProjectionY("histoDeltatTimestamppy", startBin, endBin);

    /* average time */
    histoVertexTimestamppx->GetXaxis()->SetRange(startBin, endBin);
    time[nPoints] = histoVertexTimestamppx->GetMean();
    timeerr[nPoints] = histoVertexTimestamppx->GetMeanError();

    /* fit vertex */
    if (FitPeak(fitFunc, histoVertexTimestamppy, 10., 3., 3.) != 0) {
      AliError("troubles fitting vertex, skip");
      delete histoVertexTimestamppy;
      delete histoDeltatTimestamppy;
      continue;
    }
    vertexMean[nPoints] = fitFunc->GetParameter(1);
    vertexMeanerr[nPoints] = fitFunc->GetParError(1);
    vertexSigma[nPoints] = fitFunc->GetParameter(2);
    vertexSigmaerr[nPoints] = fitFunc->GetParError(2);
    averageVertexSpread += fitFunc->GetParameter(2);

    /* fit time-zero */
    if (FitPeak(fitFunc, histoDeltatTimestamppy, 500., 2., 1.) != 0) {
      AliError("troubles fitting time-zero TRACKS, skip");
      delete histoVertexTimestamppy;
      delete histoDeltatTimestamppy;
      continue;
    }
    timeZeroMean[nPoints] = fitFunc->GetParameter(1);
    timeZeroMeanerr[nPoints] = fitFunc->GetParError(1);
    timeZeroSigma[nPoints] = fitFunc->GetParameter(2);
    timeZeroSigmaerr[nPoints] = fitFunc->GetParError(2);
    averageTimeZero += fitFunc->GetParameter(1);
    averageTimeSigma += fitFunc->GetParameter(2);

    /* delete projection-y */
    delete histoVertexTimestamppy;
    delete histoDeltatTimestamppy;

    /* increment n points */
    nPoints++;

    /* set current bin */
    ibin = endBin;

  } /* end of loop over time bins */

  /* delete projection-x */
  delete histoVertexTimestamppx;
  delete histoDeltatTimestamppx;

  /* check points */
  if (nPoints <= 0) {
    AliError("no measurement available, quit");
    fStatus = kNoMeasurement;
    return kFALSE;
  }
  AliInfo(Form("average time-zero  = %f", averageTimeZero / nPoints));
  AliInfo(Form("average time-sigma = %f", averageTimeSigma / nPoints));
  AliInfo(Form("average v-spread   = %f", averageVertexSpread / nPoints));

  /*** CREATE RUN PARAMS OBJECT ***/

#if 0
  /* get start time from GRP */
  TGrid::Connect("alien://", 0, 0, "t");
  AliCDBManager *cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("raw://");
  cdb->SetRun(runNb);
  AliGRPManager grp;
  if (!grp.ReadGRPEntry()) {
    AliError("error while reading GRP entry");
    return kFALSE;
  }
  UInt_t startTimestamp = grp.GetGRPData()->GetTimeStart();
  TTimeStamp ts;
  ts = startTimestamp;
  AliInfo(Form("got start time from GRP: %s", ts.AsString()));
#endif
    
  /* create arrays */
  UInt_t timestamp[fgkMaxNumberOfPoints];
  Float_t t0[fgkMaxNumberOfPoints];
  Float_t tofReso[fgkMaxNumberOfPoints];
  Float_t t0Spread[fgkMaxNumberOfPoints];
  for (Int_t ipoint = 0; ipoint < nPoints; ipoint++) {
    timestamp[ipoint] = (UInt_t)time[ipoint] + startTimestamp;
    t0[ipoint] = timeZeroMean[ipoint];
    tofReso[ipoint] = timeZeroSigma[ipoint];
    t0Spread[ipoint] = vertexSigma[ipoint] / 2.99792457999999984e-02;
  }
  UInt_t run[1] = {runNb};
  UInt_t runFirstPoint[1] = {0};
  UInt_t runLastPoint[1] = {nPoints - 1};
  
  /* create run params object */
  AliTOFRunParams obj(nPoints, 1);
  AliInfo(Form("create run params object for run %d with %d points", runNb, nPoints));
  obj.SetTimestamp(timestamp);
  obj.SetT0(t0);
  obj.SetTOFResolution(tofReso);
  obj.SetT0Spread(t0Spread);
  obj.SetRunNb(run);
  obj.SetRunFirstPoint(runFirstPoint);
  obj.SetRunLastPoint(runLastPoint);
  
  /*** CREATE OCDB ENTRY ***/

  if (!dbString) {
    AliError("cannot store object because of NULL string");
    fStatus = kStoreError;
    return kFALSE;
  }

  /* install run params object in OCDB */
  AliCDBManager *cdb = AliCDBManager::Instance();
  AliCDBStorage *sto = cdb->GetStorage(dbString);
  if (!sto) {
    AliError(Form("cannot get storage %s", dbString));
    fStatus = kStoreError;
    return kFALSE;
  }
  AliCDBId id("TOF/Calib/RunParams", runNb, runNb);
  AliCDBMetaData md;
  md.SetResponsible("Roberto Preghenella");
  md.SetComment("offline TOF run parameters");
  md.SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md.SetBeamPeriod(0);
  if (!sto->Put(&obj, id, &md)) {
    fStatus = kStoreError;
    AliError(Form("error while putting object in storage %s", dbString));
    return kFALSE;
  }

  /* success */
  return kTRUE;
}

//_____________________________________________________________

Int_t
AliTOFAnalysisTaskCalibPass0::FitPeak(TF1 *fitFunc, TH1D *h, Float_t startSigma, Float_t nSigmaMin, Float_t nSigmaMax)
{
  /*
   * fit peak
   */

  Double_t fitCent = h->GetBinCenter(h->GetMaximumBin());
  Double_t fitMin = fitCent - nSigmaMin * startSigma;
  Double_t fitMax = fitCent + nSigmaMax * startSigma;
  if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
  if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
  fitFunc->SetParameter(1, fitCent);
  fitFunc->SetParameter(2, startSigma);
  Int_t fitres = h->Fit(fitFunc, "WWq0", "", fitMin, fitMax);
  if (fitres != 0) return fitres;
  /* refit with better range */
  for (Int_t i = 0; i < 3; i++) {
    fitCent = fitFunc->GetParameter(1);
    fitMin = fitCent - nSigmaMin * fitFunc->GetParameter(2);
    fitMax = fitCent + nSigmaMax * fitFunc->GetParameter(2);
    if (fitMin < h->GetXaxis()->GetXmin()) fitMin = h->GetXaxis()->GetXmin();
    if (fitMax > h->GetXaxis()->GetXmax()) fitMax = h->GetXaxis()->GetXmax();
    fitres = h->Fit(fitFunc, "q0", "", fitMin, fitMax);
    if (fitres != 0) return fitres;
  }
  return fitres;
}
