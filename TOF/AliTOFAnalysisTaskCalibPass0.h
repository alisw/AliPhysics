#ifndef ALITOFANALYSISTASKCALIBPASS0_H
#define ALITOFANALYSISTASKCALIBPASS0_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides TOF pass0/passX calibration tools   //
//                                                           //
///////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"

class AliESDtrack;
class AliESDEvent;
class AliPhysicsSelection;
class AliESDtrackCuts;
class AliESDVertex;
class AliGRPManager;
class AliGRPObject;
class AliTOFcalib;
class TList;
class TH2F;
class TF1;
class TH1D;

class AliTOFAnalysisTaskCalibPass0 :
public AliAnalysisTaskSE
{

 public:

  AliTOFAnalysisTaskCalibPass0(); // default constructor
  virtual ~AliTOFAnalysisTaskCalibPass0(); // default destructor

  virtual void UserCreateOutputObjects(); // user create output objects
  virtual void UserExec(Option_t *); // user exec

  /* getters */
  AliPhysicsSelection *GetEventCuts() const {return fEventCuts;}; // getter
  AliESDtrackCuts *GetTrackCuts() const {return fTrackCuts;}; // getter
  AliTOFcalib *GetTOFcalib() const {return fTOFcalib;}; // getter
  
  /* setters */
  void SetEventSelectionFlag(Bool_t value = kTRUE) {fEventSelectionFlag = value;}; // setter
  void SetVertexSelectionFlag(Bool_t value = kTRUE) {fVertexSelectionFlag = value;}; // setter
  void SetVertexCut(Double_t value) {fVertexCut = value;}; // setter

  /* post-processing methods */
  Bool_t ProcessOutput(const Char_t *filename, const Char_t *dbString); // process output
  Bool_t DoProcessOutput(const Char_t *filename, const Char_t *dbString); // process output
  Int_t GetStatus(); // get status

  /* static setters */
  static void SetMinVertexIntegral(Double_t value) {fgMinVertexIntegral = value;}; // setter
  static void SetMinDeltatIntegal(Double_t value) {fgMinDeltatIntegral = value;}; // setter
  static void SetMinVertexIntegralSample (Double_t value) {fgMinVertexIntegralSample = value;}; // setter
  static void SetMinDeltatIntegralSample (Double_t value) {fgMinDeltatIntegralSample = value;}; // setter
  
 protected:

  AliTOFAnalysisTaskCalibPass0(const AliTOFAnalysisTaskCalibPass0 &); // copy constructor
  AliTOFAnalysisTaskCalibPass0 &operator=(const AliTOFAnalysisTaskCalibPass0 &); // operator=

  /* status codes */
  enum EStatusCode_t {
    kOk,
    kInputError, /* open file error, missing histos */
    kDataError, /* problems with histo information */
    kNotActive, /* not active in data taking and/or reconstruction */
    kLowStatistics, /* too low statistics */
    kNoMeasurement, /* no measurement performed */
    kStoreError, /* problems storing OCDB */
    kNStatusCodes
  };
  Int_t fStatus; /* status code */
  static const Char_t *fgkStatusCodeName[kNStatusCodes];

  /* methods */
  Bool_t InitRun(); // init run
  Bool_t InitEvent(); // init event
  Bool_t HasTOFMeasurement(const AliESDtrack *track) const ; // has TOF measurement

  /* post-processing methods */
  Bool_t CheckMatchingPerformance(const TH2F *histoDeltazEta, const TH2F *histoAcceptedTracksEtaPt, const TH2F *histoMatchedTracksEtaPt) const; // check matching efficiency
  Bool_t CalibrateAndStore(TH2F *histoVertexTimestamp, TH2F *histoDeltatTimestamp, const Char_t *dbString); // calibrate and store
  Int_t FitPeak(TF1 *fitFunc, TH1D *h, Float_t startSigma, Float_t nSigmaMin, Float_t nSigmaMax); // fit peak

  /* flags and cuts */
  Bool_t fInitFlag; // init flag
  Bool_t fEventSelectionFlag; // event selection flag
  Bool_t fVertexSelectionFlag; // vertex selection flag
  Double_t fVertexCut; // vertex cut

  /* ESD analysis */
  Int_t fRunNumber; // run number
  AliESDEvent *fESDEvent; // ESD event
  AliPhysicsSelection *fEventCuts; // event cuts
  AliESDtrackCuts *fTrackCuts; // track cuts
  UInt_t fStartTime; // start time
  UInt_t fEndTime; // end time
  UInt_t fEventTime; // event time
  UInt_t fElapsedTime; // event time since start
  const AliESDVertex *fkVertex; // vertex

  /* GRP related stuff */
  AliGRPManager *fGRPManager; // GRP manager
  const AliGRPObject *fkGRPObject; // GRP object
  
  /* TOF related stuff */
  AliTOFcalib *fTOFcalib; // TOF calib

  /* lists and histos */
  TList *fHistoList; // list of histograms
  TH2F *fHistoVertexTimestamp; // vertex-timestamp histo
  TH2F *fHistoDeltatTimestamp; // deltat-timestamp histo
  TH2F *fHistoDeltazEta; // deltaz-eta histo
  TH2F *fHistoDeltazCosTheta; // deltaz-costheta histo
  TH2F *fHistoAcceptedTracksEtaPt; // accepted tracks eta-pt histo
  TH2F *fHistoMatchedTracksEtaPt; // matched tracks eta-pt histo

  /* post-processing variables */
  static const Int_t fgkMaxNumberOfPoints; // max number of points
  static Double_t fgMinVertexIntegral; // min vertex integral
  static Double_t fgMinDeltatIntegral; // min vertex integral
  static Double_t fgMinVertexIntegralSample; // min vertex integral sample
  static Double_t fgMinDeltatIntegralSample; // min vertex integral sample


  ClassDef(AliTOFAnalysisTaskCalibPass0, 1);
};

#endif /* ALIANALYSISTASKEVENTTIME_H */
