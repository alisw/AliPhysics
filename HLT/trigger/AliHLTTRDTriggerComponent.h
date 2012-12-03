//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRDTRIGGERCOMPONENT_H
#define ALIHLTTRDTRIGGERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTTRDTriggerComponent.h
/// @author Felix Rettig, Stefan Kirsch
/// @date   2012-08-16
/// @brief

#include "TString.h"
#include "AliHLTTrigger.h"
#include "AliESDEvent.h"
#include "AliESDTrdTrack.h"
#include "AliExternalTrackParam.h"
#include "AliTRDgeometry.h"
#ifdef __TRDHLTDEBUG
  #include "AliTRDtrackingEventDisplay.h"
#endif

#define trd_det_lsi(det) ((det) / 6)               // convert TRD detector/chamber 0-539 index to linear stack index 0-89
#define trd_det_lyr(det) ((det) % 6)               // convert detector (=chamber) number 0-539 to local layer 0-5
#define trd_det_si(det) (((det) % 30) / 6)         // convert detector (=chamber) number 0-539 to local stack index 0-4

class TObjArray;
class TH1I;
class TH2I;
class AliHLTComponentBenchmark;
class AliESDtrack;
class AliHLTGlobalBarrelTrack;
class AliESDtrackCuts;
class AliHLTESDTrackCuts;
class AliTRDonlineTrackingDataContainer;

/**
 * @class  AliHLTTRDTriggerComponent
 */
class AliHLTTRDTriggerComponent : public AliHLTTrigger
{
 public:

  AliHLTTRDTriggerComponent();
  virtual ~AliHLTTRDTriggerComponent();

  virtual const char* GetTriggerName() const;
  virtual AliHLTComponent* Spawn();

 protected:
  int DoInit(int argc, const char** argv);
  int DoDeinit();
  int Reconfigure(const char* cdbEntry, const char* chainId);
  int ReadPreprocessorValues(const char* modules);
  int ConfigureFromCDBObject(TString cdbPath);
  int ScanConfigurationArgument(int argc, const char** argv);

 private:
  AliHLTTRDTriggerComponent (const AliHLTTRDTriggerComponent&);
  AliHLTTRDTriggerComponent& operator=(const AliHLTTRDTriggerComponent&);
  virtual int DoTrigger();

  Bool_t CheckRefTrackCuts(AliESDtrack* track);

  // TRD-trigger specific
  void ScanTriggerClasses(const char* firedTriggerClasses);
  int PrepareESDData();
  int PrepareHLTData();
  int PrepareTRDData();
  int MatchTRDTracks();
  int MatchTRDTracksESD();
  int MatchTRDTracksHLT();
  Bool_t TRDElectronTrigger(const char *ident, const Double_t minPt, const UShort_t minPID);

  Bool_t TrackPlaneIntersect(AliExternalTrackParam *trk, Double_t pnt[3], Double_t norm[3], Double_t mag);

  Int_t EstimateTrackDistance(AliExternalTrackParam *refParam,
				 const UShort_t stack,
				 const UShort_t layerMask,
				 const Float_t trklLocalY[6], const Int_t trklBinZ[6],
				 Double_t mag, Double_t *ydist, Double_t *zdist);

  Int_t EstimateTrackDistance(AliESDtrack *esd_track,
				 const UShort_t stack,
				 const UShort_t layerMask,
				 const Float_t trklLocalY[6], const Int_t trklBinZ[6],
				 Double_t mag, Double_t *ydist, Double_t *zdist);

  Double_t RateTrackMatch(Double_t distY, Double_t distZ, Double_t rpt, Double_t gpt);

  void DumpTrackingData();
  void AssignTrackInfo(TString* infoStr, const UInt_t stack, const UInt_t trackIndex, const char* flagStr = "");
#ifdef __TRDHLTDEBUG
  void RenderEvent(const Bool_t showGtuTracks = kTRUE, const Bool_t showTracklets = kTRUE, const Bool_t showRefTracks = kTRUE);
#endif
  void DbgLog(const char* prefix, ...);

  TString  fName;                           //! trigger name
  Double_t fRefTrackSelectionEtaLimit;      //! ref track preselection maximum eta
  Double_t fRefTrackSelectionVertexXYLimit; //! ref track preselection maximum distance to ip in XY plane
  Double_t fRefTrackSelectionVertexZLimit;  //! ref track preselection maximum distance to ip in Z
  Double_t fRefTrackSelectionPtThreshold;   //! pt threshold for ref track preselection in GeV/c
  Double_t fMatchRatingThreshold;           //! track match rating threshold
  Double_t fElectronTriggerPtThresholdHSE;  //! pt threshold for HSE electron trigger
  UShort_t fElectronTriggerPIDThresholdHSE; //! PID threshold for HSE electron trigger
  Double_t fElectronTriggerPtThresholdHQU;  //! pt threshold for HQU electron trigger
  UShort_t fElectronTriggerPIDThresholdHQU; //! PID threshold for HQU electron trigger
  Bool_t fApplyRefTrackCuts;                //! switch on/off ref track cuts for matching
  Bool_t fElectronTriggerOnL1TrgOnly;       //! run electron trigger only for events with L1 electron trigger
  UShort_t fHistoMode;                      //! histogramming mode, 0: single event, 1: accumulative (debugging)
  UShort_t fDebugLevel;                     //! set debug checks/output level, 0: debug off
  Bool_t fExtendedHistos;                   //! switch on/off additional histograms
  Bool_t fEventRendering;                   //! switch on/off event rendering
  Bool_t fPushHistos;                       //! switch on/off pushing of histograms event by event
  Bool_t fWriteHistos;                      //! switch on/off histogram writing on deinit

  static const char* fgkDefaultOCDBEntry;                    //! default OCDB entry
  static const char* fgkTriggerDecisionElectronHSE;          //! electron trigger flag string
  static const char* fgkTriggerDecisionElectronHQU;          //! electron trigger flag string

  static const AliHLTEventID_t fgkInvalidEventId = 0xffffffffffffffffllu;
  static const unsigned int fkTRDLayers = 6;                 //! number of layers per stack in TRD
  static const unsigned int fkTRDStacks = 90;                //! number of stacks in TRD
  static const unsigned int fkTRDSectors = 18;               //! number of sectors in TRD
  static const unsigned int fkMaxRefTracksPerStack = 25000;  //! maximum number of ref tracks per stack
  static const unsigned int fkMaxRefTracks = 25000;          //! maximum number of ref tracks

  static const unsigned int fkElectronTriggerHSE = 0x1;      //! HSE electron trigger flag
  static const unsigned int fkElectronTriggerHQU = 0x2;      //! HSE electron trigger flag

  AliHLTEventID_t fEventId;                                  //! hlt internal event id
  Int_t fRunNumber;                                          //! run number
  TString* fChunkId;                                         //! chunk identifier
  UInt_t fSectorsWithData;                                   //! data present flags for each sector
  Bool_t fIsMinBiasEvent;                                    //! indicates a minimum bias event
  Bool_t fIsTRDElectronEvent;                                //! indicates a TRD L1 electron triggered event
  Bool_t fESDtracksPresent;                                  //! indicates that ESD tracks are present
  Bool_t fHLTtracksPresent;                                  //! indicates that HLT raw tracks are present

  AliTRDgeometry* fTRDGeometry;                              //! instance of TRD geometry
  AliESDEvent* fEsdEvent;                                    //! current ESD event
  AliTRDonlineTrackingDataContainer* fTrackingData;          //! container for TRD tracking data
  vector<AliHLTGlobalBarrelTrack>* fHLTTracks;               //! HLT raw tracks
  AliESDtrackCuts* fRefTrackCuts;                            //! reference track cuts

#ifdef __TRDHLTDEBUG
  AliTRDtrackingEventDisplay* fEventDisplay;                 //! event rendering
  AliHLTComponentBenchmark* fBenchmark;                      //! benchmark instance
#endif

  TObjArray* fHistArray;
  TH1I* fHistMatchRating;                                    //! histo
  TH2I* fHistMatchRatingByPt;                                //! histo
  TH2I* fHistMatchRatingByPid;                               //! histo
  TH1I* fHistTrackPt;                                        //! histo
  TH1I* fHistTrackPtMatched;                                 //! histo
  TH2I* fHistTrackPtCorr;                                    //! histo
  TH1I* fHistTrackPid;                                       //! histo
  TH1I* fHistTrackPidMatched;                                //! histo
  TH1I* fHistElectronCandidatePt;                            //! histo
  TH1I* fHistElectronCandidateMatchedPt;                     //! histo
  TH1I* fHistElectronCandidatePid;                           //! histo
  TH1I* fHistElectronCandidateMatchedPid;                    //! histo
  TH2I* fHistRefTrackPid;                                    //! histo
  TH2I* fHistMatchedRefTrackPid;                             //! histo
  TH2I* fHistPIDvsTruncPID;                                  //! histo
  TH2I* fHistElectronFalsePIDvsTruncPID;                     //! histo
  TH2I* fHistElectronConfirmedPIDvsTruncPID;                 //! histo
  TH2I* fHistTrackMatchingCombinations;                      //! histo
  TH1I* fHistElectronTriggerBaseMinBias;                     //! histo
  TH1I* fHistElectronTriggerBaseTrdL1;                       //! histo

  ClassDef(AliHLTTRDTriggerComponent, 0)
};

#endif //ALIHLTTRDTRIGGERCOMPONENT_H
