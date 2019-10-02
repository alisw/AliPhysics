#ifndef ALITRDDIGITSTASK_H
#define ALITRDDIGITSTASK_H

// example and base class for an analysis task to analyse TRD digits
// Authors: Tom Dietel
// based on the AliAnalysisTaskPt

// This class provides basic functionality to analyse TRD digits in
// the context of an ESD event:
//
//  - Read and write TRD.Digits.root files
//  - Provide e,pi reference samples via AliESDv0KineCuts
//  - Match ESD tracks to TRD digits.
//
// The class requires the presence of an AliESDs.root and a
// TRD.Digits.root file. AODs have not been tested, because matching
// of tracks to TRD digits seemed easier with some ESD-only features.
//
// The TRD.Digits.root files are opened in the UserNotify() method,
// i.e. every time that a new ESD (or AOD) file is opened. The file
// contains digits trees, numbered according to the event number in
// the corresponding ESD file. Missing digits are allowed, to support
// filtered digits trees with only a subset of events.
//
// The class has basic support to identify electrons and pions (and in
// principle also protons) via V0 decays. It might need an option to
// tighten/tune the default cuts of the AliESDv0KineCuts class.
//
// The methods FindDigits and FindDigitsTrkl can be used to find the
// digits for a given track without the need to go via offline
// tracklets.

class AliTRDdigitsManager;
//class AliESDv0KineCuts;
class AliESDEvent;
class AliESDfriendTrack;
class AliExternalTrackParam;
class AliTRDtrackV1;
class AliTRDgeometry;

#include "AliESDv0KineCuts.h"
#include "AliAnalysisTaskSE.h"

#include <vector>

class AliTRDdigitsTask : public AliAnalysisTaskSE {
public:

  AliTRDdigitsTask(const char *name);
  virtual ~AliTRDdigitsTask() {}

  virtual void   UserCreateOutputObjects();
  virtual Bool_t UserNotify();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  // set file names (relative to directory of AliESDs.root) of
  // files for reading/writing of TRD digits
  void SetDigitsInputFilename(TString x) {fDigitsInputFileName=x;}
  void SetDigitsOutputFilename(TString x) {fDigitsOutputFileName=x;}

  // Set V0 kine cuts for generation of e/pi/... reference samples.
  // Setting this to NULL disables reference samples.
  void SetV0KineCuts(AliESDv0KineCuts *c) { fV0cuts = c; }
  AliESDv0KineCuts* GetV0KineCuts() {return fV0cuts;}

  typedef enum {
    kPidUndef = 0,
    kPidError = 1,
    kPidV0Electron = 2,
    kPidV0Pion = 3,
    kPidV0Proton = 4
  } EPID_t;

protected:

  //--------------------------------------------------------------------
  // Event bookkeeping
  Bool_t NextEvent(Bool_t preload=kFALSE);


  //--------------------------------------------------------------------
  // Interface for analysis functionality
  virtual void AnalyseEvent();

  AliESDEvent *fESDevent;     //! ESD event object
  TList       *fOutputList;   //! Output list

  //--------------------------------------------------------------------
  // Some histograms of general interest, along with functions to create
  // them

  void CreateV0Plots();
  TH2F*       fhArmenteros;                 //! 2D V0 QA Hist

  void CreateTriggerHistos();
  void FillTriggerHisto(TH1*);
  TH1F* fhTrgAll;  //! trigger classes of all events
  TH1F* fhTrgAcc;  //! trigger classes of accepted events

  //--------------------------------------------------------------------
  // TRD digits I/O
  Bool_t ReadDigits();
  Bool_t WriteDigits();

  AliTRDdigitsManager* fDigMan;     //! digits manager

  //--------------------------------------------------------------------
  // e,pi reference sample generation - to be moved to AliTRDdigitsTask
  AliESDv0KineCuts *fV0cuts; //  V0 cuts

  std::vector<EPID_t> fPidTags;  // tags for identified particles
  void FillV0PIDlist();


  //--------------------------------------------------------------------
  // Matching of ESD tracks to TRD digits
  AliTRDtrackV1* FindTRDtrackV1(AliESDfriendTrack* friendtrack);

  Int_t FindDigitsTrkl(const AliTRDtrackV1* trdTrack, Int_t layer,
		       Int_t* det, Int_t* row, Int_t* col,
		       Float_t* x, Float_t* y, Float_t* z);

  Int_t FindDigits(const AliExternalTrackParam* param,
		   Float_t bfield, Int_t layer,
		   Int_t* det, Int_t* row, Int_t* col);



private:

  TFile* OpenDigitsFile(TString inputfile, TString digfile, TString opt);

  TString fDigitsInputFileName;     //  Name of digits file for reading
  TString fDigitsOutputFileName;    //  Name of digits file for writing

  TFile*  fDigitsInputFile;         //! Digits file for reading
  TFile*  fDigitsOutputFile;        //! Digits file for writing

  AliTRDgeometry* fGeo;             //! TRD geometry

  Int_t fEventNoInFile;             //! event number in file for digits assoc
  Int_t fDigitsLoadedFlag;          //! flag to indicate if digits are avail


  AliTRDdigitsTask(const AliTRDdigitsTask&); // not implemented
  AliTRDdigitsTask& operator=(const AliTRDdigitsTask&); // not implemented

  ClassDef(AliTRDdigitsTask, 1);
};

#endif
