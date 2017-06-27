#ifndef ALITOFANALYSISTASKCALIBTREE_H
#define ALITOFANALYSISTASKCALIBTREE_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TMath.h"

class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrackCuts;
class AliESDtrack;
class AliESDVertex;
class AliESDpid;

class AliPhysicsSelection;

class AliTOFcalibHisto;
class AliTOFcalib;
class AliTOFT0maker;
class AliTOFT0v1;

class AliGRPManager;
class AliGRPObject;

class TH2F;

#define MAXHITS ((Int_t)100000)

class AliTOFAnalysisTaskCalibTree :
public AliAnalysisTaskSE
{

 public:

  AliTOFAnalysisTaskCalibTree(const Char_t* name = "TOFcalibTree"); // default constructor
  virtual ~AliTOFAnalysisTaskCalibTree(); // default destructor

  virtual void UserCreateOutputObjects(); // user create output objects
  virtual void UserExec(Option_t *); // user exec

  // getters
  AliPhysicsSelection *GetEventCuts() const {return fEventCuts;}; // getter
  AliESDtrackCuts *GetTrackCuts() const {return fTrackCuts;}; // getter
  AliTOFcalib *GetTOFcalib() const {return fTOFcalib;}; // getter
  AliTOFT0maker *GetTOFT0maker() const {return fTOFT0maker;}; // getter
  AliTOFT0v1 *GetTOFT0v1() const {return fTOFT0v1;}; // getter
  Float_t GetTimeResolution() const {return fTimeResolution;}; // time resolution
  AliESDpid *GetESDpid() const {return fESDpid;}; // get ESD PID

  // setters
  void SetEventSelectionFlag(Bool_t value = kTRUE) {fEventSelectionFlag = value;}; // setter
  void SetVertexSelectionFlag(Bool_t value = kTRUE) {fVertexSelectionFlag = value;}; // setter
  void SetVertexCut(Double_t value) {fVertexCut = value;}; // setter
  void SetDiscardPileupEventFlag(Bool_t value = kTRUE) {fDiscardPileupEventFlag = value;}; // setter
  void SetPrimaryDCASelectionFlag(Bool_t value = kTRUE) {fPrimaryDCASelectionFlag = value;}; // setter
  void SetTimeResolution(Float_t value) {fTimeResolution = value;}; // set time resolution
  void SetCalibrateTOFsignal(Bool_t value = kTRUE) {fCalibrateTOFsignal = value;}; // setter
  void SetComputeT0TOF(Bool_t value = kTRUE) {fComputeT0TOF = value;}; // setter
  void SetUseT0TOF(Bool_t value = kTRUE) {fUseT0TOF = value;}; // setter
  void SetUseLHCClockPhase(Bool_t value = kTRUE) {fUseLHCClockPhase = value;}; // setter
  void SetSpecificStorageParOffline(Char_t *value) {fSpecificStorageParOffline = value;}; // set specific storage ParOffline
  void SetSpecificStorageRunParams(Char_t *value) {fSpecificStorageRunParams = value;}; // set specific storage RunParams
  void SetSpecificStorageFineSlewing(Char_t *value) {fSpecificStorageFineSlewing = value;}; // set specific storage FineSlewing
  void SetSaveCoordinates(Bool_t value = kTRUE) {fSaveCoordinates = value;}; // set flag to save hit coordinates in tree

 protected:

  AliTOFAnalysisTaskCalibTree(const AliTOFAnalysisTaskCalibTree &); // copy constructor
  AliTOFAnalysisTaskCalibTree &operator=(const AliTOFAnalysisTaskCalibTree &); // operator=

  // methods:

  Bool_t InitRun(); // init run
  Bool_t InitEvent(); // init event
  Bool_t HasTOFMeasurement(AliESDtrack *track); // has TOF
  Bool_t HasPrimaryDCA(AliESDtrack *track); // has primary DCA

  // members:

  // flags and cuts
  Bool_t fInitFlag;                   // init flag
  Bool_t fEventSelectionFlag;         // event selection flag
  Bool_t fVertexSelectionFlag;        // vertex selection flag
  Double_t fVertexCut;                // vertex cut
  Bool_t fDiscardPileupEventFlag;     // discard pile-up event flag
  Bool_t fCalibrateTOFsignal;         // calibrate TOF signal
  Bool_t fComputeT0TOF;               // compute T0-TOF
  Bool_t fUseT0TOF;                   // use T0-TOF
  Bool_t fUseLHCClockPhase;           // use LHC-clock phase
  Bool_t fPrimaryDCASelectionFlag;    // selection of primaries with DCA cut

  // ESD related stuff
  Int_t fRunNumber;                   // run number
  AliESDEvent *fESDEvent;             //!<! ESD event
  AliPhysicsSelection *fEventCuts;    //!<! event cuts
  AliESDtrackCuts *fTrackCuts;        //!<! track cuts
  AliESDpid *fESDpid;                 //!<! ESD PID
  UInt_t fStartTime;                  // start time
  UInt_t fEndTime;                    // end time
  UInt_t fElapsedTime;                // event time since start
  Bool_t fIsCollisionCandidate;       // is collision candidate
  Bool_t fHasVertex;                  // has vertex
  const AliESDVertex *fVertex;        //!<! vertex

  // GRP related stuff
  AliGRPManager *fGRPManager;         //!<! GRP manager
  const AliGRPObject *fGRPObject;     //!<! GRP object

  // TOF related stuff
  TString fSpecificStorageParOffline; // specific storage ParOffline
  TString fSpecificStorageRunParams;  // specific storage RunParams
  TString fSpecificStorageFineSlewing; // specific storage FineSlewing
  Float_t fTimeResolution;            // time resolution
  AliTOFcalib *fTOFcalib;             //!<! TOF calib
  AliTOFT0maker *fTOFT0maker;         //!<! TOF-T0 maker
  AliTOFT0v1 *fTOFT0v1;               //!<! TOF-T0 v1

  // task related stuff
  Int_t fMaxHits;       //array parameter
  UInt_t ftimestamp;
  Float_t fVertexZ;
  Float_t ftimezero;
  Int_t fnhits;
  Float_t* fmomentum;       //[fMaxHits] momentum
  Float_t* flength;         //[fMaxHits] length
  Int_t* findex;            //[fMaxHits] index
  Float_t* ftime;           //[fMaxHits] time
  Float_t* ftot;            //[fMaxHits] time over threshold
  Float_t* ftexp;           //[fMaxHits] texp
  Float_t* fDeltax;         //[fMaxHits] delta-x
  Float_t* fDeltaz;         //[fMaxHits] delta-z
  Float_t* fDeltat;         //[fMaxHits] delta-t
  Float_t* fDeltaraw;       //[fMaxHits] delta-raw
  Bool_t fSaveCoordinates;

  TTree* fOutputTree;                 //!<! output tree

  ClassDef(AliTOFAnalysisTaskCalibTree, 4);
};

#endif /* ALIANALYSISTASKTOFCOMPACTCALIB_H */
