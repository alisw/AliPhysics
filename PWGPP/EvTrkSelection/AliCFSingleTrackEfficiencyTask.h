#ifndef AliCFSINGLETRACKEFFICIENCYTASK_H
#define AliCFSINGLETRACKEFFICIENCYTASK_H

#include "AliAnalysisTaskSE.h"
#include "AliSingleTrackEffCuts.h"

class TH1I;
class TParticle;
class AliStack;
class AliCFManager;
class AliAODTrack;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDVertex;
class AliVVertex;
class AliVParticle;


class AliCFSingleTrackEfficiencyTask : public AliAnalysisTaskSE {

 public:

  enum {
    kStepMCGenCut         = 0, // selected generated particles, event selection
    kStepMCKineCut        = 1, // generated particles after acceptance cuts
    kStepMCAccpCut        = 2, // particles with a minimum number of clusters in detector (ESD only)
    kStepReconstructed    = 3, // reconstructed tracks
    kStepRecoKineCuts     = 4, // reconstructed tracks after acceptance
    kStepReconstructedMC  = 5, // tracks passing ESD+MC trackCuts (kine properties)
    kStepRecoQualityCuts  = 6, // tracks passing ESD+MC trackCuts (reco properties)
    kStepRecoPID          = 7  // tracks after PID criteria
  };

  AliCFSingleTrackEfficiencyTask();
  AliCFSingleTrackEfficiencyTask(const Char_t* name, AliESDtrackCuts *trackcuts, AliSingleTrackEffCuts * mccuts);
  AliCFSingleTrackEfficiencyTask(const AliCFSingleTrackEfficiencyTask& c);
  AliCFSingleTrackEfficiencyTask& operator= (const AliCFSingleTrackEfficiencyTask& c);
  virtual ~AliCFSingleTrackEfficiencyTask();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  virtual void Init();
  virtual void LocalInit() { Init(); }

  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) { fCFManager = io; }   // global correction manager
  AliCFManager * GetCFManager() const { return fCFManager; }           // get corr manager
  void           SetQAList(TList* list) { fQAHistList = list; }        // set CF QA list (empty, not used, but needed)

  // Setters
  // analyze AOD:1 or ESD:0 data
  void SetReadAODData (Bool_t flag=kTRUE) { fReadAODData=flag; }
  // select track filter bit:1 or not:0 and set the bit
  void SetFilterBit (Bool_t flag=kTRUE) { fSetFilterBit=flag; }
  void SetFilterType (Int_t fbittype) { fbit=fbittype; }
  // select trigger event mask
  void SetTriggerMask(ULong64_t mask=0) { fTriggerMask=mask; }
  // set whether to evaluate centrality
  void SetUseCentrality(Bool_t flag, TString estimator=""){ fEvalCentrality=flag; fCentralityEstimator=estimator; }

  // Getters
  // analyze AOD:1 or ESD:0 data
  Bool_t IsReadAODData()   const { return fReadAODData; }
  // select trigger event mask
  ULong64_t GetTriggerMask() { return fTriggerMask; }
  // reconstructed track cuts
  AliESDtrackCuts *GetTrackCuts(){ return (AliESDtrackCuts*)fTrackCuts; } 
  // particle and event cuts 
  AliSingleTrackEffCuts *GetSingleTrackEffCuts() { return (AliSingleTrackEffCuts*)fMCCuts; }


 protected:

  // Check ESD generated particles
  void CheckESDMCParticles();
  // Check AOD generated particles
  void CheckAODMCParticles();
  // Check reconstructed particles
  void CheckReconstructedParticles();
  // Convert AOD track to an ESD track
  AliESDtrack* ConvertTrack(AliAODTrack *track);
  // Count the number of tracklets in given eta range
  Int_t GetNumberOfTrackletsInEtaRange(Double_t mineta, Double_t maxeta);
  // Evaluate the event centrality
  Double_t GetCentrality();

  Bool_t fReadAODData;       // flag for AOD/ESD input files
  AliCFManager *fCFManager;  // pointer to the CF manager slot 2
  TList *fQAHistList;        // list of QA histograms slot 3

  AliESDtrackCuts *fTrackCuts;    // track cuts (reconstructed level) slot 4
  AliSingleTrackEffCuts *fMCCuts; // particle cuts used slot 5
  ULong64_t fTriggerMask;         // event selection trigger mask

  Bool_t fSetFilterBit; // flag to decide if applying filter-bit selection to tracks
  Int_t  fbit;          // filter-bit selection to tracks

  Bool_t fEvalCentrality;        // flag to enable centrality determination
  TString fCentralityEstimator;  // centrality estimator

  TH1I  *fHistEventsProcessed;   //! histo for monitoring the number of events processed slot 1

  ClassDef(AliCFSingleTrackEfficiencyTask,2)
};

#endif
