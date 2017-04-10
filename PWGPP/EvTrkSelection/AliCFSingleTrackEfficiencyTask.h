#ifndef AliCFSINGLETRACKEFFICIENCYTASK_H
#define AliCFSINGLETRACKEFFICIENCYTASK_H

#include "AliAnalysisTaskSE.h"
#include "AliSingleTrackEffCuts.h"

class TH1I;
class TParticle;
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
    kStepMCGenCut         = 0, // selected generated particles, event selection  (kSlow only)
    kStepMCKineCut        = 1, // generated particles after acceptance cuts
    kStepMCAccpCut        = 2, // particles with a minimum number of clusters in detector (ESD only, kSlow only)
    kStepReconstructed    = 3, // reconstructed tracks (kSlow only)
    kStepRecoKineCuts     = 4, // reconstructed tracks after acceptance (kSlow only)
    kStepReconstructedMC  = 5, // tracks passing ESD+MC trackCuts (kine properties, kSlow only)
    kStepRecoQualityCuts  = 6, // tracks passing ESD+MC trackCuts (reco properties)
    kStepRecoPID          = 7  // tracks after PID criteria (kSlow only)
  };

  enum{
     kSlow = 0, // fast configuration, only a subset of variables
     kFast = 1  // slow configuration, all variables
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
  // setter for removal of fake tracks from the calculaltion (negative label at reconstructed level)
  void SetRemoveNegativeLabelTracks(Bool_t flag) { fRemoveNegativeLabelTracks=flag; }
  // flag to match or skip the matching of reconstructed to kinematic tracks
  void SetMatchToKinematicTrack(Bool_t flag) { fMatchToKinematicTrack=flag; }
  void SetUseGeneratedKine(Bool_t flag) {fUseGeneratedKine=flag;}

  // select trigger event mask
  void SetTriggerMask(ULong64_t mask=0) { fTriggerMask=mask; }
  // set whether to evaluate centrality
  void SetUseCentrality(Bool_t flag, TString estimator=""){ fEvalCentrality=flag; fCentralityEstimator=estimator; }
    
  void SetConfiguration(Int_t configuration) {
      (configuration == kSlow) ? Printf("Slow configuration chosen, all variables will be used!") : Printf("Fast configuration chosen, not all steps will be filled!");
      fConfiguration = configuration;
  }
  Int_t GetConfiguration() const { return fConfiguration; }
    
  // Getters
  // analyze AOD:1 or ESD:0 data
  Bool_t IsReadAODData()   const { return fReadAODData; }
  // select trigger event mask
  ULong64_t GetTriggerMask() { return fTriggerMask; }
  // reconstructed track cuts
  AliESDtrackCuts *GetTrackCuts(){ return (AliESDtrackCuts*)fTrackCuts; } 
  // particle and event cuts 
  AliSingleTrackEffCuts *GetSingleTrackEffCuts() { return (AliSingleTrackEffCuts*)fMCCuts; }
  // fake tracks removal flag
  Bool_t GetRemoveNegativeLabelTracks() { return fRemoveNegativeLabelTracks; }
  // flag to match or skip the matching of reconstructed to kinematic tracks
  Bool_t GetMatchToKinematicTrack() { return fMatchToKinematicTrack; }



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
  Double_t GetCentralityOldFramework();

  Bool_t fReadAODData;       // flag for AOD/ESD input files
  AliCFManager *fCFManager;  // pointer to the CF manager slot 2
  TList *fQAHistList;        // list of QA histograms slot 3

  AliESDtrackCuts *fTrackCuts;    // track cuts (reconstructed level) slot 4
  AliSingleTrackEffCuts *fMCCuts; // particle cuts used slot 5
  ULong64_t fTriggerMask;         // event selection trigger mask

  Bool_t fSetFilterBit; // flag to decide if applying filter-bit selection to tracks
  Int_t  fbit;          // filter-bit selection to tracks
  Bool_t fRemoveNegativeLabelTracks; // flag to remove fake tracks (reconstructed tracks with negative label)
  Bool_t fMatchToKinematicTrack;  // flag to check if the reconstructed track matches to a kinematic one in the good acceptance
  Bool_t fUseGeneratedKine;      // flag to use the generated pt, eta phi

  Bool_t fEvalCentrality;        // flag to enable centrality determination
  TString fCentralityEstimator;  // centrality estimator
    
  Int_t fConfiguration;          // Configuration to run: 0 = slow (all variables), 1 = fast (min variables)

  TH1I  *fHistEventsProcessed;   //! histo for monitoring the number of events processed slot 1

  ClassDef(AliCFSingleTrackEfficiencyTask,6)
};

#endif
