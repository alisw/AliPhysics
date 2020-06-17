/// \class AliMultDepSpecAnalysisTask
/// \brief Task to for pT spectra vs. multiplicity analysis


#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

#define MAX_ALLOWED_MULT_BINS 500
#define PRECISION 1e-6

#include "TSystem.h"
#include "TList.h"
#include "THnSparse.h"
#include "TH1.h"
#include "TRandom3.h"

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliVEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliEventCuts.h"
#include "AliESDZDC.h"

#include "AliVTrack.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliESDtrackCuts.h"
#include "AliVHeader.h"

#include "AliMCSpectraWeights.h"

#include <iostream>
#include "AliAnalysisHelpersHist.h"

class AliMultDepSpecAnalysisTask : public AliAnalysisTaskSE
{
public:
  
  // possible axis dimensions
  enum Dimension : int
  {
    pt_meas = 0,
    pt_true,
    eta_meas,
    eta_true,
    phi_meas,
    phi_true,
    mult_meas,
    mult_true,
    zv,
    event_cuts,
    sigma_pt,
    delta_pt,
  };
  
  AliMultDepSpecAnalysisTask();
  AliMultDepSpecAnalysisTask(const char *name);
  virtual ~AliMultDepSpecAnalysisTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t*);
  virtual void   Terminate(Option_t*);
  
  // Setters
  void SetTriggerMask(unsigned int triggermask)  {fTriggerMask = triggermask;}
  void SetIsMC(bool isMC = true){fIsMC = isMC;}
  void SetIsAOD(bool isAOD = true){fIsESD = !isAOD;}
  void SetUseDataDrivenCorrections(bool useDDC = true){fMCUseDDC = useDDC;}
  void SetUseZDCCut(bool useZDC){fUseZDCCut = useZDC;}
  void SetOverridePbPbEventCuts(bool overridePbPbEventCuts){fOverridePbPbEventCuts = overridePbPbEventCuts;}
  
  void SetAxis(Dimension dim, const std::string name, const std::string title, const std::vector<double>& binEdges, int nBins = 0);
  void SetCuts(Dimension dim, const std::pair<double, double>& cuts) {};

  
  // Acceptance cuts -> to be replaced soon
  void SetMinEta(double minEta)   {fMinEta = minEta;}
  void SetMaxEta(double maxEta)   {fMaxEta = maxEta;}
  void SetMinPt(double minPt)     {fMinPt = minPt;}
  void SetMaxPt(double maxPt)     {fMaxPt = maxPt;}
  void SetMaxZv(double maxZv)     {fMaxZv = maxZv;}
  void SetMaxCent(double maxCent) {fMaxCent = maxCent;}
  void SetMinCent(double minCent) {fMinCent = minCent;}
  
  // Configure this object for a train run
  static AliMultDepSpecAnalysisTask* AddTaskMultDepSpec(std::string dataSet, int cutModeLow = 100, int cutModeHigh = 119, TString options = "", bool isMC = false);

protected:
  // interface for derived classes
  virtual void DefineDefaultAxes(int maxMult = 100); // called in AddTask
  virtual void BookHistograms();    // called in UserCreateOutputObjects
  //virtual void AnaTracks();
  //virtual void AnaParticles();
  //virtual void AnaMCTruth();
  
  
  
private:
  
  std::vector<double> GetMultBinEdges(int maxMult);
  std::vector<double> GetMultBinEdges(std::vector<int> multSteps, std::vector<int> multBinWidth);
  
  void SaveTrainMetadata();
  bool SetupTask(std::string dataSet, TString options);
  bool InitTask(bool isMC, bool isAOD, std::string dataSet, TString options, int cutMode = 100);

  TList*              fOutputList;		        //!<! Output list
  AliEventCuts        fEventCuts;             //!<! Event cuts
  AliESDtrackCuts*    fTrackCuts;             //-> Track cuts
  TRandom3*           fRand;                  //!<! Random generator
  
  std::string       fTrainMetadata;           ///<  metadata of the train run used to generate the output
  
  bool              fIsESD;			              ///< Flag for ESD usage
  bool              fIsMC;                    ///< Flag for MC usage
  bool              fUseZDCCut;               ///< Flag for zdc cut usage
  bool              fOverridePbPbEventCuts;   ///< override centrality cut in PbPb
  bool              fMCUseDDC;                ///< Flag for data driven corrections usage
  // Cuts
  unsigned int        fTriggerMask;   ///< Trigger mask
  double              fMinEta;        ///< Minimum eta cut
  double              fMaxEta;        ///< Maximum eta cut
  double              fMinPt;			    ///< Minimum pT cut
  double              fMaxPt;			    ///< Maximum pT cut
  double              fMaxZv;			    ///< Maximum absolute z vertex cut
  double              fMinCent;       ///< Minimum centrality
  double              fMaxCent;       ///< Maximum centrality
  
  // Output Histograms
  std::map<Dimension, Hist::Axis> fAxes;           ///< Axis definitions used in the histograms
  Hist::Log<TH1D> fHistTrainInfo;                  //!<! Histogram to save train metadata string as bin lable; entries correspond to number of jobs
  Hist::Hist<TH1D> fHistEventSelection;            //!<! Histogram of event selection
  Hist::Hist<TH1D> fHistEvents;                    //!<! Histogram of measured event distribution
  Hist::Hist<THnSparseF> fHistTracks;              //!<! Histogram of measured tracks
  Hist::Hist<THnSparseF> fHistRelPtReso;           //!<! Histogram of relatvie pT resolution from covariance matrix
  
  Hist::Hist<THnSparseF> fHistMCEventEfficiency;   //!<! Histogram of selelcted events vs Nch
  Hist::Hist<THnSparseF> fHistMCRelPtReso;         //!<! Histogram of relative pt resolution from mc
  Hist::Hist<THnSparseF> fHistMCMultCorrelMatrix;  //!<! Histogram of multilicity correlation
  Hist::Hist<THnSparseF> fHistMCPtCorrelMatrix;    //!<! Histogram of pT correlation
  Hist::Hist<THnSparseF> fHistMCEtaCorrelMatrix;   //!<! Histogram of eta correlation
  Hist::Hist<THnSparseF> fHistMCPrimTrue;          //!<! Histogram of generated primaries
  Hist::Hist<THnSparseF> fHistMCPrimMeas;          //!<! Histogram of measured primaries
  Hist::Hist<THnSparseF> fHistMCSecMeas;           //!<! Histogram of measured secondaries
  
  // event related properties
  AliVEvent*          fEvent;			      //!<! Event object
  AliMCEvent*         fMCEvent;         //!<! MC event
  double              fMultMeas;        //!<! measured central barrel track multiplicity
  double              fMultTrue;        //!<! true multiplicity
  
  int                           fRunNumber;                 //!<! run number
  unsigned long                 fEventNumber;               //!<! event number
  unsigned int                  fTimeStamp;                 //!<! event time stamp
  
  // track related properties
  double                        fPt;                         //!<! track pT
  double                        fEta;                        //!<! track Eta
  double                        fSigmaPt;                    //!<! sigma(pT)/pT
  
  double                        fMCPt;                       //!<! mc pt
  double                        fMCEta;                      //!<! mc eta
  
  int                           fMCLabel;                    //!<! mc label
  bool                          fIsParticleInAcceptance;     //!<! particle in acceptance
  
  bool                          fMCIsChargedPrimary;          //!<! is charged primary?
  bool                          fMCIsChargedSecDecay;         //!<! is charged secondary from decay?
  bool                          fMCIsChargedSecMat;           //!<! is charged secondary from material?
  bool                          fMCIsChargedSecondary;        //!<! is charged secondary?
  
  double                        fMCParticleWeight;            //!<! scaling factor of particle to match data
  double                        fMCSecScaleWeight;            //!<! scaling factor of secondary to match data
  int                           fNRepetitions;                //!<! how often to repeat this particle to match data
  bool                          fUseRandomSeed;               ///<  use a random seed or a deterministic one (default)
  
  
  template<typename T>
  void BookHistogram(Hist::Hist<T>& hist, const std::string& name, const std::vector<Dimension>& axisNames, bool isFillWeigths = false);

  
  // Tracking functions
  bool AcceptTrackQuality(AliVTrack* track);
  double GetCentrality(AliVEvent* event);
  
  // Data driven correction related functions
  double GetSecScalingFactor(AliVParticle* particle);
  double GetParticleWeight(AliVParticle* particle);
  
  bool InitEvent();
  bool InitTrack(AliVTrack* track);
  
  template<typename Particle_t>
  bool InitParticle(Particle_t* particle);
  
  void LoopMeas(bool count = false);
  void LoopTrue(bool count = false);
  
  void FillEventHistos();
  void FillMeasTrackHistos();
  void FillMeasParticleHistos();
  void FillTrueParticleHistos();
  
  int GetNRepetitons(double scalingFactor);
  unsigned long GetSeed();
  
  AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&); // not implemented
  AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecAnalysisTask, 1); // example of analysis
  /// \endcond
};

#endif
