#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

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
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"

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
  enum Dimension : unsigned int {
    pt_meas = 0,
    pt_true,
    mult_meas,
    mult_true,
    sigma_pt,
    delta_pt,
    zv,
    LAST,
  };

  AliMultDepSpecAnalysisTask() : AliAnalysisTaskSE(){};
  AliMultDepSpecAnalysisTask(const char* name) : AliAnalysisTaskSE(name) { DefineOutput(1, TList::Class()); };
  virtual ~AliMultDepSpecAnalysisTask(){};

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);
  virtual void Terminate(Option_t*){};

  // Setters
  void SetTriggerMask(unsigned int triggermask) { fTriggerMask = triggermask; }
  void SetIsMC(bool isMC = true) { fIsMC = isMC; }
  void SetIsAOD(bool isAOD = true) { fIsESD = !isAOD; }
  void SetUseDataDrivenCorrections(bool useDDC = true) { fMCUseDDC = useDDC; }
  void SetUseZDCCut(bool useZDC) { fUseZDCCut = useZDC; }
  void SetIncludePeripheralEvents(bool includePeripheralEvents)
  {
    fIncludePeripheralEvents = includePeripheralEvents;
  }

  void SetAxis(unsigned int dim, const std::string name, const std::string title,
               const std::vector<double>& binEdges, int nBins = 0);

  void SetMinEta(double minEta) { fMinEta = minEta; }
  void SetMaxEta(double maxEta) { fMaxEta = maxEta; }
  void SetMinPt(double minPt) { fMinPt = minPt; }
  void SetMaxPt(double maxPt) { fMaxPt = maxPt; }

  // Configure this object for a train run
  static AliMultDepSpecAnalysisTask* AddTaskMultDepSpec(const std::string& dataSet,
                                                        int cutModeLow = 100, int cutModeHigh = 119,
                                                        TString options = "", bool isMC = false);
  void SaveTrainMetadata();
  bool SetupTask(std::string dataSet, TString options);
  virtual bool InitTask(bool isMC, bool isAOD, std::string dataSet, TString options, int cutMode = 100);

protected:
  virtual void DefineDefaultAxes(int maxMultMeas = 100, int maxMultTrue = 100); // called in AddTask
  virtual void BookHistograms();                     // called in UserCreateOutputObjects

  virtual bool InitEvent();
  virtual bool InitTrack(AliVTrack* track);

  // ugly workaround because template members cannot be virtual
  virtual bool InitParticle(AliMCParticle* particle)
  {
    return InitParticleBase(particle);
  } // called for ESDs
  virtual bool InitParticle(AliAODMCParticle* particle)
  {
    return InitParticleBase(particle);
  } // called for AODs

  virtual bool SelectTrack(bool count) { return true; }
  virtual bool SelectParticle(bool count) { return true; }

  void LoopMeas(bool count = false);
  void LoopTrue(bool count = false);

  virtual void FillEventHistos();
  virtual void FillEventHistosTrue();
  virtual void FillMeasTrackHistos();
  virtual void FillMeasParticleHistos();
  virtual void FillTrueParticleHistos();

  // this must be defined in header so it can be generated for types that are only present in derived classes
  template <typename T>
  void BookHistogram(Hist::Hist<T>& histContainer, const std::string& histName,
                     const std::vector<unsigned int>& dimensions, bool isFillWeigths = false)
  {
    for (auto& dim : dimensions) {
      if (fAxes.find(dim) == fAxes.end()) {
        AliFatal(Form("Not all axes for histogram %s were specified properly!", histName.data()));
        return;
      }
      histContainer.AddAxis(fAxes[dim]);
    }
    fOutputList->Add(histContainer.GenerateHist(histName));
  }

  bool AcceptTrackQuality(AliVTrack* track);
  bool InitCentrality();

  std::vector<double> GetMultBinEdges(int maxMult);
  std::vector<double> GetMultBinEdges(std::vector<int> multSteps, std::vector<int> multBinWidth);
  template <typename Particle_t>
  bool InitParticleBase(Particle_t* particle);
  double GetSecScalingFactor(AliVParticle* particle);
  double GetParticleWeight(AliVParticle* particle);
  unsigned long GetSeed();
  int GetNRepetitons(double scalingFactor);
  AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&);            // not implemented
  AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented

  std::unique_ptr<TList> fOutputList{};          //!<! Output list
  std::unique_ptr<AliEventCuts> fEventCuts{};    //!<! Event cuts
  std::unique_ptr<AliESDtrackCuts> fTrackCuts{}; //->  Track cuts
  std::unique_ptr<TRandom3> fRand{};             //!<! Random generator

  std::string fTrainMetadata{}; ///<  metadata of the train run used to generate the output

  bool fIsNominalSetting{};        ///< Flag to  propagate if this is the nominal cut setting
  bool fIsESD{true};               ///< Flag for ESD usage
  bool fIsMC{};                    ///< Flag for MC usage
  bool fUseZDCCut{};               ///< Flag for zdc cut usage
  bool fIncludePeripheralEvents{}; ///< include peripheral A-A events (cent>90)
  bool fMCUseDDC{};                ///< Flag for data driven corrections usage
  // Cuts
  unsigned int fTriggerMask{AliVEvent::kAnyINT}; ///< Trigger mask
  double fMinEta{-0.8};                          ///< Minimum eta cut
  double fMaxEta{0.8};                           ///< Maximum eta cut
  double fMinPt{0.15};                           ///< Minimum pT cut
  double fMaxPt{50.0};                           ///< Maximum pT cut

  // Output Histograms
  std::map<unsigned int, Hist::Axis> fAxes{}; ///< Axis definitions used in the histograms
  Hist::Log<TH1I> fHist_trainInfo{};          //!<! train metadata string as bin lable and number of jobs as fill
  Hist::Log<TH1I> fHist_runStatistics{};      //!<! number of reconstructed events per run

  Hist::Hist<TH1D> fHist_zVtx_evt_trig_gen{}; //!<! generated z vertex position
  Hist::Hist<TH1D> fHist_zVtx_evt_trig{};     //!<! measured z vertex position

  Hist::Hist<TH1D> fHist_multDist_evt_meas{};   //!<! measured event distribution
  Hist::Hist<TH2D> fHist_multPtSpec_trk_meas{}; //!<! measured tracks
  Hist::Hist<TH2D> fHist_ptReso_trk_meas{};     //!<! relatvie pT resolution from covariance matrix

  // MC-only histograms
  Hist::Hist<TH1D> fHist_multDist_evt_gen{};      //!<! multiplicity distribution of generated (triggered, z<10) events
  Hist::Hist<TH2D> fHist_ptReso_trk_true{};       //!<! relative pt resolution from mc

  Hist::Hist<THnSparseF> fHist_multCorrel_evt{};     //!<! multilicity correlation of events [sparse wins for large histos (540% of 100x100 hist), pPb (140% of 200x200), PbPb (25% of 500x500)]
  Hist::Hist<THnSparseF> fHist_multCorrel_prim{};    //!<! multiplicity correlation for reconstructed primary particles
  Hist::Hist<TH2D> fHist_ptCorrel_prim{};            //!<! pT correlation of measured primary particles
  Hist::Hist<TH2D> fHist_multPtSpec_prim_gen{};      //!<! generated primaries (that fulfil trigger condition and maybe z vertex cut)
  Hist::Hist<TH2D> fHist_multPtSpec_prim_meas{};     //!<! reconstructed primaries
  Hist::Hist<TH2D> fHist_multPtSpec_trk_prim_meas{}; //!<! measured primaries
  Hist::Hist<TH2D> fHist_multPtSpec_trk_sec_meas{};  //!<! measured secondaries
  
  // event related properties
  AliVEvent* fEvent{};    //!<! Event object
  AliMCEvent* fMCEvent{}; //!<! MC event
  double fVtxZ{};         //!<! z position of z vertex
  double fMCVtxZ{};       //!<! z position of z vertex in MC truth
  double fMultMeas{};     //!<! measured central barrel track multiplicity
  double fMultTrue{};     //!<! true multiplicity

  bool fIsFirstEventInJob{true};     //!<!
  int fRunNumber{};                  //!<! run number
  unsigned long fEventNumber{};      //!<! event number
  unsigned int fTimeStamp{};         //!<! event time stamp
  float fCent{};                     //!<! event centrality
  bool fIsAcceptedPeripheralEvent{}; //!<! event with centrality > 90% that passes the selection criteria
  bool fAcceptEvent{};               //!<! decision if current event is selected
  bool fAcceptEventMC{};             //!<! decision if current event is selected in mc truth

  // track related properties
  double fPt{};      //!<! track pt
  double fEta{};     //!<! track eta
  double fPhi{};     //!<! track phi
  double fSigmaPt{}; //!<! sigma(pt)/pt

  double fMCPt{};  //!<! mc pt
  double fMCEta{}; //!<! mc eta
  double fMCPhi{}; //!<! mc phi
  
  int fMCLabel{}; //!<! mc label

  bool fMCIsChargedPrimary{};   //!<! is charged primary?
  bool fMCIsChargedSecDecay{};  //!<! is charged secondary from decay?
  bool fMCIsChargedSecMat{};    //!<! is charged secondary from material?
  bool fMCIsChargedSecondary{}; //!<! is charged secondary?

  double fMCParticleWeight{1.}; //!<! scaling factor of particle to match data
  double fMCSecScaleWeight{1.}; //!<! scaling factor of secondary to match data
  int fNRepetitions{1};         //!<! how often to repeat this particle to match data
  bool fUseRandomSeed{};        ///<  use a random seed or a deterministic one (default)

  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecAnalysisTask, 1);
  /// \endcond
};

#endif
