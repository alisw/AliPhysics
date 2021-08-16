#ifndef AliMultDepSpecAnalysisTaskOLD_cxx
#define AliMultDepSpecAnalysisTaskOLD_cxx

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

class AliMultDepSpecAnalysisTaskOLD : public AliAnalysisTaskSE
{
public:
  // possible axis dimensions
  enum Dimension : unsigned int
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
    phi,
    sum_pt,
    LAST,
  };

  AliMultDepSpecAnalysisTaskOLD();
  AliMultDepSpecAnalysisTaskOLD(const char* name);
  virtual ~AliMultDepSpecAnalysisTaskOLD();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t*);
  virtual void Terminate(Option_t*);

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
  void SetCuts(unsigned int dim, const std::pair<double, double>& cuts){};

  // Acceptance cuts -> to be replaced soon
  void SetMinEta(double minEta) { fMinEta = minEta; }
  void SetMaxEta(double maxEta) { fMaxEta = maxEta; }
  void SetMinPt(double minPt) { fMinPt = minPt; }
  void SetMaxPt(double maxPt) { fMaxPt = maxPt; }

  // Configure this object for a train run
  static AliMultDepSpecAnalysisTaskOLD* AddTaskMultDepSpec(const std::string& dataSet,
                                                        int cutModeLow = 100, int cutModeHigh = 119,
                                                        TString options = "", bool isMC = false);
  void SaveTrainMetadata();

  bool SetupTask(std::string dataSet, TString options);
  virtual bool InitTask(bool isMC, bool isAOD, std::string dataSet, TString options,
                        int cutMode = 100);

protected:
  virtual void DefineDefaultAxes(int maxMult = 100); // called in AddTask
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
  virtual void FillMeasTrackHistos();
  virtual void FillMeasParticleHistos();
  virtual void FillTrueParticleHistos();

  // this must be defined in header so it can be generated for types that are only present in
  // derived classes
  template <typename T>
  void BookHistogram(Hist::Hist<T>& histContainer, const std::string& histName,
                     const std::vector<unsigned int>& dimensions, bool isFillWeigths = false)
  {
    for(auto& dim : dimensions)
    {
      if(fAxes.find(dim) == fAxes.end())
      {
        AliFatal(Form("Not all axes for histogram %s were specified properly!", histName.data()));
        return;
      }
      histContainer.AddAxis(fAxes[dim]);
    }
    fOutputList->Add(histContainer.GenerateHist(histName));
  }

  bool AcceptTrackQuality(AliVTrack* track);
  double GetCentrality(AliVEvent* event);

  std::vector<double> GetMultBinEdges(int maxMult);
  std::vector<double> GetMultBinEdges(std::vector<int> multSteps, std::vector<int> multBinWidth);
  template <typename Particle_t>
  bool InitParticleBase(Particle_t* particle);
  double GetSecScalingFactor(AliVParticle* particle);
  double GetParticleWeight(AliVParticle* particle);
  unsigned long GetSeed();
  int GetNRepetitons(double scalingFactor);
  AliMultDepSpecAnalysisTaskOLD(const AliMultDepSpecAnalysisTaskOLD&);            // not implemented
  AliMultDepSpecAnalysisTaskOLD& operator=(const AliMultDepSpecAnalysisTaskOLD&); // not implemented

  TList* fOutputList{};          //!<! Output list
  AliEventCuts fEventCuts{};     //!<! Event cuts
  AliESDtrackCuts* fTrackCuts{}; //-> Track cuts
  TRandom3* fRand{};             //!<! Random generator

  std::string fTrainMetadata{}; ///<  metadata of the train run used to generate the output

  bool fIsESD{ true };             ///< Flag for ESD usage
  bool fIsMC{};                    ///< Flag for MC usage
  bool fUseZDCCut{};               ///< Flag for zdc cut usage
  bool fIncludePeripheralEvents{}; ///< include peripheral A-A events (cent>90)
  bool fMCUseDDC{};                ///< Flag for data driven corrections usage
  // Cuts
  unsigned int fTriggerMask{ AliVEvent::kAnyINT }; ///< Trigger mask
  double fMinEta{ -10. };                          ///< Minimum eta cut
  double fMaxEta{ 10. };                           ///< Maximum eta cut
  double fMinPt{ 0.0 };                            ///< Minimum pT cut
  double fMaxPt{ 50.0 };                           ///< Maximum pT cut

  // Output Histograms
  std::map<unsigned int, Hist::Axis> fAxes{}; ///< Axis definitions used in the histograms
  Hist::Log<TH1D> fHistTrainInfo{};       //!<! Histogram to save train metadata string as bin lable
  Hist::Log<TH1D> fHistTrainInfoUE{};       //!<! Histogram to save train metadata string as bin lable
  Hist::Hist<TH1D> fHistEventSelection{}; //!<! event selection
  Hist::Hist<TH1D> fHistEvents{};         //!<! measured event distribution
  Hist::Hist<THnSparseF> fHistTracks{};   //!<! measured tracks
  Hist::Hist<THnSparseF> fHistRelPtReso{}; //!<! relatvie pT resolution

  Hist::Hist<THnSparseF> fHistMCEventEfficiency{};  //!<! selelcted events vs Nch
  Hist::Hist<THnSparseF> fHistMCRelPtReso{};        //!<! relative pt resolution from mc
  Hist::Hist<THnSparseF> fHistMCMultCorrelMatrix{}; //!<! multilicity correlation
  Hist::Hist<THnSparseF> fHistMCPtCorrelMatrix{};   //!<! pT correlation
  Hist::Hist<THnSparseF> fHistMCEtaCorrelMatrix{};  //!<! eta correlation
  Hist::Hist<THnSparseF> fHistMCPrimTrue{};         //!<! generated primaries
  Hist::Hist<THnSparseF> fHistMCPrimMeas{};         //!<! measured primaries
  Hist::Hist<THnSparseF> fHistMCSecMeas{};          //!<! measured secondaries

  // event related properties
  AliVEvent* fEvent{};    //!<! Event object
  AliMCEvent* fMCEvent{}; //!<! MC event
  double fMultMeas{};     //!<! measured central barrel track multiplicity
  double fMultTrue{};     //!<! true multiplicity
  double fSumPtMeas{};    //!<! measured sum pt
  double fSumPtTrue{};    //!<! true sum pt

  bool fIsFirstEventInJob{ true };   //!<!
  int fRunNumber{};                  //!<! run number
  unsigned long fEventNumber{};      //!<! event number
  unsigned int fTimeStamp{};         //!<! event time stamp
  double fCent{};                    //!<! event centrality
  bool fIsAcceptedPeripheralEvent{}; //!<! event with centrality > 90% that passes the selection
                                     //!< criteria

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

  double fMCParticleWeight{ 1. }; //!<! scaling factor of particle to match data
  double fMCSecScaleWeight{ 1. }; //!<! scaling factor of secondary to match data
  int fNRepetitions{ 1 };         //!<! how often to repeat this particle to match data
  bool fUseRandomSeed{};          ///<  use a random seed or a deterministic one (default)

  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecAnalysisTaskOLD, 1);
  /// \endcond
};

#endif
