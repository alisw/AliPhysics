#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

#define PRECISION 1e-6
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelpersHist.h"

class AliMCSpectraWeights;
class AliVTrack;
class AliVParticle;
class TRandom3;
class AliEventCuts;
class AliESDtrackCuts;
class AliESDtrack;

class AliMultDepSpecAnalysisTask : public AliAnalysisTaskSE
{
public:
  // possible axis dimensions
  enum Dimension : unsigned int {
    pt_meas = 0,
    pt_true,
    mult_meas,
    mult_true,
    LAST,
  };

  // reference ('generated') event class for the measurement [this is where multiplicity distributions and mult dep pt spectra will be corrected to]
  enum EventClass : unsigned int {
    triggered = 0, // mc events that fulfil the (experimental!) trigger condition [varies from dataset to dataset]
    fiducial,      // mc events that produce at least one charged particle within the fiducial phase space of the measurement [trusting only that trigger eff for these events is modelled properly]
    inelgt0,       // mc events that produce at least one charged particle in |eta| < 1 [trusting simulation for extrapolation to events with particles produced only in the unmeasured region]
  };

  AliMultDepSpecAnalysisTask() : AliAnalysisTaskSE(){};
  AliMultDepSpecAnalysisTask(const char* name) : AliAnalysisTaskSE(name) { DefineOutput(1, TList::Class()); };
  ~AliMultDepSpecAnalysisTask()
  {
    if (fTrackCuts) delete fTrackCuts;
  };
  void UserCreateOutputObjects();
  void UserExec(Option_t*);
  void Terminate(Option_t*){};

  // Setters
  void SetEventClass(unsigned int eventClass) { fMCEventClass = eventClass; }
  void SetTriggerMask(unsigned int triggermask) { fTriggerMask = triggermask; }
  void SetIsMC(bool isMC = true) { fIsMC = isMC; }
  void SetIsAOD(bool isAOD = true) { fIsESD = !isAOD; }
  void SetUseDataDrivenCorrections(bool useDDC = true) { fMCUseDDC = useDDC; }
  void SetIsNewReco(bool isNewReco = true) { fIsNewReco = isNewReco; }

  void SetAxis(unsigned int dim, const std::string name, const std::string title, const std::vector<double>& binEdges, int nBins = 0);

  void SetEtaRange(double minEta, double maxEta)
  {
    fMinEta = minEta;
    fMaxEta = maxEta;
  }
  void SetPtRange(double minPt, double maxPt)
  {
    fMinPt = minPt;
    fMaxPt = maxPt;
  }

  // Configure this object for a train run
  static AliMultDepSpecAnalysisTask* AddTaskMultDepSpec(const std::string& dataSet,
                                                        int cutModeLow = 100, int cutModeHigh = 119,
                                                        TString options = "", bool isMC = false);
  void SaveTrainMetadata();
  bool SetupTask(std::string dataSet, TString options);
  bool InitTask(bool isMC, bool isAOD, std::string dataSet, TString options, int cutMode = 100);

protected:
  void DefineDefaultAxes(int maxMultMeas = 100, int maxMultTrue = 100); // called in AddTask
  void BookHistograms();                                                // called in UserCreateOutputObjects
  void FillTrackQA(AliESDtrack* track);

  bool InitEvent();
  bool InitTrack(AliVTrack* track);

  void LoopMeas(bool count = false);
  void LoopTrue(bool count = false);

  template <typename T>
  void BookHistogram(Hist::Hist<T>& histContainer, const std::string& histName, const std::vector<unsigned int>& dimensions);

  bool AcceptTrackQuality(AliVTrack* track);
  bool InitCentrality();

  template <typename Particle_t>
  bool InitParticle(Particle_t* particle);
  double GetSecScalingFactor(AliVParticle* particle);
  double GetParticleWeight(AliVParticle* particle);
  unsigned long GetSeed();
  int GetNRepetitons(double scalingFactor);
  AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&);            // not implemented
  AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented

  std::unique_ptr<TList> fOutputList{};       //!<! Output list
  TList* fQAList{nullptr};                    //!<! QA list
  std::unique_ptr<AliEventCuts> fEventCuts{}; //!<! Event cuts
  AliESDtrackCuts* fTrackCuts{nullptr};       //->  Track cuts
  std::unique_ptr<TRandom3> fRand{};          //!<! Random generator

  std::string fTrainMetadata{}; ///<  metadata of the train run used to generate the output

  unsigned int fMCEventClass{EventClass::fiducial}; ///< baseline event class that this measurement should be corrected to
  bool fIsNominalSetting{};                         ///< Flag to  propagate if this is the nominal cut setting
  bool fIsESD{true};                                ///< Flag for ESD usage
  bool fIsMC{};                                     ///< Flag for MC usage
  bool fIsNewReco{};                                ///< flag for new reconstructions (after mid 2020) where tighter chi2 cut has to be used and out-of-bunch pileup is simulated in MCs
  bool fIncludePeripheralEvents{};                  ///< include peripheral A-A events (cent>90)
  bool fMCUseDDC{};                                 ///< Flag for data driven corrections usage
  bool fIsNominalPCC{true};                         ///< whether to run particle composition correction in nominal or in systematic mode
  bool fMCIsINELGT0{};                              ///< flag for INEL>0 event (at least one charged particle in abs(eta) < 1)

  // cuts
  unsigned int fTriggerMask{AliVEvent::kMB | AliVEvent::kINT7}; ///< Trigger mask
  double fMinEta{-0.8};                                         ///< Minimum eta cut
  double fMaxEta{0.8};                                          ///< Maximum eta cut
  double fMinPt{0.15};                                          ///< Minimum pT cut
  double fMaxPt{10.0};                                          ///< Maximum pT cut

  // output Histograms
  std::map<unsigned int, Hist::Axis> fAxes{}; ///< axis definitions used in the histograms

  Hist::Hist<TH1D> fHist_multDist_evt_meas{};   //!<! measured event distribution (contains contamination from events not in specified class or with wrong vertex position)
  Hist::Hist<TH2D> fHist_multPtSpec_trk_meas{}; //!<! measured tracks (contains contamination from secondary particles, particles smeared into acceptance and tracks originating from background events as defined above )

  // MC-only histograms
  Hist::Hist<TH2D> fHist_multPtSpec_trk_prim_meas{}; //!<! tracks from measured primaries (no contamination from secondaries, particles smeared into acceptance or background events)
  Hist::Hist<TH2D> fHist_multPtSpec_trk_sec_meas{};  //!<! tracks from measured secondaries (no contamination from particles smeared into acceptance or background events)  [for QA to disentangle secondaries from other contamination]

  Hist::Hist<TH2D> fHist_multPtSpec_prim_meas{};  //!<! measured primary charged particles as function of true properties (no contamination from background events)
  Hist::Hist<TH2D> fHist_multPtSpec_prim_gen{};   //!<! generated primary charged particles as function of true properties (from events within specified class and with proper vertex position)
  Hist::Hist<TH1D> fHist_multDist_evt_gen{};      //!<! generated event distribution  (from events within specified class and with proper vertex position)
  Hist::Hist<TH1D> fHist_multDist_evt_gen_trig{}; //!<! generated event distribution (from events within specified class and with proper vertex position) that in addition fulfils the trigger condition [for QA to disentangle trigger eff from reco eff ]
  Hist::Hist<THnSparseF> fHist_multCorrel_evt{};  //!<! multilicity correlation of measured events (excluding background events)
  Hist::Hist<THnSparseF> fHist_multCorrel_prim{}; //!<! multiplicity correlation of measured primary charged particles (excluding particles from background events)
  Hist::Hist<TH2D> fHist_ptCorrel_prim{};         //!<! pT correlation of measured primary charged particles  (excluding particles from background events)

  // QA histograms
  Hist::Log<TH1I> fHist_trainInfo{};       //!<! train metadata string as bin lable and number of compute jobs as bin content
  Hist::Log<TH1I> fHist_runStatistics{};   //!<! number of measured events per run (filled only for first event in job)
  Hist::Hist<TH1D> fHist_eventSelection{}; //!<! logging histogram to for event selection steps (filled only for first event in job)

  Hist::Hist<TH1D> fHist_zVtxGen{};  //!<! true z vertex position of all generated events (without z vertex position cut)
  Hist::Hist<TH1D> fHist_zVtxMeas{}; //!<! measured z vertex position of all measured events (without z vertex position cut)

  Hist::Hist<TH2D> fHist_deltaPt{}; //!<! relative track pt resolution (no contamination from particles smeared into acceptance or background events)
  Hist::Hist<TH2D> fHist_sigmaPt{}; //!<! relatvie pT resolution from covariance matrix (from all tracks, including contamination defined above)

  Hist::Hist<TH2D> fHist_dcaXY{}; //!<! dca in xy plane vs pt
  Hist::Hist<TH1D> fHist_dcaZ{};  //!<! dca along the beam line

  Hist::Hist<TH1D> fHist_signed1Pt{};         //!<!  signed 1/pt (1/(GeV/c))
  Hist::Hist<TH1D> fHist_eta{};               //!<!  pseudorapidity
  Hist::Hist<TH1D> fHist_phi{};               //!<!  azimuthal angle phi
  Hist::Hist<TH1D> fHist_itsFoundClusters{};  //!<!  found clusters ITS
  Hist::Hist<TH1D> fHist_itsHits{};           //!<!  hitmap ITS
  Hist::Hist<TH1D> fHist_itsChi2PerCluster{}; //!<!  chi2 per cluster ITS

  Hist::Hist<TH1D> fHist_tpcFindableClusters{};                //!<!  findable clusters TPC
  Hist::Hist<TH1D> fHist_tpcFoundClusters{};                   //!<!  found clusters TPC
  Hist::Hist<TH1D> fHist_tpcCrossedRows{};                     //!<!  crossed rows in TPC
  Hist::Hist<TH1D> fHist_tpcCrossedRowsOverFindableClusters{}; //!<!  rows / findable clusters TPC
  Hist::Hist<TH1D> fHist_tpcFractionSharedClusters{};          //!<!  fraction of shared clusters TPC
  Hist::Hist<TH1D> fHist_tpcChi2PerCluster{};                  //!<!  chi2 per cluster TPC
  Hist::Hist<TH1D> fHist_tpcGoldenChi2{};                      //!<! chi2 global vs tpc constrained track
  Hist::Hist<TH1D> fHist_tpcGeomLength{};                      //!<! track length in active volume of the TPC

  // event related properties
  AliVEvent* fEvent{};                      //!<! Event object
  AliMCEvent* fMCEvent{};                   //!<! MC event
  AliMCSpectraWeights* fMCSpectraWeights{}; //!<! fMCSpectraWeights for data-driven corrections (particle composition and secondary scaling)
  bool fIsTriggered{};                      //!<! if event fired trigger
  double fVtxZ{};                           //!<! z position of z vertex
  double fMCVtxZ{};                         //!<! z position of z vertex in MC truth
  double fMultMeas{};                       //!<! measured central barrel track multiplicity
  double fMultTrue{};                       //!<! true multiplicity

  bool fIsFirstEventInJob{true}; //!<!
  int fRunNumber{};              //!<! run number
  unsigned long fEventNumber{};  //!<! event number
  unsigned int fTimeStamp{};     //!<! event time stamp
  float fCent{};                 //!<! event centrality
  bool fMCIsGoodZPos{};          //!<! is mc event within acceptance (z<10)?
  bool fMCIsGoodEventClass{};    //!<! decision if current event is in specified baseline event class
  bool fAcceptEvent{};           //!<! decision if current event is selected
  bool fMCAcceptEvent{};         //!<! decision if current event is of 'generated' event class and has good vertex position

  // track related properties
  double fPt{};      //!<! track pt
  double fEta{};     //!<! track eta
  double fPhi{};     //!<! track phi
  double fSigmaPt{}; //!<! sigma(pt)/pt

  double fMCPt{};  //!<! mc pt
  double fMCEta{}; //!<! mc eta
  double fMCPhi{}; //!<! mc phi

  int fMCLabel{}; //!<! mc label

  bool fMCIsPileupParticle{};   //!<! is particle resp. track from out-of-bunch pileup?
  bool fMCIsChargedPrimary{};   //!<! is charged primary?
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
