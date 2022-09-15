#ifndef AliMultDepSpecAnalysisTask_cxx
#define AliMultDepSpecAnalysisTask_cxx

#define PRECISION 1e-6
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelpersHist.h"
#include "AliESDtrackCuts.h"

class AliMCSpectraWeights;
class AliVTrack;
class AliVParticle;
class TRandom3;
class AliEventCuts;
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
  static AliMultDepSpecAnalysisTask* AddTaskMultDepSpec(const std::string& dataSet, TString options,
                                                        int cutModeLow, int cutModeHigh, bool isMC);
  void SaveTrainMetadata();
  bool SetupTask(std::string dataSet, TString options);
  bool InitTask(bool isMC, bool isAOD, std::string dataSet, TString options, int cutMode = 100);

protected:
  void DefineDefaultAxes();
  void BookHistograms();
  void FillTrackQA(AliESDtrack* track);

  bool InitEvent();
  bool InitTrack(AliVTrack* track);

  void LoopMeas(bool count = false);
  void LoopTrue(bool count = false);

  template <typename T>
  void BookHistogram(Hist<T>& histContainer, const std::string& histName, const std::vector<unsigned int>& dimensions);

  bool InitCentrality();

  template <typename Particle_t>
  bool InitParticle(int particleID);
  unsigned long GetSeed();
  int GetNRepetitons(double scalingFactor);
  AliMultDepSpecAnalysisTask(const AliMultDepSpecAnalysisTask&);            // not implemented
  AliMultDepSpecAnalysisTask& operator=(const AliMultDepSpecAnalysisTask&); // not implemented

  std::unique_ptr<TList> fOutputList{};       //!<! output list
  TList* fQAList{nullptr};                    //!<! QA list
  std::unique_ptr<AliEventCuts> fEventCuts{}; //!<! event cuts
  AliESDtrackCuts* fTrackCuts{nullptr};       //-> track cuts
  std::unique_ptr<TRandom3> fRand{};          //!<! random generator

  std::string fTrainMetadata{};       ///<  metadata of the train run used to generate the output
  std::string fMCSelectedGenerator{}; ///<  selected generator (used in case of MCs with embedded signal to select the MB event)

  unsigned int fMCEventClass{EventClass::fiducial}; ///< baseline event class that this measurement should be corrected to
  bool fIsNominalSetting{};                         ///< flag to  propagate if this is the nominal cut setting
  bool fIsESD{true};                                ///< flag for ESD usage
  bool fIsMC{};                                     ///< flag for MC usage
  bool fIsNewReco{};                                ///< flag for new reconstructions (after mid 2020) where tighter chi2 cut has to be used and out-of-bunch pileup is simulated in MCs
  bool fMCEnableDDC{true};                          ///< flag for data driven corrections usage
  int fMCDDCMode{0};                                ///< whether to run data-driven corrections in nominal (0) or in systematic mode (-1,1)
  int fHighPtMode{0};                               ///< extend the binning and pt cuts for a high-pt analysis (mode 0: do not extend beyond 10 GeV/c, mode 1: extend to 50 GeV/c, mode 2: extend to 100 GeV/c)

  // cuts
  unsigned int fTriggerMask{AliVEvent::kMB | AliVEvent::kINT7}; ///< trigger mask
  double fMinEta{-0.8};                                         ///< minimum eta cut
  double fMaxEta{0.8};                                          ///< maximum eta cut
  double fMinPt{0.15};                                          ///< minimum pT cut
  double fMaxPt{10.0};                                          ///< maximum pT cut
  int fMaxMultMeas{100};                                        ///< maximum measured multiplicity
  int fMaxMultTrue{100};                                        ///< maximum true multiplicity

  std::map<unsigned int, Axis> fAxes{}; //!<! axis definitions used in the histograms

  Hist<TH1D> fHist_multDist_evt_meas{};   //!<! measured event distribution (contains contamination from events not belonging to specified class or with true vertex position outside fiducial z range)
  Hist<TH2D> fHist_multPtSpec_trk_meas{}; //!<! measured tracks (contains contamination from secondary particles, particles smeared into acceptance, particles from pileup events and tracks originating from background events as defined above)

  // MC-only histograms
  Hist<TH2D> fHist_multPtSpec_trk_prim_meas{};    //!<! tracks from measured primaries (no contamination from any of the sources mentioned above)
  Hist<TH2D> fHist_multPtSpec_trk_sec_meas{};     //!<! tracks from measured secondaries (no contamination from particles smeared into acceptance or background / pileup events)
  Hist<TH2D> fHist_multPtSpec_trk_meas_evtcont{}; //!<! tracks from events that are measured, but do not belong to the desired class of events or have true vertex position outside fiducial z range

  Hist<TH2D> fHist_multPtSpec_prim_meas{};        //!<! measured primary charged particles as function of true properties (no contamination from background events)
  Hist<TH2D> fHist_multPtSpec_prim_gen{};         //!<! generated primary charged particles as function of true properties (from events within specified class and with proper vertex position)
  Hist<TH2D> fHist_multPtSpec_prim_gen_evtloss{}; //!<! generated primary charged particles of events that did not pass the event selection as function of multiplicity and pt (subsample of the above)
  Hist<TH1D> fHist_multDist_evt_gen{};            //!<! generated event distribution  (from events within specified class and with proper vertex position)
  Hist<THnSparseF> fHist_multCorrel_evt{};        //!<! multilicity correlation of measured events (excluding background events)
  Hist<THnSparseF> fHist_multCorrel_prim{};       //!<! multiplicity correlation of measured primary charged particles (excluding particles from background or pileup events, secondary particles and particles smeared in from outside acceptance)
  Hist<TH2D> fHist_ptCorrel_prim{};               //!<! pT correlation of measured primary charged particles  (excluding particles from background or pileup events, secondary particles and particles smeared in from outside acceptance)

  Hist<TH1D> fHist_multDist_evt_gen_trig{};      //!<! generated event distribution (from events within specified class and with proper vertex position) that in addition fulfils the trigger condition [to disentangle trigger eff from reco eff ]
  Hist<TH2D> fHist_multPtSpec_prim_gen_notrig{}; //!<! generated primary charged particles of events that did not fulfil physics selection and trigger condition as function of multiplicity and pt

  Hist<TH2D> fHist_multPtSpec_trk_inter{}; //!<! reference for interim step of unfolding

  // QA histograms
  Log<TH1I> fHist_trainInfo{};         //!<! train metadata string as bin lable and number of compute jobs as bin content
  Log<TH1I> fHist_runStatistics{};     //!<! number of measured events per run
  Hist<TH1D> fHist_eventSelection{};   //!<! logging histogram to for event selection steps
  Hist<TH1D> fHist_eventSelectionMC{}; //!<! logging histogram to for event selection steps in mc

  Hist<TH1D> fHist_zVtxGen{};  //!<! true z vertex position of all generated events (without z vertex position cut)
  Hist<TH1D> fHist_zVtxMeas{}; //!<! measured z vertex position of all measured events (without z vertex position cut)

  Hist<TH2D> fHist_deltaPt{}; //!<! relative track pt resolution (no contamination from particles smeared into acceptance or background events)
  Hist<TH2D> fHist_sigmaPt{}; //!<! relatvie pT resolution from covariance matrix (from all tracks, including contamination defined above)

  Hist<TH2D> fHist_dcaXY{}; //!<! dca in xy plane vs pt
  Hist<TH1D> fHist_dcaZ{};  //!<! dca along the beam line

  Hist<TH1D> fHist_signed1Pt{};                     //!<! signed 1/pt (1/(GeV/c))
  Hist<TH1D> fHist_eta{};                           //!<! pseudorapidity distribution
  Hist<TH1D> fHist_phi{};                           //!<! azimuthal angle phi
  Hist<TH1D> fHist_itsNCls{};                       //!<! found clusters ITS
  Hist<TH1D> fHist_itsHits{};                       //!<! hitmap ITS
  Hist<TH1D> fHist_itsChi2NCl{};                    //!<! chi2 per cluster ITS
  Hist<TH1D> fHist_tpcNClsFindable{};               //!<! findable clusters TPC
  Hist<TH1D> fHist_tpcNClsFound{};                  //!<! found clusters TPC
  Hist<TH1D> fHist_tpcCrossedRows{};                //!<! crossed rows in TPC
  Hist<TH1D> fHist_tpcCrossedRowsOverFindableCls{}; //!<! rows / findable clusters TPC
  Hist<TH1D> fHist_tpcChi2NCl{};                    //!<! chi2 per cluster TPC
  Hist<TH1D> fHist_tpcGoldenChi2{};                 //!<! chi2 global vs tpc constrained track
  Hist<TH1D> fHist_tpcGeomLength{};                 //!<! track length in active volume of the TPC

  AliMCSpectraWeights* fMCSpectraWeights{}; //!<! fMCSpectraWeights for data-driven corrections (particle composition and secondary scaling)

  // transient event related properties
  AliVEvent* fEvent{};    //!<! reconstructed event
  AliMCEvent* fMCEvent{}; //!<! MC event
  bool fIsTriggered{};    //!<! if event fired trigger (and passes physics selection)
  double fVtxZ{};         //!<! position of z vertex
  double fMCVtxZ{};       //!<! position of MC truth z vertex
  double fMultMeas{};     //!<! measured central barrel track multiplicity after selections
  double fMultTrue{};     //!<! true multiplicity within kinematic range of measurement

  bool fIsFirstEventInJob{true}; //!<! flag if the current event is the first one in this computing job
  int fRunNumber{};              //!<! run number
  unsigned long fEventNumber{};  //!<! event number
  unsigned int fTimeStamp{};     //!<! event time stamp
  float fCent{};                 //!<! event centrality
  bool fEventPassesPhysSel{};    //!<! event passes physics selection (usually also including pileup rejection)
  bool fMCIsGoodZPos{};          //!<! is mc event within acceptance (z < 10 cm)?
  bool fMCIsGoodEventClass{};    //!<! decision if current event is in specified baseline 'generated' event class
  bool fAcceptEvent{};           //!<! decision if current event is selected
  bool fMCAcceptEvent{};         //!<! decision if current event is of 'generated' event class and has good vertex position
  bool fMCIsINELGT0{};           //!<! decision if this event is in INEL>0 event class (at least one charged particle in abs(eta) < 1)

  // transient track related properties
  double fPt{};      //!<! track pt
  double fEta{};     //!<! track eta
  double fPhi{};     //!<! track phi
  double fSigmaPt{}; //!<! sigma(pt)/pt

  double fMCPt{};  //!<! mc pt
  double fMCEta{}; //!<! mc eta
  double fMCPhi{}; //!<! mc phi

  int fMCLabel{}; //!<! mc label

  bool fMCIsInjectedSignal{};   //!<! is injected signal?
  bool fMCIsChargedPrimary{};   //!<! is charged primary?
  bool fMCIsChargedSecondary{}; //!<! is charged secondary?

  double fMCParticleWeight{1.}; //!<! scaling factor of particle to match data
  int fNRepetitions{1};         //!<! how often to repeat this particle to match data

  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecAnalysisTask, 1);
  /// \endcond
};

#endif
