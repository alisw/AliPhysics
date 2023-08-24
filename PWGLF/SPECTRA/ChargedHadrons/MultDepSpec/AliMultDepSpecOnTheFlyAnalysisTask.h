#ifndef AliMultDepSpecOnTheFlyAnalysisTask_cxx
#define AliMultDepSpecOnTheFlyAnalysisTask_cxx

#define PRECISION 1e-6
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelpersHist.h"

class AliVTrack;
class AliVParticle;
class AliESDtrack;

class AliMultDepSpecOnTheFlyAnalysisTask : public AliAnalysisTaskSE
{
public:
  // possible axis dimensions
  enum Dimension : unsigned int {
    pt_true = 0,
    mult_true,
    LAST,
  };

  AliMultDepSpecOnTheFlyAnalysisTask() : AliAnalysisTaskSE(){};
  AliMultDepSpecOnTheFlyAnalysisTask(const char* name) : AliAnalysisTaskSE(name) { DefineOutput(1, TList::Class()); };
  ~AliMultDepSpecOnTheFlyAnalysisTask(){};

  void UserCreateOutputObjects();
  void UserExec(Option_t*);
  void Terminate(Option_t*){};

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
  static AliMultDepSpecOnTheFlyAnalysisTask* AddTaskMultDepSpecOnTheFlyMC(const std::string& dataSet, TString options);
  bool SetupTask(std::string dataSet, TString options);
  bool InitTask(std::string dataSet, TString options);

protected:
  void DefineDefaultAxes();
  void BookHistograms();

  bool InitEvent();
  void LoopTrue(bool count = false);

  template <typename T>
  void BookHistogram(Hist<T>& histContainer, const std::string& histName, const std::vector<unsigned int>& dimensions);

  template <typename Particle_t>
  bool InitParticle(int particleID);

  AliMultDepSpecOnTheFlyAnalysisTask(const AliMultDepSpecOnTheFlyAnalysisTask&);            // not implemented
  AliMultDepSpecOnTheFlyAnalysisTask& operator=(const AliMultDepSpecOnTheFlyAnalysisTask&); // not implemented

  std::unique_ptr<TList> fOutputList{}; //!<! output list

  int fHighPtMode{0}; ///< extend the binning and pt cuts for a high-pt analysis (mode 0: do not extend beyond 10 GeV/c, mode 1: extend to 50 GeV/c, mode 2: extend to 100 GeV/c)

  double fMinEta{-0.8};  ///< minimum eta cut
  double fMaxEta{0.8};   ///< maximum eta cut
  double fMinPt{0.15};   ///< minimum pT cut
  double fMaxPt{10.0};   ///< maximum pT cut
  int fMaxMultTrue{100}; ///< maximum true multiplicity

  std::map<unsigned int, Axis> fAxes{}; //!<! axis definitions used in the histograms

  Hist<TH2D> fHist_multPtSpec_prim_gen{}; //!<! generated primary charged particles as function of true properties (from events within specified class and with proper vertex position)
  Hist<TH1D> fHist_multDist_evt_gen{};    //!<! generated event distribution  (from events within specified class and with proper vertex position)

  AliMCEvent* fMCEvent{}; //!<! MC event
  double fMCVtxZ{};       //!<! position of MC truth z vertex
  double fMultTrue{};     //!<! true multiplicity within kinematic range of measurement

  bool fMCIsGoodZPos{};       //!<! is mc event within acceptance (z < 10 cm)?
  bool fMCIsGoodEventClass{}; //!<! decision if current event is in specified baseline 'generated' event class
  bool fMCAcceptEvent{};      //!<! decision if current event is of 'generated' event class and has good vertex position

  double fMCPt{};  //!<! mc pt
  double fMCEta{}; //!<! mc eta

  bool fMCIsChargedPrimary{};   //!<! is charged primary?
  bool fMCIsChargedSecondary{}; //!<! is charged secondary?

  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecOnTheFlyAnalysisTask, 1);
  /// \endcond
};

#endif
