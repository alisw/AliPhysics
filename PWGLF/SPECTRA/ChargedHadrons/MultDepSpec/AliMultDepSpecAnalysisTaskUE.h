#ifndef AliMultDepSpecAnalysisTaskUE_cxx
#define AliMultDepSpecAnalysisTaskUE_cxx

#include "AliMultDepSpecAnalysisTaskOLD.h"

class AliMultDepSpecAnalysisTaskUE : public AliMultDepSpecAnalysisTaskOLD
{
public:
  // possible axis dimensions
  enum DimensionUE : unsigned int
  {
    pt_lead_meas = Dimension::LAST,
    phi_lead_meas,
    delta_phi_lead,
    delta_pt_lead,
    delta_phi
  };

  enum PhiRange : unsigned int
  {
    full = 0,
    toward,
    away,
    transverse
  };

  AliMultDepSpecAnalysisTaskUE();
  AliMultDepSpecAnalysisTaskUE(const char* name);
  virtual ~AliMultDepSpecAnalysisTaskUE();

  // Additional Setters
  void SetPhiRange(float min, float max) { fCosPhiMin = min; fCosPhiMax = max; }
  void SetPtLeadCut(float ptLeadMIN = 0.01) { fPtLeadCut = ptLeadMIN; }

  static AliMultDepSpecAnalysisTaskUE* AddTaskMultDepSpecUE(
    const std::string& dataSet, int cutModeLow = 100, int cutModeHigh = 119, TString options = "",
    bool isMC = false, AliMultDepSpecAnalysisTaskUE::PhiRange phiRange = full, double ptLeadMIN = 0.01,  const char* suffix = "");

  virtual bool InitTask(AliMultDepSpecAnalysisTaskUE::PhiRange phiRange, bool isMC, bool isAOD, std::string dataSet, TString options,
                        int cutMode = 100, float ptLeadMIN = 0.01);

protected:
  virtual void DefineDefaultAxes(int maxMult = 100); // called in AddTask
  virtual void BookHistograms();                     // called in UserCreateOutputObjects
  virtual bool InitEvent();
  virtual bool InitTrack(AliVTrack* track);
  virtual bool InitParticle(AliMCParticle* particle);    // called for ESDs
  virtual bool InitParticle(AliAODMCParticle* particle); // called for AODs

  virtual bool SelectTrack(bool count);
  virtual bool SelectParticle(bool count);

  virtual void FillEventHistos();
  virtual void FillMeasTrackHistos();
  virtual void FillMeasParticleHistos();
  virtual void FillTrueParticleHistos();

  void FindLeadingTrack();

private:
  AliMultDepSpecAnalysisTaskUE(const AliMultDepSpecAnalysisTaskUE&);            // not implemented
  AliMultDepSpecAnalysisTaskUE& operator=(const AliMultDepSpecAnalysisTaskUE&); // not implemented

  float fCosPhiMin{}; // lower bound for cosphi range
  float fCosPhiMax{}; // upper bound for cosphi range
  float fPtLeadCut{};

  // Additional UE Output histograms
  Hist::Hist<TH1D> fHistLeadPt{};        //!<! measured leading track pT
  Hist::Hist<TH1D> fHistLeadPhi{};       //!<! measured leading track phi
  Hist::Hist<TH1D> fHistDeltaPhi{};      //!<! measured phi difference between leading track and current track
  Hist::Hist<TH1D> fHistPtLeadCutLoss{}; //!<! events, discarded through pT Lead cut

  Hist::Hist<TH1D> fHistMCResoPtLead{};  //!<! Difference Leading Track Pt and corresponding Particle Pt
  Hist::Hist<TH1D> fHistMCResoPhiLead{}; //!<! Difference Leading Track Phi and corresponding Particle Phi

  Hist::Hist<TH2D> fHistPlateau{};      //!<! Plateau Histogram
  Hist::Hist<TH2D> fHistSumPt{};        //!<! Histogram for measuring pt correlations between regions
  Hist::Hist<TH2D> fHistMultCorr{};     //!<! Histogram for measuring mult correlations between regions

  Hist::Hist<TH2D> fHistMCPlateau{};      //!<! Plateau Histogram
  Hist::Hist<TH2D> fHistMCSumPt{};        //!<! MC truth corresponding to fHistSumPt
  Hist::Hist<TH2D> fHistMCMultCorr{};     //!<! MC truth corresponding to fHistMultCorr

  Hist::Hist<TH1D> fHistMeanSumPt{};    //!<! Histogram for calculating mean sum pt in region as function of sum_pt total
  Hist::Hist<TH1D> fHistMeanMult{};     //!<! Histogram for calculating mean mult in region as function of total mult

  Hist::Hist<TH1D> fHistMCMeanSumPt{}; //!<! MC Particle information corresponding to fHistMeanSumPt
  Hist::Hist<TH1D> fHistMCMeanMult{};  //!<! MC Particle information corresponding to fHistMeanMult

  // Additional event related properties
  double fPtLead{};  //!<! Leading track pT
  double fPhiLead{}; //!<! Leading track phi

  double fMCPtOfLead{};    //!<! Leading MC particle pT
  double fMCPhiOfLead{};   //!<! Leading MC particle phi

  double fMultMeasTot{};    //!<! full phi multiplicity
  double fSumPtMeasTot{};   //!<! full phi sum pt

  double fMultTrueTot{};    //!<! generated full phi multiplicity
  double fSumPtTrueTot{};   //!<! generated full phi sum pt

  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecAnalysisTaskUE, 1); // example of analysis
  /// \endcond
};

#endif
