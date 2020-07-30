/// \class AliMultDepSpecAnalysisTask
/// \brief Task to for pT spectra vs. multiplicity analysis in the underlying event


#ifndef AliMultDepSpecAnalysisTaskUE_cxx
#define AliMultDepSpecAnalysisTaskUE_cxx
#include "AliMultDepSpecAnalysisTask.h"

class AliMultDepSpecAnalysisTaskUE : public AliMultDepSpecAnalysisTask
{
public:

  // possible axis dimensions
  enum DimensionUE : unsigned int
  {
   pt_lead_meas = Dimension::END,
   phi_lead_meas,
   delta_phi_lead,
   delta_pt_lead,
  };

  AliMultDepSpecAnalysisTaskUE();
  AliMultDepSpecAnalysisTaskUE(const char *name);
  virtual ~AliMultDepSpecAnalysisTaskUE();

  // Additional Setters
  void SetIsUE(bool isUE = true){fIsUE = isUE;}

  static AliMultDepSpecAnalysisTaskUE* AddTaskMultDepSpecUE(const std::string& dataSet, int cutModeLow = 100, int cutModeHigh = 119, TString options = "", bool isMC = false, bool isUE = true);

  virtual bool InitTask(bool isUE, bool isMC, bool isAOD, std::string dataSet, TString options, int cutMode =  100);

protected:

  virtual void DefineDefaultAxes(int maxMult = 100); // called in AddTask
  virtual void BookHistograms();    // called in UserCreateOutputObjects
  virtual bool InitEvent();
  virtual bool InitTrack(AliVTrack* track);
  virtual bool InitParticle(AliMCParticle* particle);     // called for ESDs
  virtual bool InitParticle(AliAODMCParticle* particle);  // called for AODs

  virtual bool SelectTrack();
  virtual bool SelectParticle();

  virtual void FillEventHistos();
  virtual void FillMeasTrackHistos();
  virtual void FillMeasParticleHistos();
  virtual void FillTrueParticleHistos();

  void FindLeadingTrack();

private:

  AliMultDepSpecAnalysisTaskUE(const AliMultDepSpecAnalysisTaskUE&); // not implemented
  AliMultDepSpecAnalysisTaskUE& operator=(const AliMultDepSpecAnalysisTaskUE&); // not implemented

  double  fPtLeadMIN; // !<! Minimum leading track pt
  bool    fIsUE;  // flag for measuring UE or full phi range

  // Additional UE Output histograms
  Hist::Hist<TH1D> fHistLeadPt;              //!<! Histogram of measured leading track pT
  Hist::Hist<TH1D> fHistLeadPhi;             //!<! Histogram of measured leading track phi

  Hist::Hist<TH1D> fHistMCResoPtLead;          //!<! Histogram of Difference in Leading Track Pt and corresponding Particle Pt
  Hist::Hist<TH1D> fHistMCResoPhiLead;         //!<! Histogram of Difference in Leading Track Phi and corresponding Particle Phi
  Hist::Hist<TH1D> fHistDiffToMCPtLead;        //!<! Histogram of Difference in Leading Track Pt and Leading Particle Pt
  Hist::Hist<TH1D> fHistDiffToMCPhiLead;       //!<! Histogram of Difference in Leading Track Phi and Leading Particle Phi

  Hist::Hist<TH2D> fHistPlateau;         //!<! Histogram for deriving Plateau Plot

  // Additional event related properties
  double                        fPtLead;                     //!<! Leading track pT
  double                        fPhiLead;                    //!<! Leading track phi

  double                        fMCPtLead;                   //!<! Leading MC particle pT
  double                        fMCPhiLead;                  //!<! Leading MC particle phi
  double                        fMCPtOfLead;                 //!<! To leading MC track corresponding particle pT
  double                        fMCPhiOfLead;                //!<! To leading MC track corresponding particle phi

  /// \cond CLASSIMP
  ClassDef(AliMultDepSpecAnalysisTaskUE, 1); // example of analysis
  /// \endcond
};

#endif
