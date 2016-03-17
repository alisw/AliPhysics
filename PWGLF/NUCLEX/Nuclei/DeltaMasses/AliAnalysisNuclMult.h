#ifndef ALIANALYSISNUCLMULT_H
#define ALIANALYSISNUCLMULT_H

#include <AliAnalysisTaskSE.h>
#include <AliPIDResponse.h>

class TH2F;
class TH2D;
class TH1F;
class TH1I;
class TF1;
class TF2;
class TGraph;
class TH3F;
class AliESDtrackCuts;
class TProfile;
class TFile;
class TList;
class TObject;
class AliAODEvent;
class AliESDEvent;
class AliVEvent;
class AliAnalysisFilter;
class AliPPVsMultUtils;
class AliESDtrackCuts;

class AliAnalysisNuclMult : public AliAnalysisTaskSE {
 
 public:
  
  AliAnalysisNuclMult();
  AliAnalysisNuclMult(const char *name);
  
  virtual ~AliAnalysisNuclMult();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetTrackFilter(AliAnalysisFilter *filter) {fTrackFilter = filter;};
  void SetMultEstimator(Int_t tMultEstimator=1) {iMultEstimator = tMultEstimator;};
  void SetMultiplicityRange(Float_t multMin=-999, Float_t multMax=999) {//multMin=-999 -> analysis over all Minimum Bias collisions
    multiplicityMin=multMin;
    multiplicityMax=multMax;
  };
  
 private:
  
  AliAnalysisNuclMult(const AliAnalysisNuclMult &old); 
  AliAnalysisNuclMult& operator=(const AliAnalysisNuclMult &source);

  AliAODEvent* fAOD;                              //! AOD object
  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliAnalysisFilter *fTrackFilter;                //  track filter object
  AliPPVsMultUtils *fPPVsMultUtils;               //  for multiplicity estimator
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response
  
  TList *fList;                                   //! output container
  
  Int_t iMultEstimator;                           //  iMultEstimator: 0-> VZERO Amplitude Estimator; 1->Mid-psudorapidity Estimator;
  Float_t multiplicityMin;                        //  min. mult.
  Float_t multiplicityMax;                        //  max. mult.
  
  Int_t stdFlagPid[9];

  TH1I *htriggerMask;                             //! trigger mask
  TH1I *htriggerMask_noMB;                        //! trigger mask (MB excluded)
  TH1I *hNspdVertex;                              //! num. of spd vertices
  TH1F *hzvertex;                                 //! z-vertex distribution
  TH1I *hpileUp;                                  //! pileup in spd?
  TH1I *hmult;                                    //! multiplicity distribution
  TH1I *hNtrack;                                  //! number of tracks distribution
  
  TH2F *fdEdxVSp[2];                              //! dedx vs pTPC
  TProfile *hDeDxExp[9];                          //! TPC splines
  TH2F *fNsigmaTPC[18];                           //! NsigmaTPC vs. pT

  TH3F *fDca[2][18];                              //! DCAxy vs DCAz

  TH2F *fBetaTOFvspt[2];                          //! beta (TOF) vs pT
  TProfile *hBetaExp[9];                          //! TOF expected beta

  TH2F *fM2tof[2];                                //! M2 computed with TOF
  TH2F *fM2vspt[18];                              //! M2 computed with TOF (TPC cut)

  TF1 *fpTcorr[2];                                //! pT corrections for (anti-)deuteron
  
  void FillDca(Double_t DCAxy, Double_t DCAz, AliVTrack *track);
  void PtCorrection(Double_t &pt, Int_t FlagPid, Short_t charge);
  
  ClassDef(AliAnalysisNuclMult, 3);
  
};

#endif
