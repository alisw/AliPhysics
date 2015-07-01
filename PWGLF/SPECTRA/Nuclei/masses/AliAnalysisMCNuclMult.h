#ifndef ALIANALYSISMCNUCLMULT_H
#define ALIANALYSISMCNUCLMULT_H

#include <AliAnalysisTaskSE.h>
#include <AliPIDResponse.h>

class TH2F;
class TH2I;
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
class TTree;
class AliAODEvent;
class AliESDEvent;
class AliVEvent;
class AliAnalysisFilter;
class AliPPVsMultUtils;
class AliESDtrackCuts;
// for Monte Carlo:
class AliMCEventHandler;
class AliMCEvent;
class AliStack;
class AliVParticle;

class AliAnalysisMCNuclMult : public AliAnalysisTaskSE {
 
 public:
  
  AliAnalysisMCNuclMult();
  AliAnalysisMCNuclMult(const char *name);
  
  virtual ~AliAnalysisMCNuclMult();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetTrackFilter(AliAnalysisFilter *filter) {fTrackFilter = filter;};
  void SetMultEstimator(Int_t tMultEstimator=0) {iMultEstimator = tMultEstimator;};
  void SetMultiplicityRange(Float_t multMin=-999, Float_t multMax=999) {//multMin=-999 -> analysis over all Minimum Bias collisions
    multiplicityMin=multMin;
    multiplicityMax=multMax;
  };
  
 private:
  
  AliAnalysisMCNuclMult(const AliAnalysisMCNuclMult &old); 
  AliAnalysisMCNuclMult& operator=(const AliAnalysisMCNuclMult &source);

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
  
  Int_t stdPdg[9];                                //! Pdg code of e,mu,pi,K,p,d,t,He3,He4

  TH1I *htriggerMask;                             //! trigger mask
  TH1I *htriggerMask_noMB;
  TH1I *hNspdVertex;                              //! num. of spd vertices
  TH1F *hzvertex;                                 //! z-vertex distribution
  TH1I *hpileUp;                                  //! pileup in spd?
  TH2I *fNtrVsMult;
  TProfile *prNtrVsMult;
  TH1I *hmult_tot;                                //! multiplicity distribution
  TH1I *hmult;                                    //! selected multiplicity distribution
  TH1I *hNtrack;                                  //! number of tracks distribution
    
  TH1I *hpdg;                                     //! pdg label

  TH1F *hpt[7][18];                               //! pT efficiencies
  
  TH2F *fptRecoVsTrue[4][18];                     //! pT reco vs. true
  TProfile *prptRecoVsTrue[4][18];                //! pT reco vs. true
    
  TH2F *fdEdxVSp[2];                              //! dedx vs pTPC
  TProfile *hDeDxExp[9];                          //! TPC splines
  TH2F *fNsigmaTPC[2][18];                        //! NsigmaTPC vs. pT

  TH3F *fDca[2][3][18];                           //! DCAxy vs DCAz

  TH2F *fBetaTOFvspt[2];                          //! beta (TOF) vs pT
  TProfile *hBetaExp[9];                          //! TOF expected beta

  TH2F *fM2tof[2];                                //! m2 computed with TOF
  TH2F *fM2vspt[3][18];                           //! m2 computed with TOF (TPC cut)

  TF1 *fpTcorr[2];                                //! pT corrections for (anti-)deuteron

  void ForPtCorr(Double_t pt, Double_t t_pt, Int_t Pdg, Short_t t_charge, Bool_t isPrimary);
  void PtCorrection(Double_t &pt, Int_t FlagPid, Short_t charge);
  void CheckPtCorr(Double_t pt, Double_t t_pt, Int_t Pdg, Short_t t_charge);
  void FillDca(Double_t DCAxy, Double_t DCAz, Int_t Pdg, Short_t t_charge, Bool_t isPrimary, Bool_t isSecMat, Bool_t isSecWeak, Double_t t_pt, Bool_t kTOF);
    
  ClassDef(AliAnalysisMCNuclMult, 3);
};

#endif
