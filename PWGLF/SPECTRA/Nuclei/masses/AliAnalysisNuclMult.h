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
  void SetMultiplicityRange(Int_t multMin=-999, Int_t multMax=999) {//multMin=-999 -> analysis over all Minimum Bias collisions
    multiplicityMin=multMin;
    multiplicityMax=multMax;
  };
  
 private:
  
  AliAnalysisNuclMult(const AliAnalysisNuclMult &old); 
  AliAnalysisNuclMult& operator=(const AliAnalysisNuclMult &source);

  AliAODEvent* fAOD;                              //! AOD object
  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliAnalysisFilter *fTrackFilter;                //track filter object
  AliPPVsMultUtils *fPPVsMultUtils;               //for multiplicity estimators
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response
  
  TList *fList;                                   //! lists for slot
  
  Int_t multiplicityMin;
  Int_t multiplicityMax;

  TH1I *htriggerMask;                             //! trigger mask
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

  void FillDca(Double_t DCAxy, Double_t DCAz, AliVTrack *track);

  ClassDef(AliAnalysisNuclMult, 1);
};

#endif
