#ifndef ALIANALYSISNUCLMULT_H
#define ALIANALYSISNUCLMULT_H

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
class TProfile;
class TFile;
class TList;
class TObject;
class TTree;
class AliESDEvent;
class AliVEvent;
class AliESDtrackCuts;
class AliPPVsMultUtils;
class AliPIDResponse;
// for Monte Carlo:
class AliMCEventHandler;
class AliMCEvent;
class AliStack;
//class AliVParticle;

class AliAnalysisNuclMult : public AliAnalysisTaskSE {
 
 public:
  
  AliAnalysisNuclMult();
  AliAnalysisNuclMult(const char *name);
  
  virtual ~AliAnalysisNuclMult();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetIsMC(Bool_t IsMC=kFALSE) {isMC = IsMC;};

  void SetESDtrackCutsObj(AliESDtrackCuts *esdTrackCuts) {fESDtrackCuts = esdTrackCuts;};
  
  void SetPPVsMultUtilsObj(AliPPVsMultUtils *fAliPPVsMultUtils) {fPPVsMultUtils = fAliPPVsMultUtils;};
  
  void SetMultiplicityRange(Float_t multiplicityMin=0, Float_t multiplicityMax=100) {
    multMin=multiplicityMin;
    multMax=multiplicityMax;
  };
  
  void SetDCAzMax(Float_t max=1.) {DCAzMax = max;};

 private:
  
  AliAnalysisNuclMult(const AliAnalysisNuclMult &old); 
  AliAnalysisNuclMult& operator=(const AliAnalysisNuclMult &source);

  Bool_t isMC;                                    //  Are they real data or MC?

  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliMCEventHandler* eventHandler;                //! MCEventHandler
  AliMCEvent* mcEvent;                            //! MCEvent
  AliStack* fStack;                               //! Stack
  
  AliESDtrackCuts *fESDtrackCuts;                 //  ESD track cuts object
  AliPPVsMultUtils *fPPVsMultUtils;               //  for multiplicity estimator
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response

  TList *fList;                                   //! output container

  Int_t stdPdg[18];                               //! Pdg code of e,mu,pi,K,p,d,t,He3,He4
  
  Float_t multMin;                                //  min. multiplicity accepted
  Float_t multMax;                                //  max. multiplicity accepted
  
  Float_t DCAzMax;                                //  DCAz max 

  TH1F *htriggerMask[2];                          //! trigger mask
  TH1F *hzvertex;                                 //! z-vertex distribution
  TH1F *hNevent;                                  //! To check the event selection

  TH1F *hV0mult;                                  //! selected V0 multiplicity distribution
  TH1F *hTracklets;                               //! 
  TH2F *hTrackletsVsV0mult;                       //!

  TH1F *hrapidity[2];                             //!

  TH1F *hCheckTrackSel;                           //! To check the track selection

  TH1F *hnTPCclusters[2];                         //! number of TPC clusters
  TH1F *hchi2TPC[2];                              //! chi2 per TPC cluster 
  TH1F *hisTPCrefit[2];                           //! if kTPCrefit
  TH1F *hisITSrefit[2];                           //! if kITSrefit
  TH1F *hnSPD[2];                                 //! number of SPD rec. points
  TH1F *hnKinkDaughters[2];                       //! number of kink daughters
  TH1F *hsigmaToVtx[2];                           //! number of sigma to the vertex
  TH1F *hchi2ITS[2];                              //! chi2 per ITS cluster 
  TH1F *heta[2];                                  //! eta dist.
  TH1F *hisPropagatedToDca[2];                    //! kPropagatedToDca

  TF1 *fptCorr[1];                                //! pT corrections for (anti-)deuteron

  TH2F *fdEdxVSp[2];                              //! dedx vs pTPC
  TProfile *hDeDxExp[9];                          //! TPC splines
  TH3F *fNsigmaTPC[18];                           //! NsigmaTPC vs. pT

  TH3F *fDca[2][18];                              //! DCAxy, DCAz

  TH2F *fNsigmaTOF[18];                           //! NsigmaTOF vs. pT
  TH2F *fBetaTOFvspt[2];                          //! beta (TOF) vs pT
  TProfile *hBetaExp[9];                          //! TOF expected beta
  
  TH3F *fM2tof[2];                                //! M2 vs pT
  TH3F *fM2vspt[18];                              //! M2 vs pT (TPC cut)

  TH3F *fNsigmas[18];                             //!

  //Only for MC:
  TH2F *hpdg[1];                                  //! pdg label (after event selection)

  TH2F *hpt[4][3][18];                            //! pT distributions

  TH2F *fptRecoVsTrue[2][18];                     //! pT reco vs. true

  TH3F *fmcDca[3][2][18];                         //! DCAxy, DCAz

  TH2F *fmcNsigmaTOF[3][2][18];                   //!

  void EventSelectionMonitor();

  Bool_t IsInsideMultiplicityBin(Float_t multiplicity);
  
  Double_t GetRapidity(AliVTrack *track);

  Bool_t AcceptTrack(AliVTrack *track, Double_t &DCAxy, Double_t &DCAz);

  void TrackSelectionMonitor(Int_t nTPCclusters, Double_t chi2TPC, Bool_t isTPCrefit, Bool_t isITSrefit, Int_t nSPD, Int_t nKinkDaughters, Double_t chi2ITS, Bool_t isPropagatedToDca, Double_t DCAxy, Double_t DCAz, Double_t eta);

  void PtCorr(Double_t &pt, Double_t nsigmaTPC[9]);

  void GetNsigmaTPC(AliVTrack *track, Double_t nsigmaTPC[9]);

  Bool_t IsTOFmatching(AliVTrack *track);

  void GetExpTOFtime(AliVTrack *track, Double_t p, Double_t exptimes[9]);
  
  Double_t GetBeta(AliVTrack *track, Double_t pt, Double_t nsigmaTOF[9]);
  
  Double_t GetM2(Double_t p, Double_t beta);

  //Methods called only on MC analysis:
  void ForPtCorr(Double_t pt, Double_t t_pt, Int_t kSpec);
  
  Bool_t IsTOFgoodmatching(AliVTrack *track, Int_t label, Double_t nsigmaTOF[9], Int_t kSpec, Double_t t_pt, Bool_t isPrimary, Bool_t isSecMat, Bool_t isSecWeak);
  //---
  
  ClassDef(AliAnalysisNuclMult, 7);
};

#endif
