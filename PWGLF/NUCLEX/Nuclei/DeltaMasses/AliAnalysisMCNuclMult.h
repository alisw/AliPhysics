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
class TProfile;
class TFile;
class TList;
class TObject;
class TTree;
class AliAODEvent;
class AliESDEvent;
class AliVEvent;
class AliESDtrackCuts;
class AliPPVsMultUtils;
class AliPIDResponse;
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

  void SetESDtrackCutsObj(AliESDtrackCuts *esdTrackCuts) {fESDtrackCuts = esdTrackCuts;};
  void SetPPVsMultUtilsObj(AliPPVsMultUtils *fAliPPVsMultUtils) {fPPVsMultUtils = fAliPPVsMultUtils;};
  void SetMultiplicityRange(Float_t multiplicityMin=0, Float_t multiplicityMax=100) {
    multMin=multiplicityMin;
    multMax=multiplicityMax;
  };
  
 private:
  
  AliAnalysisMCNuclMult(const AliAnalysisMCNuclMult &old); 
  AliAnalysisMCNuclMult& operator=(const AliAnalysisMCNuclMult &source);

  AliAODEvent* fAOD;                              //! AOD object
  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliESDtrackCuts *fESDtrackCuts;                 //  ESD track cuts object
  AliPPVsMultUtils *fPPVsMultUtils;               //  for multiplicity estimator
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response

  TList *fList;                                   //! output container
  
  Float_t multMin;                                //  min. multiplicity accepted
  Float_t multMax;                                //  max. multiplicity accepted
  
  Int_t stdPdg[18];                               //! Pdg code of e,mu,pi,K,p,d,t,He3,He4
  
  TH1I *htriggerMask;                             //! trigger mask
  TH1I *htriggerMask_noMB;                        //! trigger mask for no MB events
  TH1F *hzvertex;                                 //! z-vertex distribution
  TH1I *hNevent;                                  //! To check the event selection

  TH1I *hV0mult;                                  //! selected V0 multiplicity distribution
  TH1I *hTrackMult;                               //! number of ITS+TPC tracklets
  
  TH2I *hpdg;                                     //! pdg label (after event selection)
  TH1I *htemp;                                    //! temporary
  
  TH1I *hCheckTrackSel;                           //! To check the track selection

  TH1I *hnTPCclusters[2];                         //! number of TPC clusters
  TH1D *hchi2TPC[2];                              //! chi2 per TPC cluster 
  TH1I *hisTPCrefit[2];                           //! if kTPCrefit
  TH1I *hisITSrefit[2];                           //! if kITSrefit
  TH1I *hnSPD[2];                                 //! number of SPD rec. points
  TH1I *hnKinkDaughters[2];                       //! number of kink daughters
  TH1D *hsigmaToVtx[2];                           //! number of sigma to the vertex
  TH1D *hchi2ITS[2];                              //! chi2 per ITS cluster 
  TH1D *heta[2];                                  //! eta dist.
  TH1I *hisPropagatedToDca[2];                    //! kPropagatedToDca
  
  TH2F *fdEdxVSp[2];                              //! dedx vs pTPC
  TProfile *hDeDxExp[9];                          //! TPC splines
  
  TH2F *fBetaTOFvspt[2];                          //! beta (TOF) vs pT
  TProfile *hBetaExp[9];                          //! TOF expected beta
  TH2F *fNsigmaTOF[3][2][18];                     //!

  TH1F *hpt[4][3][18];                            //! pT distributions
  
  TH2F *fptRecoVsTrue[2][18];                     //! pT reco vs. true
  TF1 *fptCorr[1];                                //! pT corrections for (anti-)deuteron

  TH3F *fDca[3][18];                              //! DCAxy vs DCAz
  
  void EventSelectionMonitor();
  Bool_t IsInsideMultiplicityBin(Float_t multiplicity);
  
  Bool_t AcceptTrack(AliVTrack *track, Double_t &DCAxy, Double_t &DCAz);
  void TrackSelectionMonitor(Int_t nTPCclusters, Double_t chi2TPC, Bool_t isTPCrefit, Bool_t isITSrefit, Int_t nSPD, Int_t nKinkDaughters, Double_t chi2ITS, Bool_t isPropagatedToDca, Double_t DCAxy, Double_t DCAz, Double_t eta);
  
  void FillDca(Double_t DCAxy, Double_t DCAz, Double_t t_pt, Int_t kSpec, Bool_t isPrimary, Bool_t isSecMat, Bool_t isSecWeak);
  
  void ForPtCorr(Double_t pt, Double_t t_pt, Int_t kSpec);
  void CheckPtCorr(Double_t pt, Double_t t_pt, Int_t kSpec);

  void CheckTPCsignal(AliVTrack *track);
  void CheckTOFsignal(AliVTrack *track, Double_t *NsigmaTOF);
  Bool_t IsTOFgoodmatching(AliVTrack *track, Int_t label, Double_t *NsigmaTOF, Int_t kSpec, Double_t t_pt, Bool_t isPrimary, Bool_t isSecMat, Bool_t isSecWeak);

  ClassDef(AliAnalysisMCNuclMult, 5);
};

#endif
