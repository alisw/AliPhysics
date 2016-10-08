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
//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;
//class AliVParticle;

class AliAnalysisNuclMult : public AliAnalysisTaskSE {
 
 public:
  
  AliAnalysisNuclMult();
  AliAnalysisNuclMult(const char *name);
  
  virtual ~AliAnalysisNuclMult();
  
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
  
  AliAnalysisNuclMult(const AliAnalysisNuclMult &old); 
  AliAnalysisNuclMult& operator=(const AliAnalysisNuclMult &source);

  AliESDEvent* fESD;                              //! ESD object
  AliVEvent* fEvent;                              //! general object
  AliESDtrackCuts *fESDtrackCuts;                 //  ESD track cuts object
  AliPPVsMultUtils *fPPVsMultUtils;               //  for multiplicity estimator
  AliPIDResponse *fPIDResponse;                   //! pointer to PID response

  TList *fList;                                   //! output container
  
  Float_t multMin;                                //  min. multiplicity accepted
  Float_t multMax;                                //  max. multiplicity accepted
  
  TH1I *htriggerMask;                             //! trigger mask
  TH1I *htriggerMask_noMB;                        //! trigger mask for no MB events
  TH1F *hzvertex;                                 //! z-vertex distribution
  TH1I *hNevent;                                  //! To check the event selection

  TH1I *hV0mult;                                  //! selected V0 multiplicity distribution
  TH1I *hTrackMult;                               //! number of ITS+TPC tracklets
  
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

  TF1 *fptCorr[1];                                //! pT corrections for (anti-)deuteron

  TH2F *fdEdxVSp[2];                              //! dedx vs pTPC
  TProfile *hDeDxExp[9];                          //! TPC splines
  TH2F *fNsigmaTPC[18];                           //! NsigmaTPC vs. pT

  TH3F *fDca[18];                                 //! DCAxy vs DCAz

  TH2F *fNsigmaTOF[18];                           //! NsigmaTOF vs. pT
  TH2F *fBetaTOFvspt[2];                          //! beta (TOF) vs pT
  TProfile *hBetaExp[9];                          //! TOF expected beta
  
  TH2F *fM2tof[2];                                //! m2 vs pT
  TH2F *fM2vspt[18];                              //! m2 vs pT (TPC cut)

  void EventSelectionMonitor();

  Bool_t IsInsideMultiplicityBin(Float_t multiplicity);
  
  Bool_t AcceptTrack(AliVTrack *track, Double_t &DCAxy, Double_t &DCAz);

  void TrackSelectionMonitor(Int_t nTPCclusters, Double_t chi2TPC, Bool_t isTPCrefit, Bool_t isITSrefit, Int_t nSPD, Int_t nKinkDaughters, Double_t chi2ITS, Bool_t isPropagatedToDca, Double_t DCAxy, Double_t DCAz, Double_t eta);

  void PtCorr(Double_t &pt, Double_t nsigmaTPC[9]);

  void GetNsigmaTPC(AliVTrack *track, Double_t nsigmaTPC[9]);

  Bool_t IsTOFmatching(AliVTrack *track);

  void GetExpTOFtime(AliVTrack *track, Double_t p, Double_t exptimes[9]);
  
  Double_t GetBeta(AliVTrack *track, Double_t pt);
  
  Double_t GetM2(Double_t p, Double_t beta);

  ClassDef(AliAnalysisNuclMult, 5);
};

#endif
