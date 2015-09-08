#ifndef AliAnalysisTaskSELambdac_H
#define AliAnalysisTaskSELambdac_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskSELambdac
/// \brief AliAnalysisTaskSE for the Lambdac candidates Invariant Mass Histogram and
/// comparison of heavy-flavour decay candidates
/// to MC truth (kinematics stored in the AOD)
/// \author r.romita@gsi.de
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TArrayD.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCuts.h"
#include "TClonesArray.h"
//#include "AliAODpidUtil.h"
#include "AliPIDResponse.h"
#include "AliNormalizationCounter.h"

class AliAnalysisTaskSELambdac : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSELambdac();
  AliAnalysisTaskSELambdac(const char *name, Bool_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana, AliRDHFCutsLctopKpi *lccutsprod);
  virtual ~AliAnalysisTaskSELambdac();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMCPid(){fMCPid=kTRUE;fReadMC=kTRUE;fRealPid=kFALSE;fResPid=kFALSE;return;}
  void SetRealPid(){fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetResonantPid(){fResPid=kTRUE;fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetCutsKF(Float_t cutsKF[2]){for(Int_t i=0;i<2;i++){fCutsKF[i]=cutsKF[i];}return;}
  void SetUseKF(Bool_t useKF=kTRUE){fUseKF=useKF;}
  void SetAnalysis(Bool_t analysis=kTRUE){fAnalysis=analysis;}
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetPtBinLimit(Int_t n, Float_t *limitarray);
  void SetFillVarHists(Bool_t setter) {fFillVarHists=setter;return;}
  void SetMultiplicityHists(Bool_t setter) {fMultiplicityHists=setter;return;}
  void SetPriorsHists(Bool_t setter) {fPriorsHists=setter;return;}
  void SetUseFilterBitCut(Bool_t setter)   { fLcCut = setter;	 return; }  
  void SetUseFilterBitPID(Bool_t setter)   { fLcPIDCut = setter; return; }

  Float_t GetUpperMassLimit() const {return fUpmasslimit;}
  Float_t GetLowerMassLimit() const {return fLowmasslimit;}
  Int_t GetNBinsPt() const {return fNPtBins;}
  Double_t GetPtBinLimit(Int_t ibin) const ;
  Bool_t IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  void IspiKpResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const ;
  void IspKpiResonant(AliAODRecoDecayHF3Prong *d,Double_t field,Int_t *resNumber) const ;
  Bool_t VertexingKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const ;
  Int_t MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t GetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC) const ;

  void FillMassHists(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC, AliRDHFCutsLctopKpi *cuts,Int_t *nSelectedloose,Int_t *nSelectedtight);
  void FillVarHists(AliAODRecoDecayHF3Prong *part, TClonesArray *arrMC, AliRDHFCutsLctopKpi *cuts, /*TList *listout,*/ AliAODEvent *aod);
  Bool_t Is3ProngFromPDG(AliAODRecoDecayHF3Prong *part, TClonesArray *arrMC, Int_t pdgToBeCompared=4);
  Bool_t IsTrackFromPDG(const AliAODTrack *daugh, TClonesArray *arrayMC, Int_t pdgToBeCompared);
  Bool_t IsThereAGeneratedLc(TClonesArray *arrayMC);
  Int_t NumberPrimaries(const AliAODEvent *aods);
  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSELambdac(const AliAnalysisTaskSELambdac &source);
  AliAnalysisTaskSELambdac& operator=(const AliAnalysisTaskSELambdac& source); 
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*6;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*6+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*6+2;}
  Int_t GetLbHistoIndex(Int_t iPtBin) const { return iPtBin*6+3;}
  Int_t GetcOnlyHistoIndex(Int_t iPtBin) const { return iPtBin*6+4;}
  Int_t GetNoQuarkHistoIndex(Int_t iPtBin) const { return iPtBin*6+5;}
  //  Int_t GetLSHistoIndex(Int_t iPtBin)const { return iPtBin*7;}

  Bool_t ReconstructKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const;
  void FillAPrioriConcentrations(AliAODRecoDecayHF3Prong *part, AliRDHFCutsLctopKpi *cuts,
				 AliAODEvent* aod, TClonesArray *arrMC);
  void MultiplicityStudies(AliAODRecoDecayHF3Prong *part, AliRDHFCutsLctopKpi *cuts,
			   AliAODEvent* aod, TClonesArray *arrMC,
			   Bool_t &flag1,Bool_t &flag2,Bool_t &flag3,
			   Bool_t &flag4, Bool_t &flag5, Bool_t &flag6); 
  enum {kMaxPtBins=10};

  TList   *fOutput; //!<! list send on output slot 0
  TH1F    *fHistNEvents; //!<!hist. for No. of events
  TH1F    *fhChi2; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3TC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3Kp; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3KpTC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3Lpi; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3LpiTC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3Dk; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater3DkTC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater33Pr; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater33PrTC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2TC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2Kp; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2KpTC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2Lpi; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2LpiTC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2Dk; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater2DkTC; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater23Pr; //!<!hist. for No. of events
  TH1F    *fhMassPtGreater23PrTC; //!<!hist. for No. of events
  TH1F *fMassHist[6*kMaxPtBins]; //!<!hist. for inv mass (LC)
  TH1F *fMassHistTC[6*kMaxPtBins]; //!<!hist. for inv mass (TC)
  TH1F *fMassHistLpi[3*kMaxPtBins]; //!<!hist. for inv mass (LC)
  TH1F *fMassHistLpiTC[3*kMaxPtBins]; //!<!hist. for inv mass (TC)
  TH1F *fMassHistKp[3*kMaxPtBins]; //!<!hist. for inv mass (LC)
  TH1F *fMassHistKpTC[3*kMaxPtBins]; //!<!hist. for inv mass (TC)
  TH1F *fMassHistDk[3*kMaxPtBins]; //!<!hist. for inv mass (LC)
  TH1F *fMassHistDkTC[3*kMaxPtBins]; //!<!hist. for inv mass (TC)
  TH1F *fMassHist3Pr[3*kMaxPtBins]; //!<!hist. for inv mass (LC)
  TH1F *fMassHist3PrTC[3*kMaxPtBins]; //!<!hist. for inv mass (TC)
  TH2F *fhMassLcPt; //!<!hist Mass x pt
  TH2F *fhMassLcplusPt; //!<!hist Mass x pt Lcplus
  TH2F *fhMassLcminusPt; //!<!hist Mass x pt Lcminu
  TH2F *fhEta3Prong;    //!<!hist. for 3-prong Eta
  TH2F *fhEta3ProngAcc; //!<!hist. for 3-prong Eta fiducial acc
  TH2F *fhEta3ProngProd; //!<!hist. for 3-prong Eta fiducial  Prod Cuts
  TH2F *fhEta3ProngAn; //!<!hist. for 3-prong Eta fiducial An Cuts
  TH2F *fhRap3Prong;    //!<!hist. for 3-prong y 
  TH2F *fhRap3ProngAcc; //!<!hist. for 3-prong Eta fiducial acc
  TH2F *fhRap3ProngProd; //!<!hist. for 3-prong Eta fiducial Prod cuts
  TH2F *fhRap3ProngAn; //!<!hist. for 3-prong Eta fiducial An cuts
  TH1F *fhSelectBit; //!<! hist for Filter Bit 	
  TH2F *fhProtonPtProngLcPt; //!<!hist for var_LcPt
  TH2F *fhBProtonPtProngLcPt;//!<!hist for var_LcPt
  TH2F *fhProtond0ProngLcPt;//!<!hist for var_LcPt
  TH2F *fhBProtond0ProngLcPt;//!<!hist for var_LcPt
  TH2F *fhKaonPtProngLcPt;//!<!hist for var_LcPt
  TH2F *fhBKaonPtProngLcPt;//!<!hist for var_LcPt
  TH2F *fhKaond0ProngLcPt;//!<!hist for var_LcPt
  TH2F *fhBKaond0ProngLcPt;//!<!hist for var_LcPt
  TH2F *fhPionPtProngLcPt;//!<!hist for var_LcPt
  TH2F *fhBPionPtProngLcPt;//!<!hist for var_LcPt
  TH2F *fhPiond0ProngLcPt;//!<!hist for var_LcPt
  TH2F *fhBPiond0ProngLcPt;//!<!hist for var_LcPt
  TH2F *fhDist12PrimLcPt;//!<!hist for var_LcPt
  TH2F *fhBDist12PrimLcPt;//!<!hist for var_LcPt
  TH2F *fhSigmaVertLcPt;//!<!hist for var_LcPt
  TH2F *fhBSigmaVertLcPt;//!<!hist for var_LcPt
  TH2F *fhdcasLcPt;//!<!hist for var_LcPt
  TH2F *fhBdcasLcPt;//!<!hist for var_LcPt
  TH2F *fhCosPointingAngleLcPt;//!<!hist for var_LcPt
  TH2F *fhBCosPointingAngleLcPt;//!<!hist for var_LcPt
  TH2F *fhDecayLengthLcPt;//!<!hist for var_LcPt
  TH2F *fhBDecayLengthLcPt;//!<!hist for var_LcPt
  TH2F *fhSum2LcPt;//!<!hist for var_LcPt
  TH2F *fhBSum2LcPt;//!<!hist for var_LcPt
  TH2F *fhPtMaxLcPt;//!<!hist for var_LcPt
  TH2F *fhBPtMaxLcPt;//!<!hist for var_LcPt
  TNtuple *fNtupleLambdac; //!<! output ntuple
  Float_t fUpmasslimit;  /// upper inv mass limit for histos
  Float_t fLowmasslimit; /// lower inv mass limit for histos
  Float_t fCutsKF[2]; /// cuts with KF vertexer
  Int_t fNPtBins; /// number of bins in Pt for histograms
  AliRDHFCutsLctopKpi *fRDCutsAnalysis; /// Cuts for Analysis
  AliRDHFCutsLctopKpi *fRDCutsProduction; /// Production Cuts
  TList *fListCuts; /// list of cuts
  Double_t fArrayBinLimits[kMaxPtBins+1]; /// limits for the Pt bins
  Bool_t fFillNtuple;   /// flag for filling ntuple
  Bool_t fReadMC;    /// flag for access to MC
  Bool_t fMCPid;    /// flag for access to MC
  Bool_t fRealPid;    /// flag for real PID
  Bool_t fResPid;      /// flag for PID with resonant channels
  Bool_t fUseKF;      /// flag to cut with KF vertexer
  Bool_t fAnalysis;      /// apply analysis cuts
  AliAnalysisVertexingHF *fVHF;  /// Vertexer heavy flavour (used to pass the cuts)
  Bool_t fFillVarHists;  /// flag for creation and fill of histograms with vars
  Bool_t fMultiplicityHists;  /// flag for activation of multiplcity histos
  Bool_t fPriorsHists;  /// flag for histos with priors
  Bool_t fLcCut;  /// flag for Lc filter bit cut
  Bool_t fLcPIDCut;  /// flag for Lc filter bit PID
  TH1F *fNentries;      /// histo with number of entries
  TList *fOutputMC;     /// output1
  TList *fAPriori;      /// output2
  TList *fMultiplicity; /// output3
  //AliAODpidUtil* fUtilPid;
  AliPIDResponse *fPIDResponse;     //!<! PID response object
  AliNormalizationCounter *fCounter;//!<!AliNormalizationCounter on output slot 7

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSELambdac,8); /// AliAnalysisTaskSE for the invariant mass analysis of heavy-flavour decay candidates (Lambdac)
  /// \endcond
};

#endif

