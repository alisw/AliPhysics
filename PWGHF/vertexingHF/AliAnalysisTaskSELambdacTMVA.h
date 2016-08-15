#ifndef ALIANALYSISTASKSELAMBDACTMVA_H
#define ALIANALYSISTASKSELAMBDACTMVA_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class Class AliAnalysisTaskSELambdacTMVA
/// \brief AliAnalysisTaskSE for the Lambdac candidates, for TMVA analysis,
/// and checks on MC generated and reconstructed Lambdac
///
/// Modified from AliAnalysisTaskSELambdac
/// \author Authors: Jaime Norman (jaime.norman@cern.ch)
/// \author          Marcel Figueredo (marcel.figueredo@cern.ch)
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
#include "AliVertexingHFUtils.h"

class AliAnalysisTaskSELambdacTMVA : public AliAnalysisTaskSE
{
	/// mimic steps in efficiency task for own steps in efficiency
	enum {
		kGeneratedLimAcc = 0,
		kGeneratedAll = 1,
		kGenerated = 2,
		kGeneratedAcc = 3,
		kReco3Prong = 4,
		kLcBit = 5,
		kIsSelectedTracks = 6,
		kIsInFidAcc	= 7,
		kPtRange = 8,
		kIsSelectedCandidate = 9,
		kIsSelectedPID = 10,
		kIsSelectedNtuple = 11 
	};	

 public:

  AliAnalysisTaskSELambdacTMVA();
  AliAnalysisTaskSELambdacTMVA(const char *name, Int_t fillNtuple,Int_t fillNtupleReco,AliRDHFCutsLctopKpi *lccutsana);
  virtual ~AliAnalysisTaskSELambdacTMVA();

	void SetKeepLcNotFromQuark(Bool_t keep = kTRUE) {fKeepLcNotFromQuark = keep;}
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMCPid(){fMCPid=kTRUE;fReadMC=kTRUE;fRealPid=kFALSE;fResPid=kFALSE;return;}
  void SetRealPid(){fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetResonantPid(){fResPid=kTRUE;fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetCutsKF(Float_t cutsKF[2]){for(Int_t i=0;i<2;i++){fCutsKF[i]=cutsKF[i];}return;}
  void SetUseKF(Bool_t useKF=kTRUE){fUseKF=useKF;}
	void SetIsHijing(Bool_t isHijing=kTRUE){fIsHijing=isHijing;}
	void SetKeepBkgNt(Bool_t keepBkgNt=kTRUE){fKeepBkgNt=keepBkgNt;}
  void SetAnalysis(Bool_t analysis=kTRUE){fAnalysis=analysis;}
  void SetUseFilterBitCut(Bool_t setter)   { fLcCut = setter;	 return; }  
  void SetUseFilterBitPID(Bool_t setter)   { fLcPIDCut = setter; return; }
	void SetCollisionSystem(Int_t syst) {fSyst = syst; return; }
  void SetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC, Bool_t &isInAcc) {fIsLcResonant=LambdacDaugh(part,arrayMC,isInAcc);} 
  void SetIsLcGen(AliAODMCParticle *partMC, TClonesArray *arrayMC);
  void SetIsLcReco(AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC);
	void SetUseNchWeight(Bool_t opt = kTRUE) {fUseNchWeight = opt;}
	void SetMCNchHisto(TH1F* h){
		if(fHistoMCNch) delete fHistoMCNch;
		fHistoMCNch=new TH1F(*h);
	}

  Bool_t GetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC) const {Bool_t dummy=kTRUE; return LambdacDaugh(part,arrayMC,dummy)>=1 ? kTRUE : kFALSE;} 
	Int_t GetPIDselectionMaxProb(AliAODRecoDecayHF3Prong *part); 
	Double_t GetNchWeight(Int_t nch);

  Bool_t IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Int_t MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Int_t LambdacDaugh(AliAODMCParticle *part,TClonesArray *arrayMC, Bool_t &isInAcc) const;
	void FillMassHists(AliAODEvent *aod,AliAODRecoDecayHF3Prong *d, TClonesArray *arrayMC, Int_t selection, Int_t selectionProb);
  void FillNtuple(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC, Int_t selection);
  void FillRecoNtuple(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC);
  void FillEffHists(Int_t kStep);
	void FillSelectionBits(AliAODRecoDecayHF3Prong *d, TH2F *hSelectionBits);

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSELambdacTMVA(const AliAnalysisTaskSELambdacTMVA &source);
  AliAnalysisTaskSELambdacTMVA& operator=(const AliAnalysisTaskSELambdacTMVA& source); 
  
  Bool_t ReconstructKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field) const;
  
  TList   *fOutput; //!<! list send on output slot 0
  TH1F    *fHistNEvents; //!<!hist. for No. of events
  TH1F    *fHistNEventsRejTM; //!<!hist. for Rejected events from null trigger mask 
  TH1F *fhSelectBit; //!<! hist for Filter Bit 	
	TH2F *fhSelectionBits;    //!<! hist for ALL Filter bits
	TH2F *fhSelectionBitsSigc; //!<! hist for ALL Filter bits Lc from c
	TH2F *fhSelectionBitsSigb; //!<! hist for ALL Filter bits Lc from b
	TH1F *fhSetIsLc; //!<! hist for before/after reco check MC LC
	TH2F *fhPIDmassLcPt; //!<!Lc Bkg+signal invariant mass vs pt
	TH2F *fhPIDmassLcPtSig; //!<!Lc signal invariant mass vs pt
	TH2F *fhPIDmassLcPtSigc; //!<!Lc from c signal invariant mass vs pt
	TH2F *fhPIDmassLcPtSigb; //!<!Lc from b signal invariant mass vs pt
	TH2F *fhMCmassLcPt; //!<!Lc Bkg+signal invariant mass vs pt
	TH2F *fhMCmassLcPtSig; //!<!Lc signal invariant mass vs pt
	TH2F *fhMCmassLcPtSigc; //!<!Lc from c signal invariant mass vs pt
	TH2F *fhMCmassLcPtSigb; //!<!Lc from b signal invariant mass vs pt
	TH2F *fhProbmassLcPt; //!<!Lc Bkg+signal invariant mass vs pt
	TH2F *fhProbmassLcPtSig; //!<!Lc signal invariant mass vs pt
	TH2F *fhProbmassLcPtSigc; //!<!Lc from c signal invariant mass vs pt
	TH2F *fhProbmassLcPtSigb; //!<!Lc from b signal invariant mass vs pt
	TH1F *fhIsLcResonantGen; //!<!hist for resonant flag gen
	TH1F *fhIsLcResonantReco;//!<!hist for resonant flag reco
	TH1F *fhIsLcGen;//!<!hist for resonant flag gen
	TH1F *fhIsLcReco;//!<!hist for resonant flag reco
	TH1F *fhRecoPDGmom; //!<!hist for Reco pdg
  TH1F *fhNBkgNI[12];	       //!<! hist. for n bkg, pT
  TH1F *fhNLc[12];             //!<! hist. for n Lc tot., pT
  TH1F *fhNLcc[12];            //!<! hist. for n Lc tot. from c, pT
  TH1F *fhNLcNonRc[12];        //!<! hist. for n Lc from c non resonant, pT
  TH1F *fhNLcL1520c[12];       //!<! hist. for n Lc from c L1520 + pi, pT
  TH1F *fhNLcKstarc[12];       //!<! hist. for n Lc from c K* + p, pT
  TH1F *fhNLcDeltac[12];       //!<! hist. for n Lc from c Delta++ + K, pT
  TH1F *fhNLcb[12];            //!<! hist. for n Lc tot. from b, pT
  TH1F *fhNLcNonRb[12];        //!<! hist. for n Lc from b non resonant, pT
  TH1F *fhNLcL1520b[12];       //!<! hist. for n Lc from b L1520 + pi, pT
  TH1F *fhNLcKstarb[12];       //!<! hist. for n Lc from b K* + p, pT
  TH1F *fhNLcDeltab[12];       //!<! hist. for n Lc from b Delta++ + K, pT
  TH2F *fhPtEtaBkgNI[12];      //!<! hist. for n bkg, pT vs eta                    
  TH2F *fhPtEtaLc[12];         //!<! hist. for n Lc tot., pT vs eta
  TH2F *fhPtEtaLcc[12];        //!<! hist. for n Lc tot. from c, pT vs eta
  TH2F *fhPtEtaLcNonRc[12];    //!<! hist. for n Lc from c non resonant, pT vs eta
  TH2F *fhPtEtaLcL1520c[12];   //!<! hist. for n Lc from c L1520 + pi, pT vs eta
  TH2F *fhPtEtaLcKstarc[12];   //!<! hist. for n Lc from c K* + p, pT vs eta
  TH2F *fhPtEtaLcDeltac[12];   //!<! hist. for n Lc from c Delta++ + K, pT vs eta
  TH2F *fhPtEtaLcb[12];        //!<! hist. for n Lc tot. from b, pT vs eta
  TH2F *fhPtEtaLcNonRb[12];    //!<! hist. for n Lc from b non resonant, pT vs eta
  TH2F *fhPtEtaLcL1520b[12];   //!<! hist. for n Lc from b L1520 + pi, pT vs eta
  TH2F *fhPtEtaLcKstarb[12];   //!<! hist. for n Lc from b K* + p, pT vs eta
  TH2F *fhPtEtaLcDeltab[12];   //!<! hist. for n Lc from b Delta++ + K, pT vs eta
  TH2F *fhPtYBkgNI[12];        //!<! hist. for n bkg, pT vs rapidity                   
  TH2F *fhPtYLc[12];           //!<! hist. for n Lc tot., pT vs rapidity
  TH2F *fhPtYLcc[12];          //!<! hist. for n Lc tot. from c, pT vs rapidity
  TH2F *fhPtYLcNonRc[12];      //!<! hist. for n Lc from c non resonant, pT vs rapidity
  TH2F *fhPtYLcL1520c[12];     //!<! hist. for n Lc from c L1520 + pi, pT vs rapidity
  TH2F *fhPtYLcKstarc[12];     //!<! hist. for n Lc from c K* + p, pT vs rapidity
  TH2F *fhPtYLcDeltac[12];     //!<! hist. for n Lc from c Delta++ + K, pT vs rapidity
  TH2F *fhPtYLcb[12];          //!<! hist. for n Lc tot. from b, pT vs rapidity
  TH2F *fhPtYLcNonRb[12];      //!<! hist. for n Lc from b non resonant, pT vs rapidity
  TH2F *fhPtYLcL1520b[12];     //!<! hist. for n Lc from b L1520 + pi, pT vs rapidity
  TH2F *fhPtYLcKstarb[12];     //!<! hist. for n Lc from b K* + p, pT vs rapidity
  TH2F *fhPtYLcDeltab[12];     //!<! hist. for n Lc from b Delta++ + K, pT vs rapidity
  TH2F *fhPtPhiBkgNI[12];      //!<! hist. for n bkg, pT vs phi                   
  TH2F *fhPtPhiLc[12];         //!<! hist. for n Lc tot., pT vs phi
  TH2F *fhPtPhiLcc[12];        //!<! hist. for n Lc tot. from c, pT vs phi
  TH2F *fhPtPhiLcNonRc[12];    //!<! hist. for n Lc from c non resonant, pT vs phi
  TH2F *fhPtPhiLcL1520c[12];   //!<! hist. for n Lc from c L1520 + pi, pT vs phi
  TH2F *fhPtPhiLcKstarc[12];   //!<! hist. for n Lc from c K* + p, pT vs phi
  TH2F *fhPtPhiLcDeltac[12];   //!<! hist. for n Lc from c Delta++ + K, pT vs phi
  TH2F *fhPtPhiLcb[12];        //!<! hist. for n Lc tot. from b, pT vs phi
  TH2F *fhPtPhiLcNonRb[12];    //!<! hist. for n Lc from b non resonant, pT vs phi
  TH2F *fhPtPhiLcL1520b[12];   //!<! hist. for n Lc from b L1520 + pi, pT vs phi
  TH2F *fhPtPhiLcKstarb[12];   //!<! hist. for n Lc from b K* + p, pT vs phi
  TH2F *fhPtPhiLcDeltab[12];   //!<! hist. for n Lc from b Delta++ + K, pT vs phi
	TH1F *fhPtMisIdpKpi; //!<! hist for pt pKpi signal mis id'd as piKp
	TH1F *fhPtMisIdpiKp; //!<! hist for pt pKpi signal mis id'd as piKp
	TH1F *fhPtCorrId; //!<! hist for correctly id'd pKpi
	TH2F *fhInvMassMisIdpKpi; //!<! hist for inv mass pKpi signal mis id'd as piKp
	TH2F *fhInvMassMisIdpiKp; //!<! hist for inv mass pKpi signal mis id'd as piKp
	TH1F *fhPtMisIdpKpiProb; //!<! hist for pt pKpi signal mis id'd as piKp most prob PID
	TH1F *fhPtMisIdpiKpProb; //!<! hist for pt pKpi signal mis id'd as piKp most prob PID
	TH1F *fhPtCorrIdProb; //!<! hist for correctly id'd pKpi most prob PID
	TH2F *fhInvMassMisIdpKpiProb; //!<! hist for inv mass pKpi signal mis id'd as piKp most prob PID
	TH2F *fhInvMassMisIdpiKpProb; //!<! hist for inv mass pKpi signal mis id'd as piKp most prob PID
  TNtuple *fNtupleLambdac; //!<! output ntuple
  TNtuple *fNtupleLambdacReco; //!<! output ntuple after reconstruction
	TF1 *fFuncWeightPythia; //!<! weight function for Pythia vs pPb prod.
	TF1 *fFuncWeightFONLL7overLHC10f6a; //!<! weight function for FONLL vs p prod.
	TF1 *fFuncWeightFONLL5overLHC13d3; //!<! weight function for FONLL vs pPb prod.
	TF1 *fFuncWeightFONLL5overLHC10f6a; //!<! weight function for FONLL vs p prod.
	TF1 *fFuncWeightFONLL5overLHC13d3Lc; //!<! weight function for FONLL vs pPb prod. Lc
	TF1 *fFuncWeightFONLL7overLHC11b2Lc; //!<! weight function for FONLL vs p prod. Lc
	TF1 *fFuncWeightFONLL7overLHC10f7aLc; //!<! weight function for FONLL vs p prod. Lc
	Bool_t fUseNchWeight; /// flag for using multiplicity weights
	TH1F *fHistoMCNch; /// multiplicity weight histogram
  Float_t fCutsKF[2]; /// cuts with KF vertexer
  Int_t fIsLc; /// is MC Lc - 0=not Lc, 1=Lc from c, 2=Lc from b
  Int_t fIsLcResonant; /// is Lc resonant - 1=non resonant, 2=via L1520 + pi, 3=via K* + p, 4=via Delta++ + K
  Float_t fCandidateVars[4]; /// candidate variables, 0=Pt, 1=Eta, 2=Y, 3=Phi
  Float_t fPtLc; /// pt of Lc candidate
  Float_t fUpmasslimit;  /// upper inv mass limit for histos
  Float_t fLowmasslimit; /// lower inv mass limit for histos
  AliRDHFCutsLctopKpi *fRDCutsAnalysis; /// Analysis cuts
  TList *fListCuts; /// list of cuts
  Int_t fFillNtuple;   ///  filling ntuple type
  Int_t fFillNtupleReco;   ///  filling ntuple type reco
	Bool_t fKeepLcNotFromQuark; /// flag to keep Lc not from quark
	Bool_t fKeepBkgNt; /// flag to keep background in 
	Int_t fSyst; /// flag for collision system. 0=pp, 1=PbPb, 2=pPb
  Bool_t fReadMC;    /// flag for access to MC
  Bool_t fMCPid;    /// flag for access to MC
  Bool_t fRealPid;    /// flag for real PID
  Bool_t fResPid;      /// flag for PID with resonant channels
  Bool_t fUseKF;      /// flag to cut with KF vertexer
  Bool_t fAnalysis;      /// apply analysis cuts
  AliAnalysisVertexingHF *fVHF;  ///  Vertexer heavy flavour (used to pass the cuts)
  Bool_t fLcCut;  ///  flag for Lc filter bit cut
  Bool_t fLcPIDCut;  /// flag for Lc filter bit PID
  Bool_t fIsHijing; /// flag for whether Lc is from Hijing 
  TH1F *fNentries;      /// histo with number of entries
  //AliAODpidUtil* fUtilPid;
  AliPIDResponse *fPIDResponse;     //!<! PID response object
  AliNormalizationCounter *fCounter;//!<!AliNormalizationCounter on output slot 7
	AliVertexingHFUtils *fVertUtil;         /// vertexing HF Util

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskSELambdacTMVA,11); /// AliAnalysisTaskSE for the invariant mass analysis of heavy-flavour decay candidates (Lambdac)
  /// \endcond
};

#endif

