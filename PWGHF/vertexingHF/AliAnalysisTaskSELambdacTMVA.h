#ifndef ALIANALYSISTASKSELAMBDACTMVA_H
#define ALIANALYSISTASKSELAMBDACTMVA_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// Class AliAnalysisTaskSELambdacTMVA
// AliAnalysisTaskSE for the Lambdac candidates, for TMVA analysis,
// and checks on MC generated and reconstructed Lambdac
// 
// Modified from AliAnalysisTaskSELambdac
// Authors: Jaime Norman (jaime.norman@cern.ch)
//          Marcel Figueredo (marcel.figueredo@cern.ch)
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

class AliAnalysisTaskSELambdacTMVA : public AliAnalysisTaskSE
{
	//mimic steps in efficiency task for own steps in efficiency
	enum {
		kGeneratedLimAcc = 0,
		kGenerated = 1,
		kGeneratedAcc = 2,
		kReco3Prong = 3,
		kLcBit = 4,
		kIsSelectedTracks = 5,
		kIsInFidAcc	= 6,
		kPtRange = 7,
		kIsSelectedCandidate = 8,
		kIsSelectedPID = 9,
		kIsSelectedNtuple = 10
	};	

 public:

  AliAnalysisTaskSELambdacTMVA();
  AliAnalysisTaskSELambdacTMVA(const char *name, Int_t fillNtuple,AliRDHFCutsLctopKpi *lccutsana);
  virtual ~AliAnalysisTaskSELambdacTMVA();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMCPid(){fMCPid=kTRUE;fReadMC=kTRUE;fRealPid=kFALSE;fResPid=kFALSE;return;}
  void SetRealPid(){fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetResonantPid(){fResPid=kTRUE;fRealPid=kTRUE;fMCPid=kFALSE;return;}
  void SetCutsKF(Float_t cutsKF[2]){for(Int_t i=0;i<2;i++){fCutsKF[i]=cutsKF[i];}return;}
  void SetUseKF(Bool_t useKF=kTRUE){fUseKF=useKF;}
  void SetAnalysis(Bool_t analysis=kTRUE){fAnalysis=analysis;}
  void SetUseFilterBitCut(Bool_t setter)   { fLcCut = setter;	 return; }  
  void SetUseFilterBitPID(Bool_t setter)   { fLcPIDCut = setter; return; }
  void SetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC, Bool_t &isInAcc) {fIsLcResonant=LambdacDaugh(part,arrayMC,isInAcc);} 
  void SetIsLc(AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC);

  Bool_t GetLambdacDaugh(AliAODMCParticle *part, TClonesArray *arrayMC) const {Bool_t dummy=kTRUE; return LambdacDaugh(part,arrayMC,dummy)>=1 ? kTRUE : kFALSE;} 

  Bool_t IspiKpMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Bool_t IspKpiMC(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Int_t MatchToMCLambdac(AliAODRecoDecayHF3Prong *d,TClonesArray *arrayMC) const ;
  Int_t LambdacDaugh(AliAODMCParticle *part,TClonesArray *arrayMC, Bool_t &isInAcc) const;
  void FillNtuple(AliAODEvent *aod,AliAODRecoDecayHF3Prong *part, TClonesArray *arrayMC, Int_t selection);
  void FillEffHists(Int_t kStep);

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
  
  TList   *fOutput; //! list send on output slot 0
  TH1F    *fHistNEvents; //!hist. for No. of events
  TH1F    *fHistNEventsRejTM; //!hist. for Rejected events from null trigger mask 
  TH1F *fhSelectBit; //! hist for Filter Bit 	
  TH1F *fhNBkgNI[11];	       //! hist. for n bkg, pT
  TH1F *fhNLc[11];             //! hist. for n Lc tot., pT
  TH1F *fhNLcc[11];            //! hist. for n Lc tot. from c, pT
  TH1F *fhNLcNonRc[11];        //! hist. for n Lc from c non resonant, pT
  TH1F *fhNLcL1520c[11];       //! hist. for n Lc from c L1520 + pi, pT
  TH1F *fhNLcKstarc[11];       //! hist. for n Lc from c K* + p, pT
  TH1F *fhNLcDeltac[11];       //! hist. for n Lc from c Delta++ + K, pT
  TH1F *fhNLcb[11];            //! hist. for n Lc tot. from b, pT
  TH1F *fhNLcNonRb[11];        //! hist. for n Lc from b non resonant, pT
  TH1F *fhNLcL1520b[11];       //! hist. for n Lc from b L1520 + pi, pT
  TH1F *fhNLcKstarb[11];       //! hist. for n Lc from b K* + p, pT
  TH1F *fhNLcDeltab[11];       //! hist. for n Lc from b Delta++ + K, pT
  TH2F *fhPtEtaBkgNI[11];      //! hist. for n bkg, pT vs eta                    
  TH2F *fhPtEtaLc[11];         //! hist. for n Lc tot., pT vs eta
  TH2F *fhPtEtaLcc[11];        //! hist. for n Lc tot. from c, pT vs eta
  TH2F *fhPtEtaLcNonRc[11];    //! hist. for n Lc from c non resonant, pT vs eta
  TH2F *fhPtEtaLcL1520c[11];   //! hist. for n Lc from c L1520 + pi, pT vs eta
  TH2F *fhPtEtaLcKstarc[11];   //! hist. for n Lc from c K* + p, pT vs eta
  TH2F *fhPtEtaLcDeltac[11];   //! hist. for n Lc from c Delta++ + K, pT vs eta
  TH2F *fhPtEtaLcb[11];        //! hist. for n Lc tot. from b, pT vs eta
  TH2F *fhPtEtaLcNonRb[11];    //! hist. for n Lc from b non resonant, pT vs eta
  TH2F *fhPtEtaLcL1520b[11];   //! hist. for n Lc from b L1520 + pi, pT vs eta
  TH2F *fhPtEtaLcKstarb[11];   //! hist. for n Lc from b K* + p, pT vs eta
  TH2F *fhPtEtaLcDeltab[11];   //! hist. for n Lc from b Delta++ + K, pT vs eta
  TH2F *fhPtYBkgNI[11];        //! hist. for n bkg, pT vs rapidity                   
  TH2F *fhPtYLc[11];           //! hist. for n Lc tot., pT vs rapidity
  TH2F *fhPtYLcc[11];          //! hist. for n Lc tot. from c, pT vs rapidity
  TH2F *fhPtYLcNonRc[11];      //! hist. for n Lc from c non resonant, pT vs rapidity
  TH2F *fhPtYLcL1520c[11];     //! hist. for n Lc from c L1520 + pi, pT vs rapidity
  TH2F *fhPtYLcKstarc[11];     //! hist. for n Lc from c K* + p, pT vs rapidity
  TH2F *fhPtYLcDeltac[11];     //! hist. for n Lc from c Delta++ + K, pT vs rapidity
  TH2F *fhPtYLcb[11];          //! hist. for n Lc tot. from b, pT vs rapidity
  TH2F *fhPtYLcNonRb[11];      //! hist. for n Lc from b non resonant, pT vs rapidity
  TH2F *fhPtYLcL1520b[11];     //! hist. for n Lc from b L1520 + pi, pT vs rapidity
  TH2F *fhPtYLcKstarb[11];     //! hist. for n Lc from b K* + p, pT vs rapidity
  TH2F *fhPtYLcDeltab[11];     //! hist. for n Lc from b Delta++ + K, pT vs rapidity
  TH2F *fhPtPhiBkgNI[11];      //! hist. for n bkg, pT vs phi                   
  TH2F *fhPtPhiLc[11];         //! hist. for n Lc tot., pT vs phi
  TH2F *fhPtPhiLcc[11];        //! hist. for n Lc tot. from c, pT vs phi
  TH2F *fhPtPhiLcNonRc[11];    //! hist. for n Lc from c non resonant, pT vs phi
  TH2F *fhPtPhiLcL1520c[11];   //! hist. for n Lc from c L1520 + pi, pT vs phi
  TH2F *fhPtPhiLcKstarc[11];   //! hist. for n Lc from c K* + p, pT vs phi
  TH2F *fhPtPhiLcDeltac[11];   //! hist. for n Lc from c Delta++ + K, pT vs phi
  TH2F *fhPtPhiLcb[11];        //! hist. for n Lc tot. from b, pT vs phi
  TH2F *fhPtPhiLcNonRb[11];    //! hist. for n Lc from b non resonant, pT vs phi
  TH2F *fhPtPhiLcL1520b[11];   //! hist. for n Lc from b L1520 + pi, pT vs phi
  TH2F *fhPtPhiLcKstarb[11];   //! hist. for n Lc from b K* + p, pT vs phi
  TH2F *fhPtPhiLcDeltab[11];   //! hist. for n Lc from b Delta++ + K, pT vs phi
  TNtuple *fNtupleLambdac; //! output ntuple
  Float_t fCutsKF[2]; // cuts with KF vertexer
  Int_t fIsLc; // is MC Lc - 0=not Lc, 1=Lc from c, 2=Lc from b
  Int_t fIsLcResonant; // is Lc resonant - 1=non resonant, 2=via L1520 + pi, 3=via K* + p, 4=via Delta++ + K
  Float_t fCandidateVars[4]; // candidate variables, 0=Pt, 1=Eta, 2=Y, 3=Phi
  Float_t fPtLc; // pt of Lc candidate
  AliRDHFCutsLctopKpi *fRDCutsAnalysis; // Analysis cuts
  TList *fListCuts; // list of cuts
  Int_t fFillNtuple;   //  filling ntuple type
  Bool_t fReadMC;    // flag for access to MC
  Bool_t fMCPid;    // flag for access to MC
  Bool_t fRealPid;    // flag for real PID
  Bool_t fResPid;      // flag for PID with resonant channels
  Bool_t fUseKF;      // flag to cut with KF vertexer
  Bool_t fAnalysis;      // apply analysis cuts
  AliAnalysisVertexingHF *fVHF;  //  Vertexer heavy flavour (used to pass the cuts)
  Bool_t fLcCut;  //  flag for Lc filter bit cut
  Bool_t fLcPIDCut;  // flag for Lc filter bit PID
  TH1F *fNentries;      // histo with number of entries
  //AliAODpidUtil* fUtilPid;
  AliPIDResponse *fPIDResponse;     //! PID response object
  AliNormalizationCounter *fCounter;//!AliNormalizationCounter on output slot 7


  ClassDef(AliAnalysisTaskSELambdacTMVA,2); // AliAnalysisTaskSE for the invariant mass analysis of heavy-flavour decay candidates (Lambdac)
};

#endif

