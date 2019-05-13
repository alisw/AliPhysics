#ifndef AliAnalysisSphericityTask_H
#define AliAnalysisSphericityTask_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */


// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TProfile.h>
#include <TTreeStream.h>
#include <TRandom.h>
#include <TObject.h>
#include <vector>  
#include <TH1F.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TF1.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliMCEvent.h>
#include <AliAnalysisFilter.h>
#include <AliStack.h>
#include <AliGenEventHeader.h>
#include <AliVHeader.h>
#include <AliAODMCParticle.h> 
#include <AliESDtrackCuts.h>
#include <AliTransverseEventShape.h>
#include <AliTPCPIDResponse.h>
#include <AliPPVsMultUtils.h>
#include <AliESDVertex.h>
#include <AliPhysicsSelection.h>
#include <AliPIDResponse.h>



class AliAnalysisSphericityTask : public AliAnalysisTaskSE {
 public:
  enum PIDMode { kRatio = 0, kSigma};
  
  AliAnalysisSphericityTask();
  AliAnalysisSphericityTask(const char *name);
  virtual ~AliAnalysisSphericityTask();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t*);
  virtual Float_t GetTest();

  Bool_t   GetAnalysisMC() { return fAnalysisMC; }   
  Double_t GetVtxCut() { return fVtxCut; }   
  Double_t GetEtaCut() { return fEtaCut; }   

  //for Strangeness
  
  virtual void   SetMinDaughterTpcClusters(const Int_t minDaughterTpcClusters = 80) {fMinDaughterTpcClusters = minDaughterTpcClusters;}
  virtual void   SetQualityCutTPCrefit(Bool_t qualityCutTPCrefit = kTRUE) { fkQualityCutTPCrefit =  qualityCutTPCrefit;} 
  virtual void   SetOADBPath(const char* path) {fOADBPath=path;}
  const char*    GetOADBPath() const { return fOADBPath.Data(); }
  virtual void   SetPIDMode(PIDMode pidmode, Double_t cut1, Double_t cut2, Double_t bound) 
  {fPIDMode = pidmode;
  if(fPIDMode == kSigma) { fNSigma1 = cut1; fNSigma2 = cut2; fNBoundP = bound;}
  if(fPIDMode == kRatio) { fNRatio1 = cut1; fNRatio2 = cut2; fNBoundP = bound;}
  }

  virtual void   SetExtraSelections(Bool_t extraSelections = 0) { fkExtraSelections =  extraSelections;} 
  virtual void   SetExtraSelectionsCut(Bool_t extraSelectionscut = 0) { fkExtraSelectionsCut =  extraSelectionscut;} 
  virtual void   SetRerunVertexers(Bool_t Rerun = kFALSE){fRerunVertexers = Rerun;}
  virtual void   SetV0Cuts(Double_t* v0Selections){fV0Cuts = v0Selections;}
  virtual void   SetCascadeCuts(Double_t* CascSelections){fCascadeCuts = CascSelections;}
  virtual void   SetMaxV0Rapidity(const Float_t maxV0Rapidity = 1.) {fMaxV0Rapidity = maxV0Rapidity;}
  virtual void   SetInvMassCutKaon(Float_t InvMassCutKaon = 0) { fInvMassCutKaon =  InvMassCutKaon;}
  virtual void   SetInvMassCutLambda(Float_t InvMassCutLambda = 0) { fInvMassCutLambda =  InvMassCutLambda;}
  virtual void   SetInvMassCutXi(Float_t InvMassCutXi = 0) { fInvMassCutXi =  InvMassCutXi;}
  virtual void   SetInvMassCutOmega(Float_t InvMassCutOmega = 0) { fInvMassCutOmega =  InvMassCutOmega;} 
  
  //for event shape analysis
  virtual void  SetUseHybridESA(Bool_t usehyb) {fUseHybrid = usehyb;}
  virtual void  SetTrackFilterESAHyb1(AliAnalysisFilter* trackH1F) {fTrackFilterHybrid1 = trackH1F;}
  virtual void  SetTrackFilterESAHyb2(AliAnalysisFilter* trackH2F) {fTrackFilterHybrid2 = trackH2F;}
  virtual void  SetTrackFilterESA(AliAnalysisFilter* trackF) {fTrackFilterESA = trackF;}
  virtual void  SetMinMultForESA(Int_t minnch)     {fMinMultESA = minnch;}
  virtual void  SetStepSizeESA(Float_t sizestep)   {fSizeStepESA = sizestep;}
  virtual void  SetIsEtaAbsESA(Bool_t isabseta){fIsAbsEtaESA = isabseta;}
  virtual void  SetTrackEtaMinESA(Float_t etaminF) {fEtaMinCutESA = etaminF;}
  virtual void  SetTrackEtaMaxESA(Float_t etamaxF) {fEtaMaxCutESA = etamaxF;}
  virtual void  SetTrackPtMinESA(Float_t ptminF) {fPtMinCutESA = ptminF;}
  virtual void  SetTrackPtMaxESA(Float_t ptmaxF) {fPtMaxCutESA = ptmaxF;}
  virtual void  SetTrigger(UInt_t ktriggerInt) {ftrigBit = ktriggerInt;}
  virtual void  SetCentralityEstimator(const char * centEst) {fCentEst = centEst;}
  virtual void  SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  virtual void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}
  virtual void  SetVtxCut(Double_t vtxCut){fVtxCut = vtxCut;}
  virtual void  SetEtaCut(Double_t etaCut){fEtaCut = etaCut;}
  virtual void  SetPileUpRej(Bool_t isrej) {fPileUpRej = isrej;}   
  virtual void  SetMinCent(Float_t minvalc) {fMinCent = minvalc;}
  virtual void  SetMaxCent(Float_t maxvalc) {fMaxCent = maxvalc;}
  virtual void  SetStoreMcIn(Bool_t value) {fStoreMcIn = value;}
  virtual void  SetAnalysisPbPb(Bool_t isanaPbPb) {fAnalysisPbPb = isanaPbPb;}
  virtual void  SetStrangeness(Bool_t Strange) {fStrangeness = Strange;}   

  private:
  virtual Float_t GetVertex(const AliVEvent* event) const;

  Int_t   GetV0MIndex(Float_t V0MPercentile);
 	Int_t   GetMultiplicityIndex(Int_t Mult);
  Bool_t IsAccepted (AliESDtrack *track);
  Double_t Rapidity(Double_t,Double_t,Double_t,Int_t)const;
  Bool_t IsProton (AliESDtrack *track);
  Bool_t IsPion (AliESDtrack *track);
  Bool_t IsKaon (AliESDtrack *track);

  TString fOADBPath;                                  // OADB path to use
  AliPIDResponse *fPIDResponse;                       //! PID response Handler
  Int_t   fOldRun;                                    //! current run number
  Int_t   fRecoPass;                                  //! reconstruction pass

  AliTPCPIDResponse fTPCpid;                          // Tool data member to manage the TPC Bethe-Bloch info
  PIDMode fPIDMode; 				                          //PID mode:dE/dx Ratio-Nsigma areas
  Double_t fNBoundP;				                          //Momentum bound between two pid cuts
  Double_t fNSigma1; 				                          //N-sigma cut in the dE/dx band
  Double_t fNSigma2; 				                          //N-sigma cut in the dE/dx band
  Double_t fNRatio1; 				                          //min value of the ratio of the measured dE/dx vs the expected
  Double_t fNRatio2; 				                          //min value of the ratio of the measured dE/dx vs the expected
  
  AliESDEvent*      fESD;                             //! ESD object
  AliAODEvent*      fAOD;                             //! AOD object
  AliPPVsMultUtils *fPPVsMultUtils;                   //! object for Vzero based multiplicity
  AliMCEvent*       fMC;                              //! MC object
  AliStack*         fMCStack;                         //! MC ESD stack
  TClonesArray*     fMCArray;                         //! MC array for AOD

  Bool_t            fUseHybrid;
  AliAnalysisFilter *fTrackFilterHybrid1;
  AliAnalysisFilter *fTrackFilterHybrid2;
  AliAnalysisFilter *fTrackFilterESA;                 // Track filter for Event Shapes
  Int_t             fMinMultESA;
  Float_t           fSizeStepESA;
  Bool_t            fIsAbsEtaESA;
  Float_t           fEtaMaxCutESA;
  Float_t           fEtaMinCutESA;
  Float_t           fPtMaxCutESA;
  Float_t           fPtMinCutESA;
  Int_t             fNrec;

  TString           fCentEst;                         // V0A , V0M, 
  TString           fAnalysisType;                    //  "ESD" or "AOD"
  Bool_t            fAnalysisMC;                      //  Real(kFALSE) or MC(kTRUE) flag
  Bool_t            fAnalysisPbPb;                    //  true you want to analyze PbPb data, false for pp
  UInt_t            ftrigBit;
  TRandom*          fRandom;                          //! random number generator
  Bool_t            fPileUpRej;                       // kTRUE is pile-up is rejected
  Bool_t            fStrangeness;                                                 
  Int_t             fCent;                            //minimum centrality

  Bool_t            fkExtraSelections;                // Boolean : kTRUE = apply tighter selections, before starting the analysis
  Bool_t            fkExtraSelectionsCut;                // Boolean : kTRUE = apply cuts fot PID
  Bool_t            fkQualityCutTPCrefit;             // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
  Int_t             fMinDaughterTpcClusters;          //  Minimum number of TPC clusters for the both daughter tracks of the V0
  Bool_t            fRerunVertexers;
  Double_t*         fV0Cuts;
  Double_t*         fCascadeCuts;
  Float_t           fMaxV0Rapidity;                   //  Maximum rapidity selection for the V0

  // Cuts and options

  Double_t     fVtxCut;                               // Vtx cut on z position in cm
  Double_t     fEtaCut;                               // Eta cut used to select particles
  Float_t      fMinCent;                              //minimum centrality
  Float_t      fMaxCent;                              //maximum centrality
  Bool_t       fStoreMcIn;                            // Store MC input tracks

  Float_t      fInvMassCutKaon;
  Float_t      fInvMassCutLambda;
  Float_t      fInvMassCutOmega;
  Float_t      fInvMassCutXi;

  // Help variables
  
  Short_t      fMcProcessType;                        // -1=invalid, 0=data, 1=ND, 2=SD, 3=DD
  Short_t      fTriggeredEventMB;                     // 1 = triggered, 0 = not trigged (MC only)
  Short_t      fVtxStatus;                            // -1 = no vtx, 0 = outside cut, 1 = inside cut
  Float_t      fZvtx;                                 // z vertex
  Float_t      fZvtxMC;                               // z vertex MC (truth)
  Int_t        fRun;                                  // run no
  ULong64_t    fEventId;                              // unique event id

  Bool_t       fStrangePart;
  Bool_t       fFindKaon;
  Bool_t       fFindLambda;
  Bool_t       fFindAntiLambda;
  Bool_t       fFindXiPlus;
  Bool_t       fFindXiMinus;
  Bool_t       fFindOmegaPlus;
  Bool_t       fFindOmegaMinus;

  // Output objects
  
  TList*        fListOfObjects;                       //! Output list of objects
  TH1D*         hVtxBeforeCuts;                       //! Vertex z dist before cuts
  TH1D*         hVtxAfterCuts;                        //! Vertex z dist after cuts
  TH1D*         hn1;
  TH1D*         hn2;               
  TH1D*         hso;
  TH1D*         hcent;
  TH1D*         hst;
  TH1D*         fhphiSt;        
  TH1D*         fhetaSt;         
  TH1D*         fhptSt;        
  TH1D*         HMultRef;
  TH2D*         hetaVSphi; 
  TH2D*         hetaVSphiISO; 
  TH2D*         hetaVSphiJET; 
  TH2D*         hetaVSphiMID; 
  TH2D*         hsphericityVSmulti;
  TH2D*         hsphericityVSpT; 
  TH2D*         hmultiplicityVSpT; 
  TH2D*         hmultiplicityVSpTbin; 
  TH2D*         hsphericityVScent;
  TH2D*         hsphericityVSMEANpT; 
  TH2D*         hmultiplicityVSMEANpT;
  TH2D*         hmultiplicityVSMEANpTbin;  
  TH2D*         hsphericity;     
  TH2D*         hsphericityVSmultiL;
  TH2D*         hsphericityL;
  TH2D*         hsphericityVSpTL;
  TH2D*         hsphericityVScentL;
  TH2D*         hsphericityVSMEANpTL;
  TH2D*         hmultiplicityVSMEANpTL;
  TH2D*         hmultiplicityVSMEANpTbinL;
  TH2D*         hmultiplicityVSpTL; 
  TH2D*         hmultiplicityVSpTbinL;
  TH2D*         hsphericityVSmultiAL;
  TH2D*         hsphericityAL;
  TH2D*         hsphericityVSpTAL;
  TH2D*         hsphericityVScentAL;
  TH2D*         hsphericityVSMEANpTAL;
  TH2D*         hmultiplicityVSMEANpTAL;
  TH2D*         hmultiplicityVSMEANpTbinAL;
  TH2D*         hmultiplicityVSpTAL; 
  TH2D*         hmultiplicityVSpTbinAL;
  TH2D*         hsphericityVSmultiLO;
  TH2D*         hsphericityLO;
  TH2D*         hsphericityVSpTLO;
  TH2D*         hsphericityVScentLO;
  TH2D*         hsphericityVSMEANpTLO;
  TH2D*         hmultiplicityVSMEANpTLO;
  TH2D*         hmultiplicityVSMEANpTbinLO;
  TH2D*         hmultiplicityVSpTLO; 
  TH2D*         hmultiplicityVSpTbinLO;
  TH2D*         hsphericityVSmultiALLO;
  TH2D*         hsphericityALLO;
  TH2D*         hsphericityVSpTALLO;
  TH2D*         hsphericityVScentALLO;
  TH2D*         hsphericityVSMEANpTALLO;
  TH2D*         hmultiplicityVSMEANpTALLO;
  TH2D*         hmultiplicityVSMEANpTbinALLO;
  TH2D*         hmultiplicityVSpTALLO; 
  TH2D*         hmultiplicityVSpTbinALLO;
  TH2D*         hsphericityVSmultiX;
  TH2D*         hsphericityX;
  TH2D*         hsphericityVSpTX;
  TH2D*         hsphericityVScentX;
  TH2D*         hsphericityVSMEANpTX;
  TH2D*         hmultiplicityVSMEANpTX;
  TH2D*         hmultiplicityVSMEANpTbinX;
  TH2D*         hmultiplicityVSpTX; 
  TH2D*         hmultiplicityVSpTbinX;
  TH2D*         hsphericityVSmultiXP;
  TH2D*         hsphericityXP;
  TH2D*         hsphericityVSpTXP;
  TH2D*         hsphericityVScentXP;
  TH2D*         hsphericityVSMEANpTXP;
  TH2D*         hmultiplicityVSMEANpTXP;
  TH2D*         hmultiplicityVSMEANpTbinXP;
  TH2D*         hmultiplicityVSpTXP; 
  TH2D*         hmultiplicityVSpTbinXP;
  TH2D*         hsphericityVSmultiXN;
  TH2D*         hsphericityXN;
  TH2D*         hsphericityVSpTXN;
  TH2D*         hsphericityVScentXN;
  TH2D*         hsphericityVSMEANpTXN;
  TH2D*         hmultiplicityVSMEANpTXN;
  TH2D*         hmultiplicityVSMEANpTbinXN;
  TH2D*         hmultiplicityVSpTXN; 
  TH2D*         hmultiplicityVSpTbinXN;
  TH2D*         hsphericityVSmultiXNXP;
  TH2D*         hsphericityXNXP;
  TH2D*         hsphericityVSpTXNXP;
  TH2D*         hsphericityVScentXNXP;
  TH2D*         hsphericityVSMEANpTXNXP;
  TH2D*         hmultiplicityVSMEANpTXNXP;
  TH2D*         hmultiplicityVSMEANpTbinXNXP;
  TH2D*         hmultiplicityVSpTXNXP; 
  TH2D*         hmultiplicityVSpTbinXNXP;
  TH2D*         hsphericityVSmultiO;
  TH2D*         hsphericityO;
  TH2D*         hsphericityVSpTO;
  TH2D*         hsphericityVScentO;
  TH2D*         hsphericityVSMEANpTO;
  TH2D*         hmultiplicityVSMEANpTO;
  TH2D*         hmultiplicityVSMEANpTbinO;
  TH2D*         hmultiplicityVSpTO; 
  TH2D*         hmultiplicityVSpTbinO;
  TH2D*         hsphericityVSmultiOP;
  TH2D*         hsphericityOP;
  TH2D*         hsphericityVSpTOP;
  TH2D*         hsphericityVScentOP;
  TH2D*         hsphericityVSMEANpTOP;
  TH2D*         hmultiplicityVSMEANpTOP;
  TH2D*         hmultiplicityVSMEANpTbinOP;
  TH2D*         hmultiplicityVSpTOP; 
  TH2D*         hmultiplicityVSpTbinOP;
  TH2D*         hsphericityVSmultiON;
  TH2D*         hsphericityON;
  TH2D*         hsphericityVSpTON;
  TH2D*         hsphericityVScentON;
  TH2D*         hsphericityVSMEANpTON;
  TH2D*         hmultiplicityVSMEANpTON;
  TH2D*         hmultiplicityVSMEANpTbinON;
  TH2D*         hmultiplicityVSpTON; 
  TH2D*         hmultiplicityVSpTbinON;
  TH2D*         hsphericityVSmultiONOP;
  TH2D*         hsphericityONOP;
  TH2D*         hsphericityVSpTONOP;
  TH2D*         hsphericityVScentONOP;
  TH2D*         hsphericityVSMEANpTONOP;
  TH2D*         hmultiplicityVSMEANpTONOP;
  TH2D*         hmultiplicityVSMEANpTbinONOP;
  TH2D*         hmultiplicityVSpTONOP; 
  TH2D*         hmultiplicityVSpTbinONOP;
  TH2D*         hsphericityVSmultiKA;
  TH2D*         hsphericityKA;
  TH2D*         hsphericityVSpTKA;
  TH2D*         hsphericityVScentKA;
  TH2D*         hsphericityVSMEANpTKA;
  TH2D*         hmultiplicityVSMEANpTKA;
  TH2D*         hmultiplicityVSMEANpTbinKA;
  TH2D*         hmultiplicityVSpTKA; 
  TH2D*         hmultiplicityVSpTbinKA;

  TH2D*         hsphericityVSmultiL2;
  TH2D*         hsphericityL2;
  TH2D*         hsphericityVSpTL2;
  TH2D*         hsphericityVScentL2;
  TH2D*         hsphericityVSMEANpTL2;
  TH2D*         hmultiplicityVSMEANpTL2;
  TH2D*         hmultiplicityVSMEANpTbinL2;
  TH2D*         hmultiplicityVSpTL2; 
  TH2D*         hmultiplicityVSpTbinL2;
  TH2D*         hsphericityVSmultiALO;
  TH2D*         hsphericityALO;
  TH2D*         hsphericityVSpTALO;
  TH2D*         hsphericityVScentALO;
  TH2D*         hsphericityVSMEANpTALO;
  TH2D*         hmultiplicityVSMEANpTALO;
  TH2D*         hmultiplicityVSMEANpTbinALO;
  TH2D*         hmultiplicityVSpTALO; 
  TH2D*         hmultiplicityVSpTbinALO;
  TH2D*         hsphericityVSmultiLOO;
  TH2D*         hsphericityLOO;
  TH2D*         hsphericityVSpTLOO;
  TH2D*         hsphericityVScentLOO;
  TH2D*         hsphericityVSMEANpTLOO;
  TH2D*         hmultiplicityVSMEANpTLOO;
  TH2D*         hmultiplicityVSMEANpTbinLOO;
  TH2D*         hmultiplicityVSpTLOO; 
  TH2D*         hmultiplicityVSpTbinLOO;
  TH2D*         hsphericityVSmultiALLOO;
  TH2D*         hsphericityALLOO;
  TH2D*         hsphericityVSpTALLOO;
  TH2D*         hsphericityVScentALLOO;
  TH2D*         hsphericityVSMEANpTALLOO;
  TH2D*         hmultiplicityVSMEANpTALLOO;
  TH2D*         hmultiplicityVSMEANpTbinALLOO;
  TH2D*         hmultiplicityVSpTALLOO; 
  TH2D*         hmultiplicityVSpTbinALLOO;
  TH2D*         hsphericityVSmultiXO;
  TH2D*         hsphericityXO;
  TH2D*         hsphericityVSpTXO;
  TH2D*         hsphericityVScentXO;
  TH2D*         hsphericityVSMEANpTXO;
  TH2D*         hmultiplicityVSMEANpTXO;
  TH2D*         hmultiplicityVSMEANpTbinXO;
  TH2D*         hmultiplicityVSpTXO; 
  TH2D*         hmultiplicityVSpTbinXO;
  TH2D*         hsphericityVSmultiXPO;
  TH2D*         hsphericityXPO;
  TH2D*         hsphericityVSpTXPO;
  TH2D*         hsphericityVScentXPO;
  TH2D*         hsphericityVSMEANpTXPO;
  TH2D*         hmultiplicityVSMEANpTXPO;
  TH2D*         hmultiplicityVSMEANpTbinXPO;
  TH2D*         hmultiplicityVSpTXPO; 
  TH2D*         hmultiplicityVSpTbinXPO;
  TH2D*         hsphericityVSmultiXNO;
  TH2D*         hsphericityXNO;
  TH2D*         hsphericityVSpTXNO;
  TH2D*         hsphericityVScentXNO;
  TH2D*         hsphericityVSMEANpTXNO;
  TH2D*         hmultiplicityVSMEANpTXNO;
  TH2D*         hmultiplicityVSMEANpTbinXNO;
  TH2D*         hmultiplicityVSpTXNO; 
  TH2D*         hmultiplicityVSpTbinXNO;
  TH2D*         hsphericityVSmultiXNXPO;
  TH2D*         hsphericityXNXPO;
  TH2D*         hsphericityVSpTXNXPO;
  TH2D*         hsphericityVScentXNXPO;
  TH2D*         hsphericityVSMEANpTXNXPO;
  TH2D*         hmultiplicityVSMEANpTXNXPO;
  TH2D*         hmultiplicityVSMEANpTbinXNXPO;
  TH2D*         hmultiplicityVSpTXNXPO; 
  TH2D*         hmultiplicityVSpTbinXNXPO;
  TH2D*         hsphericityVSmultiOO;
  TH2D*         hsphericityOO;
  TH2D*         hsphericityVSpTOO;
  TH2D*         hsphericityVScentOO;
  TH2D*         hsphericityVSMEANpTOO;
  TH2D*         hmultiplicityVSMEANpTOO;
  TH2D*         hmultiplicityVSMEANpTbinOO;
  TH2D*         hmultiplicityVSpTOO; 
  TH2D*         hmultiplicityVSpTbinOO;
  TH2D*         hsphericityVSmultiOPO;
  TH2D*         hsphericityOPO;
  TH2D*         hsphericityVSpTOPO;
  TH2D*         hsphericityVScentOPO;
  TH2D*         hsphericityVSMEANpTOPO;
  TH2D*         hmultiplicityVSMEANpTOPO;
  TH2D*         hmultiplicityVSMEANpTbinOPO;
  TH2D*         hmultiplicityVSpTOPO; 
  TH2D*         hmultiplicityVSpTbinOPO;
  TH2D*         hsphericityVSmultiONO;
  TH2D*         hsphericityONO;
  TH2D*         hsphericityVSpTONO;
  TH2D*         hsphericityVScentONO;
  TH2D*         hsphericityVSMEANpTONO;
  TH2D*         hmultiplicityVSMEANpTONO;
  TH2D*         hmultiplicityVSMEANpTbinONO;
  TH2D*         hmultiplicityVSpTONO; 
  TH2D*         hmultiplicityVSpTbinONO;
  TH2D*         hsphericityVSmultiONOPO;
  TH2D*         hsphericityONOPO;
  TH2D*         hsphericityVSpTONOPO;
  TH2D*         hsphericityVScentONOPO;
  TH2D*         hsphericityVSMEANpTONOPO;
  TH2D*         hmultiplicityVSMEANpTONOPO;
  TH2D*         hmultiplicityVSMEANpTbinONOPO;
  TH2D*         hmultiplicityVSpTONOPO; 
  TH2D*         hmultiplicityVSpTbinONOPO;
  TH2D*         hsphericityVSmultiKAO;
  TH2D*         hsphericityKAO;
  TH2D*         hsphericityVSpTKAO;
  TH2D*         hsphericityVScentKAO;
  TH2D*         hsphericityVSMEANpTKAO;
  TH2D*         hmultiplicityVSMEANpTKAO;
  TH2D*         hmultiplicityVSMEANpTbinKAO;
  TH2D*         hmultiplicityVSpTKAO; 
  TH2D*         hmultiplicityVSpTbinKAO;
   
  TH1F        *fHistTrackMultiplicity;                //! Track multiplicity distribution
  TH1F        *fHistMassKaon;                         //! Invariant mass of Kaon0
  TH1F        *fHistMassXiMinus;                      //! Invariant mass of XiMinus
  TH1F        *fHistMassXiPlus;                       //! Invariant mass of XiPlus
  TH1F        *fHistMassOmegaMinus;                   //! Invariant mass of OmegaMinus
  TH1F        *fHistMassOmegaPlus;                    //! Invariant mass of OmegaPlus
  TH1F        *fHistMassLambda;                       //! Invariant mass of Lambda
  TH1F        *fHistMassAntiLambda;                   //! Invariant mass of Anti-Lambda
  TH2F	      *fHistXiArmenteros;			                //! alpha(casc. cand.) Vs PtArm(casc. cand.)
  TH2F	      *fHistV0Armenteros;			                //! alpha(casc. cand.) Vs PtArm(casc. cand.)
 
  // protected:
  AliTransverseEventShape* fESASelection;             // event selection class


  AliAnalysisSphericityTask(const AliAnalysisSphericityTask&);            // not implemented
  AliAnalysisSphericityTask& operator=(const AliAnalysisSphericityTask&); // not implemented

  ClassDef(AliAnalysisSphericityTask, 1);   
};

#endif
