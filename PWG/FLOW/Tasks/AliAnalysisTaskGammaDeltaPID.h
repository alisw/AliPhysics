/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

////////////////////////////////////////////////
// AliAnalysisTaskV0nZDCGains:
// Simple Task to fill V0 and ZDC gain files.
// Works with 15o and 18q/r. Support to be added for LHC10h 
// PA: Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)
/////////////////////////////////////////////////



#ifndef ALIANALYSISTASKGAMMADELTAPID_H
#define ALIANALYSISTASKGAMMADELTAPID_H

#include "AliAnalysisTaskSE.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TProfile2D.h"

class    AliVEvent;      
class    AliVVertex;
class    AliAODv0;
class    AliESDEvent;       
class    AliAODEvent;
class    AliAODTrack;   
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;




//
class AliAnalysisTaskGammaDeltaPID : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskGammaDeltaPID();
  AliAnalysisTaskGammaDeltaPID(const char *name);
  virtual ~AliAnalysisTaskGammaDeltaPID();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
  
  //-----> User Defined Functions:
  Int_t GetCentralityScaled0to10(Double_t fCent);
  void  SetupEventAndTaskConfigInfo();
  void  SetListForTrkCorr(TList *flist)      {this->fListTRKCorr = (TList *) flist->Clone(); }
  void  SetListForNUACorr(TList *flist)      {this->fListNUACorr = (TList *) flist->Clone(); }
  void  SetListForV0MCorr(TList *flist)      {this->fListV0MCorr = (TList *) flist->Clone(); }
  void  SetParticle(Int_t part)              {this->gParticleID  =  part;}
  void  SetCumulantHarmonic(Int_t harm)      {this->gHarmonic    =  harm;}

  
  /******* Event Cut Ranges ******/
  void SetCentralityPercentileMin(Double_t cMin) {this->fCentralityMin  = cMin;}
  void SetCentralityPercentileMax(Double_t cMax) {this->fCentralityMax  = cMax;}
  void SetCentralityEstimator(TString sEstim)    {this->sCentrEstimator = sEstim;}
  void SetVzRangeMin(Double_t vzMin)             {this->fMinVzCut       = vzMin;}
  void SetVzRangeMax(Double_t vzMax)             {this->fMaxVzCut       = vzMax;}
  void SetPileUpCutParam(Float_t m,Float_t c)    {this->fPileUpSlopeParm = m; this->fPileUpConstParm = c;}
  void SetFlagSkipPileUpCuts(Bool_t b)           {this->bSkipPileUpCut  = b;}
  void SetFlagSkipAnalysis(Bool_t b)             {this->bSkipAnalysis  = b;}  


  
  /******* Track Cut Ranges ******/
  void SetNSigmaCutTPC(Double_t     nSigTPC)     {this->fNSigmaTPCCut  =  nSigTPC;}
  void SetNSigmaCutTOF(Double_t     nSigTOF)     {this->fNSigmaTOFCut  =  nSigTOF;}
  void SetTrackCutdEdxMin(Float_t  ndEdxMin)     {this->fdEdxMin    = ndEdxMin;}
  void SetTrackCutChi2Min(Double_t  chi2Min)     {this->fTrkChi2Min =  chi2Min;}
  void SetFlagUseKinkTracks(Bool_t    bKink)     {this->bUseKinkTracks = bKink;}
  
  void SetTrackCutNclusterMin(Int_t ncl)         {this->fTPCclustMin = ncl;}
  void SetFilterBit(Int_t fb)                    {this->fFilterBit   =  fb;}
  void SetEtaRangeMin(Double_t emn)              {this->fMinEtaCut   = emn;}
  void SetEtaRangeMax(Double_t emx)              {this->fMaxEtaCut   = emx;}
  void SetPtRangeMin(Double_t ptL)               {this->fMinPtCut    = ptL;}
  void SetPtRangeMax(Double_t ptH)               {this->fMaxPtCut    = ptH;}
  void SetDCAXYRangeMax(Double_t dcaxy)          {this->fDCAxyMax    = dcaxy;}
  void SetDCAZRangeMax(Double_t dcaz)            {this->fDCAzMax    =  dcaz;}
  void SetChi2Range(Double_t chi2)               {this->fChi2Max    =  chi2;}
  void SetEtaNeg(Double_t etaL)                  {this->fEtaGapNeg   = etaL;}
  void SetEtaPos(Double_t etaH)                  {this->fEtaGapPos   = etaH;}




  //V0 Chun zheng: Cuts on V0 and its daughters
  void SetFlagAnalyseLambda(Bool_t bLambda)      {this->bAnalysLambdaPairs = bLambda;}

  void SetV0PtMin(Double_t v0PtMin)                                    {this->fV0PtMin                                   = v0PtMin;}
  void SetV0CPAMin(Double_t v0CPAMin)                                  {this->fV0CPAMin                                 = v0CPAMin;}
  void SetMassMean(Double_t massMean)                                  {this->fMassMean                                 = massMean;}
  void SetLambdaMassCut(Double_t lambdaMassCut)                        {this->fLambdaMassCut                       = lambdaMassCut;}
  void SetV0RapidityMax(Double_t v0RapidityMax)                        {this->fV0RapidityMax                       = v0RapidityMax;}
  void SetV0DecayLengthMax(Double_t v0DecayLengthMax)                  {this->fV0DecayLengthMax                 = v0DecayLengthMax;}
  void SetV0DecayLengthMin(Double_t v0DecayLengthMin)                  {this->fV0DecayLengthMin                 = v0DecayLengthMin;}
  void SetV0DCAToPrimVtxMax(Double_t v0DCAToPrimVtxMax)                {this->fV0DCAToPrimVtxMax               = v0DCAToPrimVtxMax;}
  void SetV0DcaBetweenDaughtersMax(Double_t v0DcaBetweenDaughtersMax)  {this->fV0DcaBetweenDaughtersMax = v0DcaBetweenDaughtersMax;}
  //V0 Daughter Cut
  void SetDaughtersPtMax(Double_t daughtersPtMax)                      {this->fDaughtersPtMax                     = daughtersPtMax;}
  void SetDaughtersEtaMax(Double_t daughtersEtaMax)                    {this->fDaughtersEtaMax                   = daughtersEtaMax;}
  void SetDaughtersNsigma(Double_t daughtersNsigma)                    {this->fDaughtersNsigma                   = daughtersNsigma;}
  void SetDaughtersTPCNclsMin(Double_t daughtersTPCNclsMin)            {this->fDaughtersTPCNclsMin           = daughtersTPCNclsMin;}
  void SetDaughtersDCAToPrimVtxMin(Double_t daughtersDCAToPrimVtxMin)  {this->fDaughtersDCAToPrimVtxMin = daughtersDCAToPrimVtxMin;}
  //Lambda Mass Cut


  /// To be Added later <-- Rihan
  //ZDC and VZERO Calibration
  //void IsZDCEqualize(bool ZDCEqualize)                { this->bZDCEqualize                          = ZDCEqualize;}
  //void IsZDCRecenter(bool ZDCEqualize)                { this->bZDCEqualize                          = ZDCEqualize;}
  //void IsVZEROEqualize(bool VZEROEqualize)            { this->bVZEROEqualize                      = VZEROEqualize;}
  //void IsVZERORecenter(bool VZERORecenter)            { this->bVZERORecenter                      = VZERORecenter;}
  //------ End of user defined function -------

  
  
  //------ End of user defined function -------







  
 protected:

 private:


  AliVEvent             *fVevent;             //! event
  AliESDEvent           *fESD;                //! esd Event
  AliAODEvent           *fAOD;                //! aod Event
  AliPIDResponse        *fPIDResponse;        //! PID Handler
  AliMultSelection      *fMultSelection;      //! For Centrality 
  AliAnalysisUtils      *fAnalysisUtil;       //! Event Selection Options
  TList                 *fListHist;           //! OutputList

  TList                 *fListTRKCorr;        //  Supplied from Task
  TList                 *fListNUACorr;        //  Supplied from Task
  TList                 *fListV0MCorr;        //  Supplied from Task

  //Set from the AddTask:
  Float_t         fCentralityMin;  //
  Float_t         fCentralityMax;  //
  
  //Track Variables to be used:
  Int_t                 gHarmonic;  //
  Int_t               gParticleID;  //
  Int_t                fFilterBit;  //
  Int_t              fTPCclustMin;  //
  Int_t             gOldRunNumber;  //!
  Bool_t           bUseKinkTracks;  //

  
  Float_t           fNSigmaTPCCut;  //
  Float_t           fNSigmaTOFCut;  //
  Float_t               fMinPtCut;  //
  Float_t               fMaxPtCut;  //
  Float_t               fDCAxyMax;  //                                                                                                        
  Float_t                fDCAzMax;  // 
  Float_t                fChi2Max;  //
  Double_t       fPileUpSlopeParm;  //
  Double_t       fPileUpConstParm;  //
  Bool_t           bSkipPileUpCut;  //
  Double_t             fEtaGapNeg;  //
  Double_t             fEtaGapPos;  //
  Float_t              fMinEtaCut;  //
  Float_t              fMaxEtaCut;  //
  Float_t             fTrkChi2Min;  //
  Float_t                fdEdxMin;  // 
  //Event Variables to be used:
  Float_t               fMinVzCut;  //
  Float_t               fMaxVzCut;  //
  TString         sCentrEstimator;  //
  Bool_t            bSkipAnalysis;  //

  
  /// Chunzheng: V0 (Lambda) Cut parameters:
  Double_t                   fV0PtMin; //
  Double_t                  fV0CPAMin; //
  Double_t             fV0RapidityMax; //
  Double_t          fV0DecayLengthMin; //
  Double_t          fV0DecayLengthMax; //
  Double_t         fV0DCAToPrimVtxMax; //
  Double_t  fV0DcaBetweenDaughtersMax; //
  //V0 Daughter
  Double_t            fDaughtersPtMax; //
  Double_t           fDaughtersNsigma; //
  Double_t           fDaughtersEtaMax; //
  Double_t       fDaughtersTPCNclsMin; //
  Double_t  fDaughtersDCAToPrimVtxMin; //
  //Lambda Mass
  Double_t                  fMassMean; //
  Double_t             fLambdaMassCut; //
  Bool_t           bAnalysLambdaPairs; //
  //Double_t             fCurrentVtx[3]; //!  Global variable because Accessed by Multiple functions
  //-------------------------------------
  



  /////////// Default HISTOGRAMS /////////
  
  //QA histograms:
  TH1F          *fCentDistBeforCut;        //! Cent before Any Cut
  TH1F          *fCentDistAfterCut;        //! Cent After All Cuts
  TH2F          *fHistTPConlyVsCL1Before;  //!  
  TH2F          *fHistTPConlyVsV0MBefore;  //!
  TH2F          *fHistCL0VsV0MBefore;      //!       
  TH2F          *fHistTPConlyVsCL1After;   //!    
  TH2F          *fHistTPConlyVsV0MAfter;   //!   
  TH2F          *fHistCL0VsV0MAfter;       //!   
  TH2F          *fHistTPCVsESDTrkBefore;   //!  
  TH2F          *fHistTPCVsESDTrkAfter;    //!
  TH1F          *fHistPileUpCount;         //!
  TH1F          *fHistAnalysisInfo;        //!

  
  /// Functions for Pile Up Event Removal:  
  TF1           *fSPDCutPU;     //!
  TF1           *fV0CutPU;      //!
  TF1           *fMultCutPU;    //!
  TF1           *fCenCutLowPU;  //!
  TF1           *fCenCutHighPU; //!

    
  /// Store the EP <Q> vectors:
  TProfile      *fAvgCosNPsivsCentEtaPos; //! N is harmonic set By AddTask. Default is 2nd order 
  TProfile      *fAvgSinNPsivsCentEtaPos; //!
  TProfile      *fAvgCosNPsivsCentEtaNeg; //!
  TProfile      *fAvgSinNPsivsCentEtaNeg; //!

  TProfile      *fAvgCos3PsivsCentEtaPos; //!  3rd Harmonic is hardcoded.
  TProfile      *fAvgSin3PsivsCentEtaPos; //!
  TProfile      *fAvgCos3PsivsCentEtaNeg; //!
  TProfile      *fAvgSin3PsivsCentEtaNeg; //!  
  
  TH1F          *fHistVertexZcm;     //!
  TProfile      *fHistVxvsVzMinBias; //!
  TProfile      *fHistVyvsVzMinBias; //!
  
  TProfile2D    *hAvgZNACh0vsCentVz; //!
  TProfile2D    *hAvgZNCCh0vsCentVz; //!
  TProfile2D    *hAvgZNACh1vsCentVz; //!
  TProfile2D    *hAvgZNCCh1vsCentVz; //!
  TProfile2D    *hAvgZNACh2vsCentVz; //!
  TProfile2D    *hAvgZNCCh2vsCentVz; //!
  TProfile2D    *hAvgZNACh3vsCentVz; //!
  TProfile2D    *hAvgZNCCh3vsCentVz; //!
  TProfile2D    *hAvgZNACh4vsCentVz; //!
  TProfile2D    *hAvgZNCCh4vsCentVz; //!
  
  TProfile2D    *hAvgV0ChannelsvsVz; //!

  TProfile      *hAvgQNXvsCentV0C; //!
  TProfile      *hAvgQNYvsCentV0C; //!
  TProfile      *hAvgQ3XvsCentV0C; //!
  TProfile      *hAvgQ3YvsCentV0C; //!
  TProfile      *hAvgQNXvsCentV0A; //!
  TProfile      *hAvgQNYvsCentV0A; //!
  TProfile      *hAvgQ3XvsCentV0A; //!
  TProfile      *hAvgQ3YvsCentV0A; //!  


  ///Analysis Histograms:
  TH2F        *fHistTPCPsiNPosPlane; //!
  TH2F        *fHistTPCPsiNNegPlane; //!
  TH2F        *fHistTPCPsi3PosPlane; //!
  TH2F        *fHistTPCPsi3NegPlane; //!

  TH2F        *fHistV0CPsiNEventPlane; //!  
  TH2F        *fHistV0APsiNEventPlane; //!
  TH2F        *fHistV0CPsi3EventPlane; //!  
  TH2F        *fHistV0APsi3EventPlane; //!  


 
  TProfile     *hTPCPsiNCorrelation;  //!
  TProfile     *hTPCPsi3Correlation;  //!

  TProfile     *hV0CV0APsiNCorrelation;  //!
  TProfile     *hV0CTPCPsiNCorrelation;  //!
  TProfile     *hV0ATPCPsiNCorrelation;  //!
  TProfile     *hV0CV0APsi3Correlation;  //!
  TProfile     *hV0CTPCPsi3Correlation;  //!
  TProfile     *hV0ATPCPsi3Correlation;  //!  

  
  TProfile      *hAvg3pC112vsCentPP;  //!
  TProfile      *hAvg3pC112vsCentNN;  //!
  TProfile      *hAvg3pC112vsCentOS;  //!
  TProfile      *hAvg3pC123vsCentPP;  //!
  TProfile      *hAvg3pC123vsCentNN;  //!
  TProfile      *hAvg3pC123vsCentOS;  //!  
  //TProfile      *hAvg3pC224vsCentPP;  //!
  //TProfile      *hAvg3pC224vsCentNN;  //!
  //TProfile      *hAvg3pC224vsCentOS;  //!
  
  
  TProfile      *hAvgDelta1vsCentPP;  //!
  TProfile      *hAvgDelta1vsCentNN;  //!
  TProfile      *hAvgDelta1vsCentOS;  //!
  TProfile      *hAvgDelta2vsCentPP;  //!
  TProfile      *hAvgDelta2vsCentNN;  //!
  TProfile      *hAvgDelta2vsCentOS;  //!
  TProfile      *hAvgDelta3vsCentPP;  //!
  TProfile      *hAvgDelta3vsCentNN;  //!
  TProfile      *hAvgDelta3vsCentOS;  //!
  TProfile      *hAvgDelta4vsCentPP;  //!
  TProfile      *hAvgDelta4vsCentNN;  //!
  TProfile      *hAvgDelta4vsCentOS;  //!  


  //// NUA, NUE, V0 Weights for Corrections:
  TH3F         *fHCorrectNUAChrgPos; //!
  TH3F         *fHCorrectNUAChrgNeg; //!
  TH3F         *fHCorrectNUAkPIDPos; //!
  TH3F         *fHCorrectNUAkPIDNeg; //!
  
  TH1D         *fHCorrectMCPosChrg; //!
  TH1D         *fHCorrectMCPosPion; //!
  TH1D         *fHCorrectMCPosKaon; //!
  TH1D         *fHCorrectMCPosProt; //!
   
  TH1D         *fHCorrectMCNegChrg; //!
  TH1D         *fHCorrectMCNegPion; //!
  TH1D         *fHCorrectMCNegKaon; //!
  TH1D         *fHCorrectMCNegProt; //!

  TH1D         *fHCorrectTPCQxEtaPos; //!
  TH1D         *fHCorrectTPCQyEtaPos; //!
  TH1D         *fHCorrectTPCQxEtaNeg; //!
  TH1D         *fHCorrectTPCQyEtaNeg; //!



  /// chunzheng: V0 QA Histograms
  TH1D             *fHistV0Pt;              //! Raw V0s' pT
  TH1D             *fHistV0Eta;             //! Raw V0s' eta
  TH1D             *fHistV0DcatoPrimVertex; //! Raw V0s' DcatoPV
  TH1D             *fHistV0CPA;             //! Raw V0s' CPA(cosine pointing angle)
  TH1D             *fHistV0DecayLength;     //! Raw V0s' DecayLength

  ///Lambda-X correlators:
  /// Naming Style for Histograms:
  /// Lambda = L-capital, Proton = P-capital, Anti= A-capital. e.g, _AntiProton_ or AntiLambda
  
  TProfile *fProfileGammaTPC_Lambda_hPos; //!
  TProfile *fProfileGammaTPC_Lambda_hNeg; //!
  TProfile *fProfileGammaTPC_Lambda_Proton;     //!
  TProfile *fProfileGammaTPC_Lambda_AntiProton; //!

  TProfile *fProfileGammaTPC_AntiLambda_hPos; //!
  TProfile *fProfileGammaTPC_AntiLambda_hNeg; //!
  TProfile *fProfileGammaTPC_AntiLambda_Proton;     //!
  TProfile *fProfileGammaTPC_AntiLambda_AntiProton; //!





  
  // The arrays[] are always at the end of list:
 
  /// chunzheng: Lambda QA Info
  Double_t          fCurrentVtx[3];         //! Global variable: Accessed by Multiple functions
  TH1D             *fHistLambdaPt[2];       //! [0]:Before the Mass Cut [1]:After the Mass Cut
  TH1D             *fHistLambdaEta[2];      //!
  TH1D             *fHistLambdaDcaToPrimVertex[2]; //!
  TH1D             *fHistLambdaCPA[2]; //!
  TH1D             *fHistLambdaDecayLength[2]; //!
  TH1D             *fHistLambdaMass[2]; //!
  TProfile         *fProfileLambdaMassVsPt[2]; //!
  TH1D             *fHistAntiLambdaPt[2]; //!
  TH1D             *fHistAntiLambdaEta[2]; //!
  TH1D             *fHistAntiLambdaDcaToPrimVertex[2]; //!
  TH1D             *fHistAntiLambdaCPA[2]; //!
  TH1D             *fHistAntiLambdaDecayLength[2]; //!
  TH1D             *fHistAntiLambdaMass[2]; //!
  TProfile         *fProfileAntiLambdaMassVsPt[2]; //!



  ///Custom Functions:
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  //void  GetV0MCorrectionHist(Int_t run=0,Int_t kHarmonic=0);
  void  GetMCCorrectionHist();
  
  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t CheckEventIsPileUp2018(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  Bool_t CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck);

  ///// Chunzheng: V0 (Lambda)
  Bool_t IsGoodV0(AliAODv0 *aodV0);
  Bool_t IsGoodDaughterTrack(const AliAODTrack *track);

  AliAnalysisTaskGammaDeltaPID(const AliAnalysisTaskGammaDeltaPID &other);
  AliAnalysisTaskGammaDeltaPID &operator=(const AliAnalysisTaskGammaDeltaPID &other);    
  ClassDef(AliAnalysisTaskGammaDeltaPID, 1) 



};

#endif
