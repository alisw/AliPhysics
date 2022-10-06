/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */
///////////////////////////////////////////////////////////
// Task To Analyse Gamma Delta Correlator with ESE. 
// Works with 15o and 18q/r. Support to be added for LHC10h 
// PA: Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)
///////////////////////////////////////////////////////////




#ifndef ALIANALYSISTASKGAMMADELTAPID_H
#define ALIANALYSISTASKGAMMADELTAPID_H


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
#include "AliAnalysisTaskSE.h"
#include <iostream>
#include <vector>


class    AliVEvent;      
class    AliVVertex;
class    AliAODv0;
class    AliESDEvent;       
class    AliAODEvent;
class    AliAODTrack;   
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;



class AliAnalysisTaskGammaDeltaPID : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskGammaDeltaPID();
  AliAnalysisTaskGammaDeltaPID(const char *name);
  virtual ~AliAnalysisTaskGammaDeltaPID();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
  virtual void Terminate(Option_t *);
  
  //-----> Simple User Defined Functions:
  void  SetListForTrkCorr(TList *flist)      {this->fListTRKCorr = (TList *) flist->Clone(); }
  void  SetListForNUACorr(TList *flist)      {this->fListNUACorr = (TList *) flist->Clone(); }
  void  SetListForV0MCorr(TList *flist)      {this->fListV0MCorr = (TList *) flist->Clone(); }
  void  SetListForZDCCorr(TList *flist)      {this->fListZDCCorr = (TList *) flist->Clone(); }
  
  //// Event Cuts or Ranges ////
  //void SetPileUpCutParam(Float_t m,Float_t c)  {this->fPileUpSlopeParm = m; this->fPileUpConstParm = c;}
  //void SetCentralityPercentileMin(Double_t cMin) {this->fCentralityMin  = cMin;}
  //void SetCentralityPercentileMax(Double_t cMax) {this->fCentralityMax  = cMax;}
  void SetDetectorforEventPlane(TString sdetEP)  {this->sDetectorForEP  = sdetEP;}
  void SetCentralityEstimator(TString sEstim)    {this->sCentrEstimator = sEstim;}
  void SetFlagAnalyseLambda(Bool_t bL)           {this->bAnalysLambdaPairs =  bL;}
  void SetFlagSkipPileUpCuts(Bool_t bP)          {this->bSkipPileUpCut     =  bP;}
  void SetFlagSkipAnalysis(Bool_t bA)            {this->bSkipNestedLoop    =  bA;}  
  void SetVzRangeMin(Double_t vzMin)             {this->fMinVzCut       =  vzMin;}
  void SetVzRangeMax(Double_t vzMax)             {this->fMaxVzCut       =  vzMax;}

  /// Setter for ZDC
  void SetUseZDCSpectatorPlane(Bool_t bL)        {this->bUseZDCSpectatorPlane = bL;}
  void SetTypeOfRecentering(Int_t nType)        {this->gTypeOfRecentering = nType;}
  
  /// Track Cut Ranges ///
  void SetTrackCutNclusterMin(Int_t nclTPC)   {this->fTPCclustMin   = nclTPC;}
  void SetFlagUseKinkTracks(Bool_t bKink)     {this->bUseKinkTracks = bKink;}
  void SetCumulantHarmonic(Int_t harm)        {this->gHarmonic      = harm;}
  void SetTrackCutdEdxMin(Float_t mindedx)    {this->fTPCdEdxMin    = mindedx;}
  void SetTrackCutChi2Min(Double_t chi2Min)   {this->fTrkChi2Min    = chi2Min;}
  void SetTrackCutChi2Max(Double_t chi2Max)   {this->fTrkChi2Max    = chi2Max;}
  void SetNSigmaCutTPC(Double_t nSigTPC)      {this->fNSigmaTPCCut  = nSigTPC;}
  void SetNSigmaCutTOF(Double_t nSigTOF)      {this->fNSigmaTOFCut  = nSigTOF;}
  void SetDCAXYRangeMax(Double_t dcaxy)       {this->fDCAxyMax      = dcaxy;}
  void SetDCAZRangeMax(Double_t dcaz)         {this->fDCAzzMax      = dcaz;}
  void SetEtaRangeMin(Double_t etamin)        {this->fMinEtaCut     = etamin;}
  void SetEtaRangeMax(Double_t etamax)        {this->fMaxEtaCut     = etamax;}
  void SetPtRangeMin(Double_t ptLow)          {this->fMinPtCut      = ptLow;}
  void SetPtRangeMax(Double_t ptHigh)         {this->fMaxPtCut      = ptHigh;}
  void SetParticlePID(Int_t partID)           {this->gParticleID    = partID;}
  void SetTrackChi2Max(Double_t chi2)         {this->fTrkChi2Max    = chi2;}
  void SetFilterBit(Int_t fb)                 {this->fFilterBit     = fb;}
  void SetEtaGapNeg(Double_t etaL)            {this->fEtaGapNeg     = etaL;}
  void SetEtaGapPos(Double_t etaH)            {this->fEtaGapPos     = etaH;}

  //Shi use QC to calculate v2{2} and v2{4} for ZDC resolution
  void CalculateVnfromQC(Double_t cent);
  void ResetEventByEventQuantities();
  
  //V0 Chun zheng: Cuts on V0 and its daughters
  void SetV0PtMin(Double_t v0PtMin)                                    {this->fV0PtMin                                   = v0PtMin;}
  void SetV0CPAMin(Double_t v0CPAMin)                                  {this->fV0CPAMin                                 = v0CPAMin;}
  void SetV0RapidityMax(Double_t v0RapidityMax)                        {this->fV0RapidityMax                       = v0RapidityMax;}
  void SetV0DecayLengthMax(Double_t v0DecayLengthMax)                  {this->fV0DecayLengthMax                 = v0DecayLengthMax;}
  void SetV0DecayLengthMin(Double_t v0DecayLengthMin)                  {this->fV0DecayLengthMin                 = v0DecayLengthMin;}
  void SetV0DCAToPrimVtxMax(Double_t v0DCAToPrimVtxMax)                {this->fV0DCAToPrimVtxMax               = v0DCAToPrimVtxMax;}
  void SetV0DcaBetweenDaughtersMax(Double_t v0DcaBetweenDaughtersMax)  {this->fV0DcaBetweenDaughtersMax = v0DcaBetweenDaughtersMax;}
  //V0 Daughter Cut
  void SetDaughtersPIDUseTOF(Bool_t daughterPIDUseTOF)                 {this->fV0DaughterUseTOF                = daughterPIDUseTOF;}
  void SetDaughtersPtMax(Double_t daughtersPtMax)                      {this->fDaughtersPtMax                     = daughtersPtMax;}
  void SetDaughtersEtaMax(Double_t daughtersEtaMax)                    {this->fDaughtersEtaMax                   = daughtersEtaMax;}
  void SetDaughtersNsigma(Double_t daughtersNsigma)                    {this->fDaughtersNsigma                   = daughtersNsigma;}
  void SetDaughtersTPCNclsMin(Double_t daughtersTPCNclsMin)            {this->fDaughtersTPCNclsMin           = daughtersTPCNclsMin;}
  void SetDaughtersDCAToPrimVtxMin(Double_t daughtersDCAToPrimVtxMin)  {this->fDaughtersDCAToPrimVtxMin = daughtersDCAToPrimVtxMin;}
  //Lambda Mass Cut
  void SetMassMean(Double_t massMean)                                  {this->fMassMean                                 = massMean;}
  void SetLambdaMassCut(Double_t lambdaMassCut)                        {this->fLambdaMassCut                       = lambdaMassCut;}
  //Check PID Flow
  void IsCheckPIDFlow(Bool_t checkPIDFlow)                             {this->bCheckPIDFlow                         = checkPIDFlow;}          



  
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
  TList                 *fListTRKCorr;        //  Supplied from AddTask
  TList                 *fListNUACorr;        //  Supplied from AddTask
  TList                 *fListV0MCorr;        //  Supplied from AddTask  
  TList                 *fListZDCCorr;        //  Supplied from AddTask
  TList                 *fListTemp;           //! List holding temp objects
  
  /// Functions for Pile Up Event Removal:
  TF1                   *fV0CutPU;      //!
  TF1                   *fSPDCutPU;     //!
  TF1                   *fMultCutPU;    //!
  TF1                   *fCenCutLowPU;  //!
  TF1                   *fCenCutHighPU; //!

    
  
  
  //Track Variables to be used:
  Int_t                 gHarmonic;  //
  Int_t               gParticleID;  //
  Int_t                fFilterBit;  //
  Int_t              fTPCclustMin;  //
  Int_t             gOldRunNumber;  //!


  Float_t               fMinVzCut;  //
  Float_t               fMaxVzCut;  //
  Float_t               fMinPtCut;  //
  Float_t               fMaxPtCut;  //
  Float_t               fDCAxyMax;  //
  Float_t               fDCAzzMax;  //   
  Float_t              fMinEtaCut;  //
  Float_t              fMaxEtaCut;  //
  Float_t              fEtaGapNeg;  // 
  Float_t              fEtaGapPos;  //
  Float_t             fTrkChi2Min;  //
  Float_t             fTrkChi2Max;  //
  Float_t             fTPCdEdxMin;  //
  Float_t           fNSigmaTPCCut;  // 
  Float_t           fNSigmaTOFCut;  //
  Float_t           fVertexZEvent;  //!

  TString          sDetectorForEP;  //
  TString         sCentrEstimator;  //
  
  Bool_t           bUseKinkTracks;  //
  Bool_t           bSkipPileUpCut;  //   
  Bool_t          bSkipNestedLoop;  //
  Bool_t         bUseV0EventPlane;  //
  Bool_t       bAnalysLambdaPairs;  //

  /// ZDC 
  Bool_t       bUseZDCSpectatorPlane; //
  Int_t        gTypeOfRecentering;
  Bool_t       bRecenterFailOrNot;    //

  /// Chunzheng: V0 (Lambda) Cut parameters:
  Double_t                   fV0PtMin; //!
  Double_t                  fV0CPAMin; //!
  Double_t             fV0RapidityMax; //!
  Double_t          fV0DecayLengthMin; //!
  Double_t          fV0DecayLengthMax; //!
  Double_t         fV0DCAToPrimVtxMax; //!
  Double_t  fV0DcaBetweenDaughtersMax; //!
  //V0 Daughter
  Double_t            fDaughtersPtMax; //!
  Double_t           fDaughtersNsigma; //!
  Double_t           fDaughtersEtaMax; //!
  Double_t       fDaughtersTPCNclsMin; //!
  Double_t  fDaughtersDCAToPrimVtxMin; //!
  Float_t       fV0PosProtonTPCNsigma; //!
  Float_t         fV0NegPionTPCNsigma; //!
  Float_t       fV0NegProtonTPCNsigma; //!
  Float_t         fV0PosPionTPCNsigma; //!
  Bool_t            fV0DaughterUseTOF; //!
  Float_t       fV0PosProtonTOFNsigma; //!
  Float_t         fV0NegPionTOFNsigma; //!
  Float_t       fV0NegProtonTOFNsigma; //!
  Float_t         fV0PosPionTOFNsigma; //!
  //Lambda Mass
  Double_t                  fMassMean; //!
  Double_t             fLambdaMassCut; //!

  std::vector<Int_t>    vecPosEPTrkID; //!
  std::vector<Int_t>    vecNegEPTrkID; //!
  std::vector<Double_t> vecPosEPTrkPhi; //!
  std::vector<Double_t> vecNegEPTrkPhi; //!
  std::vector<Double_t> vecPosEPTrkNUAWgt; //!
  std::vector<Double_t> vecNegEPTrkNUAWgt; //!





  //// Histograms:
  TH1F          *fHistVertexZcm;     //! See Vz Distribution
  TH1F          *fHistAnalysisInfo;  //! information about cuts
  TH1F          *fCentDistBeforCut;  //! Cent before Any Cut
  TH1F          *fCentDistAfterCut;  //! Cent After All Cuts
  TH1F          *fDebugwEventCount;  //! Event after Various cuts

  /// Correction Histograms:
  TH1D         *fHCorrectQNxV0C; //!
  TH1D         *fHCorrectQNyV0C; //!
  TH1D         *fHCorrectQNxV0A; //!
  TH1D         *fHCorrectQNyV0A; //!
  TH1D         *fHCorrectQ3xV0C; //!
  TH1D         *fHCorrectQ3yV0C; //!
  TH1D         *fHCorrectQ3xV0A; //!
  TH1D         *fHCorrectQ3yV0A; //!      
  TH1D         *fHCorrectMCPosChrg; //!
  TH1D         *fHCorrectMCPosPion; //!
  TH1D         *fHCorrectMCPosKaon; //!
  TH1D         *fHCorrectMCPosProt; //!   
  TH1D         *fHCorrectMCNegChrg; //!
  TH1D         *fHCorrectMCNegPion; //!
  TH1D         *fHCorrectMCNegKaon; //!
  TH1D         *fHCorrectMCNegProt; //!
  TH1D         *fHCorrectTPCQnxEtaPos; //!
  TH1D         *fHCorrectTPCQnyEtaPos; //!
  TH1D         *fHCorrectTPCQnxEtaNeg; //!
  TH1D         *fHCorrectTPCQnyEtaNeg; //!
  TH1D         *fHCorrectTPCQ3xEtaPos; //!
  TH1D         *fHCorrectTPCQ3yEtaPos; //!
  TH1D         *fHCorrectTPCQ3xEtaNeg; //!
  TH1D         *fHCorrectTPCQ3yEtaNeg; //!  

  
  TH2F          *fHistTPCPsiNPosPlane;   //!
  TH2F          *fHistTPCPsiNNegPlane;   //!
  TH2F          *fHistTPCPsi3PosPlane;   //!
  TH2F          *fHistTPCPsi3NegPlane;   //!
  TH2F          *fHistTPCPsi4PosPlane;   //!
  TH2F          *fHistTPCPsi4NegPlane;   //!  
  TH2F          *fHistV0CPsiNEventPlane; //!  
  TH2F          *fHistV0APsiNEventPlane; //!
  TH2F          *fHistV0CPsi3EventPlane; //!  
  TH2F          *fHistV0APsi3EventPlane; //!
  TH2F          *fHistTPCPosqVectorvsCent; //!!
  TH2F          *fHistTPCNegqVectorvsCent; //!!  
  TH2F          *fHistV0CDetqVectorvsCent; //!!
  TH2F          *fHistV0ADetqVectorvsCent; //!!

  TH2F          *fHCorrectV0ChWeghts;     //!

  TH2F          *fHistTPConlyVsCL1Before; //!
  TH2F          *fHistTPConlyVsCL1After;  //!
  TH2F          *fHistTPConlyVsV0MBefore; //!
  TH2F          *fHistTPConlyVsV0MAfter;  //!
  TH2F          *fHistCentCL0VsV0MBefore; //!
  TH2F          *fHistCentCL0VsV0MAfter;  //!  
  TH2F          *fHistTPCVsESDTrkBefore;  //!
  TH2F          *fHistTPCVsESDTrkAfter;   //!

    
  TH3F          *fHCorrectNUAChrgPos;    //!
  TH3F          *fHCorrectNUAChrgNeg;    //!
  TH3F          *fHCorrectNUAkPIDPos;    //!
  TH3F          *fHCorrectNUAkPIDNeg;    //!

  /// ZDC recentering inputs
  TH1D          *fHZDCCparameters;          //!
  TH1D          *fHZDCAparameters;          //!
  
  /// ZDC QA
  TH1F      *hZDCCPsiSpectatorPlane;   //!
  TH1F      *hZDCAPsiSpectatorPlane;   //!
  TH1F      *hZDCCAPsiSpectatorPlane;  //!
  
  TH2F          *hZDCCPsivsZDCAPsi;        //!
  TH2F          *hZDCCPsivsTPCPsi2;        //!
  TH2F          *hZDCCPsivsTPCPosPsi2;     //!
  TH2F          *hZDCCPsivsTPCNegPsi2;     //!
  TH2F          *hZDCCPsivsV0CAPsi2;       //!
  TH2F          *hZDCCPsivsV0CPsi2;        //!
  TH2F          *hZDCCPsivsV0APsi2;        //!
	
  TH2F          *hZDCAPsivsTPCPsi2;        //!
  TH2F          *hZDCAPsivsTPCPosPsi2;     //!
  TH2F          *hZDCAPsivsTPCNegPsi2;     //!
  TH2F          *hZDCAPsivsV0CAPsi2;       //!
  TH2F          *hZDCAPsivsV0CPsi2;        //!
  TH2F          *hZDCAPsivsV0APsi2;        //!
	
  TH2F          *hZDCCAPsivsTPCPsi2;       //!
  TH2F          *hZDCCAPsivsTPCPosPsi2;    //!
  TH2F          *hZDCCAPsivsTPCNegPsi2;    //!
  TH2F          *hZDCCAPsivsV0CAPsi2;      //!
  TH2F          *hZDCCAPsivsV0CPsi2;       //!
  TH2F          *hZDCCAPsivsV0APsi2;       //!

  //// Store different types of v2
  TProfile      *hV2TPCvsCent;        //!
  TProfile      *hV2V0CvsCent;        //!
  TProfile      *hV2V0AvsCent;        //!
  TProfile      *hV2V0CAvsCent;       //!
  TProfile      *hV2ZDCCvsCent;       //!
  TProfile      *hV2ZDCAvsCent;       //!
  TProfile      *hV2ZDCCACombinevsCent;  //!
  TProfile      *hV2ZDCCAvsCent;      //!
  
  //// Store variables for QC v2{4} and v2{2}
  TH1D      *hPOIQRe;             //!
  TH1D      *hPOIQIm;             //!
  TH1D      *hPOIMul;          //!
  TProfile      *hFlowQC2Pro;         //!
  TProfile      *hFlowQC4Pro;         //!
  TProfile      *hFlowwCov24Pro;         //!
  
  //// Analysis Result Histograms:
  TProfile      *hAvg3pC112vsCentPP;  //!
  TProfile      *hAvg3pC112vsCentNN;  //!
  TProfile      *hAvg3pC112vsCentOS;  //!
  TProfile      *hAvg3pC123vsCentPP;  //!
  TProfile      *hAvg3pC123vsCentNN;  //!
  TProfile      *hAvg3pC123vsCentOS;  //!  
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
  
  //// Analysis Result Histograms for double ratio ZDC
  TProfile      *hAvg3pC112vsCentPPUsingZDCC;  //!
  TProfile      *hAvg3pC112vsCentNNUsingZDCC;  //!
  TProfile      *hAvg3pC112vsCentOSUsingZDCC;  //!
  TProfile      *hAvg3pC112vsCentPPUsingZDCA;  //!
  TProfile      *hAvg3pC112vsCentNNUsingZDCA;  //!
  TProfile      *hAvg3pC112vsCentOSUsingZDCA;  //!
  TProfile      *hAvg3pC112vsCentPPUsingZDCCACombine;  //!
  TProfile      *hAvg3pC112vsCentNNUsingZDCCACombine;  //!
  TProfile      *hAvg3pC112vsCentOSUsingZDCCACombine;  //!
  TProfile      *hAvg3pC112vsCentPPUsingZDCCA;  //!
  TProfile      *hAvg3pC112vsCentNNUsingZDCCA;  //!
  TProfile      *hAvg3pC112vsCentOSUsingZDCCA;  //!
  
  
  //// Store the EP <Q> vectors:
  TProfile      *hAvgQNXvsCentV0C;  //!
  TProfile      *hAvgQNYvsCentV0C;  //!
  TProfile      *hAvgQ3XvsCentV0C;  //!
  TProfile      *hAvgQ3YvsCentV0C;  //!
  TProfile      *hAvgQNXvsCentV0A;  //!
  TProfile      *hAvgQNYvsCentV0A;  //!
  TProfile      *hAvgQ3XvsCentV0A;  //!
  TProfile      *hAvgQ3YvsCentV0A;  //!

  TProfile     *hTPCPsiNCorrelation;  //!
  TProfile     *hTPCPsi3Correlation;  //!
  TProfile     *hTPCPsi4Correlation;  //!
  TProfile     *hV0CV0APsiNCorrelation;  //!
  TProfile     *hV0CTPCPsiNCorrelation;  //!
  TProfile     *hV0ATPCPsiNCorrelation;  //!
  TProfile     *hV0CV0APsi3Correlation;  //!
  TProfile     *hV0CTPCPsi3Correlation;  //!
  TProfile     *hV0ATPCPsi3Correlation;  //!  
  
  TProfile      *fAvgCos2PsivsCentEtaPos; //!  Nth Harmonic usually N=2 
  TProfile      *fAvgSin2PsivsCentEtaPos; //!
  TProfile      *fAvgCos2PsivsCentEtaNeg; //!
  TProfile      *fAvgSin2PsivsCentEtaNeg; //!
  TProfile      *fAvgCos3PsivsCentEtaPos; //!  3rd Harmonic is hardcoded.
  TProfile      *fAvgSin3PsivsCentEtaPos; //!
  TProfile      *fAvgCos3PsivsCentEtaNeg; //!
  TProfile      *fAvgSin3PsivsCentEtaNeg; //!
  TProfile      *fAvgCos4PsivsCentEtaPos; //!  4th Harmonic is hardcoded.
  TProfile      *fAvgSin4PsivsCentEtaPos; //!
  TProfile      *fAvgCos4PsivsCentEtaNeg; //!
  TProfile      *fAvgSin4PsivsCentEtaNeg; //!   

  
  //TProfile     *fHistVxvsVzMinBias;  //!
  //TProfile     *fHistVyvsVzMinBias;  //!
  TProfile2D     *hAvgV0ChannelsvsVz;   //!


  /// chunzheng: V0 QA Histograms
  TH1D           *fHistV0Pt;              //! Raw V0s' pT
  TH1D           *fHistV0Eta;             //! Raw V0s' eta
  TH1D           *fHistV0DcatoPrimVertex; //! Raw V0s' DcatoPV
  TH1D           *fHistV0CPA;             //! Raw V0s' CPA(cosine pointing angle)
  TH1D           *fHistV0DecayLength;     //! Raw V0s' DecayLength

  ///Lambda-X correlators: Naming Style for Histograms ->
  /// Lambda = L-capital, Proton = P-capital, Anti= A-capital.
  TProfile       *fProfileDelta_Lambda_hPos; //!
  TProfile       *fProfileDelta_Lambda_hNeg; //!
  TProfile       *fProfileDelta_Lambda_Proton;     //!
  TProfile       *fProfileDelta_Lambda_AntiProton; //!
  TProfile       *fProfileDelta_AntiLambda_hPos; //!
  TProfile       *fProfileDelta_AntiLambda_hNeg; //!
  TProfile       *fProfileDelta_AntiLambda_Proton;     //!
  TProfile       *fProfileDelta_AntiLambda_AntiProton; //!
  TProfile       *fProfileGammaTPC_Lambda_hPos; //!
  TProfile       *fProfileGammaTPC_Lambda_hNeg; //!
  TProfile       *fProfileGammaTPC_Lambda_Proton;     //!
  TProfile       *fProfileGammaTPC_Lambda_AntiProton; //!
  TProfile       *fProfileGammaTPC_AntiLambda_hPos; //!
  TProfile       *fProfileGammaTPC_AntiLambda_hNeg; //!
  TProfile       *fProfileGammaTPC_AntiLambda_Proton;     //!
  TProfile       *fProfileGammaTPC_AntiLambda_AntiProton; //!
  
  TProfile       *fProfileGammaV0C_Lambda_hPos; //!
  TProfile       *fProfileGammaV0C_Lambda_hNeg; //!
  TProfile       *fProfileGammaV0C_Lambda_Proton;     //!
  TProfile       *fProfileGammaV0C_Lambda_AntiProton; //!
  TProfile       *fProfileGammaV0C_AntiLambda_hPos; //!
  TProfile       *fProfileGammaV0C_AntiLambda_hNeg; //!
  TProfile       *fProfileGammaV0C_AntiLambda_Proton;     //!
  TProfile       *fProfileGammaV0C_AntiLambda_AntiProton; //!

  TProfile       *fProfileGammaV0A_Lambda_hPos; //!
  TProfile       *fProfileGammaV0A_Lambda_hNeg; //!
  TProfile       *fProfileGammaV0A_Lambda_Proton;     //!
  TProfile       *fProfileGammaV0A_Lambda_AntiProton; //!
  TProfile       *fProfileGammaV0A_AntiLambda_hPos; //!
  TProfile       *fProfileGammaV0A_AntiLambda_hNeg; //!
  TProfile       *fProfileGammaV0A_AntiLambda_Proton;     //!
  TProfile       *fProfileGammaV0A_AntiLambda_AntiProton; //!
  
  TProfile       *fProfileGammaZNC_Lambda_hPos; //!
  TProfile       *fProfileGammaZNC_Lambda_hNeg; //!
  TProfile       *fProfileGammaZNC_Lambda_Proton;     //!
  TProfile       *fProfileGammaZNC_Lambda_AntiProton; //!
  TProfile       *fProfileGammaZNC_AntiLambda_hPos; //!
  TProfile       *fProfileGammaZNC_AntiLambda_hNeg; //!
  TProfile       *fProfileGammaZNC_AntiLambda_Proton;     //!
  TProfile       *fProfileGammaZNC_AntiLambda_AntiProton; //!

  TProfile       *fProfileGammaZNA_Lambda_hPos; //!
  TProfile       *fProfileGammaZNA_Lambda_hNeg; //!
  TProfile       *fProfileGammaZNA_Lambda_Proton;     //!
  TProfile       *fProfileGammaZNA_Lambda_AntiProton; //!
  TProfile       *fProfileGammaZNA_AntiLambda_hPos; //!
  TProfile       *fProfileGammaZNA_AntiLambda_hNeg; //!
  TProfile       *fProfileGammaZNA_AntiLambda_Proton;     //!
  TProfile       *fProfileGammaZNA_AntiLambda_AntiProton; //!
  
  TH1F          *hEmptyPointerFortheList;  //!  

  Double_t           fCurrentVtx[3];//!
 
  /// chunzheng: Lambda QA Info
  TH1D             *fHistLambdaPt[2];      //! [0]:Before the Mass Cut [1]:After the Mass Cut
  TH1D             *fHistLambdaEta[2];     //!
  TH1D             *fHistLambdaPhi[2];
  TH1D             *fHistLambdaDcaToPrimVertex[2]; //!
  TH1D             *fHistLambdaCPA[2];     //!
  TH1D             *fHistLambdaDecayLength[2];     //!
  TH1D             *fHistLambdaMass[2];    //!
  TH1D             *fHistAntiLambdaPt[2];  //!
  TH1D             *fHistAntiLambdaEta[2]; //!
  TH1D             *fHistAntiLambdaPhi[2];
  TH1D             *fHistAntiLambdaDcaToPrimVertex[2]; //!
  TH1D             *fHistAntiLambdaCPA[2]; //!
  TH1D             *fHistAntiLambdaDecayLength[2]; //!
  TH1D             *fHistAntiLambdaMass[2];//!  
  TProfile         *fProfileLambdaMassVsPt[2];     //!
  TProfile         *fProfileAntiLambdaMassVsPt[2]; //!
  ///Check PID Flow
  Bool_t           bCheckPIDFlow;
  TProfile2D       *fProfile2DRawFlowCentPthPos[5];//TPC/V0C/V0A/ZNC/ZNA
  TProfile2D       *fProfile2DRawFlowCentPthNeg[5];  
  TProfile2D       *fProfile2DRawFlowCentPtProton[5];
  TProfile2D       *fProfile2DRawFlowCentPtAntiProton[5];
  TProfile2D       *fProfile2DRawFlowCentPtLambda[6];
  TProfile2D       *fProfile2DRawFlowCentPtAntiLambda[6];//TPCPos/TPCNeg/V0C/V0A/ZNC/ZNA

  
  //// Some more functions:
  void  SetupQAHistograms();
  void  SetupAnalysisHistograms();
  void  SetupPileUpRemovalFunctions();
  void  SetupEventAndTaskConfigInfo();
  void  SetupLambdaPIDFlowHistograms();

  void  GetMCCorrectionHist();
  void  GetV0MCorrectionHist(Int_t run=0);
  void  GetZDCCorrectionHist(Int_t run=0);
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C, Double_t &qnyV0C, Double_t &qnxV0A, Double_t &qnyV0A);
  void  ApplyTPCqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t& qxEtaNeg, Double_t& qyEtaNeg,Double_t& qxEtaPos,Double_t& qyEtaPos);

  ///// Chunzheng: V0 (Lambda)
  Bool_t IsGoodV0(AliAODv0 *aodV0);
  Bool_t IsGoodDaughterTrack(const AliAODTrack *track);
  
  Bool_t   CheckEventIsPileUp2018(AliAODEvent* faod);
  Bool_t   CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck);
  Bool_t   GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A); 
  Bool_t   GetTPCQvectAndRemovePileUp2018(AliAODEvent *faod,Double_t *qnxEtaNeg,Double_t *qnyEtaNeg,Double_t *qnxEtaPos,Double_t *qnyEtaPos,Double_t &multNeg,Double_t &multPos);

 
  Double_t GetNUAWeightForTrack(Double_t fVtxZ=0,Double_t fPhi=0,Double_t fEta=0,Int_t gChrg=1);
  Double_t GetNUAWeightForTrackPID(Double_t fVtxZ=0,Double_t fPhi=0,Double_t fEta=0,Int_t gChrg=1);
  Double_t GetMCEfficiencyWeightForTrack(Double_t fPt=1.0,Int_t gChrg=1,Int_t kPID=0);
  Int_t    GetLambdaCode(const AliAODTrack *pTrack, const AliAODTrack *ntrack); 


  
  AliAnalysisTaskGammaDeltaPID(const AliAnalysisTaskGammaDeltaPID &other);
  AliAnalysisTaskGammaDeltaPID &operator=(const AliAnalysisTaskGammaDeltaPID &other);    
  ClassDef(AliAnalysisTaskGammaDeltaPID, 2) 

};

#endif
