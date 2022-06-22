/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */
///////////////////////////////////////////////////////////
// Task To Analyse Gamma Delta Correlator with ESE. 
// Works with 15o and 18q/r. Support to be added for LHC10h 
// PA: Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)
///////////////////////////////////////////////////////////




#ifndef ALIANALYSISTASKGAMMADELTAPIDSAVEQVECSIMPLE_H
#define ALIANALYSISTASKGAMMADELTAPIDSAVEQVECSIMPLE_H


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

class    AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple;
class    AliVEvent;      
class    AliVVertex;
class    AliAODv0;
class    AliESDEvent;       
class    AliAODEvent;
class    AliAODTrack;   
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;



class AliAnalysisTaskGammaDeltaPIDSaveQvecSimple : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskGammaDeltaPIDSaveQvecSimple();
  AliAnalysisTaskGammaDeltaPIDSaveQvecSimple(const char *name);
  virtual ~AliAnalysisTaskGammaDeltaPIDSaveQvecSimple();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
  virtual void Terminate(Option_t *);
  
  //-----> Simple User Defined Functions:
  void  SetListForTrkCorr(TList *flist)      {this->fListTRKCorr = (TList *) flist->Clone(); }
  void  SetListForNUACorr(TList *flist)      {this->fListNUACorr = (TList *) flist->Clone(); }
  void  SetListForV0MCorr(TList *flist)      {this->fListV0MCorr = (TList *) flist->Clone(); }
  
  void  SetupQvecSavingObjects();
  void  CalculateCMESPPP();
  //// Event Cuts or Ranges ////
  //void SetPileUpCutParam(Float_t m,Float_t c)  {this->fPileUpSlopeParm = m; this->fPileUpConstParm = c;}
  //void SetCentralityPercentileMin(Double_t cMin) {this->fCentralityMin  = cMin;}
  //void SetCentralityPercentileMax(Double_t cMax) {this->fCentralityMax  = cMax;}
  void SetDetectorforEventPlane(TString sdetEP)  {this->sDetectorForEP  = sdetEP;}
  void SetCentralityEstimator(TString sEstim)    {this->sCentrEstimator = sEstim;}
  void SetFlagSkipPileUpCuts(Bool_t bP)          {this->bSkipPileUpCut     =  bP;}
  void SetFlagSkipAnalysis(Bool_t bA)            {this->bSkipNestedLoop    =  bA;}  
  void SetVzRangeMin(Double_t vzMin)             {this->fMinVzCut       =  vzMin;}
  void SetVzRangeMax(Double_t vzMax)             {this->fMaxVzCut       =  vzMax;}


  
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
  TList                 *fTempList;           //! Holding Temperary Histograms
  TList                 *fEventList;          //!
  TTree                 *treeEvent;           //!
  
  AliAnalysisTaskGammaDeltaPIDSaveQvecEventSimple* fpQvecEvent; //!
  /// Functions for Pile Up Event Removal:
  TF1                   *fV0CutPU;      //!
  TF1                   *fSPDCutPU;     //!
  TF1                   *fMultCutPU;    //!
  TF1                   *fCenCutLowPU;  //!
  TF1                   *fCenCutHighPU; //!

  // variables for Qvec
  const static Int_t fCRCnHar = 2;
  const static Int_t fCMEnEtaBin = 4;
  Double_t fCRCEtaMin = -0.8;
  Double_t fCRCEtaMax = 0.8;
  
  TH1D *fCMEQRe[4][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEQIm[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEMult[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
  TH1D *fCMEQ2Re4[2]; //! w^2*cos(4phi)
  TH1D *fCMEQ3Re2[2]; //! w^3*cos(2phi)
  TH1D *fCMEQ2Im4[2]; //! w^2*sin(4phi)
  TH1D *fCMEQ3Im2[2]; //! w^3*sin(2phi)
  TH1D *fCMEw0[2];    //! w^0
  TH1D *fCMEw1[2];    //! w^1
  TH1D *fCMEw2[2];    //! w^2
  TH1D *fCMEw3[2];    //! w^3
  TH1D *fCMEw4[2];    //! w^4
  
  
  // QA hist for ZDC-C and ZDC-A tower energy
  //const static Int_t nBinQAZDCAvgTowEnergyFraction = 18;
  //TProfile *fZDCCAvgTow1EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  //TProfile *fZDCCAvgTow2EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  //TProfile *fZDCCAvgTow3EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  //TProfile *fZDCCAvgTow4EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  //TProfile *fZDCAAvgTow1EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  //TProfile *fZDCAAvgTow2EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  //TProfile *fZDCAAvgTow3EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  //TProfile *fZDCAAvgTow4EnergyFraction[nBinQAZDCAvgTowEnergyFraction]; //0-90 centrality
  

  
  
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
  TH3F          *fHCorrectNUAkPIDPosPion;    //!
  TH3F          *fHCorrectNUAkPIDNegPion;    //!
  TH3F          *fHCorrectNUAkPIDPosKaon;    //!
  TH3F          *fHCorrectNUAkPIDNegKaon;    //!
  TH3F          *fHCorrectNUAkPIDPosProton;    //!
  TH3F          *fHCorrectNUAkPIDNegProton;    //!

					
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

  
  
  //// Some more functions:
  void  SetupQAHistograms();
  void  SetupAnalysisHistograms();
  void  SetupPileUpRemovalFunctions();
  void  SetupEventAndTaskConfigInfo();

  
  void  GetMCCorrectionHist();
  void  GetV0MCorrectionHist(Int_t run=0);
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  ApplyV0XqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t &qnxV0C, Double_t &qnyV0C, Double_t &qnxV0A, Double_t &qnyV0A);
  void  ApplyTPCqVectRecenter(Float_t fCent,Int_t gPsiN,Double_t& qxEtaNeg, Double_t& qyEtaNeg,Double_t& qxEtaPos,Double_t& qyEtaPos);
  
  
  Bool_t   CheckEventIsPileUp2018(AliAODEvent* faod);
  Bool_t   CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck);
  Bool_t   GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A); 
  Bool_t   GetGainCorrectedV0Qvector(AliAODEvent *faod,Double_t fVtxZ,Int_t gPsiN,Double_t &qnxV0C,Double_t &qnyV0C,Double_t &qnxV0A,Double_t &qnyV0A,Double_t &sumMultV0C,Double_t &sumMultV0A); 
  Bool_t   GetTPCQvectAndRemovePileUp2018(AliAODEvent *faod,Double_t *qnxEtaNeg,Double_t *qnyEtaNeg,Double_t *qnxEtaPos,Double_t *qnyEtaPos,Double_t &multNeg,Double_t &multPos);

 
  Double_t GetNUAWeightForTrack(Double_t fVtxZ=0,Double_t fPhi=0,Double_t fEta=0,Int_t gChrg=1);
  Double_t GetNUAWeightForTrackPID(Double_t fVtxZ=0,Double_t fPhi=0,Double_t fEta=0,Int_t gChrg=1);
  Double_t GetNUAWeightForTrackPID(Double_t fVtxZ=0,Double_t fPhi=0,Double_t fEta=0,Int_t gChrg=1,Int_t gPID=0);
  Double_t GetMCEfficiencyWeightForTrack(Double_t fPt=1.0,Int_t gChrg=1,Int_t kPID=0);


  
  
  AliAnalysisTaskGammaDeltaPIDSaveQvecSimple(const AliAnalysisTaskGammaDeltaPIDSaveQvecSimple &other);
  AliAnalysisTaskGammaDeltaPIDSaveQvecSimple &operator=(const AliAnalysisTaskGammaDeltaPIDSaveQvecSimple &other);    
  ClassDef(AliAnalysisTaskGammaDeltaPIDSaveQvecSimple, 1) 

};

#endif
