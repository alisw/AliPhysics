/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */
///////////////////////////////////////////////////////////
// Task To Analyse Gamma Delta Correlator with ESE. 
// Works with 15o and 18q/r. Support to be added for LHC10h 
// PA: Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)
///////////////////////////////////////////////////////////




#ifndef ALIANALYSISTASKGAMMADELTAPIDSAVEQVEC_H
#define ALIANALYSISTASKGAMMADELTAPIDSAVEQVEC_H


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

class    AliAnalysisTaskGammaDeltaPIDSaveQvecEvent;
class    AliVEvent;      
class    AliVVertex;
class    AliAODv0;
class    AliESDEvent;       
class    AliAODEvent;
class    AliAODTrack;   
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;



class AliAnalysisTaskGammaDeltaPIDSaveQvec : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskGammaDeltaPIDSaveQvec();
  AliAnalysisTaskGammaDeltaPIDSaveQvec(const char *name);
  virtual ~AliAnalysisTaskGammaDeltaPIDSaveQvec();
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
  
  AliAnalysisTaskGammaDeltaPIDSaveQvecEvent* fpQvecEvent; //!
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
  
  TH1D *fCMEQRePion[4][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEQImPion[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEMultPion[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
  TH1D *fCMEQ2Re4Pion[2]; //! w^2*cos(4phi)
  TH1D *fCMEQ3Re2Pion[2]; //! w^3*cos(2phi)
  TH1D *fCMEQ2Im4Pion[2]; //! w^2*sin(4phi)
  TH1D *fCMEQ3Im2Pion[2]; //! w^3*sin(2phi)
  TH1D *fCMEw0Pion[2];    //! w^0
  TH1D *fCMEw1Pion[2];    //! w^1
  TH1D *fCMEw2Pion[2];    //! w^2
  TH1D *fCMEw3Pion[2];    //! w^3
  TH1D *fCMEw4Pion[2];    //! w^4
  TH1D *fCMEQReKaon[4][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEQImKaon[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEMultKaon[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
  TH1D *fCMEQ2Re4Kaon[2]; //! w^2*cos(4phi)
  TH1D *fCMEQ3Re2Kaon[2]; //! w^3*cos(2phi)
  TH1D *fCMEQ2Im4Kaon[2]; //! w^2*sin(4phi)
  TH1D *fCMEQ3Im2Kaon[2]; //! w^3*sin(2phi)
  TH1D *fCMEw0Kaon[2];    //! w^0
  TH1D *fCMEw1Kaon[2];    //! w^1
  TH1D *fCMEw2Kaon[2];    //! w^2
  TH1D *fCMEw3Kaon[2];    //! w^3
  TH1D *fCMEw4Kaon[2];    //! w^4
  TH1D *fCMEQReProton[4][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEQImProton[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEMultProton[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
  TH1D *fCMEQ2Re4Proton[2]; //! w^2*cos(4phi)
  TH1D *fCMEQ3Re2Proton[2]; //! w^3*cos(2phi)
  TH1D *fCMEQ2Im4Proton[2]; //! w^2*sin(4phi)
  TH1D *fCMEQ3Im2Proton[2]; //! w^3*sin(2phi)
  TH1D *fCMEw0Proton[2];    //! w^0
  TH1D *fCMEw1Proton[2];    //! w^1
  TH1D *fCMEw2Proton[2];    //! w^2
  TH1D *fCMEw3Proton[2];    //! w^3
  TH1D *fCMEw4Proton[2];    //! w^4
  
  //@shi add Qvector for both charge (begin)
  TH1D *fCMEQReBothCharge[2][fCRCnHar]; //! real part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQImBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEMultBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQRePionBothCharge[2][fCRCnHar]; //! real part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQImPionBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEMultPionBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQReKaonBothCharge[2][fCRCnHar]; //! real part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQImKaonBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEMultKaonBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQReProtonBothCharge[2][fCRCnHar]; //! real part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQImProtonBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEMultProtonBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  //@shi add Qvector for both charge (end)
  
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
  
  //// V2 
  TProfile      *hAvgV2TPCvsCent;  //!
  TProfile      *hAvgV2TPCvsCentPion;  //!
  TProfile      *hAvgV2TPCvsCentKaon;  //!
  TProfile      *hAvgV2TPCvsCentProton;  //!
  
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

  // PID Analysis Result Histograms
  TProfile      *hAvg3pC112vsCentOSPionPion;  //!
  TProfile      *hAvg3pC123vsCentOSPionPion;  //!
  TProfile      *hAvgDelta1vsCentOSPionPion;  //!
  TProfile      *hAvgDelta2vsCentOSPionPion;  //!
  TProfile      *hAvgDelta3vsCentOSPionPion;  //!
  TProfile      *hAvgDelta4vsCentOSPionPion;  //!
  TProfile      *hAvg3pC112vsCentOSKaonKaon;  //!
  TProfile      *hAvg3pC123vsCentOSKaonKaon;  //!
  TProfile      *hAvgDelta1vsCentOSKaonKaon;  //!
  TProfile      *hAvgDelta2vsCentOSKaonKaon;  //!
  TProfile      *hAvgDelta3vsCentOSKaonKaon;  //!
  TProfile      *hAvgDelta4vsCentOSKaonKaon;  //!
  TProfile      *hAvg3pC112vsCentOSProtonProton;  //!
  TProfile      *hAvg3pC123vsCentOSProtonProton;  //!
  TProfile      *hAvgDelta1vsCentOSProtonProton;  //!
  TProfile      *hAvgDelta2vsCentOSProtonProton;  //!
  TProfile      *hAvgDelta3vsCentOSProtonProton;  //!
  TProfile      *hAvgDelta4vsCentOSProtonProton;  //!

  TProfile      *hAvg3pC112vsCentOSPionCharge;  //!
  TProfile      *hAvg3pC123vsCentOSPionCharge;  //!
  TProfile      *hAvgDelta1vsCentOSPionCharge;  //!
  TProfile      *hAvgDelta2vsCentOSPionCharge;  //!
  TProfile      *hAvgDelta3vsCentOSPionCharge;  //!
  TProfile      *hAvgDelta4vsCentOSPionCharge;  //!		
  TProfile      *hAvg3pC112vsCentOSKaonCharge;  //!
  TProfile      *hAvg3pC123vsCentOSKaonCharge;  //!
  TProfile      *hAvgDelta1vsCentOSKaonCharge;  //!
  TProfile      *hAvgDelta2vsCentOSKaonCharge;  //!
  TProfile      *hAvgDelta3vsCentOSKaonCharge;  //!
  TProfile      *hAvgDelta4vsCentOSKaonCharge;  //!
  TProfile      *hAvg3pC112vsCentOSProtonCharge;  //!
  TProfile      *hAvg3pC123vsCentOSProtonCharge;  //!
  TProfile      *hAvgDelta1vsCentOSProtonCharge;  //!
  TProfile      *hAvgDelta2vsCentOSProtonCharge;  //!
  TProfile      *hAvgDelta3vsCentOSProtonCharge;  //!
  TProfile      *hAvgDelta4vsCentOSProtonCharge;  //!	
  

  TProfile      *hAvg3pC112vsCentPPPionPion;  //!
  TProfile      *hAvg3pC123vsCentPPPionPion;  //!
  TProfile      *hAvgDelta1vsCentPPPionPion;  //!
  TProfile      *hAvgDelta2vsCentPPPionPion;  //!
  TProfile      *hAvgDelta3vsCentPPPionPion;  //!
  TProfile      *hAvgDelta4vsCentPPPionPion;  //!
  TProfile      *hAvg3pC112vsCentPPKaonKaon;  //!
  TProfile      *hAvg3pC123vsCentPPKaonKaon;  //!
  TProfile      *hAvgDelta1vsCentPPKaonKaon;  //!
  TProfile      *hAvgDelta2vsCentPPKaonKaon;  //!
  TProfile      *hAvgDelta3vsCentPPKaonKaon;  //!
  TProfile      *hAvgDelta4vsCentPPKaonKaon;  //!
  TProfile      *hAvg3pC112vsCentPPProtonProton;  //!
  TProfile      *hAvg3pC123vsCentPPProtonProton;  //!
  TProfile      *hAvgDelta1vsCentPPProtonProton;  //!
  TProfile      *hAvgDelta2vsCentPPProtonProton;  //!
  TProfile      *hAvgDelta3vsCentPPProtonProton;  //!
  TProfile      *hAvgDelta4vsCentPPProtonProton;  //!

  TProfile      *hAvg3pC112vsCentPPPionCharge;  //!
  TProfile      *hAvg3pC123vsCentPPPionCharge;  //!
  TProfile      *hAvgDelta1vsCentPPPionCharge;  //!
  TProfile      *hAvgDelta2vsCentPPPionCharge;  //!
  TProfile      *hAvgDelta3vsCentPPPionCharge;  //!
  TProfile      *hAvgDelta4vsCentPPPionCharge;  //!		
  TProfile      *hAvg3pC112vsCentPPKaonCharge;  //!
  TProfile      *hAvg3pC123vsCentPPKaonCharge;  //!
  TProfile      *hAvgDelta1vsCentPPKaonCharge;  //!
  TProfile      *hAvgDelta2vsCentPPKaonCharge;  //!
  TProfile      *hAvgDelta3vsCentPPKaonCharge;  //!
  TProfile      *hAvgDelta4vsCentPPKaonCharge;  //!
  TProfile      *hAvg3pC112vsCentPPProtonCharge;  //!
  TProfile      *hAvg3pC123vsCentPPProtonCharge;  //!
  TProfile      *hAvgDelta1vsCentPPProtonCharge;  //!
  TProfile      *hAvgDelta2vsCentPPProtonCharge;  //!
  TProfile      *hAvgDelta3vsCentPPProtonCharge;  //!
  TProfile      *hAvgDelta4vsCentPPProtonCharge;  //!	

  
  TProfile      *hAvg3pC112vsCentNNPionPion;  //!
  TProfile      *hAvg3pC123vsCentNNPionPion;  //!
  TProfile      *hAvgDelta1vsCentNNPionPion;  //!
  TProfile      *hAvgDelta2vsCentNNPionPion;  //!
  TProfile      *hAvgDelta3vsCentNNPionPion;  //!
  TProfile      *hAvgDelta4vsCentNNPionPion;  //!
  TProfile      *hAvg3pC112vsCentNNKaonKaon;  //!
  TProfile      *hAvg3pC123vsCentNNKaonKaon;  //!
  TProfile      *hAvgDelta1vsCentNNKaonKaon;  //!
  TProfile      *hAvgDelta2vsCentNNKaonKaon;  //!
  TProfile      *hAvgDelta3vsCentNNKaonKaon;  //!
  TProfile      *hAvgDelta4vsCentNNKaonKaon;  //!
  TProfile      *hAvg3pC112vsCentNNProtonProton;  //!
  TProfile      *hAvg3pC123vsCentNNProtonProton;  //!
  TProfile      *hAvgDelta1vsCentNNProtonProton;  //!
  TProfile      *hAvgDelta2vsCentNNProtonProton;  //!
  TProfile      *hAvgDelta3vsCentNNProtonProton;  //!
  TProfile      *hAvgDelta4vsCentNNProtonProton;  //!

  TProfile      *hAvg3pC112vsCentNNPionCharge;  //!
  TProfile      *hAvg3pC123vsCentNNPionCharge;  //!
  TProfile      *hAvgDelta1vsCentNNPionCharge;  //!
  TProfile      *hAvgDelta2vsCentNNPionCharge;  //!
  TProfile      *hAvgDelta3vsCentNNPionCharge;  //!
  TProfile      *hAvgDelta4vsCentNNPionCharge;  //!		
  TProfile      *hAvg3pC112vsCentNNKaonCharge;  //!
  TProfile      *hAvg3pC123vsCentNNKaonCharge;  //!
  TProfile      *hAvgDelta1vsCentNNKaonCharge;  //!
  TProfile      *hAvgDelta2vsCentNNKaonCharge;  //!
  TProfile      *hAvgDelta3vsCentNNKaonCharge;  //!
  TProfile      *hAvgDelta4vsCentNNKaonCharge;  //!
  TProfile      *hAvg3pC112vsCentNNProtonCharge;  //!
  TProfile      *hAvg3pC123vsCentNNProtonCharge;  //!
  TProfile      *hAvgDelta1vsCentNNProtonCharge;  //!
  TProfile      *hAvgDelta2vsCentNNProtonCharge;  //!
  TProfile      *hAvgDelta3vsCentNNProtonCharge;  //!
  TProfile      *hAvgDelta4vsCentNNProtonCharge;  //!	

					
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


  
  
  AliAnalysisTaskGammaDeltaPIDSaveQvec(const AliAnalysisTaskGammaDeltaPIDSaveQvec &other);
  AliAnalysisTaskGammaDeltaPIDSaveQvec &operator=(const AliAnalysisTaskGammaDeltaPIDSaveQvec &other);    
  ClassDef(AliAnalysisTaskGammaDeltaPIDSaveQvec, 2) 

};

#endif
