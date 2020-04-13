#ifndef ALIANALYSISTASKSELC2V0BACHELORTMVA_H
#define ALIANALYSISTASKSELC2V0BACHELORTMVA_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliAnalysisTaskSELc2V0bachelorTMVA.h 61835 2013-04-05 23:07:23Z fprino $ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTPCPIDResponse.h"
#include "AliRDHFCutsLctoV0.h"
#include "AliNormalizationCounter.h"
#include "AliVertexingHFUtils.h"
#include "AliAODRecoCascadeHF.h"

/// \class AliAnalysisTaskSELc2V0bachelorTMVA

class TH1F;
class TH1D;

class AliAnalysisTaskSELc2V0bachelorTMVA : public AliAnalysisTaskSE 
{
  
 public:

  enum EBachelor {
    kBachInvalid = -1,
    kBachFake = 0,
    kBachNoProton = 1,
    kBachPrimary = 2,
    kBachNoLambdaMother = 3,
    kBachDifferentLambdaMother = 4,
    kBachCorrectLambdaMother = 5 };

  enum EK0S {
    kK0SInvalid = -1,
    kK0SFake = 0,
    kK0SNoK0S = 1,
    kK0SWithoutMother = 2,
    kK0SNotFromK0 = 3,
    kK0Primary = 4, 
    kK0NoLambdaMother = 5,
    kK0DifferentLambdaMother = 6,
    kK0CorrectLambdaMother = 7 };    

  enum EAnalysisType { /// enum for setting analysis system/year (for loading profile histograms for multiplicity correction)
     kpPb2013 = 0,
     kpPb2016 = 1,
     kpp2016 = 2,
     kpp2010 = 3};
  
  AliAnalysisTaskSELc2V0bachelorTMVA();
  AliAnalysisTaskSELc2V0bachelorTMVA(const Char_t* name, AliRDHFCutsLctoV0* cutsA,
				 Bool_t useOnTheFly=kFALSE);
  virtual ~AliAnalysisTaskSELc2V0bachelorTMVA();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 
  /// histos
  void FillLc2pK0Sspectrum(AliAODRecoCascadeHF *part, Int_t isLc,
			   Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal,
			   TClonesArray *mcArray, Int_t iLctopK0s, AliAODEvent *aod);

  void MakeAnalysisForLc2prK0S(AliAODEvent *aodEvent,
             TClonesArray *arrayLctopK0s,
			       TClonesArray *mcArray,
			       Int_t &nSelectedAnal, AliRDHFCutsLctoV0 *cutsAnal, 
			       TClonesArray *array3Prong, AliAODMCHeader *aodheader);
 
  /// set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}

  void SetK0sAnalysis(Bool_t a) {fIsK0sAnalysis=a;}
  Bool_t GetK0sAnalysis() const {return fIsK0sAnalysis;}

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }

  void SetIspA(Bool_t a) { fIspA=a; }
  Bool_t GetIspA() { return fIspA; }

  void SetFillOnlySgn(Bool_t a) { fFillOnlySgn=a; }
  Bool_t GetFillOnlySgn() { return fFillOnlySgn; }

  void SetTopoConstraint(Bool_t a) { ftopoConstraint=a; }
  Bool_t GetTopoConstraint() { return ftopoConstraint; }

  void SetCallKFVertexing(Bool_t a) { fCallKFVertexing=a; }
  Bool_t GetCallKFVertexing() { return fCallKFVertexing; }

  void SetKeepingKeepingOnlyHIJINGBkg(Bool_t a) { fKeepingOnlyHIJINGBkg = a;}
  Bool_t GetKeepingOnlyHIJINGBkg() {return fKeepingOnlyHIJINGBkg;}

  void SetKFCutChi2NDF(Float_t a) {fCutKFChi2NDF = a;}
  Float_t GetKFCutChi2NDF() {return fCutKFChi2NDF;}

  void SetKFCutDeviationFromVtx(Float_t a) {fCutKFDeviationFromVtx = a;}
  Float_t GetKFCutDeviationFromVtx() {return fCutKFDeviationFromVtx;}

  void SetKFCutDeviationFromVtxV0(Float_t a) {fCutKFDeviationFromVtxV0 = a;}
  Float_t GetKFCutDeviationFromVtxV0() {return fCutKFDeviationFromVtxV0;}

  void SetKeepingKeepingOnlyPYTHIABkg(Bool_t a) { fKeepingOnlyPYTHIABkg = a;}
  Bool_t GetKeepingOnlyPYTHIABkg() {return fKeepingOnlyPYTHIABkg;}
	
  void SetTriggerMask(ULong64_t c) { fTriggerMask = c;}	
	
  void SetMCNchHisto(TH1F* h){
    if(fHistoMCNch) delete fHistoMCNch;
    fHistoMCNch = new TH1F(*h);
  }
  
  void SetAnalysisType(Int_t mode){
     fAnalysisType = mode;
  }
  
  
  void SetMultVsZProfileLHC13b(TProfile* hprof){
    if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
    fMultEstimatorAvg[0]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC13c(TProfile* hprof){
    if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
    fMultEstimatorAvg[1]=new TProfile(*hprof);
  }

  void SetMultVsZProfileLHC16qt1stBunch(TProfile* hprof){
    if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
    fMultEstimatorAvg[0]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC16qt2ndBunch(TProfile* hprof){
    if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
    fMultEstimatorAvg[1]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC16qt3rdBunch(TProfile* hprof){
    if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
    fMultEstimatorAvg[2]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC16qt4thBunch(TProfile* hprof){
    if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
    fMultEstimatorAvg[3]=new TProfile(*hprof);
  }

  void SetMultVsZProfileLHC16j(TProfile* hprof){
    if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
    fMultEstimatorAvg[0]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC16k(TProfile* hprof){
    if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
    fMultEstimatorAvg[1]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC16l(TProfile* hprof){
    if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
    fMultEstimatorAvg[2]=new TProfile(*hprof);
  }

  void SetMultVsZProfileLHC10b(TProfile* hprof){
    if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
    fMultEstimatorAvg[0]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC10c(TProfile* hprof){
    if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
    fMultEstimatorAvg[1]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC10d(TProfile* hprof){
    if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
    fMultEstimatorAvg[2]=new TProfile(*hprof);
  }
  void SetMultVsZProfileLHC10e(TProfile* hprof){
    if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
    fMultEstimatorAvg[3]=new TProfile(*hprof);
  }


  void SetReferenceMultiplicity(Double_t rmu){fRefMult=rmu;}
  
  


 private:
  
  EBachelor CheckBachelor(AliAODRecoCascadeHF *part, AliAODTrack* bachelor, TClonesArray *mcArray);
  EK0S CheckK0S(AliAODRecoCascadeHF *part, AliAODv0* v0part, TClonesArray *mcArray);
  //EK0S CheckK0S(AliAODRecoCascadeHF *part, AliAODTrack* v0part, TClonesArray *mcArray );
  Int_t FindV0Label(AliAODRecoDecay* v0part, TClonesArray *mcArray) const;
  Int_t FindLcLabel(AliAODRecoCascadeHF* cascade, TClonesArray *mcArray) const;
  Int_t CallKFVertexing(AliAODRecoCascadeHF *cascade, AliAODv0* v0part, AliAODTrack* bach, TClonesArray *mcArray, 
			Double_t* V0KF, Double_t* errV0KF, Double_t* LcKF, Double_t* errLcKF,
			Double_t* distances, Double_t* armPolKF);

  void FillMCHisto(TClonesArray *mcArray);

  AliAnalysisTaskSELc2V0bachelorTMVA(const AliAnalysisTaskSELc2V0bachelorTMVA &source);
  AliAnalysisTaskSELc2V0bachelorTMVA& operator=(const AliAnalysisTaskSELc2V0bachelorTMVA& source); 
  TProfile* GetEstimatorHistogram(const AliVEvent *event);
  
  Bool_t fUseMCInfo;          /// Use MC info
  TList *fOutput;             //!<! User output1: list of trees

  // define the histograms
  TH1F *fCEvents;                    //!<! Histogram to check selected events
  AliPIDResponse *fPIDResponse;      //!<! PID response object
  AliPIDCombined *fPIDCombined;      //!<! combined PID response object
  Bool_t fIsK0sAnalysis;             /// switch between Lpi and K0sp
  AliNormalizationCounter *fCounter; //!<! AliNormalizationCounter on output slot 4
  AliRDHFCutsLctoV0 *fAnalCuts;      /// Cuts - sent to output slot 5
  TList *fListCuts;                  //!<! list of cuts
  TList *fListWeight;                /// list of weights
  TList *fListCounters;              //!<! list of counters on output slot 2
  TList *fListProfiles;              /// list of profiles for z-vtx correction of multiplicity
  Bool_t fUseOnTheFlyV0;             /// flag to analyze also on-the-fly V0 candidates
  Bool_t fIsEventSelected;           /// flag for event selected

  TTree   *fVariablesTreeSgn;        //!<! tree of the candidate variables after track selection (Signal)
  TTree   *fVariablesTreeBkg;        //!<! tree of the candidate variables after track selection (Background)
  Float_t *fCandidateVariables;      //!<! variables to be written to the tree

  Bool_t fIspA;                       /// flag for running on pA

  TH1F* fHistoEvents;                 //!<! histogram with number of events analyzed
  TH1F* fHistoLc;                     //!<! histogram with number of Lc
  TH1F* fHistoLcOnTheFly;             //!<! histogram with number of Lc with on-the-fly V0
  Bool_t fFillOnlySgn;                /// flag to fill only signal (speeding up processing)
  TH1F* fHistoLcBeforeCuts;           //!<! histogram with number of Lc before any cut 
  TH1F* fHistoFiducialAcceptance;     //!<! histogram to check FiducialAcceptance cut
  TH2F* fHistoCodesSgn;               //!<! histogram with codes for bachelor and V0 for signal
  TH2F* fHistoCodesBkg;               //!<! histogram with codes for bachelor and V0 for background
  TH1F* fHistoLcpKpiBeforeCuts;       //!<! histogram number of true Lc-->pKpi (3 prong) before any cut
  AliAODVertex *fVtx1;                /// primary vertex

  TH1D* fHistoDistanceLcToPrimVtx;    //!<! KF: distance Lc vertex from primary vertex   
  TH1D* fHistoDistanceV0ToPrimVtx;    //!<! KF: distance V0 vertex from primary vertex   
  TH1D* fHistoDistanceV0ToLc;         //!<! KF: distance V0 vertex from Lc vertex    

  TH1D* fHistoDistanceLcToPrimVtxSgn; //!<! KF: distance of signal Lc vertex from primary vertex    
  TH1D* fHistoDistanceV0ToPrimVtxSgn; //!<! KF: distance for signal Lc of V0 vertex from primary vertex   
  TH1D* fHistoDistanceV0ToLcSgn;      //!<! KF: distance for signal Lc of V0 vertex from Lc vertex 
         
  TH1D* fHistoVtxLcResidualToPrimVtx; //!<! KF: residual wrt MC of distance Lc vertex from primary vertex (MC - KF)
  TH1D* fHistoVtxV0ResidualToPrimVtx; //!<! KF: residual wrt MC of distance V0 vertex from primary vertex (MC - KF)
  TH1D* fHistoVtxV0ResidualToLc;      //!<! KF: residual wrt MC of distance V0 vertex from Lc vertex (MC - KF)

  TH1D* fHistoMassV0All;              //!<! KF: mass for all V0 reconstructed with KF
  TH1D* fHistoDecayLengthV0All;       //!<! KF: decay length for all V0 reconstructed with KF
  TH1D* fHistoLifeTimeV0All;          //!<! KF: life time for all V0 reconstructed with KF

  TH1D* fHistoMassV0True;             //!<! KF: mass for true V0 reconstructed with KF
  TH1D* fHistoDecayLengthV0True;      //!<! KF: decay length for true V0 reconstructed with KF
  TH1D* fHistoLifeTimeV0True;         //!<! KF: life time for true V0 reconstructed with KF

  TH1D* fHistoMassV0TrueFromAOD;      //!<! KF: AOD mass for true V0 reconstructed with KF

  TH1D* fHistoMassV0TrueK0S;          //!<! KF: mass for true V0 which are really K0S reconstructed with KF
  TH1D* fHistoDecayLengthV0TrueK0S;   //!<! KF: decay length for true V0 which are really K0S reconstructed with KF
  TH1D* fHistoLifeTimeV0TrueK0S;      //!<! KF: life time for true V0 which are really K0S reconstructed with KF

  TH1D* fHistoMassV0TrueK0SFromAOD;   //!<! KF: AOD mass for true V0 which are really K0S reconstructed with KF

  TH1D* fHistoMassLcAll;              //!<! KF: mass for all Lc reconstructed with KF
  TH1D* fHistoDecayLengthLcAll;       //!<! KF: decay length for all Lc reconstructed with KF
  TH1D* fHistoLifeTimeLcAll;          //!<! KF: life time for all Lc reconstructed with KF

  TH1D* fHistoMassLcTrue;             //!<! KF: mass for true cascades reconstructed with KF
  TH1D* fHistoDecayLengthLcTrue;      //!<! KF: decay length for true cascades reconstructed with KF
  TH1D* fHistoLifeTimeLcTrue;         //!<! KF: life time for true cascades reconstructed with KF

  TH1D* fHistoMassLcTrueFromAOD;      //!<! KF: AOD mass for true cascades reconstructed with KF

  TH1D* fHistoMassV0fromLcAll;        //!<! KF: mass of V0 for all cascades reconstructed with KF
  TH1D* fHistoDecayLengthV0fromLcAll; //!<! KF: decay length of V0 for all cascades reconstructed with KF
  TH1D* fHistoLifeTimeV0fromLcAll;    //!<! KF: life time of V0 for all cascades reconstructed with KF

  TH1D* fHistoMassV0fromLcTrue;       //!<! KF: mass of V0 for true cascades reconstructed with KF
  TH1D* fHistoDecayLengthV0fromLcTrue;//!<! KF: decay length of V0 for true cascades reconstructed with KF
  TH1D* fHistoLifeTimeV0fromLcTrue;   //!<! KF: life time of V0 for true cascades reconstructed with KF

  TH1D* fHistoMassLcSgn;              //!<! KF: mass of signal Lc reconstructed with KF
  TH1D* fHistoMassLcSgnFromAOD;       //!<! KF: AOD mass of signal Lc reconstructed with KF
  TH1D* fHistoDecayLengthLcSgn;       //!<! KF: decay length of signal Lc reconstructed with KF
  TH1D* fHistoLifeTimeLcSgn;          //!<! KF: life time of signal Lc reconstructed with KF

  TH1D* fHistoMassV0fromLcSgn;        //!<! KF: mass of V0 for signal Lc reconstructed with KF
  TH1D* fHistoDecayLengthV0fromLcSgn; //!<! KF: decay length of V0 for signal Lc reconstructed with KF
  TH1D* fHistoLifeTimeV0fromLcSgn;    //!<! KF: life time of V0 for signal Lc reconstructed with KF

  TH2D* fHistoKF;                     //!<! KF: V0 code vs Lc code from KF (mass, decaylength, lifetime considered) 
  TH1D* fHistoKFV0;                   //!<! KF: V0 code from KF (mass, decaylength, lifetime considered) 
  TH1D* fHistoKFLc;                   //!<! KF: Lc code from KF (mass, decaylength, lifetime considered) 

  TH2D* fHistoMassKFV0;               //!<! KF: mass vs mass error for V0 from KF  
  TH2D* fHistoDecayLengthKFV0;        //!<! KF: decay length vs decay length error for V0 from KF
  TH2D* fHistoLifeTimeKFV0;           //!<! KF: life time vs life time error for V0 from KF

  TH2D* fHistoMassKFLc;               //!<! KF: mass vs mass error for Lc from KF
  TH2D* fHistoDecayLengthKFLc;        //!<! KF: decay length vs decay length error for Lc from KF
  TH2D* fHistoLifeTimeKFLc;           //!<! KF: life time vs life time error for Lc from KF

  TH2D* fHistoArmenterosPodolanskiV0KF;      //!<! KF: Armeteros-Podolanski plot for all V0 from KF
  TH2D* fHistoArmenterosPodolanskiV0KFSgn;   //!<! KF: Armeteros-Podolanski plot for V0 from signal Lc from KF
  TH2D* fHistoArmenterosPodolanskiV0AOD;     //!<! KF: AOD Armeteros-Podolanski plot for all V0 from KF
  TH2D* fHistoArmenterosPodolanskiV0AODSgn;  //!<! KF: AOD Armeteros-Podolanski plot for V0 from signal Lc from KF

  TList *fOutputKF;             //!<! User output1: list of histograms from KF

  Int_t fmcLabelLc;             /// label of candidate
  Bool_t ftopoConstraint;       /// flag to use topological constraints in KF
  Bool_t fCallKFVertexing;      /// flag to decide whether to call or not KF
  Bool_t fKeepingOnlyHIJINGBkg; /// flag to fill bkg with only candidates that have daughters generated by HIJING (to be used for enriched MC)
  AliVertexingHFUtils* fUtils;  /// AliVertexingHFUtils used to check the generator of a specific candidate
  TH1F* fHistoBackground;       //!<! histo to check the number of candidates with at least one daughter for the injected signal
  Float_t fCutKFChi2NDF;        /// cut for KF on chi2/NDF
  Float_t fCutKFDeviationFromVtx; /// cut for KF on distance to primary vtx
  Float_t fCutKFDeviationFromVtxV0; /// cut for KF on distance to primary vtx for V0
  Int_t fCurrentEvent;              /// current event number - for debug purposes
  Double_t fBField;                   /// magnetic field of current event
  Bool_t fKeepingOnlyPYTHIABkg;       /// flag to allow to use only PYTHIA tracks for background
  TH1F* fHistoMCLcK0SpGen;            //!<! histo with MC Lc --> K0S + p
  TH1F* fHistoMCLcK0SpGenAcc;         //!<! histo with MC Lc --> K0S + p
  TH1F* fHistoMCLcK0SpGenLimAcc;      //!<! histo with MC Lc --> K0S + p

  ULong64_t fTriggerMask;			  /// mask to the trigger word returned by the physics selection

  TF1 *fFuncWeightPythia; //!<! weight function for Pythia vs pPb prod.
  TF1 *fFuncWeightFONLL5overLHC13d3; //!<! weight function for FONLL vs pPb prod.
  TF1 *fFuncWeightFONLL5overLHC13d3Lc; //!<! weight function for FONLL vs pPb prod.
  TH1F* fHistoMCNch;  //!<! histogram with Nch distribution from MC production
 
  Int_t fAnalysisType; /// switch to change system/year in use for loading of mult. estimators
  Int_t fNTracklets; /// tracklet multiplicity in event
  TProfile* fMultEstimatorAvg[4]; /// TProfile with mult vs. Z per period
  Double_t fRefMult; /// reference multiplicity
  
  TList *fListMultiplicityHistograms; //!<! list of multiplicity-related histograms on output slot 8 
  TH2F *fHistNtrVsZvtx; //!<! hist of ntracklets vs. z_vtx
  TH2F *fHistNtrCorrVsZvtx; //!<! hist of corrected ntracklets vs. z_vtx
  
  TH1F* fHistNtrUnCorrEvSel; //!<! hist. of ntracklets for selected events
  TH1F* fHistNtrUnCorrEvWithCand; //!<! hist. of ntracklets for evnts with a candidate
  TH1F* fHistNtrCorrEvSel; //!<! hist. of ntracklets for selected events
  TH1F* fHistNtrCorrEvWithCand; //!<! hist. of ntracklets for evnts with a candidate

  AliNormalizationCounter *fCounterC;           //!<!Counter for normalization, corrected multiplicity
  AliNormalizationCounter *fCounterU;           //!<!Counter for normalization, uncorrected multiplicity
  AliNormalizationCounter *fCounterCandidates;  //!<!Counter for normalization, corrected multiplicity for candidates

  
  /// \cond CLASSIMP    
  ClassDef(AliAnalysisTaskSELc2V0bachelorTMVA, 12); /// class for Lc->p K0
  /// \endcond    
};

#endif

