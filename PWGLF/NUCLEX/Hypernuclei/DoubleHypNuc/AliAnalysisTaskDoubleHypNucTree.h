//--- Task for investigation of the DoubleHyperHydrogen4 ---
//---     Author: Janik Ditzel; janik.ditzel@cern.ch     ---

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+   4LLH decays in the first step to a 4LHe and a pion   +
//+   and in the second step the 4LHe decays to a 3He a    +
//+   proton and a pion.                                   +
//+   We implemented 3 different methods using V0 finder   +
//+   and AliVertexer to reconstruct the 4LLH              +
//+   The first method uses 3 V0s [(4LHe, pi), (3He, pi)   + 
//+   and (p, pi)] and a combination of the tracks.        +
//+   The second method only uses tracks and the           +
//+   AliVertexer Class to build a secondary and tertiary  +
//+   Vertex.                                              +
//+   The third class uses a combination of both: The      +
//+   (4LHe,pi) will be found by the V0 finder and the     +
//+   tertiary vertex will be build by finding the 3       +
//+   tracks and reconstruct the vertex with AliVertexer.  +
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifndef ALIANALYSISTASKDOUBLEHYPNUCTREE_H
#define ALIANALYSISTASKDOUBLEHYPNUCTREE_H

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDpid;
class AliESDtrackCuts;
class AliESDv0;
class AliESDVertex;
class AliESDInputHandler;
class AliESDtrack;

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliEventCuts.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TClonesArray.h>
#include "AliPID.h"
#include "AliVertexerTracks.h"


class AliAnalysisTaskDoubleHypNucTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskDoubleHypNucTree();
  AliAnalysisTaskDoubleHypNucTree(const char *name);
  virtual ~AliAnalysisTaskDoubleHypNucTree();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SelectPIDcheckOnly(Bool_t pidch = kFALSE) {fPIDCheckOnly = pidch;};
  void SetPeriod(Int_t period = 2015) {fPeriod = period;};
  void SetTriggerMask(UInt_t triggerMask = AliVEvent::kINT7) {fTriggerMask = triggerMask;};
  void SetBetheSplines(Bool_t betheSplines = kTRUE ) {fBetheSplines = betheSplines;};
  void SetParamsHe(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsHe[i] = params[i];};
  void SetParamsT(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsT[i] = params[i];};
  void SetMethod(Int_t methodnum = 0) {
    if(methodnum == 0) fV0Analysis = kTRUE;
    if(methodnum == 1) ftrackAnalysis = kTRUE;
    if(methodnum == 2) fV0Combination = kTRUE;
    if(methodnum == 3) fLi4Analysis = kTRUE;
  }
 private:
  AliESDInputHandler    *fInputHandler;        //!<! Input handler
  AliESDpid             *fPID;                 //!<! ESD pid
  AliESDEvent           *fESDevent;            //!<! ESD event
  AliStack              *fStack;               //!<! MC stack
  AliESDv0              *fV0;                  //!<! ESD v0 - He4 + pi
  AliESDv0              *fV01;                  //!<! ESD v0 - H3 + pi
  AliESDv0              *fV02;                  //!<! ESD v0 - p + pi
  TH2F                  *fHistdEdx;            //<   Histogram of Tpc dEdx for pid qa
  TH2F                  *fHistdEdxV0;          //<   Histogram of Tpc dEdx for pid qa
  TH1F                  *fHistNumEvents;       //<   Histogram of number of events
  TH1F			            *fHistTrigger;	 	//<   Histogram of trigger for all events 
  TH1F			            *fHistV0;	 	//<   Histogram of trigger for all V0s 
  TTree                 *fTree;                //<   Tree containing reduced events
  TList                 *fHistogramList;       //<   List of histograms
  TVector3              fPrimaryVertex;       //!<! Vector of primary vertex of collision
  Double_t              fMagneticField;       //!<! Magnetic field
  Int_t                 fNV0Cand;             //!<! Number of V0 candidates in a event
  Bool_t                fPIDCheckOnly;        //< Flag to reduce the task only to PID check for Hypertriton daughters
  Bool_t                fMCtrue;              //< Flag for activating MC analysis (set automatically)
  Bool_t                ftrackAnalysis;        // ..> Bool for different Analysis type --> 1 = on; 0 = 0ff
  Bool_t                fV0Analysis;           // ..> Bool for different Analysis type --> 1 = on; 0 = 0ff
  Bool_t                fV0Combination;        // ..> Bool for different Analysis type --> 1 = on; 0 = 0ff
  Bool_t                fLi4Analysis;        // ..> Bool for different Analysis type --> 1 = on; 0 = 0ff
  AliEventCuts          fEventCuts;           //< 2015 event cuts as advised by PDG (AliEventCuts)
  UInt_t		            fTriggerMask;		//< Triggermask for event cuts
  Int_t		              fPeriod;              //< Data period for centrality selector
  Bool_t                fBetheSplines;        //< Switch between built in BetheSplines and personal Fit
  Double_t              fBetheParamsHe[6];    //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Double_t              fBetheParamsT[6];     //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Float_t               fmHHe4Pi, fmDaughterSum1, fmDaughterSum2, fmLi4, fm4LH; 
  Float_t               fpHHe4Pi, fctHHe4Pi, fptHHe4Pi, fp4LH, fpt4LH; 
  Float_t               fpSum1, fctSum1, fptSum1, fct4LH; 
  Float_t               fpSum2, fctSum2, fptSum2; 
  Float_t               fpLi4, fptLi4;
  Float_t               fdcaHe4Pi, fdcaHe3Pi, fdcaPPi; 
  Float_t               fDcaHe3P, fDcaHe3Pi2, fDcaPPi1;
  Float_t               fcosHe4Pi, fcosHe3Pi, fcosPPi; 
  Float_t               fPA4LHe, fPA4LHe1, fPAHe3He4, fPAPHe4;
  Float_t               fyHHe4pi, fySum1, fySum2; 
  Float_t               fpiDca, fpi1Dca, fpi2Dca, fpiDcaSec, fpi1DcaTert, fpi2DcaTert; 
  Float_t               fhe4Dca, fhe3Dca, fpDca, fhe4DcaSec, fhe3DcaTert, fpDcaTert; 
  Float_t               fhe4Ncls, fpiNcls, fhe3Ncls, fpi1Ncls, fpNcls, fpi2Ncls, fhe4NclsITS, fpiNclsITS, fhe3NclsITS, fpi1NclsITS, fpNclsITS, fpi2NclsITS;
  Float_t               fhe4DedxSigma, fpiDedxSigma, fhe3DedxSigma; 
  Float_t               fpi1DedxSigma, fpDedxSigma, fpi2DedxSigma, ftDedxSigma; 
  Float_t               fhe4P, fpiP, fhe3P, fpi1P, fpP, fpi2P; 
  Float_t               fhe4Dedx, fpiDedx, fhe3Dedx, fpi1Dedx, fpDedx, fpi2Dedx; 
  Float_t               farmalpha, farmpt, farmalpha1, farmpt1, farmalpha2, farmpt2; 
  Float_t               ftrig, fVertDiff;
  Float_t               fthetaP, fthetaN;
  Float_t               fEtaHe4, fEtaHe3, fEtaP, fEtaPi, fEtaPi1, fEtaPi2;
  Float_t               fPhiHe4, fPhiHe3, fPhiP, fPhiPi, fPhiPi1, fPhiPi2;
  Float_t               fGeoLengthHe4, fGeoLengthHe3, fGeoLengthP, fGeoLengthPi, fGeoLengthPi1, fGeoLengthPi2;
  Float_t               fTOFSignalHe4, fTOFSignalHe3, fTOFSignalP, fTOFSignalPi, fTOFSignalPi1, fTOFSignalPi2;
  Int_t                 fz, fmc, frunnumber, fBR;
  Int_t                 fCharge;           //< anti or particle
  Int_t                 fParticleSpecies;  //< particle species
  Int_t                 fonTheFly;
  Int_t                 fonTheFly1;
  Int_t                 fonTheFly2;
  TVector3              fVertexPosition; //< position of primary vertex
  UShort_t              fNumberV0s;      //< number of v0s in event
  Int_t                 fCentrality;     //< centrality of event
  UShort_t              fTrigger;        //< array of Triggers
  TString               fTriggerClasses; //< fired trigger classes

  Bool_t TriggerSelection();
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  Double_t GeoLength(const AliESDtrack& track);
  void dEdxCheck();
  void V0Analysis(AliESDtrackCuts trackCutsV0, Double_t xthiss, Double_t xpp, Bool_t ITSClusterCut, Int_t SecClusters, Int_t TertClusters);
  void TrackAnalysis(AliESDtrackCuts trackCutsSec, AliESDtrackCuts trackCutsTert, AliESDtrackCuts trackCutsPi, AliVertexerTracks *vertexer, AliVertexerTracks *vertexer1, Double_t *dn, Double_t *dd, Double_t xthiss, Double_t xpp, Bool_t ITSClusterCut, Int_t SecClusters, Int_t TertClusters);
  void CombinedAnalysis(AliESDtrackCuts trackCutsV0, AliESDtrackCuts trackCutsTert, AliESDtrackCuts trackCutsPi, AliVertexerTracks *vertexer, AliVertexerTracks *vertexer1, Double_t *dn, Double_t *dd, Double_t xthiss, Double_t xpp, Bool_t ITSClusterCut, Int_t SecClusters, Int_t TertClusters);
  void Li4Analysis(AliESDtrackCuts trackCutsSec);
  AliAnalysisTaskDoubleHypNucTree(const AliAnalysisTaskDoubleHypNucTree&);
  AliAnalysisTaskDoubleHypNucTree &operator=(const AliAnalysisTaskDoubleHypNucTree&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDoubleHypNucTree, 1);
  /// \endcond
};



#endif
