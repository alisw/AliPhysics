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
    if(methodnum == 1) ftrackAnalysis = kTRUE;
    
  }
 private:
  AliESDInputHandler    *fInputHandler;        //!<! Input handler
  AliESDpid             *fPID;                 //!<! ESD pid
  AliESDEvent           *fESDevent;            //!<! ESD event
  AliStack              *fStack;               //!<! MC stack
  AliESDv0              *fV0;                  //!<! ESD v0 - He4 + pi
  TH2F                  *fHistdEdx;            //<   Histogram of Tpc dEdx for pid qa
  TH2F                  *fHistdEdxV0;          //<   Histogram of Tpc dEdx for pid qa
  TH1F                  *fHistNumEvents;       //<   Histogram of number of events
  TH1F			            *fHistTrigger;	 	//<   Histogram of trigger for all events 
  TH1F			            *fHistV0;	 	//<   Histogram of trigger for all V0s 
  TTree                 *aTree, *bTree, *cTree, *dTree, *eTree, *fTree, *gTree;                //<   Tree containing reduced events
  TList                 *fHistogramList;       //<   List of histograms
  TVector3              fPrimaryVertex;       //!<! Vector of primary vertex of collision
  Double_t              fMagneticField;       //!<! Magnetic field
  Int_t                 fNV0Cand;             //!<! Number of V0 candidates in a event
  Bool_t                fPIDCheckOnly;        //< Flag to reduce the task only to PID check for Hypertriton daughters
  Bool_t                fMCtrue;              //< Flag for activating MC analysis (set automatically)
  Bool_t                ftrackAnalysis;        // ..> Bool for different Analysis type --> 1 = on; 0 = 0ff
  AliEventCuts          fEventCuts;           //< 2015 event cuts as advised by PDG (AliEventCuts)
  UInt_t		            fTriggerMask;		//< Triggermask for event cuts
  Int_t		              fPeriod;              //< Data period for centrality selector
  Bool_t                fBetheSplines;        //< Switch between built in BetheSplines and personal Fit
  Double_t              fBetheParamsHe[6];    //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Double_t              fBetheParamsT[6];     //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Float_t               fmDaughterSum, fmLi4, fm4LH, fm4LHe, fm5LHe; 
  Float_t               fp4LH, fpt4LH, fct4LH, fy4LH, fz4LH; 
  Float_t               fp4LHe, fpt4LHe, fct4LHe, fy4LHe, fz4LHe;
  Float_t               fp5LHe, fpt5LHe, fct5LHe, fy5LHe, fz5LHe;
  Float_t               fpSum, fctSum, fptSum, fySum, fz4LLH; 
  Float_t               fpLi4, fptLi4, fyLi4, fzLi4; 
  Float_t               fDCAHe3P, fDCAHe3Pi, fDCAPPi, fDCAHe4Pi, fDCAHe4P, fDCAHe3d, fDCAdPi, fDCATPi, fDCATP, fDCAHe3Pi1;
  Float_t               fPA4LHe, fPA4LH, fPA5LHe, fPA4LLH, fAngle4LLH; 
  Float_t               fpiDca, fpi1Dca, fpiDcaSec, fpi1DcaSec; 
  Float_t               fhe4Dca, fhe3Dca, fpDca, fhe4DcaSec, fhe3DcaSec, fpDcaSec, fdDca, ftDca, fdDcaSec, ftDcaSec; 
  Float_t               fhe4Ncls, fpiNcls, fhe3Ncls, fpi1Ncls, fpNcls, fdNcls, ftNcls, fhe4NclsITS, fpiNclsITS, fhe3NclsITS, fpi1NclsITS, fpNclsITS, fdNclsITS, ftNclsITS;
  Float_t               fhe4DedxSigma, fpiDedxSigma, fhe3DedxSigma, fdDedxSigma, ftDedxSigma; 
  Float_t               fpi1DedxSigma, fpDedxSigma; 
  Float_t               fhe4P, fpiP, fhe3P, fpi1P, fpP, fdP, ftP; 
  Float_t               fhe4Dedx, fpiDedx, fhe3Dedx, fpi1Dedx, fpDedx, fdDedx, ftDedx; 
  Float_t               farmalpha, farmpt; 
  Float_t               ftrig;
  Float_t               fthetaP, fthetaN;
  Float_t               fEtaHe4, fEtaHe3, fEtaP, fEtaPi, fEtaPi1, fEtaD, fEtaT;
  Float_t               fPhiHe4, fPhiHe3, fPhiP, fPhiPi, fPhiPi1, fPhiD, fPhiT;
  Float_t               fGeoLengthHe4, fGeoLengthHe3, fGeoLengthP, fGeoLengthPi, fGeoLengthPi1, fGeoLengthD, fGeoLengthT;
  Float_t               fTOFSignalHe4, fTOFSignalHe3, fTOFSignalP, fTOFSignalPi, fTOFSignalPi1, fTOFSignalD, fTOFSignalT;
  Int_t                 fmc, frunnumber;
  Int_t                 fonTheFly;
  TVector3              fVertexPosition; //< position of primary vertex
  UShort_t              fNumberV0s;      //< number of v0s in event
  Int_t                 fCentrality;     //< centrality of event
  UShort_t              fTrigger;        //< array of Triggers
  TString               fTriggerClasses; //< fired trigger classes

  Bool_t TriggerSelection();
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  Double_t GeoLength(const AliESDtrack& track);
  void dEdxCheck();
  void TrackAnalysis(AliESDtrackCuts trackCutsNuc, AliESDtrackCuts trackCutsP, AliESDtrackCuts trackCutsPi, AliVertexerTracks *vertexer, Double_t *dn, Double_t *dd, Double_t xthiss, Double_t xpp);
  AliAnalysisTaskDoubleHypNucTree(const AliAnalysisTaskDoubleHypNucTree&);
  AliAnalysisTaskDoubleHypNucTree &operator=(const AliAnalysisTaskDoubleHypNucTree&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDoubleHypNucTree, 2);
  /// \endcond
};



#endif
