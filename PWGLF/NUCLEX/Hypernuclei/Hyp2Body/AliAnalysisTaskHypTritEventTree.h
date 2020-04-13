/// \class AliAnalysisTaskHypTritEventTree
/// \brief Hypertriton Analysis in two particle decay channel
///
/// Hypertriton candidates are identified using the on-the-fly V0 finder.
/// Events with Hypertriton candidates are filled in a tree
/// using the AliReducedHypTritEvent class.
///
/// \author Lukas Kreis <lukas.kreis@cern.ch>, GSI
/// \date Sep 1, 2016
#ifndef ALIANALYSISTASKHYPTRITEVENTTREE_H
#define ALIANALYSISTASKHYPTRITEVENTTREE_H

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDpid;
class AliESDtrackCuts;
class AliESDv0;
class AliESDVertex;
class AliESDInputHandler;
class AliESDtrack;
class AliReducedHypTritEvent;

#include "AliReducedHypTritEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliEventCuts.h"

class AliAnalysisTaskHypTritEventTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHypTritEventTree();
  AliAnalysisTaskHypTritEventTree(const char *name);
  virtual ~AliAnalysisTaskHypTritEventTree();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SetPidQa(Bool_t pidQa = kTRUE) {fPidQa = pidQa;};
  void SetUseAnalysisTrackSelection(Bool_t trkSel = kTRUE) {fUseAnalysisTrkSel = trkSel;};
  void SelectPIDcheckOnly(Bool_t pidch = kFALSE) {fPIDCheckOnly = pidch;};
  void SetPeriod(Int_t period = 2015) {fPeriod = period;};
  void SetTriggerMask(UInt_t triggerMask = AliVEvent::kINT7) {fTriggerMask = triggerMask;};
  void SetBetheSplines(Bool_t betheSplines = kTRUE ) {fBetheSplines = betheSplines;};
  void SetParamsHe(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsHe[i] = params[i];};
  void SetParamsT(Double_t params[6]) { for(Int_t i=0; i < 6; i++) fBetheParamsT[i] = params[i];};

 private:
  AliESDInputHandler     *fInputHandler;        //!<! Input handler
  AliESDpid              *fPID;                 //!<! ESD pid
  AliESDEvent            *fESDevent;            //!<! ESD event
  AliReducedHypTritEvent *fReducedEvent;        //<   Reduced event containing he3 and pi
  AliReducedHypTritEvent *fReducedEventMCGen;   //<   Reduced MC event containing he3 and pi
  AliStack               *fStack;               //!<! MC stack
  AliESDv0               *fV0;                  //!<! ESD v0
  TClonesArray           *fV0Array;             //<   Array of v0s in a event
  TH2F                   *fHistdEdx;            //<   Histogram of Tpc dEdx for pid qa
  TH2F                   *fHistdEdxV0;          //<   Histogram of Tpc dEdx for pid qa
  TH1F                   *fHistNumEvents;       //<   Histogram of number of events
  TH1F			 *fHistTrigger;	 	//<   Histogram of trigger for all events 
  TH1F			 *fHistV0;	 	//<   Histogram of trigger for all V0s 
  TH1F                   *fHistMcGen;           //<   Histogram of generated Hypertriton
  TH1F                   *fHistMcRec;           //<   Histogram of reconstructed Hypertriton
  TTree                  *fTree;                //<   Tree containing reduced events
  TTree                  *fTreeMCGen;           //<   Tree containing reduced MC events
  TObjArray              *fMCGenRecArray;       //<   Array used for matching reconstructed MC with generated MC particles in one event
  TList                  *fHistogramList;       //<   List of histograms
  Int_t                   fMCGenRec[40];        //!<! Array containing MC labels of generated hypertriton 40 per event
  TLorentzVector          fMomPos;              //!<! Momentum of positive decay product
  TLorentzVector          fMomNeg;              //!<! Momentum of negative decay product
  TVector3                fPrimaryVertex;       //!<! Vector of primary vertex of collision
  Double_t                fMagneticField;       //!<! Magnetic field
  Int_t                   fNV0Cand;             //!<! Number of V0 candidates in a event
  Int_t                   fMcGenRecCounter;     //!<! Number of matched particles in one MC event
  Bool_t                  fPidQa;               //< Flag for activating pid qa histogram
  Bool_t                  fUseAnalysisTrkSel;   //< Flag to select track for analysis (true) or for good plot (false)
  Bool_t                  fPIDCheckOnly;        //< Flag to reduce the task only to PID check for Hypertriton daughters
  Bool_t                  fMCtrue;              //< Flag for activating MC analysis (set automatically)
  AliEventCuts            fEventCuts;           //< 2015 event cuts as advised by PDG (AliEventCuts)
  UInt_t		  fTriggerMask;		//< Triggermask for event cuts
  Int_t		          fPeriod;              //< Data period for centrality selector
  Bool_t                  fBetheSplines;        //< Switch between built in BetheSplines and personal Fit
  Double_t                fBetheParamsHe[6];    //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t
  Double_t                fBetheParamsT[6];     //< Bethe Aleph He3 Parameter + TPC sigma: [0][i] he3 [2][i] t

  void MCStackLoop(AliStack *stack);
  void SetMomentum(Int_t charge, Bool_t v0Charge);
  void CalculateV0(const AliESDtrack& trackN, const AliESDtrack& trackP, AliPID::EParticleType typeNeg, AliPID::EParticleType typePos);
  Bool_t TriggerSelection();
  Float_t GetInvPtDevFromBC(Int_t b, Int_t c);
  void SetMultiplicity();
  Double_t Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params);
  Bool_t McCuts(const AliReducedHypTritV0& v0, const AliReducedHypTritTrack& he, const AliReducedHypTritTrack& pi);
  Double_t GeoLength(const AliESDtrack& track);
	void SetBetheBlochParams(Int_t runNumber);
  AliAnalysisTaskHypTritEventTree(const AliAnalysisTaskHypTritEventTree&);
  AliAnalysisTaskHypTritEventTree &operator=(const AliAnalysisTaskHypTritEventTree&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskHypTritEventTree, 4);
  /// \endcond
};
#endif
