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

class AliAnalysisTaskHypTritEventTree : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHypTritEventTree();
  AliAnalysisTaskHypTritEventTree(const char *name);
  virtual ~AliAnalysisTaskHypTritEventTree();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SetPidQa(Bool_t pidQa = kTRUE) {fPidQa = pidQa;};

 private:
  AliESDEvent            *fEvent;               //!<! ESD event
  AliESDInputHandler     *fInputHandler;        //!<! Input handler
  AliESDpid              *fPID;                 //!<! ESD pid
  AliReducedHypTritEvent *fReducedEvent;        //<   Reduced event containing he3 and pi
  AliReducedHypTritEvent *fReducedEventMCGen;   //<   Reduced MC event containing he3 and pi
  AliStack               *fStack;               //!<! MC stack
  AliESDv0               *fV0;                  //!<! ESD v0
  TClonesArray           *fV0Array;             //<   Array of v0s in a event
  TH2F                   *fHistdEdx;            //<   Histogram of Tpc dEdx for pid qa
  TH1F                   *fHistNumEvents;       //<   Histogram of number of events
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
  Bool_t                  fMCtrue;              //< Flag for activating MC analysis (set automatically)

  void MCStackLoop(AliStack *stack);
  void SetMomentum(Int_t charge, Bool_t v0Charge);
  void CalculateV0(const AliESDtrack& trackN, const AliESDtrack& trackP, AliPID::EParticleType typeNeg, AliPID::EParticleType typePos);
  Bool_t TriggerSelection();

  AliAnalysisTaskHypTritEventTree(const AliAnalysisTaskHypTritEventTree&);
  AliAnalysisTaskHypTritEventTree &operator=(const AliAnalysisTaskHypTritEventTree&);
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskHypTritEventTree, 2);
  /// \endcond
};
#endif
