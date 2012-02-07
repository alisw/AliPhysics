#ifndef ALIANALYSISTASKMULTIDIELECTRON_H
#define ALIANALYSISTASKMULTIDIELECTRON_H
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#####################################################
//#                                                   # 
//#        Basic Analysis task for Dielectron         #
//#          single event analysis                    #
//#                                                   #
//#  by WooJin J. Park, GSI / W.J.Park@gsi.de         #
//#     Ionut C. Arsene, GSI / I.C.Arsene@gsi.de      #
//#     Magnus Mager, CERN / Magnus.Mager@cern.ch     #
//#     Jens Wiechula, Uni HD / Jens.Wiechula@cern.ch #
//#                                                   #
//#####################################################

#include "TList.h"

#include "AliAnalysisTaskSE.h"

// #include "AliDielectronPID.h"

//class AliDielectron;
class AliDielectronTaku;
class TH1D;
class TH2D;
class AliAnalysisCuts;
class AliTriggerAnalysis;
class AliDielectronEventCuts;
class AliDielectronVarManager;
class AliESDtrackCuts;

//class AliFlowEventSimple;
//class AliFlowLYZEventPlane;


class AliAnalysisTaskMultiDielectronNewTaku : public AliAnalysisTaskSE {
  
public:
  AliAnalysisTaskMultiDielectronNewTaku();
  AliAnalysisTaskMultiDielectronNewTaku(const char *name, AliDielectronEventCuts *fCutsEvent);
  virtual ~AliAnalysisTaskMultiDielectronNewTaku(){  }

  virtual void UserExec(Option_t *option);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();
  //temporary
//   virtual void NotifyRun(){AliDielectronPID::SetCorrVal((Double_t)fCurrentRunNumber);}
  
  void UsePhysicsSelection(Bool_t phy=kTRUE) {fSelectPhysics=phy;}
  void SetTriggerMask(UInt_t mask) {fTriggerMask=mask;}
  UInt_t GetTriggerMask() const { return fTriggerMask; }

  void SetEventFilter(AliAnalysisCuts * const filter) {fEventFilter=filter;}
  void SetTriggerOnV0AND(Bool_t v0and=kTRUE)    { fTriggerOnV0AND=v0and;    }
  void SetRejectPileup(Bool_t pileup=kTRUE)     { fRejectPileup=pileup;     }
  void AddDielectron(AliDielectronTaku * const die) { fListDielectron.Add(die); }
  void SetBranches(TTree *t);
  void FillEvent(AliVEvent * const ev);
  void MomentumEnergyMatch(const AliVParticle *track, double *par);
  
protected:
  enum {kAllEvents=0, kSelectedEvents, kV0andEvents, kFilteredEvents, kPileupEvents, kNbinsEvent};
  TList fListDielectron;             // List of dielectron framework instances
  TList fListHistos;                 //! List of histogram manager lists in the framework classes
  TList fListTree;                   //! List of trees manager lists in the framework classes
  TList fListCF;                     //! List with CF Managers
  TTree *fTree;

  Bool_t fSelectPhysics;             // Whether to use physics selection
  UInt_t fTriggerMask;               // Event trigger mask
  Bool_t fTriggerOnV0AND;            // if to trigger on V0and
  Bool_t fRejectPileup;              // pileup rejection wanted
  
  AliTriggerAnalysis *fTriggerAnalysis; //! trigger analysis class

  AliAnalysisCuts *fEventFilter;     // event filter
  AliDielectronEventCuts *fCutsEvent;

  AliESDtrackCuts *fCutsMother;   
  
  TH1D *fEventStat;                  //! Histogram with event statistics
  TH1D *fEvent;
  TH2D *fdEdXvsPt;
  TH2D *fdEdXnSigmaElecvsPt;
  TH2D *fTOFbetavsPt;
  TH2D *fTOFnSigmaElecvsPt;
  TH2D *fTPCcrossedRowsvsPt;
  TH2D *fTPCchi2vsPt;

  Double_t fgValues[AliDielectronVarManager::kNMaxValues];
  Double_t fgData[AliDielectronVarManager::kNMaxValues][1000];

  AliAnalysisTaskMultiDielectronNewTaku(const AliAnalysisTaskMultiDielectronNewTaku &c);
  AliAnalysisTaskMultiDielectronNewTaku& operator= (const AliAnalysisTaskMultiDielectronNewTaku &c);

  TString fName;

  Double_t fkTriggerMask;
  Int_t fkTriggerCent;
  Int_t fNEvent;
  Double_t fkNCut;
  Double_t fkRunNumber;
  Double_t fkCentrality;
  Double_t fkXvPrim;
  Double_t fkYvPrim;
  Double_t fkZvPrim;
  Double_t fkXRes;
  Double_t fkYRes;
  Double_t fkZRes;
  Double_t fkNTrk;
  Double_t fkTracks;
  Double_t fkNacc;
  Double_t fkNaccTrcklts;
  Double_t fkNch;
  Double_t fkZDCN1E;
  Double_t fkZDCP1E;
  Double_t fkZDCN2E;
  Double_t fkZDCP2E;
  Double_t fkV0A;
  Double_t fkV0C;

  Double_t fkRP;
  Double_t fkRPQx;
  Double_t fkRPQy;
  Double_t fkRPsub1;
  Double_t fkRPsub1Qx;
  Double_t fkRPsub1Qy;
  Double_t fkRPsub2;
  Double_t fkRPsub2Qx;
  Double_t fkRPsub2Qy;
  Double_t fkQsubRes;

  TObjArray *fkTriggerInfo;
  TObjString *fkTrigName;

  Int_t   fkNPar ;

  //  AliFlowEventSimple*    fFlowEvent;  // input event
  //  AliFlowLYZEventPlane*  fLyzEp;         //LYZ EP object

  TVector2  *fQsum;                           // flow vector sum
  Double_t  fQ2sum;                           // flow vector sum squared
  Double_t  fMag;

  ClassDef(AliAnalysisTaskMultiDielectronNewTaku, 1); //Analysis Task handling multiple instances of AliDielectron
};
#endif
