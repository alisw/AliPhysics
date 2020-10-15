/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// author: I.Lofnes, ingrid.mckibben.lofnes@cern.ch
// Analysis task for determining the vdm-scan stability
//

#ifndef ALIANALYSISTASKVDMSTABILITY_H
#define ALIANALYSISTASKVDMSTABILITY_H

#include "TList.h"

#include "AliAnalysisTaskSE.h"


class AliESDEvent;
class TH1D;
class AliAnalysisCuts;
class AliTriggerAnalysis;

class AliAnalysisTaskVdmStability : public AliAnalysisTaskSE {
    
public:
    //Class constructors
    AliAnalysisTaskVdmStability();
    AliAnalysisTaskVdmStability(const char *taskname);
    //class destructor
    virtual         ~AliAnalysisTaskVdmStability();
    //event loop - called for each event
    virtual void    UserExec(Option_t *);
    //define output - called once at beginning of runtime
    virtual void    UserCreateOutputObjects();
    //terminate - called at the end of analysis
    virtual void    Terminate(Option_t *);
    
    void SetNRuns(Int_t n){ fNRuns = n;}
    void SetNCases(Int_t c){ fNSelectionCases = c;}
    void SetFillTTree(Bool_t set){ fFillTTree = set;}
    
private:
    // bits toggled in the fEventTag data member
    enum EventTagBits {kAllEvents=0,
        kSelectedEvents,
        kV0TimeEvents,
        kPileupEvents,
        kGoodZEvents,
        kV0andPUEvents,
        kV0andZEvents,
        kV0andPUandZEvents,
        kNbinsEvent};
    
    //AliESDEvent*    fESD;           //!<! Input ESD event
    TList          fOutputList;         //! output list
    
    UInt_t fRunNumber;                  ///< Run number
    Double_t fVtxZ;                      ///< Vertex Z
    Int_t  fnVtxCont;                   ///< Number of Vtx contributors
    Bool_t  fIsGoodZ;                   ///< event with good z-vertex
    Bool_t fSelectPhysics;              ///< Physics selected event
    Bool_t fIsV0ANDfired;               ///< V0and triggered event
    Bool_t fIsT0fired;                  ///< T0 triggered event
    Bool_t fPileupEvent;                ///< pileup event
    Float_t fTV0A;                      ///< V0 A time info
    Float_t fTV0C;                      ///< V0 C time info
    Bool_t fGoodTime;                   ///< event within V0 time window
    ULong64_t fEventTag;                ///< Event tags
    
    AliESDEvent* fEvent;                 //! ESD event
    TH1D *fEventStatV0;                  //! Histogram with event statistics
    TH1D *fEventStatT0;                  //! Histogram with event statistics
    TTree* fEventTree;                   //! Event tree
    
    Int_t fNRuns;						//
    Int_t fNSelectionCases;				//
    Bool_t fFillTTree;					//
    TH1D* v0_H[21];						//!
	TH1D* t0_H[21];						//!
	TH2D* v0_Timing[21];					//!
	TH2D* t0_Timing[21];					//!
    
    void AddEventTreeVariables(TTree* &tree);
	Bool_t CheckTime(Float_t timeA, Float_t timeC);
	Bool_t BadTimingRun(Int_t run);
	Bool_t CheckZVtx(Double_t zVtx, Int_t nCont, Double_t zCut, Bool_t contCut, Int_t nContCut);
    
    AliAnalysisTaskVdmStability(const AliAnalysisTaskVdmStability &c);
    AliAnalysisTaskVdmStability& operator= (const AliAnalysisTaskVdmStability &c);
    
    ClassDef(AliAnalysisTaskVdmStability,4);
    
};
#endif

