/* Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// authors: I.Lofnes, ingrid.mckibben.lofnes@cern.ch
//          H.Degenhardt, hermann.franz.degenhardt@cern.ch
// Analysis task for determining the vdm-scan stability
//

#ifndef ALIANALYSISTASKVDMSTABILITY_H
#define ALIANALYSISTASKVDMSTABILITY_H

#include "TList.h"

#include "AliAnalysisTaskSE.h"


class AliESDEvent;
class TH1D;
class AliAnalysisCuts;

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
    void SetDiffTimingCut(Float_t min, Float_t max){fV0TimeDiff_min = min; fV0TimeDiff_max = max;}
    void SetSumTimingCut(Float_t min, Float_t max){fV0TimeSum_min = min; fV0TimeSum_max = max;}
    
private:
    // bits toggled in the fEventTag data member
    enum EventTagBits {
		kAllEvents=0,
        kSelectedEvents,
        kV0TimeEvents,
        kPileupEvents,
        kGoodZEvents,
        kV0andPUEvents,
        kV0andZEvents,
        kV0andPUandZEvents,
        kNbinsEvent
	};
    
    //AliESDEvent*    fESD;           //!<! Input ESD event
    TList          fOutputList;         //! output list
    
    UInt_t fRunNumber;                  ///< Run number
    Double_t fVtxZ;                      ///< Vertex Z
    Int_t  fnVtxCont;                   ///< Number of Vtx contributors
    Bool_t  fIsGoodZ;                   ///< event with good z-vertex
    Bool_t fSelectPhysics;              ///< Physics selected event
    Bool_t fIsV0ANDfired;               ///< V0and triggered event
    Bool_t fIsT0fired;                  ///< T0 triggered event
    Bool_t fIsEMCALfired;               ///< EMCAL triggered event
    Bool_t fIsMUfired;                  ///< MU triggered event
    Bool_t fIsDIMUfired;                ///< DIMU triggered event
    Bool_t fPileupEvent;                ///< pileup event
    Float_t fTV0A;                      ///< V0 A time info
    Float_t fTV0C;                      ///< V0 C time info
    Bool_t fGoodTime;                   ///< event within V0 time window
    
    Bool_t fzCut30;              		///< event with zVertex +-30
    Bool_t fzCut10;              		///< event with zVertex +-10
    Bool_t fzCut30nCont0;               ///< event with zVertex +-30 and nCont > 0
    Bool_t fzCut10nCont0;               ///< event with zVertex +-10 and nCont > 0
    Bool_t fzCut30nCont1;               ///< event with zVertex +-30 and nCont > 1
    Bool_t fzCut10nCont1;               ///< event with zVertex +-10 and nCont > 1
    
    Float_t fV0TimeDiff_min;
    Float_t fV0TimeDiff_max;
    Float_t fV0TimeSum_min;
    Float_t fV0TimeSum_max;
    
    ULong64_t fEventTag;                ///< Event tags
    
    AliESDEvent* fEvent;                 //! ESD event
    TH1D *fEventStatV0;                  //! Histogram with event statistics
    TH1D *fEventStatT0;                  //! Histogram with event statistics
    TTree* fEventTree;                   //! Event tree
    
    Int_t fNRuns;						//
    Int_t fNSelectionCases;				//
    Bool_t fFillTTree;					//
    TH1D* v0_H[25];						//!
	TH1D* t0_H[25];						//!
	TH2D* v0_Timing[25];				//!
	TH2D* t0_Timing[25];				//!
    TH1D* emcal_H[25];					//!
	TH1D* muon_H[25];					//!
	TH1D* dimuon_H[25];					//!
	TH2D* emcal_Timing[25];				//!
	TH2D* muon_Timing[25];				//!
	TH2D* dimuon_Timing[25];			//!
    
    void AddEventTreeVariables(TTree* &tree);
	Bool_t CheckTime(Float_t timeA, Float_t timeC);
	Bool_t BadTimingRun(Int_t run);
	Bool_t CheckZVtx(Double_t zVtx, Int_t nCont, Double_t zCut, Bool_t contCut, Int_t nContCut);
	Bool_t IsSelected(Int_t nSelection);
    
    AliAnalysisTaskVdmStability(const AliAnalysisTaskVdmStability &c);
    AliAnalysisTaskVdmStability& operator= (const AliAnalysisTaskVdmStability &c);
    
    ClassDef(AliAnalysisTaskVdmStability,4);
    
};
#endif

