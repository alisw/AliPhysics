/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliMultSelectionTask_H
#define AliMultSelectionTask_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;


//#include "TString.h"
//#include "AliESDtrackCuts.h"
//#include "AliAnalysisTaskSE.h"

class AliMultSelectionTask : public AliAnalysisTaskSE {
public:
    //Estimators
    enum eEstimator {
        kV0M,
        kV0A,
        kV0C,
        kV0MEq,
        kV0AEq,
        kV0CEq,
        kV0B,
        kV0Apartial,
        kV0Cpartial,
        kV0S,
        kV0SB,
        kAverage
    };
    
    AliMultSelectionTask();
    AliMultSelectionTask(const char *name);
    virtual ~AliMultSelectionTask();
    
    Bool_t IsEventSelected( AliESDEvent* lEvent );

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //---------------------------------------------------------------------------------------

private:
    const Int_t fNEstimators = 12;
    const TString fEstimatorNames[12] = {
        "V0M", //1
        "V0A", //2
        "V0C", //3
        "V0MEq", //4
        "V0AEq", //5
        "V0CEq", //6
        "V0B", //7
        "V0Apartial", //8
        "V0Cpartial", //9
        "V0S", //10
        "V0SB", //11
        "Average" //12
    };
    
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms
    TTree  *fTreeEvent;              //! Output Tree, Events

    //Objects Controlling Task Behaviour
    Bool_t fkCalibration;              //if true, save Calibration object

    AliESDtrackCuts *fESDtrackCuts;
    
    //To perform event selection
    //AliMultSelectionCuts *fMultSelectionCuts; //To perform event selections
    
    //===========================================================================
    //   Variables for Multiplicity Determination
    //===========================================================================
    Float_t fAmplitude_V0A;   //!
    Float_t fAmplitude_V0C;   //!
    Float_t fAmplitude_V0M;   //!
    Float_t fAmplitude_V0Apartial;   //!
    Float_t fAmplitude_V0Cpartial;   //!
    // A.T.
    Float_t fAmplitude_V0A1;   //!
    Float_t fAmplitude_V0A2;   //!
    Float_t fAmplitude_V0A3;   //!
    Float_t fAmplitude_V0A4;   //!
    Float_t fAmplitude_V0C1;   //!
    Float_t fAmplitude_V0C2;   //!
    Float_t fAmplitude_V0C3;   //!
    Float_t fAmplitude_V0C4;   //!

    Float_t fAmplitude_V0AEq; //!
    Float_t fAmplitude_V0CEq; //!
    Float_t fAmplitude_V0MEq; //!
    
    Float_t fCentrality_V0A;         //!
    Float_t fCentrality_V0C;         //!
    Float_t fCentrality_V0M;         //!
    Float_t fCentrality_V0AEq;       //!
    Float_t fCentrality_V0CEq;       //!
    Float_t fCentrality_V0MEq;       //!
    Int_t fRefMultEta5;              //!
    Int_t fRefMultEta8;              //!
    Int_t fRunNumber;                //!

    //Event Characterization Variables - optional
    Bool_t fEvSel_HasAtLeastSPDVertex;      //!
    Bool_t fEvSel_VtxZCut;                  //!
    Bool_t fEvSel_IsNotPileup;              //!
    Bool_t fEvSel_IsNotPileupMV;            //!
    Bool_t fEvSel_IsNotPileupInMultBins;    //!
    Bool_t fEvSel_Triggered;                //!
    Bool_t fEvSel_INELgtZERO;               //! //done with SPD tracklets

    //Other Selections: more dedicated filtering to be studied!
    Int_t   fEvSel_nTracklets;              //!
    Int_t   fEvSel_nTrackletsEta10; //!
    Int_t   fEvSel_nSPDClusters;            //!
    // A.T.
    Int_t   fEvSel_nSPDClusters0;            //!
    Int_t   fEvSel_nSPDClusters1;            //!
    Float_t fEvSel_VtxZ; //! the actual value
    Int_t   fEvSel_nSPDPrimVertices; //! pileup vertices
    Float_t fEvSel_distZ; //! distance between largest vertices
    Int_t   fEvSel_nContributors; //!

    // A.T.
    AliESDtrackCuts* fTrackCuts;  //! optional track cuts
    Float_t  fCentrality_TRK;            // percentile centrality from tracks
    Float_t  fCentrality_CL1;            // percentile centrality from SPD cluster 1
    Float_t  fCentrality_ZNA;            // percentile centrality from ZNA 
    Float_t  fCentrality_ZNC;            // percentile centrality from ZNC 
    Float_t  fCentrality_ZPA;            // percentile centrality from ZPA 
    Float_t  fCentrality_ZPC;            // percentile centrality from ZPC 
    Float_t  fCentrality_ZEMvsZDC;       // percentile centrality from ZEM vs ZDC

    Float_t  fZncEnergy = 0.;          //!  ZNC Energy
    Float_t  fZpcEnergy = 0.;          //!  ZPC Energy
    Float_t  fZnaEnergy = 0.;          //!  ZNA Energy
    Float_t  fZpaEnergy = 0.;          //!  ZPA Energy
    Float_t  fZem1Energy = 0.;         //!  ZEM1 Energy
    Float_t  fZem2Energy = 0.;         //!  ZEM2 Energy
    Double_t fZnaTower = 0.;           //! common PMT of ZNA 
    Double_t fZncTower = 0.;           //! common PMT of ZNC 
    Double_t fZpaTower = 0.;           //! common PMT of ZPA 
    Double_t fZpcTower = 0.;           //! common PMT of ZPC 
    Bool_t   fZnaFired = kFALSE;
    Bool_t   fZncFired = kFALSE;
    Bool_t   fZpaFired = kFALSE;
    Bool_t   fZpcFired = kFALSE;
    
    Int_t    fNTracks = 0;             //!  no. tracks


    //Histograms / Anything else as needed
    TH1D *fHistEventCounter; //!

    AliMultSelectionTask(const AliMultSelectionTask&);            // not implemented
    AliMultSelectionTask& operator=(const AliMultSelectionTask&); // not implemented

    ClassDef(AliMultSelectionTask, 1);
};

#endif
