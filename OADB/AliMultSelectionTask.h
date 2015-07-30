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
    
    AliMultSelectionTask();
    AliMultSelectionTask(const char *name);
    virtual ~AliMultSelectionTask();
    
    Bool_t IsEventSelected( AliESDEvent* lEvent );
    //Bool_t SetupRun( AliVEvent* lEvent );
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //---------------------------------------------------------------------------------------

private:
    
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

    Float_t  fZncEnergy;          //!  ZNC Energy
    Float_t  fZpcEnergy;          //!  ZPC Energy
    Float_t  fZnaEnergy;          //!  ZNA Energy
    Float_t  fZpaEnergy;          //!  ZPA Energy
    Float_t  fZem1Energy;         //!  ZEM1 Energy
    Float_t  fZem2Energy;         //!  ZEM2 Energy
    Double_t fZnaTower;           //! common PMT of ZNA
    Double_t fZncTower;           //! common PMT of ZNC
    Double_t fZpaTower;           //! common PMT of ZPA
    Double_t fZpcTower;           //! common PMT of ZPC
    Bool_t   fZnaFired;
    Bool_t   fZncFired;
    Bool_t   fZpaFired;
    Bool_t   fZpcFired;
    
    Int_t    fNTracks;             //!  no. tracks


    //Histograms / Anything else as needed
    TH1D *fHistEventCounter; //!

    AliMultSelectionTask(const AliMultSelectionTask&);            // not implemented
    AliMultSelectionTask& operator=(const AliMultSelectionTask&); // not implemented

    ClassDef(AliMultSelectionTask, 1);
};

#endif
