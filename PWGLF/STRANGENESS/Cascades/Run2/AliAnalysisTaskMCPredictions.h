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

#ifndef AliAnalysisTaskMCPredictions_H
#define AliAnalysisTaskMCPredictions_H

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

class AliAnalysisTaskMCPredictions : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskMCPredictions();
    AliAnalysisTaskMCPredictions(const char *name, Int_t lNSmallBinning = 1000, Int_t lNLargeBinning = 4000, Int_t lRebinFactor = 1);
    virtual ~AliAnalysisTaskMCPredictions();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    //Check MC type
    Bool_t IsHijing()  const;
    Bool_t IsDPMJet()  const;
    Bool_t IsEPOSLHC() const;
    
    Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2) const;
    
    void SetDo2pc( Bool_t lOpt = kTRUE ) { fkDo2pc = lOpt; }
    void SetSelectINELgtZERO ( Bool_t lOpt ) { fkSelectINELgtZERO = lOpt; } 
    
//---------------------------------------------------------------------------------------

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms

    //Histograms (Desired objects in this cross-checking task) 
    TH1D *fHistEventCounter; //! histogram for event counting
    TH1D *fHistChargedEta; //! histogram for event counting
    
    Int_t fSmallMultRange;
    Int_t fLargeMultRange;
    Int_t fRebinFactor;
    Bool_t fkSelectINELgtZERO;
    
    //Basic Histograms for counting events as a function of V0M percentiles...
    TH1D *fHistV0MMult;
    TH1D *fHistSPDMult;
    TH2D *fHistNchVsV0MMult;
    TH2D *fHistNchVsSPDMult;
    TH1D *fHistNpart;
    TH2D *fHistNchVsNpart;
    TH1D *fHistB;
    TH2D *fHistNchVsB;
    TH1D *fHistNMPI;
    TH2D *fHistNchVsNMPI;
    
    TH1D *fHistPt[27];              //! for keeping track of base spectra
    TH1D *fHistEta[27];              //! for keeping track of base spectra
    TH2D *fHistPtVsV0MMult[27];     //! for keeping track of base spectra
    TH2D *fHistPtVsSPDMult[27];     //! for keeping track of base spectra
    TH2D *fHistPtVsNpart[27];       //! for keeping track of base spectra
    TH2D *fHistPtVsB[27];           //! for keeping track of base spectra
    TH2D *fHistPtVsNMPI[27];       //! for keeping track of base spectra
  
    Bool_t fkDo2pc;
    Float_t fMinPtTriggerCharged; //for charged trigger
    Float_t fMinPtTriggerXi; //for xi trigger
    Float_t fMinPtTriggerPhi; //for phi trigger
    TH1D *fEtaTriggerCharged;
    TH1D *fEtaTriggerXi;
    TH1D *fEtaTriggerPhi;
    TH3D *fHist3d2pcSE[27]; //!
    TH3D *fHist3d2pcXiSE[27]; //!
    TH3D *fHist3d2pcPhiSE[27]; //!
    
    AliAnalysisTaskMCPredictions(const AliAnalysisTaskMCPredictions&);            // not implemented
    AliAnalysisTaskMCPredictions& operator=(const AliAnalysisTaskMCPredictions&); // not implemented

    ClassDef(AliAnalysisTaskMCPredictions, 1);
};

#endif

