#ifndef AliAnalysisTaskhStrCorr_H
#define AliAnalysisTaskhStrCorr_H

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

//-----------------------------------------------------------------
//      AliAnalysisTaskhStrCorr class
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliESDEvent;
class AliAODEvent;
class AliPhysicsSelection;
class AliCFContainer;

//#include "TString.h"
//#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskhStrCorr : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskhStrCorr();
    AliAnalysisTaskhStrCorr(const char *name);
    virtual ~AliAnalysisTaskhStrCorr();
    
    virtual void	 UserCreateOutputObjects();
    virtual void	 UserExec(Option_t *option);
    virtual void	 Terminate(Option_t *);
    
    //---------------------------------------------------------------------------------------
    //Setters for the V0 Extraction
    void SetV0SelectionMaxChisquare   ( Double_t lParameter ){ fV0Sels[0] = lParameter; }
    void SetV0SelectionDCAFirstToPV   ( Double_t lParameter ){ fV0Sels[1] = lParameter; }
    void SetV0SelectionDCASecondtoPV  ( Double_t lParameter ){ fV0Sels[2] = lParameter; }
    void SetV0SelectionDCAV0Daughters ( Double_t lParameter ){ fV0Sels[3] = lParameter; }
    void SetV0SelectionCosinePA       ( Double_t lParameter ){ fV0Sels[4] = lParameter; }
    void SetV0SelectionMinRadius      ( Double_t lParameter ){ fV0Sels[5] = lParameter; }
    void SetV0SelectionMaxRadius      ( Double_t lParameter ){ fV0Sels[6] = lParameter; }
    //---------------------------------------------------------------------------------------
    void SetCascSelectionMaxChisquare         ( Double_t lParameter ){ fCascSels[0] = lParameter; }
    void SetCascSelectionMinV0ImpactParameter ( Double_t lParameter ){ fCascSels[1] = lParameter; }
    void SetCascSelectionV0MassWindow         ( Double_t lParameter ){ fCascSels[2] = lParameter; }
    void SetCascSelectionDCABachToPV          ( Double_t lParameter ){ fCascSels[3] = lParameter; }
    void SetCascSelectionDCACascadeDaughters  ( Double_t lParameter ){ fCascSels[4] = lParameter; }
    void SetCascSelectionCascadeCosinePA      ( Double_t lParameter ){ fCascSels[5] = lParameter; }
    void SetCascSelectionCascadeMinRadius     ( Double_t lParameter ){ fCascSels[6] = lParameter; }
    void SetCascSelectionCascadeMaxRadius     ( Double_t lParameter ){ fCascSels[7] = lParameter; }
    //---------------------------------------------------------------------------------------
    
private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    
    //Top Structure
    TList	*fOutput;	              //! List of output objects
    
    //Per-Event Histograms
    TH1D* fHistEvent;             //! Event Selection Histogram (no further info)
    
    //Trigger properties analysis
    TH1D *fHistTriggerPt; //overall trigger candidate pT
    
    //Particle spectra vs pT (for trigger analysis)
    TH2D *fHistK0ShortMassVsPt; //!
    TH2D *fHistLambdaMassVsPt; //!
    TH2D *fHistAntiLambdaMassVsPt; //!
    TH2D *fHistXiMinusMassVsPt; //!
    TH2D *fHistXiPlusMassVsPt; //!
    TH2D *fHistOmegaMinusMassVsPt; //!
    TH2D *fHistOmegaPlusMassVsPt; //!
    
    //Particle correlation analysis: same-event correlation functions
    TH3D *fHistCorrFuncPeakK0ShortVsPt; //!
    TH3D *fHistCorrFuncPeakLambdaVsPt; //!
    TH3D *fHistCorrFuncPeakAntiLambdaVsPt; //!
    TH3D *fHistCorrFuncPeakXiMinusVsPt; //!
    TH3D *fHistCorrFuncPeakXiPlusVsPt; //!
    TH3D *fHistCorrFuncPeakOmegaMinusVsPt; //!
    TH3D *fHistCorrFuncPeakOmegaPlusVsPt; //!
    
    TH3D *fHistCorrFuncSideBandK0ShortVsPt; //!
    TH3D *fHistCorrFuncSideBandLambdaVsPt; //!
    TH3D *fHistCorrFuncSideBandAntiLambdaVsPt; //!
    TH3D *fHistCorrFuncSideBandXiMinusVsPt; //!
    TH3D *fHistCorrFuncSideBandXiPlusVsPt; //!
    TH3D *fHistCorrFuncSideBandOmegaMinusVsPt; //!
    TH3D *fHistCorrFuncSideBandOmegaPlusVsPt; //!
    
    //Functions to calculate the
    TF1 *fParametricK0ShortMean;
    TF1 *fParametricK0ShortSigma;
    TF1 *fParametricLambdaMean;
    TF1 *fParametricLambdaSigma;
    TF1 *fParametricXiMean;
    TF1 *fParametricXiSigma;
    
    //Objects Controlling Task Behaviour
    AliPIDResponse *fPIDResponse;     // PID response object
    
    //Objects Controlling Task Behaviour: has to be streamed!
    Double_t  fV0Sels[7];   // Array to store the 7 values for the different V0-related selections
    Double_t  fCascSels[8]; // Array to store the 8 values for the different cascade-related selections
    
    AliAnalysisTaskhStrCorr(const AliAnalysisTaskhStrCorr&);            // not implemented
    AliAnalysisTaskhStrCorr& operator=(const AliAnalysisTaskhStrCorr&); // not implemented
    
    ClassDef(AliAnalysisTaskhStrCorr, 1);
};

#endif
