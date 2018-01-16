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
    AliAnalysisTaskMCPredictions(const char *name);
    virtual ~AliAnalysisTaskMCPredictions();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    //Check MC type
    Bool_t IsHijing()  const;
    Bool_t IsDPMJet()  const;
    Bool_t IsEPOSLHC() const;
    
//---------------------------------------------------------------------------------------

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms

    //Histograms (Desired objects in this cross-checking task) 
    TH1D *fHistEventCounter; //! histogram for event counting
    
    //Basic Histograms for counting events as a function of V0M percentiles...
    TH1D *fHistV0MMult;
    TH2D *fHistNchVsV0MMult;
    TH1D *fHistNpart;
    TH2D *fHistNchVsNpart;
    TH1D *fHistB;
    TH2D *fHistNchVsB;
    
    TH1D *fHistPt[9];              //! for keeping track of base spectra
    TH2D *fHistPtVsV0MMult[9];     //! for keeping track of base spectra
    TH2D *fHistPtVsNpart[9];       //! for keeping track of base spectra
    TH2D *fHistPtVsB[9];           //! for keeping track of base spectra
    
    AliAnalysisTaskMCPredictions(const AliAnalysisTaskMCPredictions&);            // not implemented
    AliAnalysisTaskMCPredictions& operator=(const AliAnalysisTaskMCPredictions&); // not implemented

    ClassDef(AliAnalysisTaskMCPredictions, 1);
};

#endif
