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
// This is a simple task designed to test AliPPVsMultUtils functionality.
// --- david.dobrigkeit.chinellato@cern.ch
// Also has been added analysis for dN/dEta and Multiplicity + Multiplicity correlations with V0 estimators+ V0 sectors
// ---- hector.bello.martinez@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisMultCorrTaskQA_H
#define AliAnalysisMultCorrTaskQA_H

class TList;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class TVector3;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliESDEvent;
class AliPhysicsSelection;

class AliESDtrack;
class AliPPVsMultUtils;
class AliAnalysisFilter;

#include "AliAnalysisTaskSE.h"

class AliAnalysisMultCorrTaskQA : public AliAnalysisTaskSE {
public:
    AliAnalysisMultCorrTaskQA();
    AliAnalysisMultCorrTaskQA(const char *name);
    virtual ~AliAnalysisMultCorrTaskQA();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    virtual void   SetCuts(AliESDtrackCuts* cuts=0){fCuts = cuts;}
    virtual void   SetTrackFilterSKE(AliAnalysisFilter* trackF) {fTrackFilterSKE = trackF;}
    
    void   LoopESD(AliESDEvent *);
private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of histograms
    AliESDtrackCuts* fCuts; //!
    AliAnalysisFilter* fTrackFilterSKE; //!    

//===========================================================================================
//   Histograms
//===========================================================================================
    
    TH1D *fHistEventCounter; //!
    TH1D *fHistRefMult08; //!
    TH1D *fHistRefMult05; //!
    TH1D *fHistV0M; //!
    TH1D *fHistV0A; //!
    TH1D *fHistV0C; //!
    TH1D *fdNdeta; //!
    TH2D *fcorrRef05Ref08; //!
    TH2D *fcorrV0ARef08; //!
    TH2D *fcorrV0CRef08; //!
    TH2D *fcorrV0MRef08; //!
    TH1D *fHistV0Aamp; //!
    TH1D *fHistV0Camp; //!
    TH1D *fHistV0Mamp; //!
    TH2D *fcorrV0AampRef08; //!
    TH2D *fcorrV0CampRef08; //!
    TH2D *fcorrV0MampRef08; //!
    TH2D *fModulesV0; //!
    
    AliPPVsMultUtils *fPPVsMultUtils;

    AliAnalysisMultCorrTaskQA(const AliAnalysisMultCorrTaskQA&);            // not implemented
    AliAnalysisMultCorrTaskQA& operator=(const AliAnalysisMultCorrTaskQA&); // not implemented

    ClassDef(AliAnalysisMultCorrTaskQA, 10); 
};

#endif

 // - david.dobrigkeit.chinellato@cern.ch
 // --hector.bello.martinez@cern.ch
