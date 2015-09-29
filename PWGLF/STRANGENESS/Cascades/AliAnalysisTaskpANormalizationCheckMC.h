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

#ifndef AliAnalysisTaskpANormalizationCheckMC_H
#define AliAnalysisTaskpANormalizationCheckMC_H

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

class AliAnalysisTaskpANormalizationCheckMC : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskpANormalizationCheckMC();
    AliAnalysisTaskpANormalizationCheckMC(const char *name);
    virtual ~AliAnalysisTaskpANormalizationCheckMC();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    Double_t MyRapidity(Double_t rE, Double_t rPz) const;
    
    //Provide Info: Is this an NSD Event?
    Bool_t IsNSD( AliMCEvent *lMCevent ) const;
    
//---------------------------------------------------------------------------------------

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;  //! List of Cascade histograms
    
    AliAnalysisUtils *fUtils; // first event in chunk check
    
    //Histograms (Desired objects in this cross-checking task)
    TH1F *fHistEventCounter; //! histogram for event counting
    
    //Basic Histograms for counting events as a function of V0M percentiles...
    TH1F *fHistV0A_DataSelection; //!
    TH1F *fHistV0A_MCSelection;   //!

    TH2F *fHistV0AVsNch_DataSelection; //!
    TH2F *fHistV0AVsNch_MCSelection;   //!
    
    //Desired Spectra: pi/K/p/K0/Lambda/Xi/Omega/Phi/K*
    //Data Selection Spectra
    //1-dimensional with transverse momentum (integrated in multiplicity)
    TH1F *fHistPt_Generated[9];     //! for keeping track of base spectra
    TH1F *fHistPt_DataSelection[9]; //!
    TH1F *fHistPt_MCSelection[9];   //!
    
    //2-dimensional with unchecked V0M percentile...
    TH2F *fHistPtVsV0A_Generated[9]; //! 9 spectra
    TH2F *fHistPtVsV0A_DataSelection[9]; //! 9 spectra
    TH2F *fHistPtVsV0A_MCSelection[9];   //! 9 spectra
    
    AliAnalysisTaskpANormalizationCheckMC(const AliAnalysisTaskpANormalizationCheckMC&);            // not implemented
    AliAnalysisTaskpANormalizationCheckMC& operator=(const AliAnalysisTaskpANormalizationCheckMC&); // not implemented

    ClassDef(AliAnalysisTaskpANormalizationCheckMC, 1);
};

#endif
