#ifndef ALIHFEPOSTANALYSIS_H
#define ALIHFEPOSTANALYSIS_H

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

/* $Id$ */ 

//
// Post analysis class
// Creating results and draw pictures
// Called in AliAnalysisTaskHFE::Terminate or in a macro
//
#ifndef ROOT_THnSparse
#include <THnSparse.h>
#endif

class AliHFEcontainer;
class TH1;
class TList;

class AliHFEpostAnalysis : public TObject{
  public:
    AliHFEpostAnalysis();
    AliHFEpostAnalysis(const AliHFEpostAnalysis &ref);
    AliHFEpostAnalysis &operator=(const AliHFEpostAnalysis &ref);
    ~AliHFEpostAnalysis();

    Int_t SetTaskResults(AliHFEcontainer *trackContainer) { fEfficiencyContainer = trackContainer; return 1; };
    Int_t SetTaskQA(const TList *qa);
    void StoreOutput(const char *filename = "HFEresults.root");

    void DrawMCSignal2Background();
    void DrawEfficiency();
    void DrawPIDperformance();
    void DrawCutEfficiency(Bool_t MC = kTRUE, Int_t source = -1);
  private:
    enum{
      kCFC,
      kPIDperf,
      kSigBackg
    };
    TH1 *CreateHistoSignalToBackgroundMC(Int_t mode, Int_t charge);
    TH1 *CreateHistoPIDperformance(Int_t mode, Int_t charge);

    TList *fResults;                          // Container for output objects
    UChar_t fAnalysisObjects;                       // S
    AliHFEcontainer *fEfficiencyContainer;     // Task Results
    THnSparseF *fPIDperformance;              // PID Performance Studies
    THnSparseF *fSignalToBackgroundMC;        // Signal To Background Studies

    ClassDef(AliHFEpostAnalysis, 1)           // Result Creator class
};

#endif
