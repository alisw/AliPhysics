/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliAnalysisTaskHistogram
///
/// This empty task is used for the analysis train to estimate the memory and CPU consumption without any user code

#ifndef ALIANALYSISTASKBASELINE_H
#define ALIANALYSISTASKBASELINE_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHistogram : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskHistogram();
    AliAnalysisTaskHistogram(const char *name);
    virtual ~AliAnalysisTaskHistogram();
    void CreateOutputObjects();
    
    virtual void     UserExec(Option_t *option);
    
    static AliAnalysisTaskHistogram* AddTask(TString suffix);
    
 private:
    AliAnalysisTaskHistogram(const AliAnalysisTaskHistogram&); // not implemented
    AliAnalysisTaskHistogram& operator=(const AliAnalysisTaskHistogram&); // not implemented
    
    TH1F* fPtHist; //!
    
    ClassDef(AliAnalysisTaskHistogram, 1); // empty analysis
};

#endif

