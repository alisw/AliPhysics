/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliAnalysisTaskBaseLine.cxx 46301 2011-01-06 14:25:27Z agheata $ */

//
// This empty task is used for the analysis train to estimate the memory and CPU consumption without any user code
//

#ifndef ALIANALYSISTASKBASELINE_H
#define ALIANALYSISTASKBASELINE_H

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskBaseLine : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskBaseLine();
    AliAnalysisTaskBaseLine(const char *name);
    virtual ~AliAnalysisTaskBaseLine();
    
    virtual void     UserExec(Option_t *option);
    
 private:
    AliAnalysisTaskBaseLine(const AliAnalysisTaskBaseLine&); // not implemented
    AliAnalysisTaskBaseLine& operator=(const AliAnalysisTaskBaseLine&); // not implemented
    
    ClassDef(AliAnalysisTaskBaseLine, 1); // empty analysis
};

#endif

