#ifndef ALIANALYSISTASKHFEPIDQA_H
#define ALIANALYSISTASKHFEPIDQA_H

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
// Task for PID QA
// Using AliHFEpidQA and AliHFEMCpidQA
// More information can be found in the source file
//
#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class TH1;
class TList;
class TFile;

class AliLog;
class AliMCEvent;

class AliHFEpidQA;

class AliAnalysisTaskHFEpidQA : public AliAnalysisTaskSE{
  public:
    AliAnalysisTaskHFEpidQA();
    AliAnalysisTaskHFEpidQA(const Char_t *name);
    ~AliAnalysisTaskHFEpidQA();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);
    virtual Bool_t UserNotify();

    Bool_t HasV0pidQA() const { return TestBit(kV0pidQA); };
    Bool_t HasRecalculateTRDpid() const { return TestBit(kRecalculateTRDpid); };
    void SetV0pidQA(Bool_t v0pidQA = kTRUE) { SetBit(kV0pidQA, v0pidQA); };
    void SetRecalculateTRDpid(Bool_t recal = kTRUE) { SetBit(kRecalculateTRDpid, recal); };

    void SetNNref(TFile *f) { fNNref = f; };

  private:
    enum{
      kV0pidQA = BIT(22),
      kRecalculateTRDpid = BIT(23)
    };
    AliAnalysisTaskHFEpidQA(const AliAnalysisTaskHFEpidQA &ref);
    AliAnalysisTaskHFEpidQA &operator=(const AliAnalysisTaskHFEpidQA &ref);
    AliHFEpidQA *fPIDqa;    //! The heart of the analysis  
    TList *fOutput;         //! Container for output histos
    TH1 *fEvents;           //! Number of Events
    TFile  *fNNref;         //  reference file for NN

    ClassDef(AliAnalysisTaskHFEpidQA, 1)
};

#endif

