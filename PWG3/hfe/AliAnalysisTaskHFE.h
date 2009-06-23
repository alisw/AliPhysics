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
#ifndef ALIANALYSISTASKHFE_H
#define ALIANALYSISTASKHFE_H

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
#endif

#ifndef ROOT_THnSparse
#include <THnSparse.h>
#endif

class AliHFEpid;
class AliHFEcuts;
class AliHFEmcQA;
class AliHFEsecVtx;
class AliCFManager;
class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class TH1I; 
class TList;

class AliAnalysisTaskHFE : public AliAnalysisTask{
  enum{
    kIsQAOn = BIT(18),
    kIsMCQAOn = BIT(19),
    kIsSecVtxOn = BIT(20),
    kIsPriVtxOn = BIT(21)
  };
  public:
    AliAnalysisTaskHFE();
    AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref);
    AliAnalysisTaskHFE& operator=(const AliAnalysisTaskHFE &ref);
    ~AliAnalysisTaskHFE();

    virtual void ConnectInputData(Option_t *);
    virtual void CreateOutputObjects();
    virtual void Exec(Option_t *);
    virtual void Terminate(Option_t *);

    Bool_t IsQAOn() const     { return TestBit(kIsQAOn); };
    Bool_t IsMCQAOn() const   { return TestBit(kIsMCQAOn); };
    Bool_t IsSecVtxOn() const { return TestBit(kIsSecVtxOn); };
    Bool_t IsPriVtxOn() const { return TestBit(kIsPriVtxOn); };
    void SetQAOn()            { SetBit(kIsQAOn, kTRUE); };
    void SetMCQAOn()          { SetBit(kIsMCQAOn, kTRUE); };
    void SetPriVtxOn()        { SetBit(kIsPriVtxOn, kTRUE); };
    void SetSecVtxOn()        { SetBit(kIsSecVtxOn, kTRUE); };
 
  protected:
    void Copy(TObject &o) const;

  private:
    void MakeParticleContainer();

    AliESDEvent *fESD;                    //! The ESD Event
    AliMCEvent *fMC;                      //! The MC Event
    AliCFManager *fCFM;                   //! Correction Framework Manager
    THnSparseF *fCorrelation;             //! response matrix for unfolding  
    THnSparseF *fFakeElectrons;           //! Contamination from Fake Electrons
    AliHFEpid *fPID;                      //! PID
    AliHFEcuts *fCuts;                    //! Cut Collection
    AliHFEsecVtx *fSecVtx;                //! Secondary Vertex Analysis
    AliHFEmcQA *fMCQA;                    //! MC QA
    TH1I *fNEvents;                       //! counter for the number of Events
    TList *fQA;                           //! QA histos for the cuts
    TList *fOutput;                       //! Container for Task Output
    TList *fHistMCQA;                     //! Output container for MC QA histograms 
    TList *fHistSECVTX;                   //! Output container for sec. vertexing results

    ClassDef(AliAnalysisTaskHFE, 1)       // The electron Analysis Task
};
#endif

