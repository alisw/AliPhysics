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
//
// Task for Efficiency studies
// Used for testing classes AliHFEcontainer and AliHFEfilter
// Creates Efficiency Histograms
//
#ifndef ALIHFEEFFICIENCY_H
#define ALIHFEEFFICIENCY_H

#include "AliAnalysisTaskSE.h"

class TList;
class TObjArray;
class AliCFAcceptanceCuts;
class AliCFEffGrid;
class AliCFContainer;
class AliHFEcollection;
class AliHFEcontainer;
class AliHFEtrackFilter;

class AliHFEefficiency : public AliAnalysisTaskSE{
  public:
    AliHFEefficiency();
    AliHFEefficiency(const Char_t *name);
    ~AliHFEefficiency();

    virtual void UserCreateOutputObjects(); 
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);

    void FilterMC();
    void CutTRD() { fCutTRD = kTRUE; }

    void SetHFEcontainer(AliHFEcontainer * const cont) { fEfficiency = cont; };
    void UnsetHFEcontainer() { fEfficiency = NULL; }

    void Load(const char* filename = "EffTask.root");
    void PostProcess();
    Bool_t IsRunningTerminate() const { return TestBit(kTerminate); }
    void SetRunTerminate(Bool_t terminate = kTRUE) { SetBit(kTerminate, terminate); }

    void CalculatePTsmearing();
    void DrawPtResolution(const TList * const l);
    void DrawSignalEfficiency(AliCFEffGrid *eff, AliCFContainer *cont, Int_t var);
    void DrawCutEfficiency(AliCFEffGrid *eff, AliCFContainer *cont, Int_t var);

  private:
    enum{
      kNTracks,
      kPt
    };
    enum{ // Bit Definition
      kTerminate = BIT(18)
    };
    AliHFEefficiency(const AliHFEefficiency &);
    AliHFEefficiency &operator=(const AliHFEefficiency &);


    AliHFEtrackFilter *fFilter;           //! Track Filter
    AliHFEcutStep *fMCcut;                //! MC Signal
    AliCFAcceptanceCuts *fAcceptanceCuts; //! MC Acceptance cuts
    AliHFEcontainer *fEfficiency;         //! Efficiency container
    AliHFEcollection *fOutput;            //! QA histo container
    Bool_t fCutTRD;                       //  Apply TRD cuts

    ClassDef(AliHFEefficiency, 1);
};
#endif
