// -*- C++ -*-
// $Id$

/*  See cxx source for full Copyright notice */
#ifndef _ALI_ANALYSIS_TASK_AD_CENT_H_
#define _ALI_ANALYSIS_TASK_AD_CENT_H_
//-----------------------------------------------------------------
//            This task is for AD calibration
//-----------------------------------------------------------------

class TTree;

#include <TBits.h>
#include <TString.h>
#include <TVectorF.h>
#include "AliAODAD.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskADCent : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskADCent(const char *name="AliAnalysisTaskADCent");
  virtual ~AliAnalysisTaskADCent();

  void SetEstimatorNames(TString n="V0M:CL0:CL1:SPDClustersCorr:SPDTracklets") {
    fEstimatorNames = n;
    fMult.ResizeTo(1+n.CountChar(':'));
    fCent.ResizeTo(1+n.CountChar(':'));
  }

  virtual void NotifyRun();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

protected:

private:
  // not implemented:
  AliAnalysisTaskADCent(const AliAnalysisTaskADCent&);
  AliAnalysisTaskADCent& operator=(const AliAnalysisTaskADCent&);

  TString fEstimatorNames;

  // TTree branches
  AliAODAD         fAD;             //!
  TBits            fIR1Map;         //!
  TBits            fIR2Map;         //!
  TVectorF         fMult;           //
  TVectorF         fCent;           //

  // task output
  TTree *fTE;                    //!

  ClassDef(AliAnalysisTaskADCent, 1);
} ;
#endif // _ALI_ANALYSIS_TASK_AD_CENT_H_
