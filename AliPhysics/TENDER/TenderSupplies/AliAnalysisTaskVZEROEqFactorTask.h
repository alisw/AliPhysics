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

#ifndef ALIANALYSISTASKVZEROEQFACTORTASK_H
#define ALIANALYSISTASKVZEROEQFACTORTASK_H

//Needed for this task
#include "AliAnalysisTaskSE.h"

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;

class AliESDpid;
class AliESDtrackCuts;
class AliAnalysisUtils;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;
class AliVZEROCalibData; 
//class AliAnalysisTaskSE;

//#include "TString.h"
//#include "AliESDtrackCuts.h"


class AliAnalysisTaskVZEROEqFactorTask : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskVZEROEqFactorTask();
  AliAnalysisTaskVZEROEqFactorTask(const char *name);
  virtual ~AliAnalysisTaskVZEROEqFactorTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  AliVZEROCalibData* GetCalibData() const;

  void SetIsAODAnalysis(Bool_t flag) {fisAOD=flag;}
  Bool_t GetIsAODAnalysis() const {return fisAOD;}
  
 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
  TList  *fListHist;  //! List of Cascade histograms

  //Re-loaded per event if required, a la AliCentrality
  TH1F*   fEqFactors;               //! Histogram with the equalization factors used in event-plane reconstruction
  AliVZEROCalibData* fCalibData;    //! calibration data
  Long_t fRunNumber;                //! Needed to make sure we haven't swapped runs
 
//===========================================================================================
//   Histograms - Event Counting only 
//===========================================================================================

  TH1D *fHistEventCounter; //! 
  Bool_t fisAOD; // flag for analysis on AODs

   AliAnalysisTaskVZEROEqFactorTask(const AliAnalysisTaskVZEROEqFactorTask&);            // not implemented
   AliAnalysisTaskVZEROEqFactorTask& operator=(const AliAnalysisTaskVZEROEqFactorTask&); // not implemented
   
   ClassDef(AliAnalysisTaskVZEROEqFactorTask, 12);
};

#endif
