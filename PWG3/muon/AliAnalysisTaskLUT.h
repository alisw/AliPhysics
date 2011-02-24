#ifndef ALIANALYSISTASKLUT_H
#define ALIANALYSISTASKLUT_H

/* $Id$ */ 

//===================================================================
//  Class AliAnalysisTaskLUT
//
//  Extract ESD muon tracks information and store in ntuple.
//  Used for offline calculation/adjustment of Look-up-Tables for the
//  trigger chambers.
//
//  Clermont-Ferrand 2008
//===================================================================

#include "AliAnalysisTask.h"

class TNtuple;
class AliESDEvent;

class AliAnalysisTaskLUT : public AliAnalysisTask {
 public:
  AliAnalysisTaskLUT(const char *name);
  virtual ~AliAnalysisTaskLUT() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:

  AliAnalysisTaskLUT(const AliAnalysisTaskLUT& atlut);
  AliAnalysisTaskLUT& operator=(const AliAnalysisTaskLUT& atlut);

  AliESDEvent  *fESDEvent;    // ESD main tree object
  TNtuple      *fTracksLUT;   // ntuple object, ESD tracks for LUT calculation

  ClassDef(AliAnalysisTaskLUT, 0); // analysis task for extracting tracks used for LUT calculation
};

#endif
