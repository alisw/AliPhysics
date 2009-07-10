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

//-----------------------------------------------------------------------
// Author : R. Vernet, INFN - Catania (it)
//-----------------------------------------------------------------------

#ifndef ALICFTASKFORUNFOLDING_H
#define ALICFTASKFORUNFOLDING_H

#include "AliAnalysisTaskSE.h"

class TH1I;
class TParticle ;
class TFile ;
class AliStack ;
class AliCFManager;
class AliESDtrack;
class AliVParticle;
class THnSparse;

class AliCFTaskForUnfolding : public AliAnalysisTaskSE {
  public:

  AliCFTaskForUnfolding();
  AliCFTaskForUnfolding(const Char_t* name);
  virtual ~AliCFTaskForUnfolding();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager() const {return fCFManager;}           // get corr manager
  void           SetCorrelationMatrix(THnSparse* h) {fCorrelation=h;}

 protected:
  
  AliCFManager   *fCFManager    ;  // pointer to the CF manager

  // Histograms
  //Number of events
  TH1I       *fHistEventsProcessed; // simple histo for monitoring the number of events processed
  THnSparse  *fCorrelation;         //response matrix for unfolding  
  ClassDef(AliCFTaskForUnfolding,0);
};

#endif
