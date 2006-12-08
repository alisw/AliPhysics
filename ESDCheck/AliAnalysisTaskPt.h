#ifndef ALIANALYSISTASKPT_H
#define ALIANALYSISTASKPT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// A test analysis task to check the pt of tracks distribution in simulated data
//
//*-- Panos
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"

class AliESD ; 
class TH1 ;

class AliAnalysisTaskPt : public AliAnalysisTask {

public:
  AliAnalysisTaskPt(const char *name);
  virtual ~AliAnalysisTaskPt() {}
  
  virtual void   Init(Option_t * opt = "");
  virtual void   Exec(Option_t * opt = "");
  virtual void   Terminate(Option_t * opt = "");
  
private:
  TTree  * fChain ; //! pointer to the analyzed TTree or TChain
  AliESD * fESD;    //! ESD object
  TH1F   * fhPt;    //! Pt spectrum
  
  TObjArray * fOutputContainer ; //! output data container

  ClassDef(AliAnalysisTaskPt, 0); // example of analysis
};
#endif // ALIANALYSISTASKPT_H

