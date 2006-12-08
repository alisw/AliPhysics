#ifndef ALIEMCALQATASK_H
#define ALIEMCALQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the EMCAL photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESD ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 

class AliEMCALQATask : public AliAnalysisTask {

public:
  AliEMCALQATask(const char *name) ;
  virtual ~AliEMCALQATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void Init(Option_t * opt = "") ; 
  virtual void Terminate(Option_t * opt = "") ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  // Histograms
  TNtuple * fhEMCALPos ;
  TNtuple * fhEMCAL ;
  TH1D    * fhEMCALEnergy ;
  TH1I    * fhEMCALDigits ;
  TH1D    * fhEMCALRecParticles ;
  TH1I    * fhEMCALPhotons ;
  TH1D    * fhEMCALInvariantMass ;
  TH1I    * fhEMCALDigitsEvent ;
   
  ClassDef(AliEMCALQATask, 0); // a EMCAL photon analysis task
};
#endif // ALIEMCALQATASK_H
