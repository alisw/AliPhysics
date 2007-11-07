#ifndef ALIEMCALQATASK_H
#define ALIEMCALQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the EMCAL photon data in simulated data
// An analysis task to check the EMCAL photon data in simulated data
// An analysis task to check the EMCAL photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class AliESD ; 
class TNtuple ;
class TH1D ; 
class TH1I ; 
class TTree ; 

class AliEMCALQATask : public AliAnalysisTask {

public:
  AliEMCALQATask(const char *name) ;
  AliEMCALQATask(const AliEMCALQATask& ap) ;   
  AliEMCALQATask& operator = (const AliEMCALQATask& ap) ;
  virtual ~AliEMCALQATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  // Histograms
  TNtuple * fhEMCALPos ; // x,y 
  TNtuple * fhEMCAL ; // all parameters
  TH1D    * fhEMCALEnergy ; // energy
  TH1I    * fhEMCALDigits ; // digits
  TH1D    * fhEMCALRecParticles ; //recparticles
  TH1I    * fhEMCALPhotons ; // photons
  TH1D    * fhEMCALInvariantMass ; // invariant mass
  TH1I    * fhEMCALDigitsEvent ; // digits per event
   
  ClassDef(AliEMCALQATask, 0); // a EMCAL photon analysis task
};
#endif // ALIEMCALQATASK_H
