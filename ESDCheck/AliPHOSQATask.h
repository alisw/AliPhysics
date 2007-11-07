#ifndef ALIPHOSQATASK_H
#define ALIPHOSQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the PHOS photon data in simulated data
// An analysis task to check the PHOS photon data in simulated data
// An analysis task to check the PHOS photon data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class AliESD ; 
class TNtuple ;
class TH1D ; 
class TH1I ;
class TTree ;  

class AliPHOSQATask : public AliAnalysisTask {

public:
  AliPHOSQATask(const char *name) ;
  AliPHOSQATask(const AliPHOSQATask& ap) ;   
  AliPHOSQATask& operator = (const AliPHOSQATask& ap) ;
 virtual ~AliPHOSQATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  // Histograms
  TNtuple * fhPHOSPos ; // x, y
  TNtuple * fhPHOS ; // all parameters
  TH1D    * fhPHOSEnergy ; // energy
  TH1I    * fhPHOSDigits ; // sdigits
  TH1D    * fhPHOSRecParticles ; // recparticles
  TH1I    * fhPHOSPhotons ; // photons 
  TH1D    * fhPHOSInvariantMass ; // invariant mass
  TH1I    * fhPHOSDigitsEvent ; // digits per event
   
  ClassDef(AliPHOSQATask, 0); // a PHOS photon analysis task 
};
#endif // ALIPHOSQATASK_H
