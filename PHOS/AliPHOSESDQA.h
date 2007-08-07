#ifndef ALIPHOSESDQA_H
#define ALIPHOSESDQA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the PHOS photon ESD data 
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class AliESD ; 
class TH1D ; 
class TH1I ; 

class AliPHOSESDQA : public AliAnalysisTask {

public:
  AliPHOSESDQA(const char *name) ;
  virtual ~AliPHOSESDQA() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container

  // Histograms
  TNtuple * fhPHOSPos ;
  TNtuple * fhPHOS ;
  TH1D    * fhPHOSEnergy ;
  TH1I    * fhPHOSDigits ;
  TH1D    * fhPHOSRecParticles ;
  TH1I    * fhPHOSPhotons ;
  TH1D    * fhPHOSInvariantMass ;
  TH1I    * fhPHOSDigitsEvent ;
   
  ClassDef(AliPHOSESDQA, 0); // a PHOS photon analysis task 
};
#endif // ALIPHOSESDQA_H
