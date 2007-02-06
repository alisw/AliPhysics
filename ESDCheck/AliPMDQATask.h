#ifndef ALIPMDQATASK_H
#define ALIPMDQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the PMD  data in simulated data
//
//*-- Yves Schutz 
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESD ; 
class TH2F ; 
class TH1F ; 

class AliPMDQATask : public AliAnalysisTask {

public:
  AliPMDQATask(const char *name) ;
  virtual ~AliPMDQATask() ;
   
  virtual void Exec(Option_t * opt = "") ;
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "") ;

private:
  void CalculateSMN( Float_t clsX, Float_t clsY, Int_t & smn) ; 
  void DrawPMDBoundary() ;
  void DrawPMDBoundarySM1() ;
  void DrawPMDBoundarySM2() ;
  void DrawPMDBoundarySM3() ;
  void DrawPMDBoundarySM4() ;

private:
  TTree   * fChain ;            //!pointer to the analyzed TTree or TChain
  AliESD  * fESD ;              //! Declaration of leave types

  TObjArray * fOutputContainer ; //! output data container


 // Histograms
  TH2F * fhPMDP1  ; //!
  TH1F * fhPMDC2  ; //!
  TH1F * fhPMDP2  ; //!
  TH1F * fhPMDC3  ; //!
  TH1F * fhPMDP3  ; //!
  TH1F * fhPMDP4  ; //!
  TH1F * fhPMDC5  ; //!
  TH1F * fhPMDP5  ; //!
  TH2F * fhPMDCP0 ; //!
  TH2F * fhPMDCP1 ; //!
  TH2F * fhPMDCP2 ; //!
  TH2F * fhPMDCP3 ; //!
  TH2F * fhPMDCP4 ; //!
  
  TH2F * fhPMDSM1  ; //!
  TH2F * fhPMDSM2  ; //!
  TH2F * fhPMDSM3  ; //!
  TH2F * fhPMDSM4  ; //!
  TH2F * fhPMDSM5  ; //!
  TH2F * fhPMDSM6  ; //!
  TH2F * fhPMDSM7  ; //!
  TH2F * fhPMDSM8  ; //!
  TH2F * fhPMDSM9  ; //!
  TH2F * fhPMDSM10 ; //!
  TH2F * fhPMDSM11 ; //!
  TH2F * fhPMDSM12 ; //!
  TH2F * fhPMDSM13 ; //!
  TH2F * fhPMDSM14 ; //!
  TH2F * fhPMDSM15 ; //!
  TH2F * fhPMDSM16 ; //!
  TH2F * fhPMDSM17 ; //!
  TH2F * fhPMDSM18 ; //!
  TH2F * fhPMDSM19 ; //!
  TH2F * fhPMDSM20 ; //!
  TH2F * fhPMDSM21 ; //!
  TH2F * fhPMDSM22 ; //!
  TH2F * fhPMDSM23 ; //!
  TH2F * fhPMDSM24 ; //!
  TH1F * fhPMDSM   ; //!
   
  ClassDef(AliPMDQATask, 0); //! // a PMD analysis task 
}; 
#endif // ALIPMDQATASK_H
