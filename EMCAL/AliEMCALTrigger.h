#ifndef ALIEMCALTrigger_H
#define ALIEMCALTrigger_H

//___________________________________________________________
//  Class for trigger analysis.
//  Digits are grouped in TRU's (384 cells? ordered fNTRUPhi x fNTRUEta). 
//  The algorithm searches all possible 4x4 cell combinations per each TRU, 
//  adding the digits amplitude and finding the maximum. Maximums are compared 
//  to triggers threshold and they are set. Thresholds need to be fixed. 
//  Last 2 modules are half size but they are treated as fullsize, then their 
//  TRU should be smaller. When this is fixed, class will be updated. 
//  Usage:
//
//  //Inside the event loop
//  AliEMCALTrigger *tr = new AliEMCALTrigger();//Init Trigger
//  tr->SetL0MBPbPbThreshold(500);
//  tr->SetL0MBppThreshold(100);
//  tr->SetL1JetLowPtThreshold(1000);
//  tr->SetL1JetMediumPtThreshold(10000);
//  tr->SetL1JetHighPtThreshold(20000);
//  tr->Trigger(); //Execute Trigger
//  tr->Print(""); //Print results
//  //are printed
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
     
// --- ROOT system ---

class TClonesArray ;


// --- AliRoot header files ---
#include "AliTriggerDetector.h"

class AliEMCALGeometry ;

class AliEMCALTrigger : public AliTriggerDetector {
  
 public:   
  AliEMCALTrigger() ; //  ctor
  AliEMCALTrigger(const AliEMCALTrigger & trig) ; // cpy ctor
  virtual ~AliEMCALTrigger() {}; //virtual dtor
  virtual void    CreateInputs();
  virtual void    Trigger();  //Make EMCAL trigger
  
  Int_t  GetL0MBPbPbThreshold() const    {return fL0MBPbPbThreshold ; } 
  Int_t  GetL0MBppThreshold() const      {return fL0MBppThreshold ; } 
  Int_t  GetL1JetLowPtThreshold() const  {return fL1JetLowPtThreshold ; }
  Int_t  GetL1JetMediumPtThreshold() const {return fL1JetMediumPtThreshold ; }
  Int_t  GetL1JetHighPtThreshold() const {return fL1JetHighPtThreshold ; }

  void         Print(const Option_t * opt ="") const ;  

  void         SetL0MBPbPbThreshold(Int_t amp)   
    {fL0MBPbPbThreshold    = amp; }
  void         SetL0MBppThreshold(Int_t amp)     
    {fL0MBppThreshold      = amp; }
  void         SetL1JetLowPtThreshold(Int_t amp) 
    {fL1JetLowPtThreshold  = amp; } 
  void         SetL1JetMediumPtThreshold(Int_t amp) 
    {fL1JetMediumPtThreshold  = amp; } 
  void         SetL1JetHighPtThreshold(Int_t amp)
    {fL1JetHighPtThreshold = amp; }
 private:
 
  void MakeSlidingCell(const TClonesArray * trus, const Int_t nTRU, 
		       const Int_t supermod, const Int_t nCellsPhi, 
		       const Int_t nCellsEta, Float_t *ampmax) ; 
  

  void SetTriggers(const Float_t * ampmax, const Int_t nTRU) ;
  
  
 private: 
  
  Int_t    fL0MBPbPbThreshold ;      //! L0 PbPb trigger energy threshold
  Int_t    fL0MBppThreshold ;        //! L0 pp trigger energy threshold
  Int_t    fL1JetLowPtThreshold ;    //! Low pT trigger energy threshold
  Int_t    fL1JetMediumPtThreshold ; //! Medium pT trigger energy threshold
  Int_t    fL1JetHighPtThreshold ;   //! High pT trigger energy threshold

  ClassDef(AliEMCALTrigger,0)
} ;


#endif //ALIEMCALTrigger_H
