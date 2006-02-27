#ifndef ALIPHOSTrigger_H
#define ALIPHOSTrigger_H

//____________________________________________________________
//  Class for trigger analysis.
//  Class for trigger analysis.
//  Digits are grouped in TRU's (16x28 ordered fNTRUPhi x fNTRUEta). 
//  The algorithm searches all possible 4x4 cell combinations per each TRU, 
//  adding the digits amplitude and finding the maximum. Maximums are compared 
//  to triggers threshold and they are set. Thresholds need to be fixed. 
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
//  tr->Print(""); //Print result, with "deb" option all data members 
//  //are printed
//
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
     
// --- ROOT system ---

class TClonesArray ;

// --- AliRoot header files ---
#include "AliTriggerDetector.h"

class AliPHOSGeometry ;

class AliPHOSTrigger : public AliTriggerDetector {
  
 public:   
  AliPHOSTrigger() ; //  ctor
  AliPHOSTrigger(const AliPHOSTrigger & trig) ; // cpy ctor
  virtual ~AliPHOSTrigger() {}; //virtual dtor
  virtual void    CreateInputs();
  virtual void    Trigger();  //Make PHOS trigger
  
  Int_t  GetNTRU() const                 {return fNTRU ; }  
  Int_t  GetNTRUZ() const                {return fNTRUZ ; }  
  Int_t  GetNTRUPhi() const              {return fNTRUPhi ; }  
  Int_t  GetL0MBPbPbThreshold() const    {return fL0MBPbPbThreshold ; } 
  Int_t  GetL0MBppThreshold() const      {return fL0MBppThreshold ; } 
  Int_t  GetL1JetLowPtThreshold() const  {return fL1JetLowPtThreshold ; }
  Int_t  GetL1JetHighPtThreshold() const {return fL1JetHighPtThreshold ; }

  void         Print(const Option_t * opt ="") const ;  

  void         SetNTRU(Int_t ntru)               {fNTRU              = ntru; }
  void         SetNTRUZ(Int_t ntru)              {fNTRUZ             = ntru; }
  void         SetNTRUPhi(Int_t ntru)            {fNTRUPhi           = ntru; }
  void         SetL0MBPbPbThreshold(Int_t amp)   {fL0MBPbPbThreshold   = amp; }
  void         SetL0MBppThreshold(Int_t amp)     {fL0MBppThreshold     = amp; }
  void         SetL1JetLowPtThreshold(Int_t amp) 
    {fL1JetLowPtThreshold  = amp; } 
  void         SetL1JetHighPtThreshold(Int_t amp)
    {fL1JetHighPtThreshold = amp; }

 private:
  TClonesArray *  FillTRU(const TClonesArray * digits, 
			  const AliPHOSGeometry * geom, const Int_t nModules, 
			  const Int_t nCrystalsPhi, const Int_t nCrystalsZ) const ;
  void MakeSlidingCell(const TClonesArray * trus, const Int_t mod,
		       const Int_t nCrystalsPhi, const Int_t nCrystalsZ, 
		       Float_t *ampmax) ;
  void SetTriggers(const Float_t * ampmax) ;


 private: 

  Int_t    fNTRU ;                 //! Number of TRUs per module
  Int_t    fNTRUZ ;                //! Number of crystal rows per Z in one TRU
  Int_t    fNTRUPhi ;              //! Number of crystal rows per Phi in one TRU
  Int_t    fL0MBPbPbThreshold ;    //! L0 PbPb trigger energy threshold
  Int_t    fL0MBppThreshold ;      //! L0 pp trigger energy threshold
  Int_t    fL1JetLowPtThreshold ;  //! Low and High pT trigger energy threshold
  Int_t    fL1JetHighPtThreshold ; //! 

  ClassDef(AliPHOSTrigger,3)
} ;


#endif //ALIPHOSTrigger_H
