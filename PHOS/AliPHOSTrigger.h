#ifndef ALIPHOSTrigger_H
#define ALIPHOSTrigger_H

//____________________________________________________________
//  Class for trigger analysis.
//  Digits are grouped in TRU's (16x28 crystals). The algorithm searches 
//  all possible 4x4 crystal combinations and per each TRU, adding the 
//  digits amplitude and finding the maximum. Maximums are compared to 
//  triggers threshold and they are set.
//  FIRST ATTEMPT TO MAKE A TRIGGER CLASS. IT WILL CHANGE WHEN CENTRAL TRIGGER CLASS FIXES 
//  THINGS
//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
     
// --- ROOT system ---
#include "TTask.h"
#include "TClonesArray.h"
#include "TMatrixD.h"

class AliPHOSGeometry ;

// --- AliRoot header files ---

class AliPHOSTrigger : public TObject {
  
 public:   
  AliPHOSTrigger() ; //  ctor
  AliPHOSTrigger(const AliPHOSTrigger & trig) ; // cpy ctor
  virtual ~AliPHOSTrigger() {}; //virtual dtor
  
  void         MakeTrigger() ; //Make PHOS trigger

  const Bool_t IsL0Set() const       {return fL0 ;}  // Is L0 trigger set?
  const Bool_t IsL1LowSet() const    {return fL1Low ;} // Is L1 trigger set?
  const Bool_t IsL1MediumSet() const {return fL1Medium ;} 
  const Bool_t IsL1HighSet() const   {return fL1High ;} 

  const Int_t  GetL0Threshold() const       {return fL0Threshold ; }  
  const Int_t  GetL1LowThreshold() const    {return fL1LowThreshold ; }
  const Int_t  GetL1MediumThreshold() const {return fL1MediumThreshold ; }
  const Int_t  GetL1HighThreshold() const   {return fL1HighThreshold ; }

  void         Print(const Option_t * opt ="") const ;  

  void         SetL0Threshold(Int_t amp)      {fL0Threshold       = amp; }
  void         SetL1LowThreshold(Int_t amp)   {fL1LowThreshold    = amp; } 
  void         SetL1MediumThreshold(Int_t amp){fL1MediumThreshold = amp; } 
  void         SetL1HighThreshold(Int_t amp)  {fL1HighThreshold   = amp; }

 private:
  TClonesArray *  FillTRU(const TClonesArray * digits) ;
  void MakeSlidingCell(const TClonesArray * trus, const Int_t mod, 
		       Float_t *ampmax) ;
  void SetTriggers(const Float_t * ampmax) ;
  void InitTriggers() ;

  void SetL0()       { fL0       = kTRUE ; } 
  void SetL1Low()    { fL1Low    = kTRUE ; }  
  void SetL1Medium() { fL1Medium = kTRUE ; }  
  void SetL1High()   { fL1High   = kTRUE ; } 


 private: 
  
  Bool_t    fL0 ;       //! Minimum Bias Trigger
  Bool_t    fL1Low ;    //! High pT triggers
  Bool_t    fL1Medium ; //!
  Bool_t    fL1High ;   //!

  Int_t    fL0Threshold ;       //! L0 trigger energy threshold
  Int_t    fL1LowThreshold ;    //! High pT trigger energy threshold
  Int_t    fL1MediumThreshold ; //! 
  Int_t    fL1HighThreshold ;   //! 

  ClassDef(AliPHOSTrigger,0)
} ;


#endif //ALIPHOSTrigger_H
