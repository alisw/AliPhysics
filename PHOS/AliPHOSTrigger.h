#ifndef ALIPHOSTrigger_H
#define ALIPHOSTrigger_H

//____________________________________________________________
//  Class for trigger analysis.
//  Digits are grouped in TRU's (16x28 crystals). The algorithm searches 
//  all possible 4x4 crystal combinations and per each TRU, adding the 
//  digits amplitude and finding the maximum. Maximums are compared to 
//  triggers threshold and they are set.

//*-- Author: Gustavo Conesa & Yves Schutz (IFIC, SUBATECH, CERN)
     
// --- ROOT system ---

class TMatrixD ;
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
  
  const Int_t  GetNTRU() const              {return fNTRU ; }  
  const Int_t  GetNTRUZ() const             {return fNTRUZ ; }  
  const Int_t  GetNTRUPhi() const           {return fNTRUPhi ; }  
  const Int_t  GetL0Threshold() const       {return fL0Threshold ; }  
  const Int_t  GetL1LowThreshold() const    {return fL1LowThreshold ; }
  const Int_t  GetL1MediumThreshold() const {return fL1MediumThreshold ; }
  const Int_t  GetL1HighThreshold() const   {return fL1HighThreshold ; }

  void         Print(const Option_t * opt ="") const ;  

  void         SetNTRU(Int_t ntru)            {fNTRU              = ntru; }
  void         SetNTRUZ(Int_t ntru)           {fNTRUZ             = ntru; }
  void         SetNTRUPhi(Int_t ntru)         {fNTRUPhi           = ntru; }
  void         SetL0Threshold(Int_t amp)      {fL0Threshold       = amp; }
  void         SetL1LowThreshold(Int_t amp)   {fL1LowThreshold    = amp; } 
  void         SetL1MediumThreshold(Int_t amp){fL1MediumThreshold = amp; } 
  void         SetL1HighThreshold(Int_t amp)  {fL1HighThreshold   = amp; }

 private:
  TClonesArray *  FillTRU(const TClonesArray * digits, 
			  const AliPHOSGeometry * geom, const Int_t nModules, 
			  const Int_t nCrystalsPhi, const Int_t nCrystalsZ) const ;
  void MakeSlidingCell(const TClonesArray * trus, const Int_t mod,
		       const Int_t nCrystalsPhi, const Int_t nCrystalsZ, 
		       Float_t *ampmax) ;
  void SetTriggers(const Float_t * ampmax) ;


 private: 

  Int_t    fNTRU ;              //! Number of TRUs per module
  Int_t    fNTRUZ ;             //! Number of crystal rows per Z in one TRU
  Int_t    fNTRUPhi ;           //! Number of crystal rows per Phi in one TRU
  Int_t    fL0Threshold ;       //! L0 trigger energy threshold
  Int_t    fL1LowThreshold ;    //! High pT trigger energy threshold
  Int_t    fL1MediumThreshold ; //! 
  Int_t    fL1HighThreshold ;   //! 

  ClassDef(AliPHOSTrigger,2)
} ;


#endif //ALIPHOSTrigger_H
