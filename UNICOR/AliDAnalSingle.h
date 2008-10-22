// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// single particle analyzer
//=============================================================================

#ifndef ALIDANALSINGLE_H
#define ALIDANALSINGLE_H

#include "AliDAnal.h"
class AliDEvent;
class AliDHN;

//=============================================================================
class AliDAnalSingle : public AliDAnal {
   
 public:
  AliDAnalSingle(Char_t *nam="single", 
	      Double_t emi=-1, Double_t ema=1, 
	      Int_t pid=0);                       // constructor
  virtual ~AliDAnalSingle(){}                        // destructor
  void Process(AliDEvent *ev);                       // fill histograms

 protected:
  Int_t    fPid;                                  // pid; 0 means all
  Double_t fMass;                                 // mass (if pid!=0)

  ClassDef(AliDAnalSingle,1)
};
//=============================================================================
#endif
