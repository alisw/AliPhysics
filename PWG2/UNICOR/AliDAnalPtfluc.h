// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2005

//=============================================================================
// pt-fluctuations analyzer
//=============================================================================

#ifndef ALIDANALPTFLUC_H
#define ALIDANALPTFLUC_H

#include "AliDAnal.h"
class AliDEvent;
class AliDHN;

//=============================================================================
class AliDAnalPtfluc : public AliDAnal {
   
 public:
  AliDAnalPtfluc(Char_t *nam="correl", Int_t pid0=0, Int_t pid1=0);  // constructor
  virtual ~AliDAnalPtfluc(){}                                        // destructor
  void Process(Int_t tmr, AliDEvent *ev0, AliDEvent *ev1);              // process event(s)

 protected:
  Int_t    fPid0;                       // particle species 0
  Int_t    fPid1;                       // particle species 1
  ClassDef(AliDAnalPtfluc,1)
};
//=============================================================================
#endif
