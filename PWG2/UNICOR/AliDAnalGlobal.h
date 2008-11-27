// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// event global variable analyzer
//=============================================================================

#ifndef ALIDANALGLOBAL_H
#define ALIDANALGLOBAL_H

#include "AliDAnal.h"
class AliDEvent;

//=============================================================================
class AliDAnalGlobal : public AliDAnal {
   
 public:
  AliDAnalGlobal(Char_t *nam="global"); // constructor
  virtual ~AliDAnalGlobal(){}           // destructor
  void Process(AliDEvent *ev);          // fill histograms

  ClassDef(AliDAnalGlobal,1)
};
//=============================================================================
#endif
