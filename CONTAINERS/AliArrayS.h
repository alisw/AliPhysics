#ifndef ALIARRAYS_H
#define ALIARRAYS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////
//   Added additional functionality  to original TArrayS              //
//   multiple inheritance from TObject to be possible use automatic   //
//   branch mechanism for tree
//   function Expand to be possible expand array without deleting     //
//   array contents                                                  //
//                                                                   //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////
#include "TObject.h"
#include "TArrayS.h"


class AliArrayS:  public TObject,public TArrayS {
public:
  void Expand(Int_t n);
  ClassDef(AliArrayS,1) // Array handling
};
#endif //ALIARRAYS_H
