#ifndef ALIARRAYI_H
#define ALIARRAYI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////
//   Added additional functionality  to original TArrayI              //
//   multiple inheritance from TObject to be possible use automatic   //
//   branch mechanism for tree
//   function Expand to be possible expand array without deleting     //
//   array contents                                                  //
//                                                                   //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "TArrayI.h"

class AliArrayI: public TObject ,public TArrayI {
public:
  void Expand(Int_t n);
  ClassDef(AliArrayI,1) // Array handling
};

#endif //ALIARRAY_I

