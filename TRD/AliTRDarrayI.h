#ifndef ALITRDARRAYI_H
#define ALITRDARRAYI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

#include <TObject.h>
#include <TArrayI.h>
 
/////////////////////////////////////////////////////////////
//                                                         //
//  Array of integers                                      //
//  Origin M.Ivanov                                        //
//                                                         //
/////////////////////////////////////////////////////////////                   

class AliTRDarrayI: public TObject, public TArrayI {

public:

  AliTRDarrayI();
  virtual ~AliTRDarrayI();
  void Copy(TObject &a) const;
  void Expand(Int_t n);  

  ClassDef(AliTRDarrayI,1)  // An array of integers

};

#endif 

