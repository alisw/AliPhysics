#ifndef ALITRDARRAYF_H
#define ALITRDARRAYF_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$*/

#include <TObject.h>
#include <TArrayF.h>
 
/////////////////////////////////////////////////////////////
//                                                         //
//  Array of floats                                        //
//  Origin M.Ivanov                                        //
//                                                         //
/////////////////////////////////////////////////////////////                   

class AliTRDarrayF: public TObject, public TArrayF {

public:

  AliTRDarrayF();
  virtual ~AliTRDarrayF();
  void Copy(TObject &a) const;
  void Expand(Int_t n);  

  ClassDef(AliTRDarrayF,1)  // An array of floats

};

#endif 

