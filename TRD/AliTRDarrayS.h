#ifndef ALITRDARRAYS_H
#define ALITRDARRAYS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDarrayS.h,v Exp $*/

#include <TObject.h>
#include <TArrayS.h>
 
/////////////////////////////////////////////////////////////
//                                                         //
//  Array of shorts                                        //
//  Origin M.Ivanov                                        //
//                                                         //
/////////////////////////////////////////////////////////////                   

class AliTRDarrayS: public TObject, public TArrayS {

public:

  AliTRDarrayS();
  virtual ~AliTRDarrayS();
  void Copy(TObject &a) const;
  void Expand(Int_t n);  

  ClassDef(AliTRDarrayS,1)  // An array of shorts

};

#endif 

