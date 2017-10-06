/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: AliHLTObjArray.h 64008 2016-08-25 13:09:59Z mkrzewic $ */

//-------------------------------------------------------------------------
//                          Class AliHLTObjArray
//essentially a TObjArray, becomes owner after deserializing
// Origin: Mikolaj Krzewicki, mkrzewic@cern.ch
//-------------------------------------------------------------------------

#ifndef ALIHLTOBJARRAY_H
#define ALIHLTOBJARRAY_H
#include "TObjArray.h"
class AliHLTObjArray : public TObjArray {
public:
AliHLTObjArray(Int_t s = TCollection::kInitCapacity, Int_t lowerBound = 0):
   TObjArray(s, lowerBound) {}
AliHLTObjArray(const AliHLTObjArray &a):TObjArray(a) {}
virtual ~AliHLTObjArray() {}
ClassDef(AliHLTObjArray,1)
};
#endif

