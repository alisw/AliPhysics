#ifndef ALIDATATYPE_H
#define ALIDATATYPE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  AliDataType                                                              //
//                                                                           //
// Defined to make unified interface to primitive types (in Root described by//
// TDataType) and TClass.                                                    //
//                                                                           //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk         // 
///////////////////////////////////////////////////////////////////////////////

#include "AliClassInfo.h"


class AliDataType : public AliClassInfo {
public:
  AliDataType(const char *name);
  const char * GetClassName(); 
  void StreamBuffer(TBuffer& b, const void *object, UInt_t size);
  void ObjectDump(void *p);
  TDataType * GetDataType(){return fDataType;}  
protected:
  AliDataType(const AliDataType & type){;}
  AliDataType &operator = (const AliDataType & type){return *this;} //assignment operator
  TDataType * fDataType;  //root type information
  ClassDefT(AliDataType,0)
};

#endif

