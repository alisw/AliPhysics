#ifndef ALIOBJECTARRAY_H
#define ALIOBJECTARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/////////////////////////////////////////////////////////////////////////
//  AliObjectArray                                                     //
//AliObjectArray is an array of clone (identical) objects.             //     
//In comparison with the TClonesArray objects in this array don't need //
//to derive from TObject. They also don't need RTTI - type information.// 
//Objects type information is stored in object fClassInfo (instance of //
//the AliClassInfo).                                                   //
//                                                                     //
/////////////////////////////////////////////////////////////////////////

#include "AliClassInfo.h"
#include "AliMemArray.h"
class TClass;

class AliObjectArray: public AliMemArray {
public:
  AliObjectArray();
  AliObjectArray(const char * classname, Int_t buffersize=0);
  AliObjectArray(const AliObjectArray &arr); //copy constructor 
  AliObjectArray & operator = (const AliObjectArray &arr);
  ~AliObjectArray();
  Bool_t SetClass(const char * classname);
  //
  TClass * GetClass() {return fClassInfo->GetClass();}
  AliClassInfo  * GetClassInfo() const {return fClassInfo;} 
  virtual void     Dump(Int_t i);
  // 
protected :
 void  CTORBuffer(void * buffer, UInt_t size)
    {fClassInfo->CTORBuffer(buffer,size);} // buffer constructor   
 void  DTORBuffer(void * buffer, UInt_t size)
    {fClassInfo->DTORBuffer(buffer,size);} // buffer constructor
private:     
  AliClassInfo      *fClassInfo;        //pointer to containg class info  
  // 
  ClassDef(AliObjectArray,0) 
};







#endif //ALIMEMARRAY_I

