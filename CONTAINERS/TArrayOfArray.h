#ifndef TARRAYOFARRAY_H
#define TARRAYOFARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "TObject.h"
#include "AliObjectArray.h"

class TArrayOfArray: public TObject { 
public:    
  virtual void * At(UInt_t index0, UInt_t index1)=0;
  //get pointer to the object
  virtual void  Dump(UInt_t index0, UInt_t index1)=0;
  virtual Int_t Resize(Int_t index, UInt_t newsize)=0; 
  //expand array with index index to newsize 
  virtual UInt_t Push(UInt_t size)=0;
  //make new array with size  - return starting index
  virtual Int_t ArraySize(UInt_t index)=0;
  virtual Int_t GetSize()=0;
  ClassDef(TArrayOfArray,1) 
};


class TArrayOfArray_vStack: public TArrayOfArray{
public:  
  TArrayOfArray_vStack();   
  TArrayOfArray_vStack(const char *classname);
  ~TArrayOfArray_vStack(); 
  Bool_t SetClass(const char * classname);
  virtual void Clear(Option_t * opt="");
  void * Unchecked1DArray(UInt_t index){return fIndex->Unchecked1DAt(index);}
  void * Unchecked1DAt(UInt_t index0, UInt_t index1);
  //
  virtual void * At(UInt_t index0, UInt_t index1);
  //get pointer to the object 
  virtual void  Dump(UInt_t index0, UInt_t index1);
  virtual Int_t Resize(Int_t index, UInt_t newsize); 
  //expand array with index index to newsize 
  virtual UInt_t Push(UInt_t size);
  //make new array with size  - return starting index    
  virtual Int_t ArraySize(UInt_t index);
  Int_t ArraySize(){ return fArray->GetSize();}
  virtual Int_t GetSize(){return (fIndex) ?(Int_t)fIndex->GetSize()-1:0;}
private:
  AliObjectArray * fIndex;
  AliObjectArray * fArray;
  ClassDef(TArrayOfArray_vStack,1) 
};

class TArrayOfArray_vList: public TArrayOfArray{
protected:
  AliObjectArray  fIndex;
  AliObjectArray  fSecondaryIndexes;
  AliObjectArray  fArray;
  ClassDef(TArrayOfArray_vList,1) 
};

inline void * TArrayOfArray_vStack::Unchecked1DAt(UInt_t index0, UInt_t index1)
{
  // unchecked return
  return fArray->Unchecked1DAt(((UInt_t*)fIndex->GetArray())[index0]+index1);
}

#endif
  

