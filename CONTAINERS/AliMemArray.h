#ifndef ALIMEMARRAY_H
#define ALIMEMARRAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  AliMemArray                                                              //                      
//  (Pseudo)Container optimised for fast  random access operator[ ].         //      
//  Operator []  time doesn’t depend on the number of elements in container O(1) - 
//    like in standard array or like in STL vector
///////////////////////////////////////////////////////////////////////////////
#include "TObject.h"

class AliMemArray: public TObject { 
public:
  AliMemArray();
  AliMemArray(Int_t objectSize, Int_t buffersize=0);
  AliMemArray(const AliMemArray & arr); 
  AliMemArray & operator = (const AliMemArray &arr);
  virtual ~AliMemArray();
  void Swap(AliMemArray & arr);
  //
  inline void *  UncheckedAt(UInt_t index) const; //
  inline void *  At(UInt_t i) const ;     //controled return pointer to object 
  inline void *  Unchecked1DAt(UInt_t i) const;  //return pointer to the object  
  inline void *  Unchecked2DAt(UInt_t i) const ;  //return pointer to the object
  inline void *  Unchecked2DAt(UInt_t i, UInt_t j) const ;  //return pointer to the object
    
  Bool_t   Is1D(){return (fBufferSize==0);}  //true if 1D array
  Bool_t   Is2D(){return (fBufferSize!=0);}  //true if 2D array
  //
  const   void * GetArray()const {return fCont;} //return pointer to the buffer
  void * GetRow(UInt_t row)const {return fBufferSize? ((void**)fCont)[row]: fCont;} 

                                                        // return pointer to the object
  //    
  const UInt_t  GetCapacity() const {return fCapacity;}     //return number of stored objects 
  const UInt_t  GetSize() const {return fSize;}     //return number of stored objects        
  //
  void   Delete(Option_t *option=""); 
  //delete memory space occupated by the array   
  void   Clear(Option_t *option="");  
  //clear memory space occupied by the array - doesn't call destructor
  void   Resize(UInt_t n);
  void   Reserve(UInt_t n);  
  //
  const    UInt_t GetBufferSize() const {return fBufferSize;}
  const    UInt_t GetObjectSize() const {return fObjectSize;}
  const UInt_t  GetNBuffers() const {return fBufferSize? fCapacity/fBufferSize :1;} 
  void     SetBufferSize(UInt_t bufsize); 
protected :  
  void     SetObjectSize(UInt_t size);
  void   Delete2D(); //delete memory space occupated by the array  
  void   Clear2D();  //clear  memory space occupated by the array  
  void   Resize2D(UInt_t n);   
  void   Reserve2D(UInt_t n);     
  //  
protected:      
  virtual void     CTORBuffer(void *buffer, UInt_t size){;} //array placement constructor
  virtual void     DTORBuffer(void *buffer, UInt_t size){;} //array placement destructor
  virtual void     CopyBuffer(void *src, void *dest,  UInt_t size); //array placement copy constructor
  UInt_t          fSize;             //total number of valid  objects  
  UInt_t          fCapacity;         //capacity of array 
  UInt_t          fObjectSize;       //           object size
  UInt_t          fBufferSize;       //number of object in one full  buffer (0 means array in one buffer)
  // 
private:  
  void *          fCont;             //!data buffer      
  ClassDef(AliMemArray,0) 
};

void *  AliMemArray::Unchecked1DAt(UInt_t i) const 
{
  return  &(((char*)fCont)[fObjectSize*i]);
}

void *  AliMemArray::Unchecked2DAt(UInt_t i, UInt_t j) const 
{
   return &(  ((char**)fCont)[i] [j*fObjectSize]); 
}

//********************************************************************
void  *  AliMemArray::Unchecked2DAt(UInt_t i) const 
{  
  //
  //  called in GePointer if we have more then one buffer 
  return &(  ((char**)fCont)[i/fBufferSize] [(i%fBufferSize)*fObjectSize])  ;
}

void *  AliMemArray::UncheckedAt(UInt_t i) const
{
 
  return (!fBufferSize) ?Unchecked1DAt(i): Unchecked2DAt(i); 
}

void * AliMemArray::At(UInt_t i) const 
{
  return  (i<fSize) ? UncheckedAt(i):0;     //controled return pointer to object 
}

#endif
