/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
*/
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  AliMemArray                                                              //                      
//  (Pseudo)Container optimised for fast  random access operator[ ]. 
//  Operator []  time doesn’t depend on the number of elements in container O(1) - 
//    like in standard array or like in STL vector
//   To achieve maximally-fast indexing and iteration in one buffer mode  
//   the vector maintains its storage as a single contiguous array of objects (one buffer mode)

//   When a vector runs out of pre-allocated storage, in order to maintain it 
//   contiguous array it must allocate a whole new (larger) chunk of storage  
//   elsewhere and copy the objects to the new storage.
//       void *  AliMemArray::Unchecked1DAt(UInt_t i) const                    
//             return  &(((char*)fCont)[fObjectSize*i]);

//   In multi buffer mode (like two dimensional array) when a vector runs out of 
//   pre-allocated storage we don’t need to copy whole array only small buffer but operator [] is slower
//       void *  AliMemArray::Unchecked2DAt(UInt_t i, UInt_t j) const                    
//       return &(  ((char**)fCont)[i] [j*fObjectSize]); 
//       void  *  AliMemArray::Unchecked2DAt(UInt_t i) const 
//       return &(  ((char**)fCont)[i/fBufferSize] [(i%fBufferSize)*fObjectSize])  ;


//Begin_Html
//<img src="../gif/AliMemArray.gif">
//End_Html

   
//  Streamer CTORBuffer and DTORBuffer are  virtual - should be implemented in derived
//  classes. For example AliObjectArray derived from AliMemArray is general array 
//  for objects with defined AliClassInfo information. 
//                                                                           //
//                                                                          //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
///////////////////////////////////////////////////////////////////////////////



#include "AliMemArray.h"
#include "iostream.h"
#include "TMath.h"
#include "TError.h"


ClassImp(AliMemArray)
 

AliMemArray::AliMemArray()
{ 
  //
  //default constructor
  //
  fCont = 0;
  fSize = 0;
  fCapacity=0;
  fObjectSize = 0;
  fBufferSize = 0;
}

AliMemArray::AliMemArray(Int_t objectSize, Int_t buffersize)
{
  //
  // AliMemArray constructor 
  fCont = 0;
  fSize = 0;
  fCapacity=0;
  fObjectSize =objectSize;
  fBufferSize  =buffersize;
}

AliMemArray::AliMemArray(const AliMemArray & arr) 
{
  //
  //copy constructor
  fCont = arr.fCont;
  fSize = arr.fSize;
  fCapacity = arr.fCapacity;
  fObjectSize = arr.fObjectSize;
  fBufferSize =arr.fBufferSize;
  
  if (arr.fBufferSize==0) {
    fCont = new char[fCapacity*fObjectSize];
    CopyBuffer(fCont, arr.fCont, fSize);
  }
  else{ 
    Int_t buffers = fCapacity/fBufferSize;
    if (fCapacity%fBufferSize) buffers++;     
    void ** table  = (void**) new char*[buffers];   
    for (Int_t i=0; i<buffers;i++){
      table[i] = new char[fBufferSize*fObjectSize];
      Int_t size = fSize - i*fBufferSize;
      if (size > (Int_t)fBufferSize) size = fBufferSize;
      if (size >0) CopyBuffer(table[i], ((void**)arr.fCont)[i],size);
    }
    fCont = (void*)table;
  }
}
    

void AliMemArray::Swap(  AliMemArray &arr)
{
  //swap contents of array
  UInt_t size = arr.fSize;
  arr.fSize = fSize;
  fSize = size;
  UInt_t capacity = arr.fCapacity;
  arr.fCapacity = fCapacity;
  fCapacity = capacity;
  UInt_t objectSize = arr.fObjectSize;
  arr.fObjectSize = fObjectSize;
  fObjectSize = objectSize;
  UInt_t bufferSize = arr.fBufferSize;
  arr.fBufferSize = fBufferSize;
  fBufferSize = bufferSize;
  void * cont = arr.fCont;
  arr.fCont = fCont;
  fCont = cont;
}

void     AliMemArray::CopyBuffer(void *dest, void *src,  UInt_t size)
{
  //
  //array placement copy constructor
  memcpy(dest, src,size*fObjectSize);
}

AliMemArray::~AliMemArray()  
{
  //
  //default destructor
  Delete();
}


void   AliMemArray::SetObjectSize(UInt_t  bufsize)
{
  //
  //set memory size for one object - it can be changed only when the array is empty 
  if (fCont){
    ::Error("AliMemArray::SetObjectSize", "forbidden to resize just allocated objects");
    return;
  };
  fBufferSize=bufsize;
}

AliMemArray & AliMemArray::operator = (const AliMemArray &arr)
{
  //
  //
  AliMemArray  tmparr(arr);
  Swap(tmparr);
  return *this;
}

void   AliMemArray::SetBufferSize(UInt_t  bufsize)
{
  //
  //set buffer size - it can be changed only when the array is empty 
  if (fCont==0) {
    fBufferSize = bufsize;
    return;
  }
  if (fBufferSize == bufsize) return;
  
  if (bufsize==0){
    char *table = new char[fObjectSize*fCapacity];
    char * p = table;
    for (UInt_t i=0; i<fSize; i++,p+=fObjectSize) 
      memcpy(p, At(i), fObjectSize);
    //delete [](char*)fCont;
    Delete();
    fCont = table;
    fBufferSize = bufsize;
  }
  else{
    Int_t buffers = fCapacity/bufsize;
    if (fCapacity%bufsize) buffers++;     
    char ** table  =  new char*[buffers];   
    for (Int_t ibuf=0; ibuf<buffers;ibuf++){
      table[ibuf] = new char[bufsize*fObjectSize];
      Int_t size = fSize - ibuf*bufsize;
      if (size > (Int_t)bufsize) size = bufsize;      
      if (size >0) for ( Int_t ip=0;ip<size;ip++)
	memcpy(&table[ibuf][ip*fObjectSize], At(ibuf*bufsize+ip), fObjectSize);
    }
    //    delete [](char**)fCont;
    Delete();
    fCont = (void*)table;
    fBufferSize = bufsize;
  }

}



void AliMemArray::Delete(Option_t *)
{
  //
  //delete memory space occupied by the array  - 
  //Use this routine when your objects allocate
  //memory (e.g. objects inheriting from TNamed or containing TStrings
  //allocate memory). If not you better use Clear() since if is faster.  
  if (fCont){
    if (fBufferSize) {
      Delete2D();
      fSize = 0;
      fCapacity = 0;
      return;
    }
    DTORBuffer(Unchecked1DAt(0),fSize);
    delete [] (char*)fCont;
    fCont = 0;
    fSize = 0;
    fCapacity = 0;
  }
}

void AliMemArray::Clear(Option_t *)
{
  //
  //clear array   - 
  // Only use this routine when your objects don't
  // allocate memory since it will not call the object dtors.
  if (fBufferSize){
    Clear2D();
    return;
  }
  if (fCont){
    memset(fCont, 0, fSize*fObjectSize);
    fSize = 0;
  }
}

void AliMemArray::Reserve(UInt_t  n)
{
  //
  //reserve arrays space
  //  
  if (fObjectSize<=0) {
    cout<<"Object length not defined\n";
    return;
  }
  if (n==fCapacity) return;
  
  if (fBufferSize>0) {
    Reserve2D(n); //if 2D buffer
    return;
  }
  //
  if (fCapacity){
    if (fSize>n) {
      DTORBuffer(Unchecked1DAt(n),fSize-n);
      memset(&((char*)fCont)[n*fObjectSize], 0, (fSize-n)*fObjectSize); 
      fSize =n;
    }
    fCont = (char*)TStorage::ReAlloc(fCont, n*fObjectSize,fCapacity*fObjectSize); 
  }
  else  fCont = new char[n*fObjectSize];

  if (!fCont) fCapacity = 0;
  else fCapacity = n;
}


void AliMemArray::Resize(UInt_t  n)
{
  //
  //resize buffer
  //   
  if (fObjectSize<=0) {
    cout<<"Object length not defined\n";
    return;
  }
  if (fBufferSize>0) {
     Resize2D(n); //if 2D buffer
     return;
  }
  //
  if (n>fCapacity) Reserve(n);   //reserve automaticaly space if sie >capacity
  if (fSize>n){ 
    DTORBuffer(Unchecked1DAt(n),fSize-n);
    memset(&((char*)fCont)[n*fObjectSize], 0, (fSize-n)*fObjectSize);
  }
  if (fSize<n)    CTORBuffer(Unchecked1DAt(fSize),n-fSize);     
  fSize = n;
  return; 
}

void AliMemArray::Delete2D()
{
  //
  //delete memory space occupied by the array 
  if (!fBufferSize) return;

  Int_t  nbuff = (fCapacity/fBufferSize);
  if ( (fCapacity%fBufferSize)!=0) nbuff++;
  for (Int_t  i=0;i<nbuff;i++) {    
    Int_t size = fSize-i*fBufferSize;
    if  (size>0)
      DTORBuffer(GetRow(i),UInt_t(size)<fBufferSize? size:fBufferSize);    
    delete [] (char*)GetRow(i);
  }
  delete [] (void**)fCont;
  fCont =0;
  fSize = 0;  
  fCapacity = 0;
}

void AliMemArray::Clear2D()
{
  //
  //clear memory space occupied by the array  - doesn't call DTOR
  Int_t  nbuff = (fCapacity/fBufferSize);
  if ( (fCapacity%fBufferSize)!=0) nbuff++;
  for (Int_t  i=0;i<nbuff;i++) memset(GetRow(i), 0, fSize*fObjectSize);
  fSize = 0;  
}




void AliMemArray::Reserve2D(UInt_t  n) 
{
  //
  // 
  Int_t buffers = n/fBufferSize;
  if (n%fBufferSize) buffers++;
  UInt_t  nobjects=buffers*fBufferSize;
  Int_t oldbuffers = GetNBuffers() ;
  if (buffers==oldbuffers) return;
  //
  void ** table  = (void**) new char*[buffers];

  Int_t max = buffers>oldbuffers ? buffers: oldbuffers;
  for (Int_t i = 0;i<max;i++) {
    if ( (i<oldbuffers) && (i<buffers))  table[i] = GetRow(i);
    if ( (i<oldbuffers)&&(i>=buffers) ){
      Int_t dsize = TMath::Min(Int_t(fSize) - i*Int_t(fBufferSize),Int_t(fBufferSize));
      if (dsize>0) DTORBuffer(GetRow(i),dsize);
      delete [] (char*)GetRow(i);
    }
    if (i>=oldbuffers)
      table[i] = new char[fBufferSize*fObjectSize];
  }
  if (fSize>nobjects) fSize=nobjects;

  fCapacity = nobjects ;    
  delete [] (void**)fCont;
  fCont = (void*)table;
}



void AliMemArray::Resize2D(UInt_t  n) 
{
  //
  //  
  if (n>fCapacity) Reserve2D(n);   //reserve automaticaly space 
  Int_t buffers = n/fBufferSize;
  if (n%fBufferSize) buffers++;

  if (fSize>n){   //call destructor if we decrease the size of array
    Int_t oldbuffer = fSize/fBufferSize;
    
    for (Int_t i=buffers;i<oldbuffer; i++){
      Int_t iold= fSize-i*fBufferSize;
      if (iold>(Int_t)fBufferSize) iold= fBufferSize; 
      Int_t inew= n -i*fBufferSize;
      if (inew<0) inew =0;
      DTORBuffer(Unchecked2DAt(i,inew),iold-inew);
    }
  }
  if (fSize<n){   //call constructor if we increase the size of array
    Int_t oldbuffer = fSize/fBufferSize;
    for (Int_t i=oldbuffer;i<buffers; i++){
      Int_t iold = fSize-i*fBufferSize;
      if (iold<0) iold = 0;
      Int_t inew =  n -i*fBufferSize;
      if (inew>(Int_t)fBufferSize) inew = fBufferSize;
      CTORBuffer(Unchecked2DAt(i,iold),inew-iold);
    }
  }
  fSize = n;  
}






