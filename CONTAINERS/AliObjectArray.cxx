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
//  AliObjectArray                                                           //
//                                                                           // 
//AliObjectArray is an array of clone (identical) objects.                   //
//In comparison with the TClonesArray objects in this array don't need       //
//to derive from TObject. They also don't need RTTI - type information.      // 
//                                                                           //
//Objects type information is stored in object fClassInfo (instance of       //
//the AliClassInfo).                                                         //
//Objects in array are stored sequentialy in buffers. Buffer is standart C++ //
//array of objects. Buffers size is equal for all of the buffers. If we specify
//fBufferSize==0 objects are stored in one big standart C++ array.
//
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//
//Begin_Html 
/*
<img src="../gif/AliObjectArray.gif">
*/
///////////////////////////////////////////////////////////////////////////////


#include "AliObjectArray.h"

ClassImp(AliObjectArray) 

AliObjectArray::AliObjectArray():AliMemArray()
{ 
  //
  //default constructor
  //
  fClassInfo = 0;
}

AliObjectArray::AliObjectArray(const char * classname, Int_t buffersize):AliMemArray()
{
  //
  // AliObject array constructor 
  // Set class information  fClassInfo according specified type and set the array size -size
  // Buufer size fBufferSize 
  fClassInfo = 0;
  fBufferSize = buffersize;
  SetClass(classname);
}

AliObjectArray::AliObjectArray(const AliObjectArray &arr)
{  
  //
  //
  *((AliMemArray*)this) = *((AliMemArray*)&arr);
  fClassInfo = arr.GetClassInfo();
}


AliObjectArray & AliObjectArray::operator=(const AliObjectArray &arr)
{  
  //
  //
  *((AliMemArray*)this) = *((AliMemArray*)&arr);
  fClassInfo = arr.GetClassInfo();
  return (*this);
}

AliObjectArray::~AliObjectArray()  
{
  //
  //default destructor
  Delete();  
}

Bool_t AliObjectArray::SetClass(const char * classname) 
{
  //
  // Set class information fClassInfo  according class name
  //  
  Delete();
  
  fClassInfo = AliClassInfo::GenerClassInfo(classname);
  fObjectSize = fClassInfo->Size();
  return (fClassInfo!=0);  
}

void   AliObjectArray::Dump(Int_t i)
{
  //dump object at position i 
  if (At(i)) fClassInfo->ObjectDump(At(i));
  else printf("index %d - out of range\n",i);
}

void AliObjectArray::Streamer(TBuffer& R__b) 
{
  //
  //Stream of the AliVector2D
  //
  TString s; 
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TObject::Streamer(R__b); 
    s.Streamer(R__b);  //read class info
    if (s!=" ") SetClass(s.Data());
    else fClassInfo=0;
    R__b >> fBufferSize;
    R__b >> fObjectSize;
    Int_t size;
    R__b >> size;
    Resize(size); 
    if (fSize>0){
      if (fBufferSize==0) fClassInfo->StreamBuffer(R__b, GetArray(), fSize);
      else
	for (UInt_t i=0;i<GetNBuffers();i++)
	  fClassInfo->StreamBuffer(R__b, GetRow(i), fBufferSize);
    }   
  } else {
    R__b.WriteVersion(AliObjectArray::IsA());
    TObject::Streamer(R__b);     
    if (fClassInfo) s = fClassInfo->GetClassName();
    else s=" ";
    s.Streamer(R__b);
    R__b << fBufferSize;
    R__b << fObjectSize;
    R__b << fSize;
   
    if (fSize>0){
      if (fBufferSize==0) fClassInfo->StreamBuffer(R__b, GetArray(), fSize);
      else
	for (UInt_t i=0;i<GetNBuffers();i++)
	  fClassInfo->StreamBuffer(R__b, GetRow(i), fBufferSize);
    }
  }
}





