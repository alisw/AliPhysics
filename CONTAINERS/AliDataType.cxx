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
//  AliDataType                                                              //
//                                                                           //
// Defined to make unified interface to primitive types (in Root described by//
// TDataType) and TClass.                                                    //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////


#include "AliDataType.h"
#include "AliCTypes.h"

#include "TMath.h"
#include "TClass.h"
#include "TDataType.h"   

#include "TROOT.h"
#include "iostream.h"


ClassImp(AliDataType)



AliDataType::AliDataType(const char *name)
{
  //
  // AliData type constructor
  fDataType = new TDataType(name);
  if (fDataType->GetType()<0)  fDataType= (gROOT->GetType(name,kTRUE));
  if ((fDataType) && (fDataType->GetType())){ 
    fSize = fDataType->Size();
    //    fType = (EDataType) fDataType->GetType();     
    SetTitle(name);
    SetName(name);
    fgList.Add(this);
  }
}


const char * AliDataType::GetClassName() 
{ 
  //class name of the object   
  return (fDataType) ? fDataType->GetName():0;
} 
  
void AliDataType::StreamBuffer(TBuffer& b, const void *object, UInt_t size)
{
  //streamer for buffer of objects  
  char * last = &((char*)object)[size*fSize];
  char * pfirst = (char*)object;
  char *p;
  if (b.IsWriting()){
    switch ((EDataType) fDataType->GetType()){
    case kChar_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((Char_t*)p);
      break;
    case kShort_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((Short_t*)p);
      break;
    case kInt_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((Int_t*)p);
    break;
    case kLong_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((Long_t*)p);
      break;
    case kUChar_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((UChar_t*)p);
      break;
    case kUShort_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((UShort_t*)p);
      break;
    case kUInt_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((UInt_t*)p);
      break;
    case kFloat_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((Float_t*)p);
      break;
    case kDouble_t:
      for (p= pfirst; p<last;p+=fSize) b<<*((Double_t*)p);
      break;
    default:
      break;
    }
  }
  else
  switch ((EDataType) fDataType->GetType()){
    case kChar_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((Char_t*)p);
      break;
    case kShort_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((Short_t*)p);
      break;
    case kInt_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((Int_t*)p);
    break;
    case kLong_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((Long_t*)p);
      break;
    case kUChar_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((UChar_t*)p);
      break;
    case kUShort_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((UShort_t*)p);
      break;
    case kUInt_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((UInt_t*)p);
      break;
    case kFloat_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((Float_t*)p);
      break;
    case kDouble_t:
      for (p= pfirst; p<last;p+=fSize) b>>*((Double_t*)p);
      break;
    default:
      break;
    }
  
}

void AliDataType::ObjectDump(void *p) 
{
  //
  // dump object information
  // assume that object p has type described by AliDataTYPE
  switch ((EDataType) fDataType->GetType()){
  case kChar_t:
    cout<<*((Char_t*)p)<<"\n";
    break;
  case kShort_t:
    cout<<*((Short_t*)p)<<"\n";
    break;
  case kInt_t:
    cout<<*((Int_t*)p)<<"\n";
    break;
  case kLong_t:
    cout<<*((Long_t*)p)<<"\n";
    break;
  case kUChar_t:
    cout<<*((UChar_t*)p)<<"\n";
    break;
  case kUShort_t:
    cout<<*((UShort_t*)p)<<"\n";
    break;
  case kUInt_t:
    cout<<*((UInt_t*)p)<<"\n";
    break;
  case kFloat_t:
    cout<<*((Float_t*)p)<<"\n";
    break;
  case kDouble_t:
    cout<<*((Double_t*)p)<<"\n";
    break;
  default:
    break;
  }
}
