#ifndef ALICLASSINFO_H
#define ALICLASSINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  AliClassInfo                                                             //
//                                                                           //
// Defined to make unified interface to primitive types (in Root described by//
// TDataType) and TClass.                                                    //
//                                                                           //
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk         // 
///////////////////////////////////////////////////////////////////////////////


#include "TNamed.h"
class TDataType;



class AliClassInfo : public TNamed { 
public:   
  AliClassInfo(){;}
  virtual ~AliClassInfo(){;}
  virtual void CTORBuffer(void * pointer, UInt_t size=1){;}
    //  {return (*ctorbuffer)(p,size);}
  virtual void DTORBuffer(void * pointer, UInt_t size=1){;}
    //{return (*dtorbuffer)(p,size);}
  virtual void StreamBuffer(TBuffer& b, const void *object, UInt_t size){;}
    //{return (*streamb)(b,object,size);}
  virtual void ObjectDump(void *p){;}
  virtual const char  * GetClassName(){ return 0;}
  virtual TClass *    GetClass(){return 0;} 
  virtual TDataType * GetDataType(){return 0;}
  const UInt_t Size(){return fSize;}
  static AliClassInfo * FindClassInfo(const char * name);
  static AliClassInfo * GenerClassInfo(const char * clname);
  static void  GenerClassInfoCode(const char * clname, Bool_t load,
				    const char *incpath, const char *outfile);
  const TList  &  GetListOfClass(){return fgList;} 
protected:
  static Bool_t GenerClassInterface(const char * clname, FILE * fout);
  static TList fgList; // list of loaded class
  UInt_t   fSize;        //size of object
  ClassDef(AliClassInfo,1) 
};


#endif //ALICLASSINFO_H

