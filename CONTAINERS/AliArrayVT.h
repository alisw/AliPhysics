#ifndef ALIARRAYVT_H
#define ALIARRAYVT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


class AliArrayVT: public TObject{
  //
  //function object
  //for AliObjectArray
  //
public: 
  AliArrayVT(){cmp=0; mkbuffer=0; delbuffer=0; stream=0; dump=0;fClass=0;}

  int      Compare(void* p1, void* p2) { return (*cmp)(p1,p2);}
  char *   MakeBuffer(UInt_t size) { return (*mkbuffer)(size);}
  void     DeleteBuffer(void * p) { return (*delbuffer)(p);}  
  void     CTORBuffer(void *p,  const UInt_t size) { return (*ctorbuffer)(p,size);}
  void     DTORBuffer(void * p, const UInt_t size) { return (*dtorbuffer)(p,size);}
  void     ObjectCTOR(void *p) { return (*objectctor)(p);}
  void     ObjectDTOR(void * p) { return (*objectdtor)(p);}  
  void     StreamObject(TBuffer& b, void * object) {return (*stream)(b, object);} 
  void     StreamBuffer(TBuffer& b, const void *object, UInt_t size) {return (*streamb)(b,object,size);}
  void     ClassDump(void *p) { return (*dump)(p);}
public: 
  int     (*cmp)(void*, void*);
  char*   (*mkbuffer)(UInt_t );
  void    (*delbuffer)(void*);
  void    (*objectctor)(void*);
  void    (*objectdtor)(void*);    
  void    (*ctorbuffer)(void*, UInt_t size);
  void    (*dtorbuffer)(void*, UInt_t size);  
 
  void    (*stream)(TBuffer &, void *);
  void    (*streamb)(TBuffer &, const void *,UInt_t);
  void    (*dump)(void*);  
  TString  fClassName;   //class name of the object
  TClass * fClass;       //class type of the object
  UInt_t   fSize;        //size of object
  ClassDef(AliArrayVT,0) 
};

#endif //ALIARRAYVT

