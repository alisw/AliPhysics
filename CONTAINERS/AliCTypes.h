#ifndef ALI_CTypes
#define ALI_CTypes
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Ali C types                                                          //
//                     
// MI                                                                   //
// Macros for defining Containers                                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TMemberInspector.h"
#include "TClass.h"
#include "TStorage.h"

class TDumpMembers : public TMemberInspector {
public:
   TDumpMembers() { }
   void Inspect(TClass *cl, const char *parent, const char *name, void *addr);
};

//////////////////////////////////////////////////////////////////////////



#define LClassDef(name,id) \
private: \
   static TClass *fgIsA; \
public: \
   static TClass *Class(); \
   static const char *Class_Name(); \
   static Version_t Class_Version() { return id; } \
   static void Dictionary(); \
   TClass *IsA() const { return name::Class(); } \
   void ShowMembers(TMemberInspector &insp, char *parent); \
   void Streamer(TBuffer &b); \
   friend TBuffer &operator>>(TBuffer &buf, name *&obj); \
   _ClassInit_(name) \
   static const char *DeclFileName() { return __FILE__; } \
   static int DeclFileLine() { return __LINE__; } \
   static const char *ImplFileName(); \
   static int ImplFileLine();  \
   void    *operator new(size_t sz) { return TStorage::ObjectAlloc(sz); } \
   void    *operator new(size_t sz, void *vp) { return vp; } \
   void  Dump();

#define LClassImp(name)  \
  void name::Dump() { \
     char parent[256]; \
     parent[0] = 0;  \
     TDumpMembers dm;  \
     ShowMembers(dm, parent); \
  } \
  ClassImp(name)


/*
#define ClassArrayDef(name) \
private: \
  static AliArrayVT * fgArrayInfo;  \
public: \
  static name * Unchecked1DAt(const AliObjectArray & arr, UInt_t index)\
    { return ((name*)arr.Unchecked1DAt(index));} \
  static name * Unchecked2DAt(const AliObjectArray & arr, UInt_t index)  \
    { return ((name*)arr.Unchecked2DAt(index));} \
  static name * UncheckedAt(const AliObjectArray & arr, UInt_t index)  \
    { return (name*)arr.UncheckedAt(index);}  \
  static name * At(const AliObjectArray & arr, UInt_t index)  \
    { return (arr.GetClassInfo()==fgArrayInfo) ? (name*)arr.At(index): 0;}  \
  static name * CastAt(const AliObjectArray & arr, UInt_t index)  \
    { return (name*) ((arr.GetClass())->DynamicCast(Class(), arr.At(index))) ;}  \
  static AliArrayVT *  GetArrayInfo(); \



//This should be in automatic generated code
#define ClassArrayImp(name) \
AliArrayVT *  name::fgArrayInfo = 0; \
char *name ##__MakeBuffer(UInt_t size) { return (char*)new  name[size];} \
void  name ##__DeleteBuffer(void *p) {delete [] (name*)p;}  \
int   name ##__Cmp(char *p1, char * p2){return ((name*)p1)->Compare((name*)p2);}  \
void  name ##__StreamObject(TBuffer&   b, char * object){ ((name*)object)->Streamer(b);} \
void  name ##__StreamBuffer(TBuffer&   b, char * object, UInt_t size ); \
void  name ##__DumpObject(char *p) {return ((name*)p)->Dump();}  \
void  name ##__StreamBuffer(TBuffer&   b, char * object, UInt_t size ) \
{ \
  for (UInt_t i=0;i<size;i++) ((name*)object)[i].Streamer(b); \
} \
void name ##__DTORBuffer(void *p) { \
  ((name*)p)->~name(); \
} \
void name ##__CTORBuffer(void *p) { \
   new (p)name; \
} \
void name ##__DTORBuffer(void *p, const UInt_t size) { \
  name * pp = (name*)p; \
  name * pmax  = pp+size; \
  while (pp<pmax) (pp++)->~name(); \
} \
void name ##__CTORBuffer(void *p, const UInt_t size) { \
   name * pp = (name*)p; \
   name * pmax  = pp+size; \
  while (pp<pmax) new (pp++)name; \
} \
AliArrayVT *name::GetArrayInfo() \
{  \
  if (!fgArrayInfo) {  \
     fgArrayInfo = new AliArrayVT;    \
     fgArrayInfo->cmp = &(name ##__Cmp); \
     fgArrayInfo->delbuffer = &(name ##__DeleteBuffer); \
     fgArrayInfo->mkbuffer = &(name ##__MakeBuffer);  \
     fgArrayInfo->ctorbuffer = &(name ##__CTORBuffer); \
     fgArrayInfo->dtorbuffer = &(name ##__DTORBuffer);  \
     fgArrayInfo->stream  =  &(name ##__StreamObject); \
     fgArrayInfo->streamb =  &(name ##__StreamBuffer); \
     fgArrayInfo->dump =   &(name ##__DumpObject); \
     fgArrayInfo->fClass = name::Class(); \
     fgArrayInfo->fClassName = "name"; \
     fgArrayInfo->fSize = sizeof(name); \
  }  \
  return fgArrayInfo; \
} 
*/


#endif  // ALI_CTypes
