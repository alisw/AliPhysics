#ifndef ALIL3ROOTTYPES_H
#define ALIL3ROOTTYPES_H


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Basic types used by level3                                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>

//#define no_root //switch on the root
#define use_aliroot //use AliL3FileHandler

#ifndef no_root
#define use_root
#include <TObject.h>
#include <Rtypes.h>

#else

#ifndef ROOT_Rtypes
//---- types -------------------------------------------------------------------

typedef char           Char_t;      //Signed Character 1 byte
typedef unsigned char  UChar_t;     //Unsigned Character 1 byte
typedef short          Short_t;     //Signed Short integer 2 bytes
typedef unsigned short UShort_t;    //Unsigned Short integer 2 bytes
#ifdef R__INT16
typedef long           Int_t;       //Signed integer 4 bytes
typedef unsigned long  UInt_t;      //Unsigned integer 4 bytes
#else
typedef int            Int_t;       //Signed integer 4 bytes
typedef unsigned int   UInt_t;      //Unsigned integer 4 bytes
#endif
#ifdef R__B64
typedef int            Seek_t;      //File pointer
typedef long           Long_t;      //Signed long integer 4 bytes
typedef unsigned long  ULong_t;     //Unsigned long integer 4 bytes
#else
typedef int            Seek_t;      //File pointer
typedef long           Long_t;      //Signed long integer 8 bytes
typedef unsigned long  ULong_t;     //Unsigned long integer 8 bytes
#endif
typedef float          Float_t;     //Float 4 bytes
typedef double         Double_t;    //Float 8 bytes
typedef char           Text_t;      //General string
typedef unsigned char  Bool_t;      //Boolean (0=false, 1=true)
typedef unsigned char  Byte_t;      //Byte (8 bits)
typedef short          Version_t;   //Class version identifier
typedef const char     Option_t;    //Option string
typedef int            Ssiz_t;      //String size
typedef float          Real_t;      //TVector and TMatrix element type

typedef void         (*VoidFuncPtr_t)();  //pointer to void function


//---- constants ---------------------------------------------------------------

#ifndef NULL
#define NULL 0
#endif

const Bool_t kTRUE   = 1;
const Bool_t kFALSE  = 0;

const Int_t  kMaxInt      = 2147483647;
const Int_t  kMaxShort    = 32767;
const size_t kBitsPerByte = 8;
const Ssiz_t kNPOS        = ~(Ssiz_t)0;


//---- ClassDef macros ---------------------------------------------------------


#define ClassDef(name,id) 

#define ClassImp(name) 
#endif  //end of Rtypes 

#endif  //end of root selection

//---- Timms AliL3EventDataType  from AliL3EventDataType.h

union AliL3EventDataTypeRoot{
  ULong_t      fID;
  unsigned char      fDescr[8];
};

typedef union AliL3EventDataTypeRoot AliL3EventDataTypeRoot;

#define ROOT_UNKNOWN_DATAID               (((ULong_t)'UNKN')<<32 | 'OWN ')
#define ROOT_COMPOSITE_DATAID             (((ULong_t)'COMP')<<32 | 'OSIT')
#define ROOT_ADCCOUNTS_DATAID             (((ULong_t)'ADCC')<<32 | 'OUNT')
#define ROOT_ADCCOUNTS_UNPACKED_DATAID    (((ULong_t)'ADCC')<<32 | 'NTUP')
#define ROOT_CLUSTERS_DATAID              (((ULong_t)'CLUS')<<32 | 'TERS')
#define ROOT_SPACEPOINTS_DATAID           (((ULong_t)'SPAC')<<32 | 'EPTS')
#define ROOT_VERTEXDATA_DATAID            (((ULong_t)'VRTX')<<32 | 'DATA')
#define ROOT_TRACKSEGMENTS_DATAID         (((ULong_t)'TRAC')<<32 | 'SEGS')
#define ROOT_SLICETRACKS_DATAID           (((ULong_t)'SLCT')<<32 | 'RCKS')
#define ROOT_TRACKS_DATAID                (((ULong_t)'TRAC')<<32 | 'KS  ')
#define ROOT_NODELIST_DATAID              (((ULong_t)'NODE')<<32 | 'LIST')
#define ROOT_EVENTTRIGGER_DATAID          (((ULong_t)'EVTT')<<32 | 'RG  ')

#endif
