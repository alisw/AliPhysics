// XEMacs -*-C++-*-
#ifndef ALIHLTHOMERDATA_H
#define ALIHLTHOMERDATA_H

/************************************************************************
**
**
** This file is property of and copyright by the Technical Computer
** Science Group, Kirchhoff Institute for Physics, Ruprecht-Karls-
** University, Heidelberg, Germany, 2001
** This file has been written by Timm Morten Steinbeck, 
** timm@kip.uni-heidelberg.de
**
**
** See the file license.txt for details regarding usage, modification,
** distribution and warranty.
** Important: This file is provided without any warranty, including
** fitness for any particular purpose.
**
**
** Newer versions of this file's package will be made available from 
** http://web.kip.uni-heidelberg.de/Hardwinf/L3/ 
** or the corresponding page of the Heidelberg Alice Level 3 group.
**
*************************************************************************/

/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/

#ifdef USE_ROOT
#include "Rtypes.h"
#endif
#include <limits.h>


// Determine the sizes of the different integer type
// homer_uint32, homer_uint64
#if !defined(USE_ROOT) && !defined(__CINT__)
// First homer_uint32
#if USHRT_MAX==4294967295
typedef unsigned short homer_uint32;
#else // USHRT_MAX==4294967295

#if UINT_MAX==4294967295
typedef unsigned int homer_uint32;
#else // UINT_MAX==4294967295

#if ULONG_MAX==4294967295l
typedef unsigned long homer_uint32;
#else // ULONG_MAX==4294967295l

#error Could not typedef homer_uint32

#endif // ULONG_MAX==4294967295l

#endif // UINT_MAX==4294967295

#endif // USHRT_MAX==4294967295

// Then homer_uint16
#if USHRT_MAX==65535
typedef unsigned short homer_uint16;
#else // USHRT_MAX==65535

#if UINT_MAX==65535
typedef unsigned int homer_uint16;
#else // UINT_MAX==65535

#if ULONG_MAX==65535
typedef unsigned long homer_uint16;
#else // ULONG_MAX==65535

#error Could not typedef homer_uint16

#endif // ULONG_MAX==65535

#endif // UINT_MAX==65535

#endif // USHRT_MAX==65535

// Then homer_uint64
#if ULONG_MAX==18446744073709551615ULL
typedef unsigned long homer_uint64;
#else // ULONG_MAX==18446744073709551615UL

#if defined __GNUC__
typedef unsigned long long homer_uint64;
#else // defined __GNUC__

#if defined __SUNPRO_CC
typedef unsigned long long homer_uint64;
#else // defined __SUNPRO_CC

#error Could not typedef homer_uint64

#endif // defined __SUNPRO_CC
#endif // defined __GNUC__

#endif // ULONG_MAX==18446744073709551615UL

typedef unsigned char homer_uint8;

#else // !USE_ROOT && !CINT


typedef UShort_t homer_uint16;
typedef UInt_t homer_uint32;
typedef ULong64_t homer_uint64;
typedef Byte_t homer_uint8;

#ifdef __CINT__
typedef int key_t;
#endif

#endif // USE_ROOT

//typedef homer_uint64 AliEventID_t;



const unsigned kAttribute_8b_StartOffset = 0;
const unsigned kByteOrderAttribute_8b_Offset = kAttribute_8b_StartOffset+0;
const unsigned kVersionAttribute_8b_Offset = kAttribute_8b_StartOffset+1;
const unsigned kAlignment_8b_StartOffset = 24;
const unsigned kUInt64Alignment_8b_Offset = kAlignment_8b_StartOffset+0;
const unsigned kUInt32Alignment_8b_Offset = kAlignment_8b_StartOffset+1;
const unsigned kUInt16Alignment_8b_Offset = kAlignment_8b_StartOffset+2;
const unsigned kUInt8Alignment_8b_Offset = kAlignment_8b_StartOffset+3;
const unsigned kDoubleAlignment_8b_Offset = kAlignment_8b_StartOffset+4;
const unsigned kFloatAlignment_8b_Offset = kAlignment_8b_StartOffset+5;


const unsigned kID_64b_Offset = 1;
const unsigned kLength_64b_Offset = 2;
const unsigned kType_64b_Offset = 4;
const unsigned kSubType1_64b_Offset = 5;
const unsigned kSubType2_64b_Offset = 6;
const unsigned kBirth_s_64b_Offset = 7;
const unsigned kBirth_us_64b_Offset = 8;
const unsigned kProducerNode_64b_Offset = 9;
const unsigned kOffset_64b_Offset = 10;
const unsigned kSize_64b_Offset = 11;
const unsigned kStatusFlags_64b_Offset = 12;
const unsigned kEnd_64b_Offset = 13;
const unsigned kCount_64b_Words = kEnd_64b_Offset;


// Possible values for fAttributes[kByteOrderAttribute]
/* Keep this consistent with BCLNetworkData.hpp kLittleEndian/kBigEndian and AliHLTSubEventDataDescriptor.hpp */
const homer_uint8 kHOMERUnknownByteOrder      = 0;
const homer_uint8 kHOMERLittleEndianByteOrder = 1;
const homer_uint8 kHOMERBigEndianByteOrder    = 2;
#ifdef __i386__
    const homer_uint8 kHOMERNativeByteOrder = kHOMERLittleEndianByteOrder;
#else
#ifdef __arm__
    const homer_uint8 kHOMERNativeByteOrder = kHOMERLittleEndianByteOrder;
#else
#ifdef __x86_64__
    const homer_uint8 kHOMERNativeByteOrder = kHOMERLittleEndianByteOrder;
#else
#ifdef __ia64__
    const homer_uint8 kHOMERNativeByteOrder = kHOMERLittleEndianByteOrder;
#else
#if defined(__powerpc__)
    const homer_uint8 kHOMERNativeByteOrder = kHOMERBigEndianByteOrder;
#else
#ifdef __CINT__
    const homer_uint8 kHOMERNativeByteOrder = kHOMERLittleEndianByteOrder;
#warning Assuming little endian format for __CINT__
#else
#error Byte format (little/big endian) currently not defined for platforms other than intel i386 compatible, x86-64 (AMD64) and arm...
#endif
#endif
#endif
#endif
#endif
#endif


//#define HOMER_BLOCK_DESCRIPTOR_TYPEID             (((homer_uint64)'HOBL')<<32 | 'KDES')
#define HOMER_BLOCK_DESCRIPTOR_TYPEID               ( (((homer_uint64)'H')<<56)|(((homer_uint64)'O')<<48)|(((homer_uint64)'B')<<40)|(((homer_uint64)'L')<<32)|(((homer_uint64)'K')<<24)|(((homer_uint64)'D')<<16)|(((homer_uint64)'E')<<8)|(((homer_uint64)'S')<<0) )
//#define HOMER_BLOCK_DESCRIPTOR_TYPEID            ( (((homer_uint64)'H')<<56)|(((homer_uint64)'O')<<48) )
#define HOMER_HEADER_CURRENT_VERSION              2


class HOMERBlockDescriptor
    {
    public:

	static unsigned GetHOMERBlockDescriptorSize()
		{
		return sizeof(homer_uint64)*kCount_64b_Words;
		}

	HOMERBlockDescriptor( void* header = 0 )
		{
		fHeader = header;
		}
	void UseHeader( void* header )
		{
		fHeader = header;
		}
	void Initialize()
		{
		if ( fHeader )
		    {
		    for ( unsigned ii=0; ii<kCount_64b_Words; ii++ )
			((homer_uint64*)fHeader)[ ii ] = (homer_uint64)0;
		    ((homer_uint64*)fHeader)[ kID_64b_Offset ] = HOMER_BLOCK_DESCRIPTOR_TYPEID;
		    ((homer_uint64*)fHeader)[ kLength_64b_Offset ] = GetHOMERBlockDescriptorSize();
		    ((homer_uint8*)fHeader)[ kByteOrderAttribute_8b_Offset ] = kHOMERNativeByteOrder;
		    ((homer_uint8*)fHeader)[ kVersionAttribute_8b_Offset ] = HOMER_HEADER_CURRENT_VERSION;
		    }
		}

	void SetByteOrder( homer_uint8 bo )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kByteOrderAttribute_8b_Offset ] = bo;
		}
	homer_uint8 GetByteOrder() const 
		{
		if ( fHeader )
		    return ((homer_uint8*)fHeader)[ kByteOrderAttribute_8b_Offset ];
		else
		    return 0xFF;
		}
	void SetVersion( homer_uint8 v )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kVersionAttribute_8b_Offset ] = v;
		}
	void SetID( homer_uint64 id  )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kID_64b_Offset ] = id;
		}
	void SetHeaderLength( homer_uint64 l  )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kLength_64b_Offset ] = l;
		}
	homer_uint64 GetHeaderLength()
		{
		if ( fHeader )
		    return ((homer_uint64*)fHeader)[ kLength_64b_Offset ];
		else
		    return 0;
		}
	void SetAlignment( homer_uint8 type, homer_uint8 align )
		{
		if ( fHeader && type<6 )
		    ((homer_uint8*)fHeader)[ kAlignment_8b_StartOffset+type ] = align;
		}
	void SetUInt64Alignment( homer_uint8 align )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kUInt64Alignment_8b_Offset ] = align;
		}
	void SetUInt32Alignment( homer_uint8 align )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kUInt32Alignment_8b_Offset ] = align;
		}
	void SetUInt16Alignment( homer_uint8 align )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kUInt16Alignment_8b_Offset ] = align;
		}
	void SetUInt8Alignment( homer_uint8 align )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kUInt8Alignment_8b_Offset ] = align;
		}
	void SetDoubleAlignment( homer_uint8 align )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kDoubleAlignment_8b_Offset ] = align;
		}
	void SetFloatAlignment( homer_uint8 align )
		{
		if ( fHeader )
		    ((homer_uint8*)fHeader)[ kFloatAlignment_8b_Offset ] = align;
		}
	void SetType( homer_uint64 t )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kType_64b_Offset ] = t;
		}
	void SetSubType1( homer_uint64 st1 )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kSubType1_64b_Offset ] = st1;
		}
	void SetSubType2( homer_uint64 st2 )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kSubType2_64b_Offset ] = st2;
		}
	void SetBirth_s( homer_uint64 bs )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kBirth_s_64b_Offset ] = bs;
		}
	void SetBirth_us( homer_uint64 bus )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kBirth_us_64b_Offset ] = bus;
		}
	void SetProducerNode( homer_uint64 pn )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kProducerNode_64b_Offset ] = pn;
		}
	void SetBlockOffset( homer_uint64 bo )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kOffset_64b_Offset ] = bo;
		}
	homer_uint64 GetBlockOffset()
		{
		if ( fHeader )
		    return ((homer_uint64*)fHeader)[ kOffset_64b_Offset ];
		else
		    return 0;
		}
	void SetBlockSize( homer_uint64 bs )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kSize_64b_Offset ] = bs;
		}
	homer_uint64 GetBlockSize()
		{
		if ( fHeader )
		    return ((homer_uint64*)fHeader)[ kSize_64b_Offset ];
		else
		    return 0;
		}
	void SetStatusFlags( homer_uint64 bs )
		{
		if ( fHeader )
		    ((homer_uint64*)fHeader)[ kStatusFlags_64b_Offset ] = bs;
		}
	homer_uint64 GetStatusFlags()
		{
		if ( fHeader )
		    return ((homer_uint64*)fHeader)[ kStatusFlags_64b_Offset ];
		else
		    return 0;
		}

	void* GetHeader() const
		{
		return fHeader;
		}
		
    protected:
	void* fHeader;
	

    };



/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/

#endif // ALIHLTHOMERDATA_H
