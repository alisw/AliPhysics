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

/** @file   AliHLTHOMERWriter.cxx
    @author Timm Steinbeck
    @date   Sep 14 2007
    @brief  HLT Online Monitoring Environment including ROOT - Writer   
    @note   migrated from PubSub HLT-stable-20070905.141318 (rev 2375)    */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTHOMERWriter.h"
#include <sys/time.h>
#include <time.h>


HOMERWriter::HOMERWriter()
  :
  fDataOffset(0),
  fBlocks()
    {
    Clear();
    }

HOMERWriter::~HOMERWriter()
    {
    }

void HOMERWriter::Clear()
    {
    fDataOffset = 0;
    fBlocks.clear();
    }

void HOMERWriter::AddBlock( const void* descriptor, const void* data )
    {
    TBlockData bd;
    memcpy( bd.fDescriptor, descriptor, HOMERBlockDescriptor::GetHOMERBlockDescriptorSize() );
    bd.fData = data;
    HOMERBlockDescriptor hbd( &bd.fDescriptor );
    hbd.SetBlockOffset( fDataOffset );
    fDataOffset += hbd.GetBlockSize();
    fBlocks.push_back( bd );
    }

homer_uint32 HOMERWriter::GetTotalMemorySize( bool includeData )
    {
    if ( includeData )
	return fDataOffset + HOMERBlockDescriptor::GetHOMERBlockDescriptorSize()*(fBlocks.size()+1);
    else
	return HOMERBlockDescriptor::GetHOMERBlockDescriptorSize()*(fBlocks.size()+1);
    }

void HOMERWriter::Copy( void* destination, homer_uint64 eventType, homer_uint64 eventNr, homer_uint64 statusFlags, homer_uint64 nodeID, bool includeData )
    {
    HOMERBlockDescriptor homerBlock;
    homer_uint8* bd = reinterpret_cast<homer_uint8*>( destination );
    struct timeval now;
    gettimeofday( &now, NULL );
    homerBlock.UseHeader( bd );
    homerBlock.Initialize();
    homerBlock.SetUInt64Alignment( HOMERWriter::DetermineUInt64Alignment() );
    homerBlock.SetUInt32Alignment( HOMERWriter::DetermineUInt32Alignment() );
    homerBlock.SetUInt16Alignment( HOMERWriter::DetermineUInt16Alignment() );
    homerBlock.SetUInt8Alignment( HOMERWriter::DetermineUInt8Alignment() );
    homerBlock.SetDoubleAlignment( HOMERWriter::DetermineDoubleAlignment() );
    homerBlock.SetFloatAlignment( HOMERWriter::DetermineFloatAlignment() );
    homerBlock.SetType( eventType );
    homerBlock.SetSubType1( eventNr );
    homerBlock.SetSubType2( fBlocks.size() );
    homerBlock.SetBirth_s( now.tv_sec );
    homerBlock.SetBirth_us( now.tv_usec );
    homerBlock.SetProducerNode( nodeID );
    homerBlock.SetBlockOffset( homerBlock.GetHeaderLength() );
    homerBlock.SetBlockSize( HOMERBlockDescriptor::GetHOMERBlockDescriptorSize()*fBlocks.size() );
    homerBlock.SetStatusFlags( statusFlags );
    bd += homerBlock.GetHeaderLength();

    //unsigned long dataOffset = HOMERBlockDescriptor::GetHOMERBlockDescriptorSize()*(fBlocks.size()+1);
    unsigned long dataOffset = homerBlock.GetBlockOffset() + homerBlock.GetBlockSize();
    std::vector<TBlockData>::iterator iter, end;
    iter = fBlocks.begin();
    end = fBlocks.end();
    while ( iter != end )
	{
	homerBlock.UseHeader( iter->fDescriptor );
	homerBlock.SetBlockOffset( homerBlock.GetBlockOffset()+dataOffset );
	memcpy( bd, iter->fDescriptor, homerBlock.GetHeaderLength() );
	bd += homerBlock.GetHeaderLength();
	if ( includeData )
	    {
	    memcpy( ((homer_uint8*)destination)+homerBlock.GetBlockOffset(), iter->fData, homerBlock.GetBlockSize() );
	    }
	iter++;
	}
    }
	



    struct HOMERWriterAlignment64TestStructure
    {
	homer_uint64 f64Fill;
	homer_uint64 f64Test64;
	homer_uint32 f32Fill;
	homer_uint64 f64Test32;
	homer_uint16 f16Fill;
	homer_uint64 f64Test16;
	homer_uint8  f8Fill;
	homer_uint64 f64Test8;
    };
    struct HOMERWriterAlignment32TestStructure
    {
	homer_uint64 f64Fill;
	homer_uint32 f32Test64;
	homer_uint32 f32Fill;
	homer_uint32 f32Test32;
	homer_uint16 f16Fill;
	homer_uint32 f32Test16;
	homer_uint8  f8Fill;
	homer_uint32 f32Test8;
    };
    struct HOMERWriterAlignment16TestStructure
    {
	homer_uint64 f64Fill;
	homer_uint16 f16Test64;
	homer_uint32 f32Fill;
	homer_uint16 f16Test32;
	homer_uint16 f16Fill;
	homer_uint16 f16Test16;
	homer_uint8  f8Fill;
	homer_uint16 f16Test8;
    };
    struct HOMERWriterAlignment8TestStructure
    {
	homer_uint64 f64Fill;
	homer_uint8 f8Test64;
	homer_uint32 f32Fill;
	homer_uint8 f8Test32;
	homer_uint16 f16Fill;
	homer_uint8 f8Test16;
	homer_uint8  f8Fill;
	homer_uint8 f8Test8;
    };
    struct HOMERWriterAlignmentDoubleTestStructure
    {
	homer_uint64 f64Fill;
	double fDoubleTest64;
	homer_uint32 f32Fill;
	double fDoubleTest32;
	homer_uint16 f16Fill;
	double fDoubleTest16;
	homer_uint8  f8Fill;
	double fDoubleTest8;
    };
    struct HOMERWriterAlignmentFloatTestStructure
    {
	homer_uint64 f64Fill;
	float fFloatTest64;
	homer_uint32 f32Fill;
	float fFloatTest32;
	homer_uint16 f16Fill;
	float fFloatTest16;
	homer_uint8  f8Fill;
	float fFloatTest8;
    };

homer_uint8 HOMERWriter::DetermineUInt64Alignment()
    {
    HOMERWriterAlignment64TestStructure test;
    if ( (unsigned long)(&test.f64Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f64Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 64 bit element does not immedately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f64Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 64 bit element does not immedately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f64Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 64 bit element does not immedately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 HOMERWriter::DetermineUInt32Alignment()
    {
    HOMERWriterAlignment32TestStructure test;
    if ( (unsigned long)(&test.f32Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f32Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 32 bit element does not immedately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f32Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 32 bit element does not immedately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f32Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 32 bit element does not immedately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 HOMERWriter::DetermineUInt16Alignment()
    {
    HOMERWriterAlignment16TestStructure test;
    if ( (unsigned long)(&test.f16Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f16Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 16 bit element does not immedately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f16Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 16 bit element does not immedately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f16Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 16 bit element does not immedately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 HOMERWriter::DetermineUInt8Alignment()
    {
    HOMERWriterAlignment8TestStructure test;
    if ( (unsigned long)(&test.f8Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f8Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 8 bit element does not immedately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f8Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 8 bit element does not immedately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f8Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 8 bit element does not immedately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 HOMERWriter::DetermineDoubleAlignment()
    {
    HOMERWriterAlignmentDoubleTestStructure test;
    if ( (unsigned long)(&test.fDoubleTest64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.fDoubleTest32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The double element does not immedately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.fDoubleTest16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The double element does not immedately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.fDoubleTest8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The double element does not immedately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 HOMERWriter::DetermineFloatAlignment()
    {
    HOMERWriterAlignmentFloatTestStructure test;
    if ( (unsigned long)(&test.fFloatTest64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.fFloatTest32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The float element does not immedately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.fFloatTest16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The float element does not immedately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.fFloatTest8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The float element does not immedately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }



/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/
