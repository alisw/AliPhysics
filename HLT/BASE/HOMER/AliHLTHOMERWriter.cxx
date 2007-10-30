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

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTHOMERWriter.h"
#include <sys/time.h>
#include <time.h>


AliHLTHOMERWriter::AliHLTHOMERWriter()
  :
  fDataOffset(0),
  fBlocks()
    {
    // Writer implementation of the HOMER interface.
    // The HLT Monitoring Environment including ROOT is
    // a native interface to ship out data from the HLT chain.
    // See pdf document shiped with the package
    // for class documentation and tutorial.
    Clear();
    }

AliHLTHOMERWriter::~AliHLTHOMERWriter()
    {
    // see header file for class documentation
    }

void AliHLTHOMERWriter::Clear()
    {
    // see header file for class documentation
    fDataOffset = 0;
    fBlocks.clear();
    }

void AliHLTHOMERWriter::AddBlock( const void* descriptor, const void* data )
    {
    // see header file for class documentation
    TBlockData bd;
    memcpy( bd.fDescriptor, descriptor, HOMERBlockDescriptor::GetHOMERBlockDescriptorSize() );
    bd.fData = data;
    HOMERBlockDescriptor hbd( &bd.fDescriptor );
    hbd.SetBlockOffset( fDataOffset );
    fDataOffset += hbd.GetBlockSize();
    fBlocks.push_back( bd );
    }

homer_uint32 AliHLTHOMERWriter::GetTotalMemorySize( bool includeData )
    {
    // see header file for class documentation
    if ( includeData )
	return fDataOffset + HOMERBlockDescriptor::GetHOMERBlockDescriptorSize()*(fBlocks.size()+1);
    else
	return HOMERBlockDescriptor::GetHOMERBlockDescriptorSize()*(fBlocks.size()+1);
    }

void AliHLTHOMERWriter::Copy( void* destination, homer_uint64 eventType, homer_uint64 eventNr, homer_uint64 statusFlags, homer_uint64 nodeID, bool includeData )
    {
    // see header file for class documentation
    HOMERBlockDescriptor homerBlock;
    homer_uint8* bd = reinterpret_cast<homer_uint8*>( destination );
    struct timeval now;
    gettimeofday( &now, NULL );
    homerBlock.UseHeader( bd );
    homerBlock.Initialize();
    homerBlock.SetUInt64Alignment( AliHLTHOMERWriter::DetermineUInt64Alignment() );
    homerBlock.SetUInt32Alignment( AliHLTHOMERWriter::DetermineUInt32Alignment() );
    homerBlock.SetUInt16Alignment( AliHLTHOMERWriter::DetermineUInt16Alignment() );
    homerBlock.SetUInt8Alignment( AliHLTHOMERWriter::DetermineUInt8Alignment() );
    homerBlock.SetDoubleAlignment( AliHLTHOMERWriter::DetermineDoubleAlignment() );
    homerBlock.SetFloatAlignment( AliHLTHOMERWriter::DetermineFloatAlignment() );
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

homer_uint8 AliHLTHOMERWriter::DetermineUInt64Alignment()
    {
    // see header file for class documentation
    HOMERWriterAlignment64TestStructure test;
    if ( (unsigned long)(&test.f64Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f64Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 64 bit element does not immediately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f64Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 64 bit element does not immediately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f64Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 64 bit element does not immediately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 AliHLTHOMERWriter::DetermineUInt32Alignment()
    {
    // see header file for class documentation
    HOMERWriterAlignment32TestStructure test;
    if ( (unsigned long)(&test.f32Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f32Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 32 bit element does not immediately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f32Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 32 bit element does not immediately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f32Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 32 bit element does not immediately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 AliHLTHOMERWriter::DetermineUInt16Alignment()
    {
    // see header file for class documentation
    HOMERWriterAlignment16TestStructure test;
    if ( (unsigned long)(&test.f16Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f16Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 16 bit element does not immediately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f16Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 16 bit element does not immediately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f16Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 16 bit element does not immediately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 AliHLTHOMERWriter::DetermineUInt8Alignment()
    {
    // see header file for class documentation
    HOMERWriterAlignment8TestStructure test;
    if ( (unsigned long)(&test.f8Test64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.f8Test32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The 8 bit element does not immediately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.f8Test16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The 8 bit element does not immediately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.f8Test8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The 8 bit element does not immediately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 AliHLTHOMERWriter::DetermineDoubleAlignment()
    {
    // see header file for class documentation
    HOMERWriterAlignmentDoubleTestStructure test;
    if ( (unsigned long)(&test.fDoubleTest64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.fDoubleTest32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The double element does not immediately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.fDoubleTest16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The double element does not immediately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.fDoubleTest8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The double element does not immediately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

homer_uint8 AliHLTHOMERWriter::DetermineFloatAlignment()
    {
    // see header file for class documentation
    HOMERWriterAlignmentFloatTestStructure test;
    if ( (unsigned long)(&test.fFloatTest64) != ((unsigned long)(&test.f64Fill))+sizeof(test.f64Fill) )
	{
	// Alignment is beyond 64 bit, this is, to the best of my knowledge, currently unheard of.
	return ~(homer_uint8)0;
	}
    if ( (unsigned long)(&test.fFloatTest32) != ((unsigned long)(&test.f32Fill))+sizeof(test.f32Fill) )
	{
	// The float element does not immediately follow the 32 bit element, 
	// therefore the alignment has to be greater than 4.
	return (homer_uint8)8;
	}
    if ( (unsigned long)(&test.fFloatTest16) != ((unsigned long)(&test.f16Fill))+sizeof(test.f16Fill) )
	{
	// The float element does not immediately follow the 16 bit element, 
	// therefore the alignment has to be greater than 2.
	return (homer_uint8)4;
	}
    if ( (unsigned long)(&test.fFloatTest8) != ((unsigned long)(&test.f8Fill))+sizeof(test.f8Fill) )
	{
	// The float element does not immediately follow the 8 bit element, 
	// therefore the alignment has to be greater than 1.
	return (homer_uint8)2;
	}
    return 1;
    }

AliHLTHOMERWriter* AliHLTHOMERWriterCreate()
    {
    // see header file for function documentation
    return new AliHLTHOMERWriter();
    }

void AliHLTHOMERWriterDelete(AliHLTHOMERWriter* pInstance)
    {
    // see header file for function documentation
    if (pInstance) delete pInstance;
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
