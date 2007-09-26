// XEMacs -*-C++-*-
#ifndef _HOMERWRITER_HPP_
#define _HOMERWRITER_HPP_
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

/** @file   AliHLTHOMERWriter.h
    @author Timm Steinbeck
    @date   Sep 14 2007
    @brief  HLT Online Monitoring Environment including ROOT - Writer   
    @note   migrated from PubSub HLT-stable-20070905.141318 (rev 2375)    */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTHOMERData.h"
#include <vector>



class HOMERWriter
    {
    public:

	HOMERWriter();
	virtual ~HOMERWriter();

	void Clear();

	void AddBlock( const void* descriptor, const void* data );
	void AddBlock( const HOMERBlockDescriptor* descriptor, const void* data )
		{
		AddBlock( descriptor->GetHeader(), data );
		}

	homer_uint32 GetTotalMemorySize( bool includeData = true );
	void Copy( void* destination, homer_uint64 eventType, homer_uint64 eventNr, homer_uint64 statusFlags, homer_uint64 nodeID, bool includeData = true );

	static homer_uint8 DetermineUInt64Alignment();
	static homer_uint8 DetermineUInt32Alignment();
	static homer_uint8 DetermineUInt16Alignment();
	static homer_uint8 DetermineUInt8Alignment();
	static homer_uint8 DetermineDoubleAlignment();
	static homer_uint8 DetermineFloatAlignment();


    protected:



	struct TBlockData
	    {
	      homer_uint64 fDescriptor[kCount_64b_Words]; //!transient
	      const void* fData; //!transient
	    };

        unsigned long fDataOffset; //!transient

        std::vector<TBlockData> fBlocks; //!transient
#ifdef USE_ROOT
      ClassDef(HOMERWriter,0);
#endif
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

#endif // _HOMERWRITER_HPP_
