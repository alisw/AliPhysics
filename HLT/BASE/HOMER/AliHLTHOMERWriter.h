// XEMacs -*-C++-*-
#ifndef ALIHLTHOMERWRITER_H
#define ALIHLTHOMERWRITER_H
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



class AliHLTHOMERWriter
    {
    public:

	AliHLTHOMERWriter();
	virtual ~AliHLTHOMERWriter();

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


        struct HOMERWriterAlignment64TestStructure
        {
        	homer_uint64 f64Fill;   // !
        	homer_uint64 f64Test64; // !
        	homer_uint32 f32Fill;   // !
        	homer_uint64 f64Test32; // !
        	homer_uint16 f16Fill;   // !
        	homer_uint64 f64Test16; // !
        	homer_uint8  f8Fill;    // !
        	homer_uint64 f64Test8;  // !
        };
        struct HOMERWriterAlignment32TestStructure
        {
        	homer_uint64 f64Fill;   // !
        	homer_uint32 f32Test64; // !
        	homer_uint32 f32Fill;   // !
        	homer_uint32 f32Test32; // !
        	homer_uint16 f16Fill;   // !
        	homer_uint32 f32Test16; // !
        	homer_uint8  f8Fill;    // !
        	homer_uint32 f32Test8;  // !
        };
        struct HOMERWriterAlignment16TestStructure
        {
        	homer_uint64 f64Fill;   // !
            	homer_uint16 f16Test64; // !
        	homer_uint32 f32Fill;   // !
        	homer_uint16 f16Test32; // !
        	homer_uint16 f16Fill;   // !
        	homer_uint16 f16Test16; // !
        	homer_uint8  f8Fill;    // !
        	homer_uint16 f16Test8;  // !
        };
        struct HOMERWriterAlignment8TestStructure
        {
        	homer_uint64 f64Fill; // !
        	homer_uint8 f8Test64; // !
        	homer_uint32 f32Fill; // !
        	homer_uint8 f8Test32; // !
        	homer_uint16 f16Fill; // !
        	homer_uint8 f8Test16; // !
        	homer_uint8  f8Fill;  // !
        	homer_uint8 f8Test8;  // !
        };
        struct HOMERWriterAlignmentDoubleTestStructure
        {
        	homer_uint64 f64Fill; // !
        	double fDoubleTest64; // !
        	homer_uint32 f32Fill; // !
        	double fDoubleTest32; // !
        	homer_uint16 f16Fill; // !
        	double fDoubleTest16; // !
        	homer_uint8  f8Fill;  // !
        	double fDoubleTest8;  // !
        };
        struct HOMERWriterAlignmentFloatTestStructure
        {
        	homer_uint64 f64Fill; // !
        	float fFloatTest64;   // !
        	homer_uint32 f32Fill; // !
        	float fFloatTest32;   // !
        	homer_uint16 f16Fill; // !
        	float fFloatTest16;   // !
        	homer_uint8  f8Fill;  // !
        	float fFloatTest8;    // !
        };
    protected:



	struct TBlockData
	    {
	      homer_uint64 fDescriptor[kCount_64b_Words]; //!transient
	      const void* fData; //!transient
	    };

        unsigned long fDataOffset; //!transient

        std::vector<TBlockData> fBlocks; //!transient
#ifdef USE_ROOT
      ClassDef(AliHLTHOMERWriter,0);
#endif
    };


/** defined for backward compatibility */
typedef AliHLTHOMERWriter HOMERWriter;

// external interface of the HOMER writer
#define ALIHLTHOMERWRITER_CREATE "AliHLTHOMERWriterCreate"
#define ALIHLTHOMERWRITER_DELETE "AliHLTHOMERWriterDelete"

#ifdef __cplusplus
extern "C" {
#endif

  typedef AliHLTHOMERWriter* (*AliHLTHOMERWriterCreate_t)(const void* pBuffer, int size);
  typedef void (*AliHLTHOMERWriterDelete_t)(AliHLTHOMERWriter* pInstance);
  /**
   * Create instance of HOMER writer.
   */
  AliHLTHOMERWriter* AliHLTHOMERWriterCreate();

  /**
   * Delete instance of HOMER writer.
   */
  void AliHLTHOMERWriterDelete(AliHLTHOMERWriter* pInstance);
#ifdef __cplusplus
}
#endif



/*
***************************************************************************
**
** $Author$ - Initial Version by Timm Morten Steinbeck
**
** $Id$ 
**
***************************************************************************
*/

#endif // ALIHLTHOMERWRITER_H
