 // -*- Mode: C++ -*-

#ifndef ALIHLTTPCHWCLUSTERDECODERCOMPONENT_H
#define ALIHLTTPCHWCLUSTERDECODERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

//  @file   AliHLTTPCHWClusterDecoderComponent.h
//
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Matthias Richter <matthias.richter@cern.ch>
//  @author Torsten Alt <talt@cern.ch>
// 
//  @brief  This HLT Component converts output of FPGA clusterfinder 
//  @brief  to AliHLTTPCRawClusterData format. It also merges splitted clusters 
//  @brief  at branch borders and between the patches.
//  @note

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"

class AliHLTTPCHWCFData;
class AliHLTTPCHWClusterMergerV1;

/**
 * @class AliHLTTPCHWClusterDecoderComponent
 *
 * The format of the incoming HW clusters is as follows:
 * 
 * WORD 0 : CDH0
 * WORD 1 : CDH1
 * WORD 2 : CDH2
 * WORD 3 : CDH3
 * WORD 4 : CDH4   8 32-bit words for the common data header
 * WORD 5 : CDH5
 * WORD 6 : CDH6
 * WORD 7 : CDH7
 *
 * WORD 8  : contains cluster identification (bits 30-31), Row (6 bits (bit 24-29)) and Charge (24 bits (bit 0-23))
 * WORD 9  : contain the pad (32 bits)
 * WORD 10 : contain the time (32 bits)
 * WORD 11 : contain the pad error (32 bits)
 * WORD 12 : contain the time error (32 bits)
 *
 * WORD 13 : contains cluster identification (bits 30-31), Row (6 bits (bit 24-29)) and Charge (24 bits (bit 0-23))
 * WORD 14 : contains the pad (32 bits)
 * WORD 15 : contains the time (32 bits)
 * WORD 16 : contains the pad error (32 bits)
 * WORD 17 : contains the time error (32 bits)
  
 * WORD 18 : RCU TRAILER WORD 0
 * WORD 19 : RCU TRAILER WORD 1
 * ...
 * WORD 18+N : RCU TRAILER WORD N
 *
 * The cluster is signified with bits 30 and 31 being 11 (0x3).
 * The RCU trailer is signified with bits 30 and 31 being 10 (0x2).
 *
 * A buffer reads the incoming data block. We skip the first 8 words of the CDH.
 * Then we shift 30 bits to the right to get the last 2 bits of the next word, 30 and 31.
 * 
 * If bit3031 = 0x3 (11), we are in the beginning of a cluster. 
 * We apply a 24-bit mask to get the first 24 bits that represent the charge of the cluster.
 * We cast the word buffer to an 8-bit pointer that will point at bit 0 (for little endian),
 * then increment it by 3 bytes, which takes us to bit 24. The last byte of the word
 * contains row info and the cluster tag (11).
 * 
 * With a & bit operation the row info is retrieved (6 bits). 
 * The rest of the cluster information is contained in the following 4 words.
 *
 * If bit3031 = 0x2 (10), we have reached the RCU trailer. 
 *
 * -------------------------------------------------------
 *
 * Short note about little endian:
 *
 * DE AD BE EF (a 4 byte word) 
 * 
 * The byte order in memory to represent this word will be:
 * 
 * EF
 * BE
 * AD
 * DE
 * 
 * The least significant byte value is stored at the lowest address.
 * http://en.wikipedia.org/wiki/Endianness#Little-endian
 * 
 * @ingroup alihlt_tpc_components
 */

class AliHLTTPCHWClusterDecoderComponent : public AliHLTProcessor {
    
public:

  /** standard constructor */    
  AliHLTTPCHWClusterDecoderComponent();           
  /** destructor */
  virtual ~AliHLTTPCHWClusterDecoderComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process
      
  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();							     
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);			     
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();					     
  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);			   
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ); 
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn(); 							   
  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);

  void PrintDebug(AliHLTUInt32_t * buffer, Int_t size);
  
protected:
	
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing capabilities of the component. 

  int DoInit( int argc, const char** argv );
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		     AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		     AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  int DoDeinit();
  
  int Reconfigure(const char* cdbEntry, const char* chainId);  

  using AliHLTProcessor::DoEvent;
  
private:
   
  int ScanConfigurationArgument(int argc, const char** argv);
          
  /** copy constructor prohibited */
  AliHLTTPCHWClusterDecoderComponent(const AliHLTTPCHWClusterDecoderComponent&);

  /** assignment operator prohibited */
  AliHLTTPCHWClusterDecoderComponent& operator=(const AliHLTTPCHWClusterDecoderComponent&);

  int InitClusterMerger();

  static const char* fgkOCDBEntry; //!transient

  AliHLTTPCHWCFData* fpDecoder; // decoder
  AliHLTTPCHWClusterMergerV1* fpClusterMerger; //! merger instance
  Bool_t fDoMerge; // flag whether to merge clusters at branch borders
  Bool_t fAlreadyMerged; // flag whether the incoming clusters are already merged at branch borders
  bool fTPCPresent; //flag whether TPC is present in detector list
  Bool_t fProcessingRCU2Data; // processing of RCU2 data - clusters are not split into two input branches, only merge at patch borders
  Bool_t fCheckEdgeFlag; //Consider only clusters with edge flag set
  int fAddRandomClusters; // Add random clusters for stress tests
  int fSignificantBitsCharge; //Number of significant bits to be stored for cluster charge
  int fSignificantBitsWidth; //Same for sigma

  AliHLTComponentBenchmark fBenchmark; // benchmarks
};

#endif
