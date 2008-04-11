// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPDUMPCOMPONENT_H
#define ALIHLTTPCCOMPDUMPCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCompDumpComponent.h
    @author Timm Steinbeck
    @date   
    @brief  Declaration of a copy component. */


#include "AliHLTProcessor.h"

/**
 * @class AliHLTTPCCompDumpComponent
 * @author Timm Steinbeck
 * @brief A dummy HLT processing component. 
 * @date 05-03-2008
 *
 * An implementiation of a copy component that just copies its input data
 * to debug a components input data
 * @ingroup alihlt_tpc
 */
class AliHLTTPCCompDumpComponent : public AliHLTProcessor
    {
    public:
      
      /** standard constructor */
      AliHLTTPCCompDumpComponent();

      /** standard destructor */
      virtual ~AliHLTTPCCompDumpComponent();

	// Public functions to implement AliHLTComponent's interface.
	// These functions are required for the registration process

      /** function to get component id 
       * @return const char* pointer to componentid
       */
      const char* GetComponentID();

      /** function to get input data types
       * @param list vecotr of AliHLTComponent_DataType
       */
      void GetInputDataTypes( vector<AliHLTComponent_DataType>& list);

      /** function to get output data type
       * @return AliHLTComponent_DataType
       */
      AliHLTComponent_DataType GetOutputDataType();

      /** function to get output data size
       * @param constBase address of an unsigned long
       * @param inputMultiplier address of a double
       */
      virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

      /** spawn function
       * @return AliHLTComponent* pointer to instance
       */
      AliHLTComponent* Spawn();
	
    protected:
	
	// Protected functions to implement AliHLTComponent's interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component. 

      /** initialisation function
       * @param argc integer counting number of input arguments
       * @param argv const char** for parameter values
       * @return zero upon success
       */
	int DoInit( int argc, const char** argv );

      /** deinitialisation function
       * @return zero upon success
       */
	int DoDeinit();

      /** do event function
       * @param evtData      const AliHLTComponent_EventData& to event data
       * @param blocks       const AliHLTComponent_BlockData* to blocks of event data
       * @param trigData     AliHLTComponent_TriggerData& of trigger data
       * @param outputPtr    AliHLTUInt8_t* pointer to output data
       * @param size         AliHLTUInt32_t& of output size
       * @param outputBlocks vector<AliHLTComponent_BlockData>& of output block data
       * @return zero upon success
      */
      int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
		   AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		   AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks );
      
      
      /** function to initialise bit input
       * @param input     AliHLTUInt8_t
       * @param inputSize UInt_t
       */
      void InitBitDataInput( AliHLTUInt8_t* input, UInt_t inputSize );

      /** function to close bit data input */	
      void CloseBitDataInput()
      {
      }

      /** function to get current byte input position
       * @return unsigned long value of current byte input position
       */
      unsigned long GetCurrentByteInputPosition() const
      {
	return (unsigned long)( fBitDataCurrentInput - fBitDataCurrentInputStart );
      }
      
      /** function to get current bit input position
       * @return fBitDataCurrentPosInWord unsigned value
       */
      unsigned GetCurrentBitInputPosition() const
      {
	return fBitDataCurrentPosInWord;
      }

      /** function to get current input byte
       * @return fBitDataCurrentWord AliHLTUInt8_t value
       */
      AliHLTUInt8_t GetCurrentInputByte() const
      {
	return fBitDataCurrentWord;
      }
      
      /** function to get end of bit input
       * @return boolean if input is at end or not
       */
      bool EndOfBitInput() const
      {
	return (fBitDataCurrentInput>=fBitDataCurrentInputEnd);
      }
      
      /** function got get bit data input size bytes
       * @return UInt_t value
       */
      UInt_t GetBitDataInputSizeBytes() const
      {
	return fBitDataCurrentInput-fBitDataCurrentInputStart;
      }
      
      /** function to determine input bit
       * @param value  AliHLTUInt8_t &
       * @return boolean (input bit = 1 or = 0)
       */
      bool InputBit( AliHLTUInt8_t & value );

      /** function to determine input bits below 8 bits
       * @param value     AliHLTUInt8_t &
       * @param bitCount  UInt_t const &
       * @return boolean
       */
      bool InputBits( AliHLTUInt8_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 8 and 16 bits
       * @param value     AliHLTUInt16_t &
       * @param bitCount  UInt_t const &
       * @return boolean
       */
      bool InputBits( AliHLTUInt16_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 16 and 32 bits
       * @param value     AliHLTUInt32_t &
       * @param bitCount  UInt_t const &
       * @return boolean
       */
      bool InputBits( AliHLTUInt32_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 16 and 32 bits II
       * @param value     Int_t &
       * @param bitCount  UInt_t const &
       * @return boolean
       */
      bool InputBits( Int_t & value, UInt_t const & bitCount );

      /** function to determine input bits between 32 and 64 bits
       * @param value     AliHLTUInt64_t &
       * @param bitCount  UInt_t const &
       * @return boolean
       */
      bool InputBits( AliHLTUInt64_t & value, UInt_t const & bitCount );

      /** pad function for 8 bits */
      void Pad8Bits();

      /** function for input bytes 
       * @param data      AliHLTUInt8_t*
       * @param byteCount UInt_t const &
       */
      bool InputBytes( AliHLTUInt8_t* data, UInt_t const & byteCount );

    private:
      /** copy constructor prohibited */
      AliHLTTPCCompDumpComponent(const AliHLTTPCCompDumpComponent&);
      /** assignment operator prohibited */
      AliHLTTPCCompDumpComponent& operator=(const AliHLTTPCCompDumpComponent&);
      
      /** member varibable for current bit data word */
      AliHLTUInt8_t fBitDataCurrentWord; // member varibable for current bit data word
      /** member variable for current bit data position in word */
      UInt_t fBitDataCurrentPosInWord; // member variable for current bit data position in word
      /** member variable for current input bit data */
      AliHLTUInt8_t *fBitDataCurrentInput; // member variable for current input bit data
      /** member variable for current bit data input start */
      AliHLTUInt8_t *fBitDataCurrentInputStart; // member variable for current bit data input start
      /** member variable for current bit data input end */
      AliHLTUInt8_t *fBitDataCurrentInputEnd; // member variable for current bit data input end

      ClassDef(AliHLTTPCCompDumpComponent, 0)

    };
#endif
