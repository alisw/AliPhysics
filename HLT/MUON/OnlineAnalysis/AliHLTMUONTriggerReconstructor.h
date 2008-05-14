#ifndef AliHLTMUONTRIGGERRECONSTRUCTOR_H
#define AliHLTMUONTRIGGERRECONSTRUCTOR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/**********************************************************************
 Created on : 16/05/2007
 Purpose    : This class is supposed to read the trigger DDL files and 
              give the output AliHLTMUONTriggerRecordStruct
 Author     : Indranil Das, HEP Division, SINP
 Email      : indra.das@saha.ac.in | indra.ehep@gmail.com

 Artur Szostak <artursz@iafrica.com>:
  Completely reimplemented the lookup table to a simplified format.
**********************************************************************/

///
///  @file   AliHLTMUONTriggerReconstructor.h
///  @author Indranil Das <indra.das@saha.ac.in>,
///          Artur Szostak <artursz@iafrica.com>
///  @date   16 May 2007
///  @brief  Declaration of the AliHLTMUONTriggerReconstructor class for processing trigger DDL data.
///

#include "AliHLTLogging.h"
#include "AliHLTMUONDataTypes.h"
#include "AliMUONTriggerDDLDecoder.h"
#include "AliMUONTriggerDDLDecoderEventHandler.h"

extern "C" struct AliHLTMUONTriggerRecordStruct;


class AliHLTMUONTriggerReconstructor : public AliHLTLogging
{
public:

	AliHLTMUONTriggerReconstructor();
	virtual ~AliHLTMUONTriggerReconstructor();

	bool Run(
			const AliHLTUInt8_t* rawData,
			AliHLTUInt32_t rawDataSize,
			AliHLTMUONTriggerRecordStruct* trigRecord,
			AliHLTUInt32_t& nofTrigRec,
			bool suppressPartialTrigs = false
		);

	/**
	 * Returns a pointer to the lookup table used to map between channel
	 * addresses and geometrical strip positions.
	 */
	AliHLTMUONTriggerRecoLookupTable* LookupTableBuffer() { return fDecoder.GetHandler().LookupTableBuffer(); }
	
	/// Returns the size of the lookup table.
	size_t LookupTableSize() const { return fDecoder.GetHandler().LookupTableSize(); }

private:

	class AliDecoderHandler : public AliMUONTriggerDDLDecoderEventHandler, public AliHLTLogging
	{
	public:
	
		AliDecoderHandler();
		
		/// Default destructor.
		virtual ~AliDecoderHandler() {}
		
		/// Returns a pointer to the lookup table.
		AliHLTMUONTriggerRecoLookupTable* LookupTableBuffer() { return &fLookupTable; }
		
		/// Returns the size of the lookup table.
		size_t LookupTableSize() const { return sizeof(fLookupTable); }
		
		/**
		 * Returns the maximum number of trigger records that can be
		 * written to the output buffer.
		 */
		AliHLTUInt32_t MaxOutputTrigRecs() const { return fMaxOutputTrigRecs; }
		
		/**
		 * Sets the maximum number of trigger records that can be
		 * written to the output buffer.
		 */
		void MaxOutputTrigRecs(AliHLTUInt32_t n) { fMaxOutputTrigRecs = n; }
		
		/**
		 * Returns the number of reconstructed trigger records actually
		 * stored in the output buffer.
		 */
		AliHLTUInt32_t OutputTrigRecsCount() const { return fOutputTrigRecsCount; }
		
		/**
		 * Returns the pointer to the output buffer which stores reconstructed
		 * trigger records.
		 */
		AliHLTMUONTriggerRecordStruct* OutputTrigRecs() const { return fOutputTrigRecs; }
		
		/**
		 * Sets the pointer to the output buffer which stores reconstructed
		 * trigger records. Also resets the number of trigger records stored.
		 */
		void OutputTrigRecs(AliHLTMUONTriggerRecordStruct* buf)
		{
			fOutputTrigRecs = buf;
			fOutputTrigRecsCount = 0;
			fOverflowed = false;
		}
		
		/**
		 * Returns the flag specifying if we should suppress partial trigger records.
		 */
		bool SuppressPartialTriggers() const { return fSuppressPartialTriggers; }
		
		/**
		 * Sets the flag specifying if we should suppress partial trigger records.
		 */
		void SuppressPartialTriggers(bool s) { fSuppressPartialTriggers = s; }
		
		/**
		 * Returns true if the output buffer was overflowed.
		 */
		bool OverflowedOutputBuffer() const { return fOverflowed; }
		
		
		// Methods inherited from AliMUONTriggerDDLDecoderEventHandler:
		
		/// Called for each new buffer.
		void OnNewBuffer(const void* buffer, UInt_t /*bufferSize*/);
		
		/// Increments a counter for every new regional structure.
		void OnNewRegionalStruct(
				const AliMUONRegionalHeaderStruct* /*regionalStruct*/,
				const AliMUONRegionalScalarsStruct* /*scalars*/,
				const void* /*data*/
			)
		{
			fCurrentRegional++;
		}
		
		/// Reset the counter for every new regional structure.
		void OnEndOfRegionalStruct(
				const AliMUONRegionalHeaderStruct* /*regionalStruct*/,
				const AliMUONRegionalScalarsStruct* /*scalars*/,
				const void* /*data*/
			)
		{
			// Start from -1 since we increment immediately in OnLocalStruct.
			fCurrentLocal = -1;
		}
		
		/// Converts a local trigger structure from the L0 into a trigger record.
		void OnLocalStruct(
				const AliMUONLocalInfoStruct* localStruct,
				const AliMUONLocalScalarsStruct* scalars
			);
		
		/// Logs an error message if there was a decoding problem with the DDL payload.
		void OnError(ErrorCode code, const void* location);
	
	private:
		// Do not allow copying of this class.
		/// Not implemented
		AliDecoderHandler(const AliDecoderHandler& rhs); // copy constructor
		/// Not implemented
		AliDecoderHandler& operator = (const AliDecoderHandler& rhs); // assignment operator
		
		AliHLTMUONTriggerRecoLookupTable fLookupTable;  ///< The lookup table used for mapping between channel addresses and geometrical information.
		const void* fBufferStart; ///< Pointer to the start of the current DDL payload buffer.
		AliHLTUInt32_t fMaxOutputTrigRecs;  ///< Maximum number of reconstructed trigger records that can be stored in fOutputTrigRecs.
		AliHLTUInt32_t fOutputTrigRecsCount;  ///< The number of reconstructed trigger records actually stored in fOutputTrigRecs.
		AliHLTMUONTriggerRecordStruct* fOutputTrigRecs;  ///< Pointer to the output buffer of trigger records structures.
		AliHLTInt32_t fTrigRecId;  ///< A running counter for the trigger record ID.
		AliHLTInt32_t fCurrentRegional;  ///< The current regional trigger structure number.
		AliHLTInt32_t fCurrentLocal;  ///< The current local trigger structure number.
		bool fSuppressPartialTriggers;  ///< Flag to indicate if we should suppres partial triggers.
		bool fOverflowed;  ///< Flag to indicate if we overflowed the output buffer.
	};

	/// Not implemented
	AliHLTMUONTriggerReconstructor(const AliHLTMUONTriggerReconstructor& rhs); // copy constructor
	/// Not implemented
	AliHLTMUONTriggerReconstructor& operator = (const AliHLTMUONTriggerReconstructor& rhs); // assignment operator

	AliMUONTriggerDDLDecoder<AliDecoderHandler> fDecoder; ///< Raw DDL data decoder.
};

#endif // AliHLTMUONTRIGGERRECONSTRUCTOR_H
