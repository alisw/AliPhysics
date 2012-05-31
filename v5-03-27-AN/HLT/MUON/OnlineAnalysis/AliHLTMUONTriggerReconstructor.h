#ifndef AliHLTMUONTRIGGERRECONSTRUCTOR_H
#define AliHLTMUONTRIGGERRECONSTRUCTOR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

// $Id$

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

extern "C" struct AliHLTMUONRecHitStruct;
extern "C" struct AliHLTMUONTriggerRecordStruct;
extern "C" struct AliHLTMUONTrigRecInfoStruct;
extern "C" struct AliHLTMUONTriggerRecoLutRow;


class AliHLTMUONTriggerReconstructor : public AliHLTLogging
{
public:

	AliHLTMUONTriggerReconstructor();
	virtual ~AliHLTMUONTriggerReconstructor();

	bool Run(
			const AliHLTUInt8_t* rawData,
			AliHLTUInt32_t rawDataSize,
			bool scalarEvent,
			AliHLTMUONTriggerRecordStruct* trigRecord,
			AliHLTUInt32_t& nofTrigRec
		);

	/**
	 * Returns a pointer to the lookup table used to map between channel
	 * addresses and geometrical strip positions.
	 */
	AliHLTMUONTriggerRecoLookupTable* LookupTableBuffer() { return fDecoder.GetHandler().LookupTableBuffer(); }
	
	/// Returns the size of the lookup table.
	size_t LookupTableSize() const { return fDecoder.GetHandler().LookupTableSize(); }
	
	/// Returns true if the decoder is set to enable recovery logic if
	/// raw data errors are found.
	bool TryRecover() const { return fDecoder.TryRecover(); };
	
	/// Sets if the decoder should enable the error recovery logic.
	void TryRecover(bool value);
	
	/**
	 * Sets the DDL bit according to the DDL value given.
	 * It is used to keep the trigger record IDs unique.
	 */
	void SetDDL(AliHLTInt32_t ddl) { fDecoder.GetHandler().SetDDL(ddl); }
	
	/**
	 * Returns true if the output buffer was overflowed in the last call to Run.
	 */
	bool OverflowedOutputBuffer() const { return fDecoder.GetHandler().OverflowedOutputBuffer(); }
	
	/**
	 * Returns the flag specifying if we should suppress partial trigger records.
	 * Partial triggers do not pass the 3/4'ths coincidence requirement.
	 */
	bool SuppressPartialTriggers() const { return fDecoder.GetHandler().SuppressPartialTriggers(); }
	
	/**
	 * Sets the flag specifying if we should suppress partial trigger records.
	 * Partial triggers do not pass the 3/4'ths coincidence requirement.
	 */
	void SuppressPartialTriggers(bool s) { fDecoder.GetHandler().SuppressPartialTriggers(s); }
	
	/**
	 * Returns true if the crate ID as found in the regional header
	 * will be used for lookups in the LUT, rather than the sequencial
	 * index number of the header.
	 */
	bool UseCrateId() const { return fDecoder.GetHandler().UseCrateId(); }
	
	/**
	 * Sets the flag indicating if the crate ID as found in the regional
	 * header should be used for lookups in the LUT, rather than the
	 * sequencial index number of the header.
	 */
	void UseCrateId(bool value) { fDecoder.GetHandler().UseCrateId(value); }
	
	/**
	 * Returns true if the local board ID as found in the local structure
	 * will be used for lookups in the LUT, rather than the sequencial
	 * index number of the structure.
	 */
	bool UseLocalId() const { return fDecoder.GetHandler().UseLocalId(); }
	
	/**
	 * Sets the flag indicating if the local board ID as found in the local
	 * structure should be used for lookups in the LUT, rather than the
	 * sequencial index number of the structure.
	 */
	void UseLocalId(bool value) { fDecoder.GetHandler().UseLocalId(value); }
	
	/// Return the flag indicating if the debug information is stored during decoding.
	bool StoreDebugInfo() const { return fDecoder.GetHandler().StoreDebugInfo(); }
	
	/// Sets the flag indicating if the debug information should be stored.
	void StoreDebugInfo(bool value) { fDecoder.GetHandler().StoreDebugInfo(value); }
	
	/// Returns the number of elements in the debug information buffer.
	AliHLTUInt32_t InfoBufferCount() const { return fDecoder.GetHandler().InfoBufferCount(); }
	
	/// Returns the debug information buffer.
	const AliHLTMUONTrigRecInfoStruct* InfoBuffer() const { return fDecoder.GetHandler().InfoBuffer(); }
	
	/// Empty the info buffer.
	void ZeroInfoBuffer() { fDecoder.GetHandler().ZeroInfoBuffer(); }
	
	/**
	 * Returns the flag indicating if the error message for a wrong event type found
	 * in the DARC header should be suppressed.
	 */
	bool DontPrintWrongEventError() const { return fDecoder.GetHandler().DontPrintWrongEventError(); }
	
	/**
	 * Sets the flag indicating if the error message for a wrong event type found
	 * in the DARC header should be suppressed.
	 */
	void DontPrintWrongEventError(bool value) { fDecoder.GetHandler().DontPrintWrongEventError(value); }
	
private:

	class AliDecoderHandler : public AliMUONTriggerDDLDecoderEventHandler, public AliHLTLogging
	{
	public:
	
		AliDecoderHandler();
		
		/// Default destructor.
		virtual ~AliDecoderHandler();
		
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
		
		/**
		 * Returns true if the OnError handler method will only generate warning
		 * messages and rather than error messages.
		 */
		bool WarnOnly() const { return fWarnOnly; }
		
		/**
		 * Sets the flag indicating if the OnError method should only generate
		 * warnings rather than error messages.
		 */
		void WarnOnly(bool value) { fWarnOnly = value; }
		
		/**
		 * Returns true if the crate ID as found in the regional header
		 * will be used for lookups in the LUT, rather than the sequencial
		 * index number of the header.
		 */
		bool UseCrateId() const { return fUseCrateId; }
		
		/**
		 * Sets the flag indicating if the crate ID as found in the regional
		 * header should be used for lookups in the LUT, rather than the
		 * sequencial index number of the header.
		 */
		void UseCrateId(bool value) { fUseCrateId = value; }
		
		/**
		 * Returns true if the local board ID as found in the local structure
		 * will be used for lookups in the LUT, rather than the sequencial
		 * index number of the structure.
		 */
		bool UseLocalId() const { return fUseLocalId; }
		
		/**
		 * Sets the flag indicating if the local board ID as found in the local
		 * structure should be used for lookups in the LUT, rather than the
		 * sequencial index number of the structure.
		 */
		void UseLocalId(bool value) { fUseLocalId = value; }
		
		/**
		 * Sets the DDL bit according to the DDL value given.
		 */
		void SetDDL(AliHLTInt32_t ddl) { fDDLBit = (ddl == 20 ? 0x00 : 0x80); }
		
		// Methods inherited from AliMUONTriggerDDLDecoderEventHandler:
		
		/// Called for each new buffer. Just remember the start location of the buffer.
		void OnNewBuffer(const void* buffer, UInt_t /*bufferSize*/)
		{
			assert( buffer != NULL );
			fBufferStart = buffer;
			fHadWrongEventTypeError = fHadNonWrongEventTypeError = false;
		}
		
		/**
		 * Sets the regional structure sequencial number and decodes the crate ID.
		 * Also zero the local structure pointers.
		 * If fUseCrateId is false then we use the sequencial number instead. This
		 * might be necessary for incorrectly generated or buggy raw data.
		 */
		void OnNewRegionalStructV2(
				UInt_t num,
				const AliMUONRegionalHeaderStruct* regionalStruct,
				const AliMUONRegionalScalarsStruct* /*scalars*/,
				const void* /*data*/
			);
		
		/**
		 * Updates the local trigger structure pointers and processes the
		 * the last local trigger.
		 */
		void OnEndOfRegionalStructV2(
				UInt_t num,
				const AliMUONRegionalHeaderStruct* regionalStruct,
				const AliMUONRegionalScalarsStruct* scalars,
				const void* data
			);
		
		/**
		 * Updates the local trigger structure pointers and processes the
		 * current local trigger.
		 */
		void OnLocalStructV2(
				UInt_t num,
				const AliMUONLocalInfoStruct* localStruct,
				const AliMUONLocalScalarsStruct* scalars
			);
		
		/// Logs an error message if there was a decoding problem with the DDL payload.
		void OnError(ErrorCode code, const void* location);
		
		/// Return the flag indicating if the debug information is stored during decoding.
		bool StoreDebugInfo() const { return fStoreInfo; }
		
		/// Sets the flag indicating if the debug information should be stored.
		void StoreDebugInfo(bool value) { fStoreInfo = value; }
		
		/// Returns the number of elements in the debug information buffer.
		AliHLTUInt32_t InfoBufferCount() const { return fInfoBufferCount; }
		
		/// Returns the debug information buffer.
		const AliHLTMUONTrigRecInfoStruct* InfoBuffer() const { return fInfoBuffer; }
		
		/// Empty the info buffer.
		void ZeroInfoBuffer() { fInfoBufferCount = 0; }
		
		/**
		 * Returns the flag indicating if the error message for a wrong event type found
		 * in the DARC header should be suppressed.
		 */
		bool DontPrintWrongEventError() const { return fDontPrintWrongEventError; }
		
		/**
		 * Sets the flag indicating if the error message for a wrong event type found
		 * in the DARC header should be suppressed.
		 */
		void DontPrintWrongEventError(bool value) { fDontPrintWrongEventError = value; }
		
		/// Returns true if the last decoded DDL had a wrong event type error in the DARC header.
		bool HadWrongEventTypeError() const { return fHadWrongEventTypeError; }
		
		/**
		 * Returns true if the last decoded DDL had a different error than just a
		 * wrong event type error in the DARC header.
		 */
		bool HadNonWrongEventTypeError() const { return fHadNonWrongEventTypeError; }
		
	private:
		// Do not allow copying of this class.
		/// Not implemented
		AliDecoderHandler(const AliDecoderHandler& rhs); // copy constructor
		/// Not implemented
		AliDecoderHandler& operator = (const AliDecoderHandler& rhs); // assignment operator
		
		/**
		 * Finds the strip bits / positions on MT1 that were fired in
		 * the current local trigger structure decision.
		 */
		bool FindStripsOnMT1(AliHLTInt32_t& xPos, AliHLTInt32_t& yPos);
		
		/**
		 * Selects the correct X strip patterns to use in FindXStrips.
		 * \param [out] strips  Resulting array of X strip patterns to use
		 *    for chambers 11 to 14.
		 */
		void SelectXPatterns(AliHLTUInt64_t strips[4]);
		
		/**
		 * Selects the correct Y strip patterns to use in FindYStrips and local IDs for
		 * finding the correct row in the lookup table.
		 * \param [in] xpos Array of X strip positions generated by FindXStrips.
		 *    Values are in the range [0..47].
		 * \param [out] strips  Resulting array of Y strip patterns to use.
		 * \param [out] locId  Resulting array of local IDs to use for the lookup table.
		 */
		void SelectYPatterns(AliHLTInt32_t xpos[4], AliHLTUInt32_t strips[4], AliHLTUInt8_t locId[4]);
		
		/**
		 * Tries to find the fired X strips for chambers 11 to 14.
		 */
		void FindXStrips(AliHLTInt32_t startPos, AliHLTUInt64_t strips[4], AliHLTInt32_t pos[4]);
		
		/**
		 * Tries to find the fired Y strips for chambers 11 to 14.
		 */
		void FindYStrips(AliHLTInt32_t startPos, AliHLTUInt32_t strips[4], AliHLTInt32_t pos[4]);
		
		/**
		 * Fetches the appropriate LUT row for a given X strip position.
		 */
		const AliHLTMUONTriggerRecoLutRow& GetLutRowX(AliHLTInt32_t xPos, AliHLTUInt8_t chamber);
		
		/**
		 * Reconstructs a hit with global position coordinates from strip
		 * information on a given chamber.
		 */
		void ReconstructHit(
				AliHLTUInt64_t xStrips, AliHLTUInt32_t yStrips,
				AliHLTInt32_t xPos, AliHLTInt32_t yPos, AliHLTUInt8_t yLocId,
				AliHLTUInt8_t chamber, AliHLTMUONRecHitStruct& hit
			);
		
		/**
		 * Converts the fCurrentStruct local trigger structure from the L0
		 * into a trigger record.
		 */
		void ProcessLocalStruct();
		
		AliHLTMUONTriggerRecoLookupTable fLookupTable;  ///< The lookup table used for mapping between channel addresses and geometrical information.
		const void* fBufferStart; ///< Pointer to the start of the current DDL payload buffer.
		AliHLTUInt32_t fMaxOutputTrigRecs;  ///< Maximum number of reconstructed trigger records that can be stored in fOutputTrigRecs.
		AliHLTUInt32_t fOutputTrigRecsCount;  ///< The number of reconstructed trigger records actually stored in fOutputTrigRecs.
		AliHLTMUONTriggerRecordStruct* fOutputTrigRecs;  ///< Pointer to the output buffer of trigger records structures.
		AliHLTInt32_t fTrigRecId;  ///< A running counter for the trigger record ID.
		AliHLTInt32_t fDDLBit;  ///< The DDL bit used to generate unique trigger record IDs.
		bool fSuppressPartialTriggers;  ///< Flag to indicate if we should suppres partial triggers.
		bool fOverflowed;  ///< Flag to indicate if we overflowed the output buffer.
		bool fWarnOnly;  ///< Flag indicating if the OnError method should generate warnings rather than error messages.
		bool fUseLocalId;  ///< Flag to indicate if the local structure ID as found in the local structures should be used or not.
		bool fUseCrateId;  ///< Flag to indicate if the crate ID as found in the regional header structures should be used or not.
		AliHLTInt8_t fCurrentCrateId;  ///< The current trigger crate ID number from the regional header.
		UInt_t fCurrentRegional;  ///< Index number of current regional structure being decoded.
		UInt_t fNextLocalIndex;  ///< Index number of fNextStruct local structure being decoded.
		const AliMUONLocalInfoStruct* fPrevStruct;  ///< Previous local trigger structure.
		const AliMUONLocalInfoStruct* fCurrentStruct;  ///< Current local trigger structure.
		const AliMUONLocalInfoStruct* fNextStruct;  ///< Next local trigger structure.
		bool fStoreInfo;  ///< Store debug information in fInfoBuffer.
		AliHLTUInt32_t fInfoBufferSize;  ///< Number of elements storable in fInfoBuffer.
		AliHLTUInt32_t fInfoBufferCount;  ///< Number of elements stored in the fInfoBuffer.
		AliHLTMUONTrigRecInfoStruct* fInfoBuffer;  ///< Buffer for storing the debug information.
		bool fDontPrintWrongEventError;    ///< Flag indicating if the error message for kWrongEventType is suppressed or not.
		bool fHadWrongEventTypeError;  ///< Flag indicating if a kWrongEventType error was found in the last decoded DDL.
		bool fHadNonWrongEventTypeError;  ///< Flag indicating if a different error than kWrongEventType was found in the last decoded DDL.

		static const AliMUONLocalInfoStruct fgkNullStruct; ///< Empty structure marker.
	};

	/// Not implemented
	AliHLTMUONTriggerReconstructor(const AliHLTMUONTriggerReconstructor& rhs); // copy constructor
	/// Not implemented
	AliHLTMUONTriggerReconstructor& operator = (const AliHLTMUONTriggerReconstructor& rhs); // assignment operator

	AliMUONTriggerDDLDecoder<AliDecoderHandler> fDecoder; ///< Raw DDL data decoder.
};

#endif // AliHLTMUONTRIGGERRECONSTRUCTOR_H
