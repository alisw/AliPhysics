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

extern "C" struct AliHLTMUONRecHitStruct;
extern "C" struct AliHLTMUONTriggerRecordStruct;


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
		 * Sets the DDL bit according to the DDL value given.
		 */
		void SetDDL(AliHLTInt32_t ddl) { fDDLBit = (ddl == 20 ? 0x00 : 0x80); }
		
		/**
		 * Generates reconstructed hits from strip information.
		 */
		bool GenerateHits(
				AliHLTMUONRecHitStruct* outputBuffer,
				AliHLTUInt32_t& maxEntries
			);
		
		// Methods inherited from AliMUONTriggerDDLDecoderEventHandler:
		
		/// Called for each new buffer. Just remember the start location of the buffer.
		void OnNewBuffer(const void* buffer, UInt_t /*bufferSize*/)
		{
			assert( buffer != NULL );
			fBufferStart = buffer;
		}
		
		/**
		 * Sets the regional structure sequencial number and decodes the crate ID.
		 * if fUseCrateId is false then we use the sequencial number instead. This
		 * might be necessary for for incorrectly generated or buggy raw data.
		 */
		void OnNewRegionalStructV2(
				UInt_t num,
				const AliMUONRegionalHeaderStruct* regionalStruct,
				const AliMUONRegionalScalarsStruct* /*scalars*/,
				const void* /*data*/
			)
		{
			fCurrentRegional = num;
			fCurrentCrateId = (fUseCrateId ? GetRegionalId(regionalStruct) : num);
		}
		
		/// Converts a local trigger structure from the L0 into a trigger record.
		void OnLocalStructV2(
				UInt_t num,
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
		
		/**
		 * Finds the strip bits / positions on MT1 that were fired given
		 * the local trigger structure decision.
		 */
		bool FindStripsOnMT1(
				const AliMUONLocalInfoStruct* localStruct,
				AliHLTInt32_t& xPos, AliHLTInt32_t& yPos
			);
		
		/**
		 * Tries to find the fired X strips for chambers 11 to 14.
		 */
		void FindXStrips(
				const AliMUONLocalInfoStruct* localStruct,
				AliHLTInt32_t startPos, AliHLTInt32_t pos[4]
			);
		
		/**
		 * Tries to find the fired Y strips for chambers 11 to 14.
		 */
		void FindYStrips(
				const AliMUONLocalInfoStruct* localStruct,
				AliHLTInt32_t startPos, AliHLTInt32_t pos[4]
			);
		
		/**
		 * Reconstructs a hit with global position coordinates from strip
		 * information on a given chamber.
		 */
		void ReconstructHit(
				AliHLTUInt32_t xStrips, AliHLTUInt32_t yStrips,
				AliHLTInt32_t xPos, AliHLTInt32_t yPos,
				AliHLTUInt8_t crateId, AliHLTUInt8_t locId,
				AliHLTUInt8_t chamber, AliHLTMUONRecHitStruct& hit
			);
		
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
		bool fUseCrateId;  ///< Flag to indicate if the crate ID as found in the regional header structures should be used or not.
		AliHLTInt8_t fCurrentCrateId;  ///< The current trigger crate ID number from the regional header.
		UInt_t fCurrentRegional;  ///< Index number of current regional structure being decoded.
	};

	/// Not implemented
	AliHLTMUONTriggerReconstructor(const AliHLTMUONTriggerReconstructor& rhs); // copy constructor
	/// Not implemented
	AliHLTMUONTriggerReconstructor& operator = (const AliHLTMUONTriggerReconstructor& rhs); // assignment operator

	AliMUONTriggerDDLDecoder<AliDecoderHandler> fDecoder; ///< Raw DDL data decoder.
};

#endif // AliHLTMUONTRIGGERRECONSTRUCTOR_H
