#ifndef ALIHLTMUONHITRECONSTRUCTOR_H
#define ALIHLTMUONHITRECONSTRUCTOR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////
//Author : Indranil Das, SINP, INDIA
//         Sukalyan Chattopadhyay, SINP, INDIA
//         
//Email :  indra.das@saha.ac.in
//         sukalyan.chattopadhyay@saha.ac.in 
///////////////////////////////////////////////

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <cmath>
#include <map>
#include <cassert>

#include "TString.h"
#include "AliHLTLogging.h"
#include "AliMUONTrackerDDLDecoder.h"
#include "AliMUONTrackerDDLDecoderEventHandler.h"
#include "AliHLTMUONDataTypes.h"

#include "AliRawDataHeader.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif


extern "C" struct AliHLTMUONRecHitStruct;

//TODO: Change code to not use std::map to avoid too many AliRoot coding rule violations.
typedef std::map<AliHLTInt32_t, AliHLTInt32_t> IdManuChannelToEntry;


class AliHLTMUONHitReconstructor : public AliHLTLogging
{
public:
	
	AliHLTMUONHitReconstructor();
	virtual ~AliHLTMUONHitReconstructor(void);
	
	void SetLookUpTable(
			const AliHLTMUONHitRecoLutRow* lookupTable,
			const IdManuChannelToEntry* idToEntry
		);
	
	void SetDCCut(AliHLTInt32_t dcCut) { fDCCut = dcCut; }
	AliHLTInt32_t GetDCCut() const { return fDCCut; }
	
	bool Run(
			const AliHLTUInt32_t* rawData,
			AliHLTUInt32_t rawDataSize,
			AliHLTMUONRecHitStruct* recHit,
			AliHLTUInt32_t& nofHit
		);
	
	static AliHLTInt32_t GetkDetectorId() { return fgkDetectorId; }
	static AliHLTInt32_t GetkDDLOffSet() { return fgkDDLOffSet; }
	static AliHLTInt32_t GetkNofDDL() { return fgkNofDDL; }
	static AliHLTInt32_t GetkDDLHeaderSize() { return fgkDDLHeaderSize; }
	
	/// Returns true if the decoder is set to enable recovery logic if
	/// raw data errors are found.
	bool TryRecover() const { return fHLTMUONDecoder.TryRecover(); };
	
	/// Sets if the decoder should enable the error recovery logic.
	void TryRecover(bool value);

private:

	static const AliHLTInt32_t fgkDetectorId;      // DDL Offset
	static const AliHLTInt32_t fgkDDLOffSet;       // DDL Offset
	static const AliHLTInt32_t fgkNofDDL;          // Number of DDL 
	static const AliHLTInt32_t fgkDDLHeaderSize;   // DDL header size
	static const AliHLTInt32_t fgkLutLine;         // nof Line in LookupTable

protected:

	AliHLTMUONHitReconstructor(const AliHLTMUONHitReconstructor& rhs); // copy constructor
	AliHLTMUONHitReconstructor& operator=(const AliHLTMUONHitReconstructor& rhs); // assignment operator

private:

	struct AliHLTMUONPad
	{
		AliHLTInt32_t fDetElemId;  // The detector element ID of the pad.
		AliHLTInt32_t fIX, fIY;  // The X,Y number of the pad.
		AliHLTFloat32_t fRealX, fRealY, fRealZ;   // The real coordinate of the pad.
		AliHLTFloat32_t fHalfPadSize; // half padsize in X and Y
		AliHLTInt32_t fPlane;   // The plane and PCB zone ID numbers.
		AliHLTFloat32_t fCharge;  // The charge measured on the pad.
	};
	

	class AliHLTMUONRawDecoder : public AliMUONTrackerDDLDecoderEventHandler, public AliHLTLogging
	{
	public:
	
		AliHLTMUONRawDecoder();
		virtual ~AliHLTMUONRawDecoder();
		
		void OnData(UInt_t dataWord, bool /*parityError*/);
		void OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* /*data*/) 
		{
			fBusPatchId = int(header->fBusPatchId);
		};
		
		void OnNewBuffer(const void* buffer, UInt_t bufferSize);
		void OnError(ErrorCode code, const void* location);
		
		void SetDCCut(AliHLTInt32_t dcCut) {fDCCut = dcCut;}
		void SetPadData(AliHLTMUONPad* padData) {fPadData = padData;}
		void SetLookUpTable(const AliHLTMUONHitRecoLutRow* lookUpTableData) {fLookUpTableData = lookUpTableData;}
		void SetNofFiredDetElemId(AliHLTInt32_t& nofFiredDetElem) {fNofFiredDetElem = &nofFiredDetElem;}
		void SetIdManuChannelToEntry(const IdManuChannelToEntry* idToEntry) {fIdToEntry = idToEntry;}
		void SetMaxFiredPerDetElem(AliHLTInt32_t* maxFiredPerDetElem) {fMaxFiredPerDetElem = maxFiredPerDetElem;}
		
		AliHLTInt32_t GetDataCount() {return fDataCount;}
		
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
	
	private:
		// Do not allow copying of this class.
                /// Not implemented
		AliHLTMUONRawDecoder(const AliHLTMUONRawDecoder& rhs); // copy constructor
                /// Not implemented
		AliHLTMUONRawDecoder& operator=(const AliHLTMUONRawDecoder& rhs); // assignment operator
	
		const void* fBufferStart;   // Pointer to the start of the current DDL payload buffer.
		AliHLTInt32_t fBusPatchId;  // buspatchId
		AliHLTInt32_t fDCCut;       // DC Cut value
		AliHLTMUONPad* fPadData;    // pointer to the array containing the information of each padhits
		const AliHLTMUONHitRecoLutRow* fLookUpTableData;   // pointer to the array of Lookuptable data
		AliHLTInt32_t* fNofFiredDetElem;         // counter for detector elements that are fired 
		AliHLTInt32_t* fMaxFiredPerDetElem;      // array for detector elements that are fired 
		const IdManuChannelToEntry* fIdToEntry;  // Mapping between Linenumber to IdManuChannel;
		
		AliHLTInt32_t fDataCount;           // Data Counter
		AliHLTInt32_t fPrevDetElemId;       // previous detection elementId
		AliHLTInt32_t fPadCharge;           // pedestal subtracted pad charge
		AliHLTFloat32_t fCharge;            //calibrated pad charge 
		AliHLTInt32_t fIdManuChannel;       // id manu channel
		AliHLTInt32_t fLutEntry;            // i-th entry in lookuptable
		
		bool fWarnOnly;  ///< Flag indicating if the OnError method should generate warnings rather than error messages.
	};

	AliMUONTrackerDDLDecoder<AliHLTMUONRawDecoder> fHLTMUONDecoder; // robust HLTMUON Decoder
	
	AliHLTInt32_t fkBlockHeaderSize;      // Block header size
	AliHLTInt32_t fkDspHeaderSize;        // DSP header size
	AliHLTInt32_t fkBuspatchHeaderSize;   // buspatch header size
	
	AliHLTInt32_t fDCCut;                 // DC Cut value
	
	AliHLTMUONPad* fPadData;  // pointer to the array containing the information of each padhits
	const AliHLTMUONHitRecoLutRow* fLookUpTableData;  // pointer to the array of Lookuptable data (The memory is not owned by this component).
	
	AliHLTMUONRecHitStruct* fRecPoints;      // Reconstructed hits
	AliHLTUInt32_t *fRecPointsCount;         // nof reconstructed hit.
	AliHLTUInt32_t fMaxRecPointsCount;       // max nof reconstructed hit.
	
	AliHLTInt32_t fCentralCountB, fCentralCountNB;   // centeral hits.
	AliHLTInt32_t fDigitPerDDL;                      // Total nof Digits perDDL.
	
	AliHLTInt32_t *fCentralChargeB, *fCentralChargeNB;       // pointer to an array of central hit
	AliHLTFloat32_t *fRecX, *fRecY;                          // pointer to an array of reconstructed hit
	AliHLTFloat32_t *fAvgChargeX, *fAvgChargeY;              // average charge on central pad found using CG method
	AliHLTInt32_t *fNofBChannel, *fNofNBChannel;             // number of channels bending and non-bending.
	AliHLTInt32_t fGetIdTotalData[336][237][2];              // an array of idManuChannel with argument of centralX, centralY and planeType.
	AliHLTInt32_t fNofFiredDetElem,fMaxFiredPerDetElem[130];  // counter for detector elements that are fired
	const IdManuChannelToEntry* fIdToEntry;   // Mapping between Linenumber to IdManuChannel (The object is not owned by this component).
	
	bool DecodeDDL(const AliHLTUInt32_t* rawData, AliHLTUInt32_t rawDataSize);
	void FindCentralHits(AliHLTInt32_t minPadId, AliHLTInt32_t maxPadId);
	bool FindRecHits();
	void RecXRecY();
	bool MergeRecHits();
	void Clear();

};

#endif // ALIHLTMUONHITRECONSTRUCTOR_H
