#ifndef ALIHLTMUONHITRECONSTRUCTOR_H
#define ALIHLTMUONHITRECONSTRUCTOR_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

// $Id$

///////////////////////////////////////////////
//Author : Indranil Das, SINP, INDIA
//         Sukalyan Chattopadhyay, SINP, INDIA
//         
//Email :  indra.das@saha.ac.in
//         sukalyan.chattopadhyay@saha.ac.in 
///////////////////////////////////////////////

#include <map>

#include "AliHLTLogging.h"
#include "AliMUONTrackerDDLDecoder.h"
#include "AliMUONTrackerDDLDecoderEventHandler.h"
#include "AliHLTMUONDataTypes.h"

#if __GNUC__ && __GNUC__ < 3
#define std
#endif


extern "C" struct AliHLTMUONRecHitStruct;
extern "C" struct AliHLTMUONClusterStruct;
extern "C" struct AliHLTMUONChannelStruct;
typedef std::map<AliHLTInt32_t, AliHLTInt32_t> IdManuChannelToEntry;
typedef IdManuChannelToEntry MaxEntryPerBusPatch;

class AliHLTMUONHitReconstructor : public AliHLTLogging
{
public:
	
	AliHLTMUONHitReconstructor();
	virtual ~AliHLTMUONHitReconstructor(void);
	
	void SetLookUpTable(
			const AliHLTMUONHitRecoLutRow* lookupTable,
			const IdManuChannelToEntry* idToEntry,
			const MaxEntryPerBusPatch* maxEntryPerBP
		);
	
	void SetDCCut(AliHLTInt32_t dcCut) { fDCCut = dcCut; }
	AliHLTInt32_t GetDCCut() const { return fDCCut; }
	
	bool Run(
			AliHLTUInt32_t* rawData,
			AliHLTUInt32_t rawDataSize,
			AliHLTMUONRecHitStruct* const recHit,
			AliHLTUInt32_t& nofHit
		);
	
	/**
	 * Fills the output clusters array with the extra cluster information generated.
	 * If the method GenerateClusterInfo(true) was not called, then no cluster information
	 * is generated and this method will not do anything.
	 * [out] \param clusters  This is the output array that will be filled.
	 * [in/out] \param nofClusters Initially this contains the maximum number of elements
	 *     that can be stored in the clusters array. The method will fill this with
	 *     the actual number of elements stored.
	 * \returns true if all elements were copied and false if there is not enough space in
	 *     the output array.
	 */
	bool FillClusterData(AliHLTMUONClusterStruct* clusters, AliHLTUInt32_t& nofClusters);
	
	/**
	 * Fills the output channels array with the extra channel information generated.
	 * If the method GenerateChannelInfo(true) was not called, then no extra channel
	 * information is generated and this method will not do anything.
	 * [out] \param channels  This is the output array that will be filled.
	 * [in/out] \param nofChannels Initially this contains the maximum number of elements
	 *     that can be stored in the channels array. The method will fill this with
	 *     the actual number of elements stored.
	 * \returns true if all elements were copied and false if there is not enough space in
	 *     the output array.
	 */
	bool FillChannelData(AliHLTMUONChannelStruct* channels, AliHLTUInt32_t& nofChannels);
	
	static AliHLTInt32_t GetkDetectorId() { return fgkDetectorId; }
	static AliHLTInt32_t GetkDDLOffSet() { return fgkDDLOffSet; }
	static AliHLTInt32_t GetkNofDDL() { return fgkNofDDL; }
	static AliHLTInt32_t GetkDDLHeaderSize() { return fgkDDLHeaderSize; }
	static AliHLTInt32_t GetkNofDetElemInDDL(Int_t iDDL);
	static AliHLTInt32_t GetkMinDetElemIdInDDL(Int_t iDDL);
	
	/// The error recovery mode used for TryRecover.
	enum ERecoveryMode
	{
		kDontTryRecover = 0,  /// Will not try recover from errors.
		kRecoverFull,  /// Try recover from all errors.
		kRecoverJustSkip,  /// Just skip any corrupt structures.
		kRecoverFromParityErrorsOnly  /// Recover only from parity errors.
	};
	
	/// Returns the recovery mode used for the TryRecover option.
	/// This controls if the decoder is set to enable recovery logic if
	/// raw data errors are found.
	ERecoveryMode TryRecover() const { return fRecoveryMode; }
	
	/// Sets if the decoder should enable the error recovery logic and how.
	void TryRecover(ERecoveryMode mode);
	
	/// Returns true if ADC digits with parity errors are skipped.
	bool SkipParityErrors() const { return fHLTMUONDecoder.GetHandler().SkipParityErrors(); }
	
	/// Sets the flag indicating if ADC digits with parity errors are skipped.
	void SkipParityErrors(bool value) { fHLTMUONDecoder.GetHandler().SkipParityErrors(value); }
	
	/// Returns true if messages about parity errors are not printed.
	bool DontPrintParityErrors() const { return fHLTMUONDecoder.GetHandler().DontPrintParityErrors(); }
	
	/// Sets the flag indicating if messages about parity errors are not printed.
	void DontPrintParityErrors(bool value) { fHLTMUONDecoder.GetHandler().DontPrintParityErrors(value); }
	
	/// Returns true if the extra cluster information should be generated.
	bool GenerateClusterInfo() const { return fGenerateClusterInfo; }
	
	/// Sets the flag to indicate if the extra cluster information should be generated.
	void GenerateClusterInfo(bool value) { fGenerateClusterInfo = value; }
	
	/// Returns true if the extra channel information should be generated.
	bool GenerateChannelInfo() const { return fGenerateChannelInfo; }
	
	/// Sets the flag to indicate if the extra channel information should be generated.
	void GenerateChannelInfo(bool value) { fGenerateChannelInfo = value; }
	
	/// Returns the maximum channel multiplicity allowed per cluster.
	bool MaxChannelMultiplicity() const { return fMaxChannelMult; }
	
	/// Sets the maximum channel multiplicity allowed per cluster.
	/// \note This effects the memory allocation required. Generally M*N number of
	///   channel structures will be allocated, where M = fMaxChannelMult and
	///   N = the maximum number of hits that can be filled into the output buffers.
	void MaxChannelMultiplicity(bool value) { fMaxChannelMult = value; }
	
	/// Returns the DDL number the component expects data to be received from.
	AliHLTInt32_t DDLNumber() const { return fDDL; }
	
	/// Sets the DDL number the component expects data to be received from.
	void DDLNumber(AliHLTInt32_t value) { fDDL = (value & 0x1F); }  // 0x1F forces value into our required range.

	bool InitDetElemInDDLArray();
	bool DeInitDetElemInDDLArray();

private:

	static const AliHLTInt32_t fgkDetectorId;      // DDL Offset
	static const AliHLTInt32_t fgkDDLOffSet;       // DDL Offset
	static const AliHLTInt32_t fgkNofDDL;          // Number of DDL 
	static const AliHLTInt32_t fgkDDLHeaderSize;   // DDL header size
	static const AliHLTInt32_t fgkLutLine;         // nof Line in LookupTable
	static const AliHLTInt32_t fgkMaxNofDataPerDetElem;         // maximum allowed data points per detlem
	static const AliHLTInt32_t fgkNofDetElemInDDL[20] ;         // nof Detelem in a given ddl
	static const AliHLTInt32_t fgkMinDetElemIdInDDL[20] ;       // the detelem which has minimum value in ddl

	AliHLTMUONHitReconstructor(const AliHLTMUONHitReconstructor& rhs); // copy constructor
	AliHLTMUONHitReconstructor& operator=(const AliHLTMUONHitReconstructor& rhs); // assignment operator

	struct AliHLTMUONPad
	{
		AliHLTInt32_t fDetElemId;  // The detector element ID of the pad.
		AliHLTInt32_t fIX, fIY;  // The X,Y number of the pad.
		AliHLTFloat32_t fRealX, fRealY, fRealZ;   // The real coordinate of the pad.
		AliHLTFloat32_t fHalfPadSize; // half padsize in X and Y
	        AliHLTFloat32_t fPadSizeXY; // padsize in Y for bending plane and X for nonbending
		AliHLTInt32_t fPlane;   // The plane and PCB zone ID numbers.
		AliHLTFloat32_t fCharge;  // The charge measured on the pad.
		AliHLTInt32_t fBusPatch;  // The bus patch of the raw data word from the DDL stream.
		AliHLTUInt32_t fRawData;  // The raw data word from the DDL stream.
	};
	

	class AliHLTMUONRawDecoder : public AliMUONTrackerDDLDecoderEventHandler, public AliHLTLogging
	{
	public:
	
		AliHLTMUONRawDecoder();
		virtual ~AliHLTMUONRawDecoder();
		
		void OnData(UInt_t dataWord, bool parityError);
		void OnNewBusPatch(const AliMUONBusPatchHeaderStruct* header, const void* data) ;
		
		void OnNewBuffer(const void* buffer, UInt_t bufferSize);
		void OnError(ErrorCode code, const void* location);
		
		void SetDCCut(AliHLTInt32_t dcCut) {fDCCut = dcCut;}
		void SetPadData(AliHLTMUONPad* const padData) {fPadData = padData;}
		void SetLookUpTable(const AliHLTMUONHitRecoLutRow* lookUpTableData) {fkLookUpTableData = lookUpTableData;}

		void SetIdManuChannelToEntry(const IdManuChannelToEntry* idToEntry) {fkIdToEntry = idToEntry;}
		void SetMaxEntryPerBusPatch(const MaxEntryPerBusPatch* maxEntryPerBP) {fkMaxEntryPerBusPatch = maxEntryPerBP;}

		AliHLTInt32_t DDLNumber() const { return fDDL; }
		void DDLNumber(AliHLTInt32_t value) { fDDL = (value & 0x1F); }  // 0x1F forces value into our required range.
		// The following two methods have to called after set the ddl
		void SetNofFiredDetElemId(AliHLTUInt16_t* const nofDataInDetElem) {fNofDataInDetElem = nofDataInDetElem;}
		void SetMaxFiredPerDetElem(AliHLTUInt16_t** const dataCountListPerDetElem) {fDataCountListPerDetElem = dataCountListPerDetElem;}		


		AliHLTInt32_t GetDataCount() const {return fDataCount;}
		
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
		 * Returns true if ADC digits with parity errors are skipped.
		 */
		bool SkipParityErrors() const { return fSkipParityErrors; }
		
		/**
		 * Sets the flag indicating if ADC digits with parity errors are skipped.
		 */
		void SkipParityErrors(bool value) { fSkipParityErrors = value; }
		
		/**
		 * Returns true if messages about parity errors are not printed.
		 */
		bool DontPrintParityErrors() const { return fDontPrintParityErrors; }
		
		/**
		 * Sets the flag indicating if messages about parity errors are not printed.
		 */
		void DontPrintParityErrors(bool value) { fDontPrintParityErrors = value; }
		
		/**
		 * Returns true if parity error messages are printed as warnings.
		 */
		bool PrintParityErrorAsWarning() const { return fPrintParityErrorAsWarning; }
		
		/**
		 * Sets the flag indicating if parity error messages are printed as warnings.
		 */
		void PrintParityErrorAsWarning(bool value) { fPrintParityErrorAsWarning = value; }
		
		/**
		 * Returns true if a parity error was found during the last call to Decode.
		 */
		bool ParityErrorFound() const { return fParityErrorFound; }
		
		/**
		 * Returns true if a non-parity error was found during the last call to Decode.
		 */
		bool NonParityErrorFound() const { return fNonParityErrorFound; }
	
	private:
		// Do not allow copying of this class.
                /// Not implemented
		AliHLTMUONRawDecoder(const AliHLTMUONRawDecoder& rhs); // copy constructor
                /// Not implemented
		AliHLTMUONRawDecoder& operator=(const AliHLTMUONRawDecoder& rhs); // assignment operator
	
		const void* fkBufferStart;   // Pointer to the start of the current DDL payload buffer.
		AliHLTInt32_t fBusPatchId;  // buspatchId
		AliHLTInt32_t fDCCut;       // DC Cut value
		AliHLTMUONPad* fPadData;    // pointer to the array containing the information of each padhits
		const AliHLTMUONHitRecoLutRow* fkLookUpTableData;   // pointer to the array of Lookuptable data
		/* AliHLTInt32_t* fNofFiredDetElem;         // counter for detector elements that are fired  */
		/* AliHLTInt32_t* fMaxFiredPerDetElem;      // array for detector elements that are fired  */
		const IdManuChannelToEntry* fkIdToEntry;  // Mapping between Linenumber to IdManuChannel;
		const MaxEntryPerBusPatch* fkMaxEntryPerBusPatch;   // Maximum allowed entry per Buspatch.
		
		AliHLTInt32_t fDDL;                 // DDL number
		AliHLTInt32_t fDataCount;           // Data Counter
		AliHLTInt32_t fPrevDetElemId;       // previous detection elementId
		AliHLTInt32_t fPadCharge;           // pedestal subtracted pad charge
		AliHLTFloat32_t fCharge;            //calibrated pad charge 
		AliHLTInt32_t fIdManuChannel;       // id manu channel
		AliHLTInt32_t fLutEntry;            // i-th entry in lookuptable
		
		AliHLTUInt16_t **fDataCountListPerDetElem;	         ///< List of datacounts associated with given ddl
		AliHLTUInt16_t *fNofDataInDetElem;	                 ///< Nof datacount in a ddl
		bool fWarnOnly;  ///< Flag indicating if the OnError method should generate warnings rather than error messages.
		bool fSkipParityErrors;  ///< Flag indicating if ADC digits with parity errors should be skipped.
		bool fDontPrintParityErrors;  ///< Flag for controlling if messages about parity errors should be printed.
		bool fPrintParityErrorAsWarning;  ///< Flag for controlling if parity error messages are printed as warnings.
		bool fParityErrorFound;  ///< Flag if a parity error was found after a decoding pass.
		bool fNonParityErrorFound;  ///< Flag which indicates if a non parity error code was found after a decoding pass.
		bool fIsMuchNoisy;  ///< tag for noisy buspatch.
	};

	AliMUONTrackerDDLDecoder<AliHLTMUONRawDecoder> fHLTMUONDecoder; // robust HLTMUON Decoder
	
	AliHLTInt32_t fkBlockHeaderSize;      // Block header size
	AliHLTInt32_t fkDspHeaderSize;        // DSP header size
	AliHLTInt32_t fkBuspatchHeaderSize;   // buspatch header size
	
	AliHLTInt32_t fDCCut;                 // DC Cut value
	
	AliHLTMUONPad* fPadData;  // pointer to the array containing the information of each padhits
	const AliHLTMUONHitRecoLutRow* fkLookUpTableData;  // pointer to the array of Lookuptable data (The memory is not owned by this component).
	
	AliHLTMUONRecHitStruct* fRecPoints;      // Reconstructed hits
	AliHLTUInt32_t *fRecPointsCount;         // nof reconstructed hit.
	AliHLTUInt32_t fMaxRecPointsCount;       // max nof reconstructed hit.
	
	AliHLTMUONClusterStruct* fClusters;  // Array of output cluster infromation.
	AliHLTUInt32_t fClusterCount;       // Number of elemenets in fClusters.
	AliHLTUInt32_t fMaxClusters;        // Maximum number of clusters in fClusters.
	bool fGenerateClusterInfo;   // Flag indicating if extra cluster information should be generated.
	AliHLTInt32_t fNewClusterId;  // The ID to use for a new cluster structure.
	AliHLTInt32_t fDDL;   // The source DDL number of the raw data.
	
	AliHLTMUONChannelStruct* fChannels;  // Array of output channel infromation.
	AliHLTUInt32_t fChannelCount;       // Number of elemenets in fChannels.
	AliHLTUInt32_t fMaxChannels;        // Maximum number of channels in fChannels.
	bool fGenerateChannelInfo;   // Flag indicating if extra channel information should be generated.
	AliHLTUInt32_t fMaxChannelMult;  // Indicates the maximum channel multiplicity per cluster allowed.
	
	AliHLTInt32_t fCentralCountB, fCentralCountNB;   // centeral hits.
	AliHLTInt32_t fDigitPerDDL;                      // Total nof Digits perDDL.

	AliHLTUInt16_t **fDataCountListPerDetElem;	         // List of datacounts associated with given ddl
	AliHLTUInt16_t *fNofDataInDetElem;	                 // Nof datacount in a ddl
	AliHLTInt32_t *fCentralChargeB, *fCentralChargeNB;       // pointer to an array of central hit
	AliHLTFloat32_t *fRecX, *fRecY;                          // pointer to an array of reconstructed hit
	AliHLTFloat32_t *fAvgChargeX, *fAvgChargeY;              // average charge on central pad found using CG method
	AliHLTFloat32_t *fTotChargeX, *fTotChargeY;              // Total charge in bending and nonbending direction
	AliHLTInt32_t *fNofBChannel, *fNofNBChannel;             // number of channels bending and non-bending.
	AliHLTInt32_t *fNofYNeighbour;                           // number of neighbour pad in y direction, needed for y-resolution correction
	AliHLTInt32_t fGetIdTotalData[336][237][2];              // an array of idManuChannel with argument of centralX, centralY and planeType.
	//AliHLTInt32_t fNofFiredDetElem,fMaxFiredPerDetElem[130];  // counter for detector elements that are fired
	const IdManuChannelToEntry* fkIdToEntry;   // Mapping between Linenumber to IdManuChannel (The object is not owned by this component).
	const MaxEntryPerBusPatch* fkMaxEntryPerBusPatch;   // Maximum allowed entry per Buspatch.
	
	ERecoveryMode fRecoveryMode;  ///< The recovery mode for the decoder.
	
	bool DecodeDDL(AliHLTUInt32_t* rawData, AliHLTUInt32_t rawDataSize);
	void FindCentralHits(AliHLTInt32_t iDet);
	bool FindRecHits();
	void RecXRecY();
	bool MergeSlatRecHits();
	bool MergeQuadRecHits();
	void Clear();

};

#endif // ALIHLTMUONHITRECONSTRUCTOR_H
