#ifndef AliHLTMUONRAWDATAHISTOCOMPONENT_H
#define AliHLTMUONRAWDATAHISTOCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

///
///  @file   AliHLTMUONRawDataHistoComponent.h
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   30 April 2008
///  @brief  Declaration of a component to generate basic monitoring histograms of raw data.
///

#include "AliHLTMUONProcessor.h"
#include "AliHLTMUONDataTypes.h"
#include "AliMUONTrackerDDLDecoder.h"
#include "AliMUONTriggerDDLDecoder.h"
#include "TH1D.h"

/**
 * @class AliHLTMUONRawDataHistoComponent
 * @brief Dimuon HLT component for generating basic monitoring histograms for raw data.
 */
class AliHLTMUONRawDataHistoComponent : public AliHLTMUONProcessor
{
public:
	AliHLTMUONRawDataHistoComponent();
	virtual ~AliHLTMUONRawDataHistoComponent();

	// Public functions to implement the AliHLTProcessor interface.
	// These functions are required for the registration process.
	virtual const char* GetComponentID();
	virtual void GetInputDataTypes(AliHLTComponentDataTypeList& list);
	virtual AliHLTComponentDataType GetOutputDataType();
	virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
	virtual AliHLTComponent* Spawn();

protected:

	// Protected functions to implement the AliHLTProcessor interface.
	// These functions provide initialization as well as the actual processing
	// capabilities of the component.
	virtual int DoInit(int argc, const char** argv);
	virtual bool IgnoreArgument(const char* arg) const;
	virtual int DoDeinit();
	virtual int DoEvent(
			const AliHLTComponentEventData& evtData,
			AliHLTComponentTriggerData& trigData
		);
	
	using AliHLTProcessor::DoEvent;

private:

	class AliDecoderHandler : public AliHLTLogging
	{
	public:
	
		AliDecoderHandler() : AliHLTLogging(), fErrorHist(NULL) {}
		virtual ~AliDecoderHandler() {}
		
		/// Returns the error codes histogram.
		TH1D* ErrorHist() const { return fErrorHist; }
		
		/// Sets the error codes histogram.
		void ErrorHist(TH1D* hist) { fErrorHist = hist; }
	
	protected:
	
		/// Fills the error histogram with the given error code.
		void FillErrorHist(Int_t code);
	
		TH1D* fErrorHist;   /// Histogram of error codes found.
		
	private:
	
		// Do not allow copying of this object.
		/// Not implemented.
		AliDecoderHandler(const AliDecoderHandler& obj);
		/// Not implemented.
		AliDecoderHandler& operator = (const AliDecoderHandler& obj);
	};

	class AliTrackerDecoderHandler :
		public AliMUONTrackerDDLDecoderEventHandler, public AliDecoderHandler
	{
	public:
		AliTrackerDecoderHandler() :
			AliMUONTrackerDDLDecoderEventHandler(),
			AliDecoderHandler(),
			fManuHist(NULL),
			fSignalHist(NULL)
		{}
		
		virtual ~AliTrackerDecoderHandler() {}
		
		/// Returns the signals per MANU histogram.
		TH1D* ManuHist() const { return fManuHist; }
		
		/// Sets the signals per MANU histogram.
		void ManuHist(TH1D* hist) { fManuHist = hist; }
		
		/// Returns the signals histogram.
		TH1D* SignalHist() const { return fSignalHist; }
		
		/// Sets the signals histogram.
		void SignalHist(TH1D* hist) { fSignalHist = hist; }
		
		// Methods inherited from AliMUONTrackerDDLDecoderEventHandler:
		
		/// Called for each new data word found.
		void OnData(UInt_t data, bool /*parityError*/);
		
		/// Fills the fErrorHist histogram with the error code received.
		void OnError(ErrorCode code, const void* /*location*/) { FillErrorHist(Int_t(code)); }
	
	private:
	
		// Do not allow copying of this object.
		/// Not implemented.
		AliTrackerDecoderHandler(const AliTrackerDecoderHandler& obj);
		/// Not implemented.
		AliTrackerDecoderHandler& operator = (const AliTrackerDecoderHandler& obj);
		
		TH1D* fManuHist;    /// Histogram of signal distribution per MANU.
		TH1D* fSignalHist;  /// Histogram of the ADC signal distribution.
	};
	
	class AliTriggerDecoderHandler :
		public AliMUONTriggerDDLDecoderEventHandler, public AliDecoderHandler
	{
	public:
		AliTriggerDecoderHandler() :
			AliMUONTriggerDDLDecoderEventHandler(),
			AliDecoderHandler()
		{}
		
		virtual ~AliTriggerDecoderHandler() {}
		
		// Methods inherited from AliMUONTriggerDDLDecoderEventHandler:
		
		/// Fills the fErrorHist histogram with the error code received.
		void OnError(ErrorCode code, const void* /*location*/) { FillErrorHist(Int_t(code)); }
	
	private:
	
		// Do not allow copying of this object.
		/// Not implemented.
		AliTriggerDecoderHandler(const AliTriggerDecoderHandler& obj);
		/// Not implemented.
		AliTriggerDecoderHandler& operator = (const AliTriggerDecoderHandler& obj);
	};

	// Do not allow copying of this class.
	AliHLTMUONRawDataHistoComponent(const AliHLTMUONRawDataHistoComponent& /*obj*/);
	AliHLTMUONRawDataHistoComponent& operator = (const AliHLTMUONRawDataHistoComponent& /*obj*/);

	void ProcessTrackerDDL(const AliHLTComponentBlockData* block);
	void ProcessTriggerDDL(const AliHLTComponentBlockData* block);
	void FreeObjects();
	
	AliMUONTrackerDDLDecoder<AliTrackerDecoderHandler> fTrackerDecoder;  // Raw data decoder for the tracker data.
	AliMUONTriggerDDLDecoder<AliTriggerDecoderHandler> fTriggerDecoder;  // Raw data decoder for the trigger data.
	
	double fLastPublishTime;  /// Timestamp for the last time we published data (seconds).
	double fCurrentEventTime;  /// Timestamp for the current event being processed (seconds).
	double fPublishDelay;  /// Delay in second to wait between publishing data.
	TH1D* fErrorHist[22]; /// Histograms for error codes per DDL.
	TH1D* fManuHist[20]; /// Histograms for MANU distributions per DDL.
	TH1D* fSignalHist[20]; /// Histograms for signal distributions per DDL.
	bool fSuppressEmptyHists;  /// Flag indicating if empty histograms should be published or not.
	bool fProcessDataEventsOnly;  /// Flag indicating if only data events should be processed.
	
	ClassDef(AliHLTMUONRawDataHistoComponent, 0);  // Trigger decision component for the dimuon HLT.
};

//-----------------------------------------------------------------------------

inline void AliHLTMUONRawDataHistoComponent::AliDecoderHandler::FillErrorHist(Int_t code)
{
	/// Fills the error code into the error code histogram.
	
	assert(fErrorHist != NULL);
	Int_t mincode = Int_t( fErrorHist->GetXaxis()->GetBinCenter(1) );
	Int_t maxcode = Int_t( fErrorHist->GetXaxis()->GetBinCenter(fErrorHist->GetNbinsX()) );
	if (code < mincode or maxcode < code)
	{
		HLTError("Filling an error code which is out of range."
			" Received code %d, but expected it to be in the range [%d..%d]",
			int(code), mincode, maxcode
		);
	}
	fErrorHist->Fill(code);
}


inline void AliHLTMUONRawDataHistoComponent::AliTrackerDecoderHandler::OnData(
		UInt_t data, bool /*parityError*/
	)
{
	/// Fills the signals histogram and also the signal per MANU histogram.

	UInt_t minadc = UInt_t( fSignalHist->GetXaxis()->GetBinCenter(1) );
	UInt_t maxadc = UInt_t( fSignalHist->GetXaxis()->GetBinCenter(fSignalHist->GetNbinsX()) );
	UInt_t minmanu = UInt_t( fManuHist->GetXaxis()->GetBinCenter(1) );
	UInt_t maxmanu = UInt_t( fManuHist->GetXaxis()->GetBinCenter(fManuHist->GetNbinsX()) );
	
	UShort_t manuId; UChar_t channelId; UShort_t adc;
	UnpackADC(data, manuId, channelId, adc);
	
	if (adc < minadc or maxadc < adc)
	{
		HLTError("Filling a signal value which is out of range. Received ADC value"
			" of %d channels, but expected it to be in the range [%d..%d]",
			int(adc), minadc, maxadc
		);
	}
	fSignalHist->Fill(adc);
	
	if (manuId < minmanu or maxmanu < manuId)
	{
		HLTError("Filling a MANU ID value which is out of range. Received"
			" value of %d, but expected it to be in the range [%d..%d]",
			int(manuId), minmanu, maxmanu
		);
	}
	fManuHist->Fill(manuId);
}

#endif // AliHLTMUONRAWDATAHISTOCOMPONENT_H

