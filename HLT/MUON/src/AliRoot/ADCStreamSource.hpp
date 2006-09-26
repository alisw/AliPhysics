////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

/* AliHLTMUONADCStreamSource is used to store the raw DDL data generated
   from an AliRoot simulation.
   This class is used as a storage class for the input dataset for
   AliHLTMUONMicrodHLT.
 */

#ifndef ALIHLTMUONADCSTREAMSOURCE_H
#define ALIHLTMUONADCSTREAMSOURCE_H

#include "TROOT.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TString.h"
#include "AliRoot/ADCStream.hpp"


class AliHLTMUONADCStreamSource : public TObject
{
public:  // Unfortunately ROOT requires the following to be public.
	
	class AliDataBlock : public TObject
	{
	public:
		AliDataBlock() : fEventNumber(0), fStream() {}
		virtual ~AliDataBlock() {};
		
		Int_t& EventNumber() { return fEventNumber; }
		AliHLTMUONADCStream& Stream() { return fStream; }

	private:
	
		Int_t fEventNumber;  // Event number of the stream.
		AliHLTMUONADCStream fStream;  // The ADC stream block.
		
		ClassDef(AliDataBlock, 1)  // Data per event.
	};

public:

	AliHLTMUONADCStreamSource();
	virtual ~AliHLTMUONADCStreamSource();

	/* Fills the internal data structures from the specified file
	 */
	void FillFromFile(const TString& filename, Int_t eventnumber);
	
	/* Fills the internal data structures from the specified directory.
	   FillFromFile is called for every file in the directory that is
	   prefixed with MUON_ and ends in .ddl
	 */
	void FillFrom(const TString& directory, Int_t eventnumber);
	
	/* Same as the methods above except the directory name is created as
	   dirprefix + eventnumber, where eventnumber is looped from firstevent
	   to lastevent.
	 */
	void FillFrom(const TString& dirprefix, UInt_t firstevent, UInt_t lastevent);
	
	/* Clears all the internal arrays.
	 */
	virtual void Clear(Option_t* option = "");
	
	// Get methods.
	Int_t CurrentStream()  const { return fCurrentStream; };
	
	/* Returns the number of ADC streams stored.
	 */
	Int_t NumberOfStreams() const;
	
	/* Fetches the index'th ADC stream stored.
	   kTRUE is returned if the stream was found, kFALSE otherwise.
	 */
	Bool_t GetStream(Int_t index) const;
	
	/* Fetches the first ADC stream stored.
	   kTRUE is returned if the stream was found, kFALSE otherwise.
	 */
	Bool_t FirstStream() const;
	
	/* Fetches the next ADC stream stored following the currently selected one.
	   kTRUE is returned if the stream was found, kFALSE otherwise.
	 */
	Bool_t NextStream() const;
	
	/* Returns the corresponding AliRoot event number for the current stream.
	   -1 is returned if no stream is selected.
	 */
	Int_t EventNumber() const;
	
	/* Returns the current ADC stream selected.
	   kFALSE is returned if there is no stream selected.
	 */
	Bool_t FetchStream(AliHLTMUONADCStream& stream) const;
	
	/* Returns the index'th ADC stream.
	   kTRUE is returned if the stream was found, kFALSE otherwise.
	 */
	Bool_t FetchStream(Int_t index, AliHLTMUONADCStream& stream) const;
	
	/* Returns the current ADC stream selected.
	   A NULL pointer is returned if no ADC stream is selected.
	 */
	const AliHLTMUONADCStream* FetchStream() const;

private:

	/* Adds a new AliHLTMUONADCStream object to the internal arrays.
	 */
	void AddStream(AliHLTMUONADCStream& stream, UInt_t eventnumber);

	mutable Int_t fCurrentStream;  //! The currently selected stream index.

	TClonesArray fList;  // List of ADC streams.

	ClassDef(AliHLTMUONADCStreamSource, 1)  // The source of ADC stream data for dHLT.
};


#endif // ALIHLTMUONADCSTREAMSOURCE_H
