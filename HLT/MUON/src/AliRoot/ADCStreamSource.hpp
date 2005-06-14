////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_ADC_STREAM_SOURCE_HPP
#define dHLT_ALIROOT_ADC_STREAM_SOURCE_HPP

#include "TROOT.h"
#include "TObject.h"
#include "TString.h"
#include "AliRoot/ADCStream.hpp"

#include <vector>

namespace AliMUONHLT
{


class ADCStreamSource : public TObject
{
public:

	ADCStreamSource();
	virtual ~ADCStreamSource();

	/* Fills the internal data structures from the specified file
	 */
	void FillFromFile(const TString& filename, const Int_t eventnumber);
	
	/* Fills the internal data structures from the specified directory.
	   FillFromFile is called for every file in the directory that is
	   prefixed with MUON_ and ends in .ddl
	 */
	void FillFrom(const TString& directory, const Int_t eventnumber);
	
	/* Same as the methods above except the directory name is created as
	   dirprefix + eventnumber, where eventnumber is looped from firstevent
	   to lastevent.
	 */
	void FillFrom(const TString& dirprefix, const UInt_t firstevent, const UInt_t lastevent);
	
	/* Clears all the internal arrays.
	 */
	void Clear();
	
	// Get methods.
	Int_t CurrentStream()  const { return fCurrentStream; };
	
	/* Returns the number of ADC streams stored.
	 */
	Int_t NumberOfStreams() const;
	
	/* Fetches the index'th ADC stream stored.
	   kTRUE is returned if the stream was found, kFALSE otherwise.
	 */
	Bool_t GetStream(const Int_t index) const;
	
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
	Bool_t FetchStream(ADCStream& stream) const;
	
	/* Returns the index'th ADC stream.
	   kTRUE is returned if the stream was found, kFALSE otherwise.
	 */
	Bool_t FetchStream(const Int_t index, ADCStream& stream) const;
	
	/* Returns the current ADC stream selected.
	   A NULL pointer is returned if no ADC stream is selected.
	 */
	const ADCStream* FetchStream() const;

private:

	/* Adds a new ADCStream object to the internal arrays.
	 */
	void AddStream(ADCStream& stream, const UInt_t eventnumber);

	mutable Int_t fCurrentStream;  //! The currently selected stream index.
	

public:  // Unfortunately ROOT requires the following to be public.
	
	struct DataBlock
	{
		virtual ~DataBlock() {};

		Int_t fEventNumber;  // Event number of the stream.
		ADCStream fStream;  // The ADC stream block.
		
		ClassDef(DataBlock, 1);  // Data per event.
	};

private:

	std::vector<DataBlock> fList;  // List of ADC streams.

	ClassDef(ADCStreamSource, 1);  // The source of ADC stream data for dHLT.
};


}; // AliMUONHLT

#endif // dHLT_ALIROOT_ADC_STREAM_SOURCE_HPP
