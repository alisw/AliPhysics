////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIHLTMUONADCSTREAM_H
#define ALIHLTMUONADCSTREAM_H

#include <TObject.h>
#include <TArrayI.h>
#include <Riostream.h>


class AliHLTMUONADCStream : public TObject
{
	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const AliHLTMUONADCStream& s);

public:

	/* Default constructor initialises everything to zero.
	 */
	AliHLTMUONADCStream();
	
	AliHLTMUONADCStream(const UInt_t* data, UInt_t size);

	virtual ~AliHLTMUONADCStream();
	
	UInt_t Size();
	void Size(UInt_t size);
	void Fill(const UInt_t* data, UInt_t size);
	
	UInt_t* Data() { return (UInt_t*) fData.GetArray(); }
	const UInt_t* Data() const { return (UInt_t*) fData.GetArray(); }
	
// 	UInt_t& operator [] (const UInt_t index);
	UInt_t operator [] (UInt_t index) const;

private:

	// TODO: complete the ADC stream specification.
	TArrayI fData;  // The DDL raw data.

	ClassDef(AliHLTMUONADCStream, 1)  // ADC stream data.
};


#endif // ALIHLTMUONADCSTREAM_H
