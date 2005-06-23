////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#ifndef dHLT_ALIROOT_ADC_STREAM_HPP
#define dHLT_ALIROOT_ADC_STREAM_HPP

#include <TObject.h>
#include <TArrayI.h>
#include <Riostream.h>


namespace AliMUONHLT
{


class ADCStream : public TObject
{
public:

	/* Default constructor initialises everything to zero.
	 */
	ADCStream();
	
	ADCStream(const UInt_t* data, UInt_t size);

	virtual ~ADCStream();
	
	UInt_t Size();
	void Size(UInt_t size);
	void Fill(const UInt_t* data, UInt_t size);
	
	UInt_t* Data() { return (UInt_t*) fData.GetArray(); }
	const UInt_t* Data() const { return (UInt_t*) fData.GetArray(); }
	
// 	UInt_t& operator [] (const UInt_t index);
	UInt_t operator [] (UInt_t index) const;

	// ostream operator usefull for text output.
	friend ostream& operator << (ostream& os, const ADCStream& s);

private:

	// TODO: complete the ADC stream specification.
	TArrayI fData;

	ClassDef(ADCStream, 1)  // ADC stream data.
};


} // AliMUONHLT

#endif // dHLT_ALIROOT_ADC_STREAM_HPP
