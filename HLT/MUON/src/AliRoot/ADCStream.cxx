////////////////////////////////////////////////////////////////////////////////
//
// Author: Artur Szostak
// Email:  artur@alice.phy.uct.ac.za | artursz@iafrica.com
//
////////////////////////////////////////////////////////////////////////////////

#include "AliRoot/ADCStream.hpp"
#include <TMath.h>
#include "AliHLTMUONUtils.h"


ClassImp(AliHLTMUONADCStream)


AliHLTMUONADCStream::AliHLTMUONADCStream() :
	TObject(), fData()
{
	fData.Set(0);
}


AliHLTMUONADCStream::AliHLTMUONADCStream(const UInt_t* data, UInt_t size) :
	TObject(), fData()
{
	fData.Set(size, (Int_t*)data);
}


AliHLTMUONADCStream::~AliHLTMUONADCStream()
{
	fData.Reset();
}


UInt_t AliHLTMUONADCStream::Size()
{
	return fData.GetSize();
}


void AliHLTMUONADCStream::Size(UInt_t size)
{
	fData.Set(size);
}


void AliHLTMUONADCStream::Fill(const UInt_t* data, UInt_t size)
{
	fData.Set(size, (Int_t*)data);
}


// UInt_t& AliHLTMUONADCStream::operator [] (const UInt_t index)
// {
// 	Assert( index < (UInt_t) fData.GetSize() );
// 	return (UInt_t) fData[index];
// };


UInt_t AliHLTMUONADCStream::operator [] (UInt_t index) const
{
	Assert( index < (UInt_t) fData.GetSize() );
	return fData[index];
}


ostream& operator << (ostream& os, const AliHLTMUONADCStream& s)
{
	os << "{AliHLTMUONADCStream: " << (void*)s.Data() << "}";
	os << endl;
	for (Int_t i = 0; i < s.fData.GetSize(); i++)
	{
		char buffer[32];
		char* str = (char*)&buffer[0];
		sprintf(str, "0x%X", s.fData[i]);
		os << i << "\t" << str << endl;
	}
	return os;
}

