/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for TOF Online calibration: defining channel status                 //
// New object created, to use an array instead of a TObjArray.               //
// Storing all the info coming from HW FEE map, pulser runs, and noise       //
// runs in a single object (char).                                           // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <AliTOFChannelOnlineStatusArray.h>
#include <AliLog.h>

ClassImp(AliTOFChannelOnlineStatusArray)

//________________________________________________________________
AliTOFChannelOnlineStatusArray::AliTOFChannelOnlineStatusArray():
	TObject(),
	fSize(0),
	fArray(0x0)
{
	//default constructor
}
//________________________________________________________________
AliTOFChannelOnlineStatusArray::~AliTOFChannelOnlineStatusArray()
{
	//distructor
	delete [] fArray;
}
//________________________________________________________________
AliTOFChannelOnlineStatusArray::AliTOFChannelOnlineStatusArray(Int_t size):
	TObject(),
	fSize(size),
	fArray(new UChar_t[size])
{
	// ctor with size
	for (Int_t ich = 0; ich<size; ich ++){
	  SetStatus(ich,kTOFOnlineUnknown);
	}
}
//________________________________________________________________
AliTOFChannelOnlineStatusArray::AliTOFChannelOnlineStatusArray(const AliTOFChannelOnlineStatusArray & source):
      TObject(),
      fSize(source.fSize),
      fArray(source.fArray)
{ 
	// copy constructor
}
//________________________________________________________________
AliTOFChannelOnlineStatusArray &AliTOFChannelOnlineStatusArray::operator=(const AliTOFChannelOnlineStatusArray & source) 
{ 
	// assignment operator

  if (this == &source)
    return *this;

  TObject::operator=(source);
  fSize= source.fSize;
  fArray= source.fArray;
  return *this;
}
//________________________________________________________________
void AliTOFChannelOnlineStatusArray::SetStatus(Int_t pos, UChar_t parr)
{
	// setting status for channel at position = pos
	AliDebug(2,Form("status = %d",(UInt_t)parr));
	if (pos>-1 && pos < fSize)fArray[pos] = parr;
	AliDebug(2,Form("fArray[%d] = %d",pos,(UInt_t)fArray[pos]));
}
//________________________________________________________________
void AliTOFChannelOnlineStatusArray::SetHWStatus(Int_t pos, UChar_t parr)
{
	// setting status for channel at position = pos
	AliDebug(2,Form("HW status = %d",(UInt_t)parr));
	if (pos>-1 && pos < fSize) {
		fArray[pos] &= kTOFHWReset;
		fArray[pos] |= parr;
	}
	AliDebug(2,Form("fArray[%d] = %d",pos,(UInt_t)fArray[pos]));
}
//________________________________________________________________
void AliTOFChannelOnlineStatusArray::SetPulserStatus(Int_t pos, UChar_t parr)
{
	// setting status for channel at position = pos
	AliDebug(2,Form("Pulser status = %d",(UInt_t)parr));
	if (pos>-1 && pos < fSize){
		fArray[pos] &= kTOFPulserReset;
		fArray[pos] |= parr;
	}
	AliDebug(2,Form("fArray[%d] = %d",pos,(UInt_t)fArray[pos]));
}
//________________________________________________________________
void AliTOFChannelOnlineStatusArray::SetNoiseStatus(Int_t pos, UChar_t parr)
{
	// setting status for channel at position = pos
	AliDebug(2,Form("Noise status = %d",(UInt_t)parr));
	if (pos>-1 && pos < fSize){
		fArray[pos] &= kTOFNoiseReset;
		fArray[pos] |= parr;
	}
	AliDebug(2,Form("fArray[%d] = %d",pos,(UInt_t)fArray[pos]));
}
//________________________________________________________________
UChar_t AliTOFChannelOnlineStatusArray::GetStatus(Int_t pos) const 
{
	// getting the status for channel at position = pos 
	UChar_t parr = 0x0; 
	if  (pos>-1 && pos < fSize)parr = fArray[pos];
	return parr;
}
//________________________________________________________________
UChar_t AliTOFChannelOnlineStatusArray::GetHWStatus(Int_t pos) const 
{
	// getting the HW status for channel at position = pos 
	UChar_t parr = 0x0; 
	if  (pos>-1 && pos < fSize)parr = fArray[pos];
	AliDebug(2,Form("parr = %d ",(UInt_t)parr));
	UChar_t hwSt = parr & kTOFHW;
	//UChar_t hwSt = parr & 0x3;
	return hwSt;
}
//________________________________________________________________
UChar_t AliTOFChannelOnlineStatusArray::GetPulserStatus(Int_t pos) const 
{
	// getting the Pulser status for channel at position = pos 
	UChar_t parr = 0x0; 
	if  (pos>-1 && pos < fSize)parr = fArray[pos];
	AliDebug(2,Form("parr = %d ",(UInt_t)parr));
	UChar_t pulserSt = parr & kTOFPulser;
	//UChar_t pulserSt = parr & 0xc;
	return pulserSt;
    }
//________________________________________________________________
UChar_t AliTOFChannelOnlineStatusArray::GetNoiseStatus(Int_t pos) const 
{
	// getting the noise status for channel at position = pos 
	UChar_t parr = 0x0; 
	if  (pos>-1 && pos < fSize)parr = fArray[pos];
	AliDebug(2,Form("parr = %d ",(UInt_t)parr));
	UChar_t noiseSt = parr & kTOFNoise;
	//	UChar_t noiseSt = parr & 0x30;
	return noiseSt; 
}
