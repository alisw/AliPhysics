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
// class for TOF Online calibration: defining channel delay                  //
// using an array instead of a TObjArray                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <AliTOFChannelOnlineArray.h>
#include <AliLog.h>

ClassImp(AliTOFChannelOnlineArray)

//________________________________________________________________
AliTOFChannelOnlineArray::AliTOFChannelOnlineArray():
	TObject(),
	fSize(0),
	fArray(0x0)
{
	//default constructor
}
//________________________________________________________________
AliTOFChannelOnlineArray::~AliTOFChannelOnlineArray()
{
	//distructor
	delete [] fArray;
}
//________________________________________________________________
AliTOFChannelOnlineArray::AliTOFChannelOnlineArray(Int_t size):
	TObject(),
	fSize(size),
	fArray(new Float_t[size])
{
	// ctor with size
	for (Int_t ich = 0; ich<size; ich ++){
	  SetDelay(ich,0);
	}
}
//________________________________________________________________
AliTOFChannelOnlineArray::AliTOFChannelOnlineArray(const AliTOFChannelOnlineArray & source):
  TObject(source),
  fSize(source.fSize),
  fArray(source.fArray)
{ 
	// copy constructor
}
//________________________________________________________________
AliTOFChannelOnlineArray &AliTOFChannelOnlineArray::operator=(const AliTOFChannelOnlineArray & source) 
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
void AliTOFChannelOnlineArray::SetDelay(Int_t pos, Float_t parr)
{
	// setting status for channel at position = pos
	AliDebug(2,Form("status = %d",(Float_t)parr));
	if (pos>-1 && pos < fSize)fArray[pos] = parr;
	AliDebug(2,Form("fArray[%d] = %d",pos,(Float_t)fArray[pos]));
}
//________________________________________________________________
Float_t AliTOFChannelOnlineArray::GetDelay(Int_t pos) const 
{
	// getting the status for channel at position = pos 
	Float_t parr = 0x0; 
	if  (pos>-1 && pos < fSize)parr = fArray[pos];
	return parr;
}
