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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This class represents the value(s) of a the LHC DPs at a given timestamp   //
// The variuos measurement that occurred at the same timestamp are associated //
// to the same timestamp.                                                     //  
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "AliDCSArray.h"

#include "TTimeStamp.h"
#include <TString.h>

ClassImp(AliDCSArray)
	
AliDCSArray::AliDCSArray() :
		TObject(),
		fType(kInvalid),
		fnentries(0),
		fBool(0x0),
		fChar(0x0),
		fInt(0x0),
		fUInt(0x0),
		fFloat(0x0),
		fString(0x0),
		fTimeStamp(0)
{
	//
	// default constructor
	//
}

//--------------------------------------------------------------------------
AliDCSArray::AliDCSArray(Int_t nentries, Bool_t* value, TTimeStamp* timeStamp) : 
	TObject(),
	fType(kBool),
	fnentries(nentries),
	fBool(new Bool_t[fnentries]),
	fChar(0x0),
	fInt(0x0),
	fUInt(0x0),
	fFloat(0x0),
	fString(0x0),
	fTimeStamp(timeStamp)
{
	//
	// constructor for Bool
	//
	
	for (Int_t i = 0; i<fnentries; i++){
		fBool[i] = value[i];
	}
}

//--------------------------------------------------------------------------
AliDCSArray::AliDCSArray(Int_t nentries, Char_t* value, TTimeStamp* timeStamp) :
	TObject(),
	fType(kChar),
	fnentries(nentries),
	fBool(0x0),
	fChar(new Char_t[fnentries]),
	fInt(0x0),
	fUInt(0x0),
	fFloat(0x0),
	fString(0x0),
	fTimeStamp(timeStamp)
{
	//
	// constructor for Char
	//
	for (Int_t i = 0; i<fnentries; i++){
		fChar[i] = value[i];
	}
}

//-------------------------------------------------------------------------
AliDCSArray::AliDCSArray(Int_t nentries, Int_t* value, TTimeStamp* timeStamp) :
	TObject(),
	fType(kInt),
	fnentries(nentries),
	fBool(0x0),
	fChar(0x0),
	fInt(new Int_t[fnentries]),
	fUInt(0x0),
	fFloat(0x0),
	fString(0x0),
	fTimeStamp(timeStamp)
{
	//
	// constructor for Int
	//
	for (Int_t i = 0; i<fnentries; i++){
		fInt[i] = value[i];
	}
}

//-------------------------------------------------------------------------
AliDCSArray::AliDCSArray(Int_t nentries, UInt_t* value, TTimeStamp* timeStamp) :
	TObject(),
	fType(kUInt),
	fnentries(nentries),
	fBool(0x0),
	fChar(0x0),
	fInt(0x0),
	fUInt(new UInt_t[fnentries]),
	fFloat(0x0),
	fString(0x0),
	fTimeStamp(timeStamp)
{
	//
	// constructor for UInt
	//

	for (Int_t i = 0; i<fnentries; i++){
		fUInt[i] = value[i];
	}
}

//-------------------------------------------------------------------------
AliDCSArray::AliDCSArray(Int_t nentries, Float_t* value, TTimeStamp* timeStamp) :
	TObject(),
	fType(kFloat),
	fnentries(nentries),
	fBool(0x0),
	fChar(0x0),
	fInt(0x0),
	fUInt(0x0),
	fFloat(new Float_t[fnentries]),
	fString(0x0),
	fTimeStamp(timeStamp)
{
	//
	// constructor for Float
	//

	for (Int_t i = 0; i<fnentries; i++){
		fFloat[i] = value[i];
	}
}

//------------------------------------------------------------------------
AliDCSArray::AliDCSArray(Int_t nentries, TString* value, TTimeStamp* timeStamp) :
	TObject(),
	fType(kString),
	fnentries(nentries),
	fBool(0x0),
	fChar(0x0),
	fInt(0x0),
	fUInt(0x0),
	fFloat(0x0),
	fString(new TString[fnentries]),
	fTimeStamp(timeStamp)
{
	//
	// constructor for String
	//

	for (Int_t i = 0; i<fnentries; i++){
		fString[i] = value[i];
	}
}

//-----------------------------------------------------------------------
AliDCSArray::AliDCSArray(const AliDCSArray& c) :
	TObject(c),
	fType(c.fType),
	fnentries(c.fnentries),
	fBool(0x0),
	fChar(0x0),
	fInt(0x0),
	fUInt(0x0),
	fFloat(0x0),
	fString(0x0),
	fTimeStamp(c.fTimeStamp)
{
	//
	// copy constructor
	//
	
	if (fType == kBool && c.fBool){
		fBool = new Bool_t[fnentries];
		memcpy(fBool,c.fBool,fnentries*sizeof(Bool_t));
	} 
	if (fType == kChar && c.fChar){
		fChar = new Char_t[fnentries];
		memcpy(fChar,c.fChar,fnentries*sizeof(Char_t));
	} 
 	if (fType == kUInt && c.fUInt){
		fUInt = new UInt_t[fnentries];
		memcpy(fUInt,c.fUInt,fnentries*sizeof(UInt_t));
	} 
 	if (fType == kInt && c.fInt){
		fInt = new Int_t[fnentries];
		memcpy(fInt,c.fInt,fnentries*sizeof(Int_t));
	} 
 	if (fType == kFloat && c.fFloat){
		fFloat = new Float_t[fnentries];
		memcpy(fFloat,c.fFloat,fnentries*sizeof(Float_t));
	} 
 	if (fType == kString && c.fString){
		fString = new TString[fnentries];
		memcpy(fString,c.fString,fnentries*sizeof(TString));
	} 

}

//-----------------------------------------------------------------------
AliDCSArray::~AliDCSArray()
{
	//
	// destructor
	//

	if (fBool){
		fBool = 0x0;
		delete fBool;
	}
	if (fChar){
		fChar = 0x0;
		delete fChar;
	}
	if (fUInt){
		fUInt = 0x0;
		delete fUInt;
	}
	if (fInt){
		fInt = 0x0;
		delete fInt;
	}
	if (fFloat){
		fFloat = 0x0;
		delete fFloat;
	}
	if (fString){
		fString = 0x0;
		delete fString;
	}
	if (fTimeStamp){
		fTimeStamp = 0x0;
		delete fTimeStamp;
	}
}

//-----------------------------------------------------------------------
AliDCSArray &AliDCSArray::operator=(const AliDCSArray &c)
{
	//
	// operator = 
	//
	
	if(&c == this) return *this;
	if (fType == kBool && c.fBool){
		fBool = new Bool_t[fnentries];
		memcpy(fBool,c.fBool,fnentries*sizeof(Bool_t));
	} 
	if (fType == kChar && c.fChar){
		fChar = new Char_t[fnentries];
		memcpy(fChar,c.fChar,fnentries*sizeof(Char_t));
	} 
 	if (fType == kUInt && c.fUInt){
		fUInt = new UInt_t[fnentries];
		memcpy(fUInt,c.fUInt,fnentries*sizeof(UInt_t));
	} 
 	if (fType == kInt && c.fInt){
		fInt = new Int_t[fnentries];
		memcpy(fInt,c.fInt,fnentries*sizeof(Int_t));
	} 
 	if (fType == kFloat && c.fFloat){
		fFloat = new Float_t[fnentries];
		memcpy(fFloat,c.fFloat,fnentries*sizeof(Float_t));
	} 
 	if (fType == kString && c.fString){
		fString = new TString[fnentries];
		memcpy(fString,c.fString,fnentries*sizeof(TString));
	} 
	
	return *this;
}

//-----------------------------------------------------------------------
void AliDCSArray::Init()
{
	//
	// init helper, that initializes everything to 0
	//
	
	fType = kInvalid;
	
	fnentries = 0;
	fBool = 0x0;
	fChar = 0x0;
	fInt = 0x0;
	fUInt = 0x0;
	fFloat = 0x0;
	fString = 0x0;
	
	fTimeStamp = 0x0;
}
