#ifndef ALI_DCS_ARRAY_H
#define ALI_DCS_ARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// This class represents the value(s) of a the LHC DPs at a given timestamp   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TObjArray.h>
//#include <TTimeStamp.h>
class TObjString;

class AliDCSArray : public TObject {
 public:
	enum Type {
		kInvalid = 0,
		kBool = 1,
		kChar = 2,
		kInt = 3,
		kUInt = 4,
		kFloat = 5,
		kString = 6,
		kDouble = 7
	};
	
	AliDCSArray();
	AliDCSArray(const AliDCSArray& c);
	
	virtual ~AliDCSArray();
	
	AliDCSArray& operator=(const AliDCSArray& c);
	
	AliDCSArray(Int_t nentries, Bool_t* value, Double_t timeStamp);
	AliDCSArray(Int_t nentries, Char_t* value, Double_t timeStamp);
	AliDCSArray(Int_t nentries, Int_t* value, Double_t timeStamp);
	AliDCSArray(Int_t nentries, UInt_t* value, Double_t timeStamp);
	AliDCSArray(Int_t nentries, Float_t* value, Double_t timeStamp);
	AliDCSArray(Int_t nentries, Double_t* value, Double_t timeStamp);
	AliDCSArray(Int_t nentries, TObjArray* value, Double_t timeStamp);
	
	Int_t GetNEntries() const { return fnentries;}
	Bool_t* GetBool() const { return fBool; }
	Char_t* GetChar() const { return fChar; }
	Int_t* GetInt() const { return fInt; }
	UInt_t* GetUInt() const { return fUInt; }
	Float_t* GetFloat() const { return fFloat; }
	Double_t* GetDouble() const { return fDouble; }
	TObjArray* GetStringArray() const { return fStringArray; }

	Bool_t GetBool(Int_t index) const { return fBool[index]; }
	Char_t GetChar(Int_t index) const { return fChar[index]; }
	Int_t GetInt(Int_t index) const { return fInt[index]; }
	UInt_t GetUInt(Int_t index) const { return fUInt[index]; }
	Float_t GetFloat(Int_t index) const { return fFloat[index]; }
	Double_t GetDouble(Int_t index) const { return fDouble[index]; }
	TObjString* GetStringArray(Int_t index) const { return (TObjString*)fStringArray->At(index); }
	
	Type GetType() const { return fType; }
	
	Double_t GetTimeStamp() const { return fTimeStamp; }
	void SetTimeStamp(Double_t timeStamp) { fTimeStamp = timeStamp; }
	
 protected:
	
	void Init();
	
	Type fType;                // type of the value stored
	
	Int_t fnentries;           // n. of entries at the same timestamp
	Bool_t* fBool;             //[fnentries] bool value
	Char_t* fChar;             //[fnentries] char value
	Int_t* fInt;               //[fnentries] int value
	UInt_t* fUInt;             //[fnentries] uint value
	Float_t* fFloat;           //[fnentries] float value
	//	TString* fString;   //[fnentries] string value
	TObjArray* fStringArray;    //TObjArray to store TObjStrinf for string value
	
	Double_t fTimeStamp;    // timestamp of this value

	Double_t* fDouble;         //[fnentries] double value
	
	ClassDef(AliDCSArray, 2);
};

#endif
