#ifndef ALI_DCS_VALUE_H
#define ALI_DCS_VALUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class represents the value(s) of a DCS data point at a given timestamp
//

#include <TObject.h>

class AliDCSValue : public TObject {
public:
  enum Type {
    kInvalid = 0,

    kBool = 1,
    kChar = 2,
    kInt = 3,
    kUInt = 4,
    kFloat = 5,

    kDynBool = 11,
    kDynChar = 12,
    kDynInt = 13,
    kDynUInt = 14,
    kDynFloat = 15
  };

  AliDCSValue();
  AliDCSValue(const AliDCSValue& c);

  ~AliDCSValue();

  AliDCSValue& operator=(const AliDCSValue& c);
  virtual void Copy(TObject& c) const;

  AliDCSValue(Bool_t value, UInt_t timeStamp);
  AliDCSValue(Char_t value, UInt_t timeStamp);
  AliDCSValue(Int_t value, UInt_t timeStamp);
  AliDCSValue(UInt_t value, UInt_t timeStamp);
  AliDCSValue(Float_t value, UInt_t timeStamp);

  AliDCSValue(Int_t size, Bool_t* vals, UInt_t timeStamp);
  AliDCSValue(Int_t size, Char_t* vals, UInt_t timeStamp);
  AliDCSValue(Int_t size, Int_t* vals, UInt_t timeStamp);
  AliDCSValue(Int_t size, UInt_t* vals, UInt_t timeStamp);
  AliDCSValue(Int_t size, Float_t* vals, UInt_t timeStamp);

  Bool_t GetBool() const { return fBool; }
  Char_t GetChar() const { return fChar; }
  Int_t GetInt() const { return fInt; }
  UInt_t GetUInt() const { return fUInt; }
  Float_t GetFloat() const { return fFloat; }

  Bool_t GetDynBool(Int_t n) const { return fBoolPtr[n]; }
  Char_t GetDynChar(Int_t n) const { return fCharPtr[n]; }
  Int_t GetDynInt(Int_t n) const { return fIntPtr[n]; }
  UInt_t GetDynUInt(Int_t n) const { return fUIntPtr[n]; }
  Float_t GetDynFloat(Int_t n) const { return fFloatPtr[n]; }

  Type GetType() const { return fType; }
  Int_t GetDynamicSize() const { return fLength; }

  UInt_t GetTimeStamp() const { return fTimeStamp; }
  void SetTimeStamp(UInt_t timeStamp) { fTimeStamp = timeStamp; }

  Int_t GetSize() const;
  static Bool_t IsDynamic(Type type);

  const Char_t* ToString() const;

protected:
  void Init();

  Type fType;           // type of the value stored

  Bool_t fBool;         // bool value
  Char_t fChar;         // char value
  Int_t fInt;           // int value
  UInt_t fUInt;         // uint value
  Float_t fFloat;       // float value

  UInt_t fLength;       // length of the following arrays, the ones that are != 0 are expected to contain fLength entries

  Bool_t* fBoolPtr;     //[fLength] bool array
  Char_t* fCharPtr;     //[fLength] char array
  Int_t* fIntPtr;       //[fLength] int array
  UInt_t* fUIntPtr;     //[fLength] uint array
  Float_t* fFloatPtr;   //[fLength] float array

  UInt_t fTimeStamp;    // timestamp of this value

	ClassDef(AliDCSValue, 2);
};

#endif
