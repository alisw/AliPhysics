#ifndef ALITOFARRAY_H
#define ALITOFARRAY_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// ========================================================================
// Class to store variable size arrays of Float_t
// ========================================================================

class TObject;
class TArrayF;
class TCollection;

class AliTOFArray : public TObject {
  public:
	AliTOFArray(): TObject(),fSize(0),fArray(0x0){}
	AliTOFArray(Int_t size); 
	AliTOFArray(const AliTOFArray & source);
	AliTOFArray& operator=(const AliTOFArray & source);
	Int_t GetSize() const {return fSize;}
	void SetArray(Int_t pos, Int_t size=0);
	void SetAt(Int_t pos, Int_t nelements, Float_t* content);
	void SetAt(Int_t pos, Int_t ielement, Float_t content);
	void RemoveArray(Int_t pos);
	Float_t* GetArray(Int_t pos);
	Float_t GetArrayAt(Int_t pos, Int_t ielement);
	Int_t GetArraySize(Int_t pos);
	void ReSetArraySize(Int_t pos, Int_t size);
	virtual Long64_t Merge(TCollection *list);
	virtual ~AliTOFArray();

 private:
    Int_t fSize;       // Size of the array of TArrayFs
    TArrayF ** fArray; //[fSize]

    ClassDef(AliTOFArray,1)
};
#endif
