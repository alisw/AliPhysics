#ifndef ALI_SIMPLE_VALUE_H
#define ALI_SIMPLE_VALUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This class is a simple wrapper of
// all primitive types used in PVSS SCADA system.
//

#include <TObject.h>
#include <TString.h>

class AliSimpleValue: public TObject {
public:
	enum Type {
		kInvalid = 0,
		kBool = 1,
		kByte = 2,
		kInt = 3,
		kUInt = 4,
		kFloat = 5,
		kDynBool = 11,
		kDynByte = 12,
		kDynInt = 13,
		kDynUInt = 14,
		kDynFloat = 15
	};	


	AliSimpleValue();

        AliSimpleValue(const AliSimpleValue& other);

        AliSimpleValue(Type type, Int_t size = 0);

        AliSimpleValue(Bool_t val);

        AliSimpleValue(Char_t val);

        AliSimpleValue(Int_t val);

        AliSimpleValue(UInt_t val);

        AliSimpleValue(Float_t val);

        AliSimpleValue(Int_t size, const Bool_t* vals);

        AliSimpleValue(Int_t size, const Char_t* vals);

        AliSimpleValue(Int_t size, const Int_t* vals);

        AliSimpleValue(Int_t size, const UInt_t* vals);

        AliSimpleValue(Int_t size, const Float_t* vals);

        ~AliSimpleValue();


        AliSimpleValue& operator=(const AliSimpleValue& other);

        Bool_t operator==(const AliSimpleValue& other) const;

        void SetBool(Bool_t val);

        void SetByte(Char_t val);

        void SetInt(Int_t val);

        void SetUInt(UInt_t val);

        void SetFloat(Float_t val);


        Bool_t GetBool() const;

        Char_t GetByte() const;

        Int_t GetInt() const;

        UInt_t GetUInt() const;

        Float_t GetFloat() const;


	void SetDynBool(Int_t n, Bool_t val);

        void SetDynByte(Int_t n, Char_t val);

        void SetDynInt(Int_t n, Int_t val);

        void SetDynUInt(Int_t n, UInt_t val);

        void SetDynFloat(Int_t n, Float_t val);


        Bool_t GetDynBool(Int_t n) const;

        Char_t GetDynByte(Int_t n) const;

        Int_t GetDynInt(Int_t n) const;

        UInt_t GetDynUInt(Int_t n) const;

        Float_t GetDynFloat(Int_t n) const;


        Type GetType() const {return fType;};

        Int_t GetSize() const;

        Int_t GetDynamicSize() const;

        TString ToString() const;


        static Bool_t IsDynamic(Type type);

        static Int_t GetPrimitiveSize(Type type);

        static const char* GetTypeString(Type type);

private:

	class BoolHolder: public TObject {
	public:
		Bool_t fValue;

		BoolHolder() {}

		BoolHolder(Bool_t val):fValue(val) {}
		
		virtual TObject* Clone(const char* name) const; 
	
		virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(BoolHolder, 1);
	};

	class ByteHolder: public TObject {
	public:
		Char_t fValue;

		ByteHolder() {};

		ByteHolder(Char_t val):fValue(val) {}
		
		virtual TObject* Clone(const char* name) const;

		virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(ByteHolder, 1);
	};

	class IntHolder: public TObject {
	public:
		Int_t fValue;

		IntHolder() {}

		IntHolder(Int_t val):fValue(val) {}
		
		virtual TObject* Clone(const char* name) const; 

		virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(IntHolder, 1);
	};

	class UIntHolder: public TObject {
	public:
		UInt_t fValue;

		UIntHolder() {}

		UIntHolder(UInt_t val):fValue(val) {}
		
		virtual TObject* Clone(const char* name) const;

		virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(UIntHolder, 1);
	};

	class FloatHolder: public TObject {
	public:
		Float_t fValue;

		FloatHolder() {}
	
		FloatHolder(Float_t val):fValue(val) {}
		
		virtual TObject* Clone(const char* name) const;

		virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(FloatHolder, 1);
	};

	class DynHolder: public TObject {
	public:
		Int_t fSize;

		DynHolder(): fSize(0) {}
		DynHolder(Int_t size): fSize(size){}

		ClassDef(DynHolder, 0);
	};

	class DynBoolHolder: public DynHolder {
	public:
		Bool_t* fValues; //[fSize]

		DynBoolHolder(): fValues(NULL) {}

		DynBoolHolder(Int_t size, const Bool_t* buf = NULL):
			DynHolder(size) {
			fValues = new Bool_t[size];
			if (buf) memcpy(fValues, buf, size * sizeof(Bool_t));
		}

		virtual ~DynBoolHolder() {if (fValues) delete[] fValues;}

		virtual TObject* Clone(const char* name) const;

                virtual Bool_t IsEqual(const TObject* object) const; 

		ClassDef(DynBoolHolder, 1);
	};

	class DynByteHolder: public DynHolder {
	public:
		Char_t* fValues; //[fSize]
	
		DynByteHolder(): fValues(NULL) {}

		DynByteHolder(Int_t size, const Char_t* buf = NULL):
			DynHolder(size) {
			fValues = new Char_t[size];
			if (buf) memcpy(fValues, buf, size * sizeof(Char_t));
		}

		virtual ~DynByteHolder() {if (fValues) delete[] fValues;}

		virtual TObject* Clone(const char* name) const;

                virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(DynByteHolder, 1);
	};

	class DynIntHolder: public DynHolder {
	public:
		Int_t* fValues; //[fSize]

		DynIntHolder(): fValues(NULL) {}

		DynIntHolder(Int_t size, const Int_t* buf = NULL):
			DynHolder(size) {
			fValues = new Int_t[size];
			if (buf) memcpy(fValues, buf, size * sizeof(Int_t));
		}

		virtual ~DynIntHolder() {if (fValues) delete[] fValues;}

		virtual TObject* Clone(const char* name) const;

                virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(DynIntHolder, 1);
	};

	class DynUIntHolder: public DynHolder {
	public:
		UInt_t* fValues; //[fSize]

		DynUIntHolder(): fValues(NULL) {}

		DynUIntHolder(Int_t size, const UInt_t* buf = NULL):
			DynHolder(size) {
			fValues = new UInt_t[size];
			if (buf) memcpy(fValues, buf, size * sizeof(UInt_t));	
		} 

		virtual ~DynUIntHolder() {if (fValues) delete[] fValues;}

		virtual TObject* Clone(const char* name) const;

                virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(DynUIntHolder, 1);
	};

	class DynFloatHolder: public DynHolder {
	public:
		Float_t* fValues; //[fSize]

		DynFloatHolder(): fValues(NULL) {}

		DynFloatHolder(Int_t size, const Float_t* buf = NULL):
			DynHolder(size) {
			fValues = new Float_t[size];
			if (buf) memcpy(fValues, buf, size * sizeof(Float_t));
		}

		virtual ~DynFloatHolder() {if (fValues) delete[] fValues;}

		virtual TObject* Clone(const char* name) const;

                virtual Bool_t IsEqual(const TObject* object) const;

		ClassDef(DynFloatHolder, 1);
	};


	TObject* fHolder;

	Type fType;


	Bool_t TypeOk(Type type) const;

	Bool_t BoundsOk(Int_t n) const;


	ClassDef(AliSimpleValue, 1);
};

#endif
