#ifndef ALI_SIMPLE_VALUE_H

#define ALI_SIMPLE_VALUE_H



/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *

 * See cxx source for full Copyright notice                               */



/* $Id$ */



//

// This class is a simple wrapper of

// all primitive types used in PVSS SCADA system.

// bool, char, byte, int ,uint, float, and arrays of these values

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



	class AliBoolHolder: public TObject {

	

	public:

		AliBoolHolder() {}

		virtual ~AliBoolHolder() {}



		AliBoolHolder(Bool_t val):fValue(val) {}

		

		Bool_t GetValue() const {return fValue;}

		void SetValue(const Bool_t val) {fValue = val;}

		

		virtual TObject* Clone(const char* name) const; 

	

		virtual Bool_t IsEqual(const TObject* object) const;



	private:

		Bool_t fValue;		// The value



		ClassDef(AliBoolHolder, 1);

	};



	class AliByteHolder: public TObject {

	

	public:



		AliByteHolder() {};

		virtual ~AliByteHolder() {}



		AliByteHolder(Char_t val):fValue(val) {}

		

		Char_t GetValue() const {return fValue;}

		void SetValue(const Char_t val) {fValue = val;}

		

		virtual TObject* Clone(const char* name) const;



		virtual Bool_t IsEqual(const TObject* object) const;



	private:

		Char_t fValue;		// The value



		ClassDef(AliByteHolder, 1);

	};



	class AliIntHolder: public TObject {

	

	public:



		AliIntHolder() {}

		virtual ~AliIntHolder() {}



		AliIntHolder(Int_t val):fValue(val) {}

		

		Int_t GetValue() const {return fValue;}

		void SetValue(const Int_t val) {fValue = val;}

		

		virtual TObject* Clone(const char* name) const; 



		virtual Bool_t IsEqual(const TObject* object) const;





	private:

		Int_t fValue;		// The value



		ClassDef(AliIntHolder, 1);

	};



	class AliUIntHolder: public TObject {

	

	public:



		AliUIntHolder() {}

		virtual ~AliUIntHolder() {}



		AliUIntHolder(UInt_t val):fValue(val) {}

		

		UInt_t GetValue() const {return fValue;}

		void SetValue(const UInt_t val) {fValue = val;}

		

		virtual TObject* Clone(const char* name) const;



		virtual Bool_t IsEqual(const TObject* object) const;





	private:

		UInt_t fValue;		// The value



		ClassDef(AliUIntHolder, 1);

	};



	class AliFloatHolder: public TObject {

	

	public:



		AliFloatHolder() {}

		virtual ~AliFloatHolder() {}

	

		AliFloatHolder(Float_t val):fValue(val) {}

		

		Float_t GetValue() const {return fValue;}

		void SetValue(const Float_t val) {fValue = val;}

		

		virtual TObject* Clone(const char* name) const;



		virtual Bool_t IsEqual(const TObject* object) const;



	private:

		Float_t fValue;		// The value



		ClassDef(AliFloatHolder, 1);

	};



	class AliDynHolder: public TObject {

	

	public:



		AliDynHolder(): fSize(0) {}

		AliDynHolder(Int_t size): fSize(size){}

		

		Int_t GetSize() const {return fSize;}

		void  SetSize(const Int_t n) {fSize = n;}

		

	protected:

		Int_t fSize;		// The size



		ClassDef(AliDynHolder, 0);

	};



	class AliDynBoolHolder: public AliDynHolder {

	

	public:



		AliDynBoolHolder(): fValues(NULL) {}



		AliDynBoolHolder(Int_t size, const Bool_t* buf = NULL):

			AliDynHolder(size) {

			fValues = new Bool_t[size];

			if (buf) memcpy(fValues, buf, size * sizeof(Bool_t));

		}



		virtual ~AliDynBoolHolder() {if (fValues) delete[] fValues;}

		

		Bool_t GetValue(Int_t n) const {return fValues[n];}

		Bool_t* GetValues() const {return fValues;}

		void SetValue(const Int_t pos, const Bool_t val) {fValues[pos] = val;}

		

		virtual TObject* Clone(const char* name) const;



                virtual Bool_t IsEqual(const TObject* object) const; 





	private:

		AliDynBoolHolder(const AliDynBoolHolder& /*other*/): AliDynHolder() { }

		AliDynBoolHolder& operator= (const AliDynBoolHolder& /*other*/) {return *this;}

		Bool_t* fValues; //[fSize] 	 The value



		ClassDef(AliDynBoolHolder, 1);

	};



	class AliDynByteHolder: public AliDynHolder {

	

	public:

	

		AliDynByteHolder(): fValues(NULL) {}



		AliDynByteHolder(Int_t size, const Char_t* buf = NULL):

			AliDynHolder(size) {

			fValues = new Char_t[size];

			if (buf) memcpy(fValues, buf, size * sizeof(Char_t));

		}



		virtual ~AliDynByteHolder() {if (fValues) delete[] fValues;}

		

		Char_t GetValue(Int_t n) const {return fValues[n];}

		Char_t* GetValues() const {return fValues;}

		void SetValue(const Int_t pos, const Char_t val) {fValues[pos] = val;}

		

		virtual TObject* Clone(const char* name) const;



                virtual Bool_t IsEqual(const TObject* object) const;





	private:

		AliDynByteHolder(const AliDynByteHolder& /*other*/): AliDynHolder() { }

		AliDynByteHolder& operator= (const AliDynByteHolder& /*other*/) {return *this;}

		Char_t* fValues; //[fSize] 	 The value



		ClassDef(AliDynByteHolder, 1);

	};



	class AliDynIntHolder: public AliDynHolder {

	

	public:



		AliDynIntHolder(): fValues(NULL) {}



		AliDynIntHolder(Int_t size, const Int_t* buf = NULL):

			AliDynHolder(size) {

			fValues = new Int_t[size];

			if (buf) memcpy(fValues, buf, size * sizeof(Int_t));

		}



		virtual ~AliDynIntHolder() {if (fValues) delete[] fValues;}

		

		Int_t GetValue(Int_t n) const {return fValues[n];}

		Int_t* GetValues() const {return fValues;}

		void SetValue(const Int_t pos, const Int_t val) {fValues[pos] = val;}

		

		virtual TObject* Clone(const char* name) const;



                virtual Bool_t IsEqual(const TObject* object) const;





	private:

		AliDynIntHolder(const AliDynIntHolder& /*other*/): AliDynHolder() { }

		AliDynIntHolder& operator= (const AliDynIntHolder& /*other*/) {return *this;}

		Int_t* fValues; //[fSize] 	 The value



		ClassDef(AliDynIntHolder, 1);

	};



	class AliDynUIntHolder: public AliDynHolder {

	

	public:



		AliDynUIntHolder(): fValues(NULL) {}



		AliDynUIntHolder(Int_t size, const UInt_t* buf = NULL):

			AliDynHolder(size) {

			fValues = new UInt_t[size];

			if (buf) memcpy(fValues, buf, size * sizeof(UInt_t));

		} 



		virtual ~AliDynUIntHolder() {if (fValues) delete[] fValues;}

		

		UInt_t GetValue(Int_t n) const {return fValues[n];}

		UInt_t* GetValues() const {return fValues;}

		void SetValue(const Int_t pos, const UInt_t val) {fValues[pos] = val;}

		

		virtual TObject* Clone(const char* name) const;



                virtual Bool_t IsEqual(const TObject* object) const;





	private:

		AliDynUIntHolder(const AliDynUIntHolder& /*other*/): AliDynHolder() { }

		AliDynUIntHolder& operator= (const AliDynUIntHolder& /*other*/) {return *this;}

		UInt_t* fValues; //[fSize] 	 The value



		ClassDef(AliDynUIntHolder, 1);

	};



	class AliDynFloatHolder: public AliDynHolder {

	

	public:



		AliDynFloatHolder(): fValues(NULL) {}



		AliDynFloatHolder(Int_t size, const Float_t* buf = NULL):

			AliDynHolder(size) {

			fValues = new Float_t[size];

			if (buf) memcpy(fValues, buf, size * sizeof(Float_t));

				

		}



		virtual ~AliDynFloatHolder() {if (fValues) delete[] fValues;}

		

		Float_t GetValue(Int_t n) const {return fValues[n];}

		Float_t* GetValues() const {return fValues;}

		void SetValue(const Int_t pos, const Float_t val) {fValues[pos] = val;}

		

		virtual TObject* Clone(const char* name) const;



                virtual Bool_t IsEqual(const TObject* object) const;





	private:

		AliDynFloatHolder(const AliDynFloatHolder& /*other*/): AliDynHolder() { }

		AliDynFloatHolder& operator= (const AliDynFloatHolder& /*other*/) {return *this;}

		Float_t* fValues; //[fSize] 	 The value



		ClassDef(AliDynFloatHolder, 1);

	};





	TObject* fHolder;  	// Holder of the basic type value



	Type fType;  	// Type of the basic type value





	Bool_t TypeOk(Type type) const;



	Bool_t BoundsOk(Int_t n) const;





	ClassDef(AliSimpleValue, 1);

};



#endif

