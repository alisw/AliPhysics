#ifndef ALITRDTRAPALU_H
#define ALITRDTRAPALU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDtrapAlu.h 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRAP-ALU implementation                                               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliTRDtrapAlu:public TObject {
 
 public:
   
               AliTRDtrapAlu();
	       //AliTRDtrapAlu(AliTRDtrapAlu& bin); //copy constructor
  virtual      ~AliTRDtrapAlu();


  void Init(const Int_t& precom=10, const Int_t& postcom=2
          , const Int_t& lRestriction = -1, const Int_t& uRestriction = -1);

  Int_t   GetValue () const { 
    // return the value 
    return fValue;
  }

  Int_t   GetSignedValue ()const{  
    // return the value with its sign
    if(fSigned == kFALSE) return fValue;
    else return fValue*(-1);
  }

  Int_t   GetValuePre ()const{
    // return value of pre-comma part as integer
    Int_t valPre = fValue>>fPostCom;
    return valPre;   
  }

  Double_t GetValueWhole() const;  

  Int_t   GetPre()const{
    // return nr of pre-comma bits
    return fPreCom;
  }

  Int_t   GetPost()const{
    // return nr of past-comma bits
    return fPostCom;
  }

  Bool_t  GetSign()const{
    // return true if signed
    if(fSigned == kTRUE) return kTRUE;
    return kFALSE;
  }
   
  Bool_t  CheckUSize(const Int_t& val)const{
    // compare value to the upper restriction
    if(val>fuRestriction) return kFALSE;
    return kTRUE;
  }

  Bool_t  CheckLSize(const Int_t& val)const{
    // compare value to the lower restriction
    if(val<flRestriction) return kFALSE;
    return kTRUE;
  }

  void AssignFormatted(const Int_t& formVal){ 
    // assign a value with proper format; assigns formVal directly to fValue; better not use explicitely
    fValue  = formVal;
    //fValue  = fValue & (LUT(fPreCom + fPostCom) - 1); // no cut-off wanted
  } 

  void SetSign(const Int_t& s){
    // sets the sign
    if(s >= 0) fSigned = kFALSE;
    if(s <  0) fSigned = kTRUE;
  }

  void   WriteWord();  

  AliTRDtrapAlu& AssignInt(const  Int_t& first);     // in case a decimal integer is assigned to a binary; 
  AliTRDtrapAlu& AssignDouble(const  Double_t& first);  // change "Double_t" into "Float_t"
  AliTRDtrapAlu& operator=(const AliTRDtrapAlu& binary);
       
  AliTRDtrapAlu operator+(const AliTRDtrapAlu& binary); //binary is not const, because in a+(b*c) binary is reference to the object, to which Mem() is also a reference and this object is changed
  AliTRDtrapAlu operator-(const AliTRDtrapAlu& binary);
  AliTRDtrapAlu operator*(const AliTRDtrapAlu& binary);
  AliTRDtrapAlu operator/(const AliTRDtrapAlu& binary);

 protected:

  // void FastInit(const Int_t& precom = 10, const Int_t&  postcom = 2, const Int_t& formVal = 0); //meant to combine definition of format with integer-value assignment; not to apply by user  

  //the following two functions encapsulate global static members; can only be changed by member functions (visibility only inside class)

  Int_t  MakePower(const Int_t& base=1,const  Int_t& exponent=1)const;

  /*static AliTRDtrapAlu& Mem() { 
    // a global instance of the class, which is only defined once
    static AliTRDtrapAlu fAuxiliary;
    return fAuxiliary;
  }*/

  static Int_t LUT(const Int_t& index);
       
  const Int_t&  Min(const Int_t& comp1, const Int_t& comp2)const{
    // return the minimum
    if (comp1 <= comp2) return comp1;
    return comp2;
  }

  const Int_t&  Max(const Int_t& comp1, const Int_t& comp2)const{
    // return the maximum
    if (comp1 >= comp2) return comp1;
    return comp2;
  }

  //static AliTRDtrapAlu fAlu;

  Int_t    fValue;          // the value in integers
  Int_t    fPreCom;         // number of pre-comma bits
  Int_t    fPostCom;        // number of past-comma bits
  Int_t    fuRestriction;   // the upper restriction for the value
  Int_t    flRestriction;   // the lower restriction for the value
  Bool_t   fSigned;         // signed value? 

  ClassDef(AliTRDtrapAlu,1) // TRAP-ALU

};
#endif


