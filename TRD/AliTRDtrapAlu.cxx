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

/* $Id: AliTRDtrapAlu.cxx 25891 2008-05-19 14:58:18Z fca $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRAP-ALU implementation                                                  //
//                                                                           //
//  Author:                                                                  //
//    Clemens Haltebourg <halteb@physi.uni-heidelberg.de>                    //
//                                                                           //
//  Usage of the class:                                                      //
//    Declaration of class instances: AliTRDtrapAlu a,b,c;                   //
//    Initialization:                 a.Init(2,11); b.Init(4,4); c.Init(5,4);//
//    Assigning values:               a.AssignDouble(5.7); b.AssignInt(3);   //
//    (you can also do b.AssignDouble(3) with same effect);                  //
//    Calculation:                     c = a*b;                              //
//    Test if c has right value:       c.WriteWord();                        //
//    Don't declare pointers; operators not overridden for pointer types;    //
//    You have to dereference yourself;                                      //
//    Use operators +,-,*,/ only with instances of the class; don't do       //
//    things like c=a*2 but rather b.AssignInt(2); c=a*b;                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTRDtrapAlu.h"

ClassImp(AliTRDtrapAlu)

//_____________________________________________________________________________  
AliTRDtrapAlu::AliTRDtrapAlu():TObject()

  ,fValue(0)
  ,fPreCom(0)
  ,fPostCom(0)
  ,fuRestriction(0)
  ,flRestriction(0)
  ,fSigned(kFALSE)

{
  
  // default constructor
  
}

//_____________________________________________________________________________  
AliTRDtrapAlu::~AliTRDtrapAlu(){
  //destructor
}

//_____________________________________________________________________________  
void AliTRDtrapAlu::Init(const Int_t& precom, const Int_t& postcom, const Int_t& lRestriction, const Int_t& uRestriction){
   // initialization: characterizes the bit-word (nr of pre- and post-comma bits, boundaries)
   fPostCom = postcom;
   fPreCom  = precom;
   fValue   = 0;       //currently, re-initialization kills the value
   fSigned  = kFALSE;

   if (fPreCom + fPostCom > 32 || fPreCom > 31) {fPreCom = 1; fPostCom = 0;return;} // prevent pre-comma part exceeding 31 spaces
   if (fPreCom  <= 0) {fPreCom  = 1;}
   if (fPostCom <  0) {fPostCom = 0;}
   
   Int_t lut = LUT(fPreCom + fPostCom)-1;
   if (uRestriction <= -1 || uRestriction > lut) {fuRestriction = lut;}
   else {fuRestriction = uRestriction;}
   if (lRestriction <= -1 || lRestriction > fuRestriction) {flRestriction = -lut;}
   else {flRestriction = lRestriction;}
   // up to now you can only choose a non-negative lower restriction (e.g. if you want your values to be >=0) ; can't deal with asymmetric borders; have to be implemented if needed
}

//_____________________________________________________________________________  
Double_t AliTRDtrapAlu::GetValueWhole() const { 
   // get the actual value (respecting pre- and post-comma parts) in integer-description
   Double_t valPre = (Double_t)(fValue>>fPostCom);
   Double_t valPost = 0.0;
   for(Int_t i = 0; i<=fPostCom-1; i++){
     Double_t num = (fValue>>i)&1;
     Double_t denom = LUT(fPostCom-i);
     valPost = valPost + num/denom;
   }
   Double_t val = valPre + valPost;
   return val;
 }

//_____________________________________________________________________________  
void AliTRDtrapAlu::WriteWord(){
  // for debugging purposes
  printf("bit-word: ");
  if (fSigned == true) printf("-");
  for(Int_t i = fPostCom + fPreCom - 1; i >= fPostCom; i--){  //read from behind in order to write the word from left to right
    printf("%d",(fValue>>i) & 1);
  }
  printf(".");
  for (Int_t j = fPostCom - 1; j >= 0; j--){
    printf("%d",(fValue>>j) & 1);
  }
  printf("\n");
         
}

//_____________________________________________________________________________  
AliTRDtrapAlu& AliTRDtrapAlu::AssignInt(const Int_t& first){  
  // assign an integer

  // parameter "first" is an integer for the pre-comma part (not UInt in order to match the error case first<0)
  fSigned = kFALSE;
  Int_t exponent = fPreCom + fPostCom;

    
  if (first<0) {
    fValue  = 0;                      //setting fValue to 0; first should not be negative
    fValue  = fValue & 0;
    return *this;
  }

  if (CheckUSize(first<<fPostCom) == kFALSE){
    
    //setting fValue to maximum; first was to big
    fValue  = fuRestriction;
    fValue  = fValue & (LUT(exponent)-1);
    return *this;
  }

  if (CheckLSize(first<<fPostCom) == kFALSE){
    
    //setting fValue to minimum; first was to small
    fValue  = flRestriction;
    fValue  = fValue & (LUT(exponent)-1);
    return *this;
  }

  
  fValue  = first;
  fValue  = fValue<<fPostCom; 
  fValue  = fValue & (LUT(exponent)-1);
 
  return *this;
    
}

//_____________________________________________________________________________  
AliTRDtrapAlu& AliTRDtrapAlu::AssignDouble(const  Double_t& first){
  // assign a double
 
  fSigned = kFALSE;
  Int_t exponent           = fPreCom + fPostCom;
  Int_t firstPre          = 0;  //integer part of first
  Int_t firstPost         = 0;  //comma part of first (cut off with enough accuracy
  Int_t c                  = 0;
  Double_t firstPreFloat = 0;
  
  
  Int_t power1 = LUT(exponent);
  
  firstPre       = (Int_t)first;
  firstPreFloat = firstPre;
  
  if(firstPre < 0){
    fValue  = 0;
    fValue = fValue & 0;
    return *this;
  }
  
  if(CheckUSize((Int_t)(first*LUT(fPostCom))) == kFALSE){
    
    //fValue  = MakePower(2,fPreCom) - 1;
    fValue  = fuRestriction;
    fValue  = fValue & (power1 - 1);
    return *this;
  }
  
  if(CheckLSize((Int_t)(first*LUT(fPostCom))) == kFALSE){
    
    //fValue  = MakePower(2,fPreCom) - 1;
    fValue  = flRestriction;
    fValue  = fValue & (power1 - 1);
    return *this;
  }
  

  fValue = firstPre;
  
  //get post comma part with adequate accuracy
  firstPost = (Int_t)((first - firstPreFloat)*LUT(fPostCom));
  for(Int_t i = 1; i <= fPostCom; i++) {
    c = (firstPost>>(fPostCom - i)) & 1;
    fValue  = fValue<<1;
    fValue  = fValue | c;
  }

  fValue = fValue & (power1 - 1);
  return *this;
}

//_____________________________________________________________________________  
AliTRDtrapAlu& AliTRDtrapAlu::operator=(const AliTRDtrapAlu& binary){
  // assign an object of type AliTRDtrapAlu

  Int_t c    = 0;
  //Int_t exponent = fPreCom + fPostCom;
  
  
  Int_t power1 = LUT(fPreCom + fPostCom);

  fValue          = binary.GetValue();         // in case this==&binary : binary's values are overwritten
  Int_t diffPost = binary.GetPost()-fPostCom;
  Int_t check     = 0;
  if(diffPost<0) check = fValue<<(-diffPost);
  else check = fValue>>(diffPost);
  if (CheckUSize(check)==kFALSE){    //checking size of pre-comma part
    
    //setting fValue to maximum
      
           
    fValue  = fuRestriction;         // fuRestriction >= 0 
    fValue  = fValue & (power1 - 1);
    fSigned = kFALSE;
    return *this;
  }

  Int_t val = (binary.GetSign()==kFALSE) ? check : -check; 
  if (CheckLSize(val)==kFALSE){    //checking size of pre-comma part
    
    //setting fValue to minimum
      
           
    if (flRestriction < 0) {
      fValue  = -flRestriction;
      fSigned = kTRUE;
    }
    else {
      fValue  = flRestriction;
      fSigned = kFALSE;
    }
    fValue  = fValue & (power1 - 1);
    return *this;
  }
  
  if (this == & binary) return *this;
  
  fSigned = kFALSE;
  Int_t iValue = fValue;
  fValue = fValue>>(binary.GetPost());           //only keep the valid pre-comma bits
  
  //append existing post-comma bits to fValue; cut off or add 0 if post-comma numbers don`t match
  for(Int_t i = 1; i <= fPostCom; i++){
    if(i <= (binary.GetPost())){
      c = ((iValue)>>(binary.GetPost()-i)) & 1;
    }
    else{
      c = 0;
    }
    fValue  = fValue<<1;
    fValue  = fValue | c;
  }
  
  fValue = fValue & (power1 - 1);
  fSigned = binary.GetSign();
  return *this;
}

//_____________________________________________________________________________  
AliTRDtrapAlu AliTRDtrapAlu::operator+(const AliTRDtrapAlu& binary){ 
  // + operator

  //no const parameter because referenced object will be changed
     
    AliTRDtrapAlu alu;
  
  Int_t binPre     = binary.GetPre();
  Int_t binPost    = binary.GetPost();
  Int_t binVal     = binary.GetValue();
  
  Int_t min         = Min(binPost,fPostCom);
  Int_t max         = Max(binPre,fPreCom);
  
  Int_t shift       = binPost - min;
  Int_t add1        = (binVal)>>(shift);    //for addition: cut off at minimum accuracy
  shift             = fPostCom - min;
  Int_t add2        = fValue>>(shift);
  if(binary.GetSign() == kTRUE) add1 = -add1;
  if(fSigned == kTRUE) add2 = -add2;
  Int_t add = add1 + add2;
  
  /*
  //because the parameter "binary" could be a reference to the object to which Mem() is a reference, do not change Mem() until you have extracted all information from "binary"; otherwise you change the information you would like to read
  Mem().Init(max + 1,min);      //buffer: enough space for pre-comma,post-comma according to accuracy
  Mem().AssignFormatted(Max(add,-add));
  Mem().SetSign(add);
  

  //Mem().FastInit(max+1,min,add);
  return Mem();*/
 
 alu.Init(max + 1,min);      //buffer: enough space for pre-comma,post-comma according to accuracy
 alu.AssignFormatted(Max(add,-add));
 alu.SetSign(add);
 
 return alu;

}

//_____________________________________________________________________________  
AliTRDtrapAlu AliTRDtrapAlu::operator-(const AliTRDtrapAlu& binary){
  // - operator    

    AliTRDtrapAlu alu;

  Int_t binPre    = binary.GetPre();
  Int_t binPost   = binary.GetPost();
  Int_t binVal    = binary.GetValue();


  Int_t min      = Min(binPost,fPostCom);
  Int_t max      = Max(binPre,fPreCom);

  Int_t shift    = binPost - min;
  Int_t sub1 = (binVal)>>(shift); //for addition: cut off at minimum accuracy
  shift = fPostCom - min;
  Int_t sub2 = fValue>>(shift);
  if(binary.GetSign() == kTRUE) sub1 = -sub1;
  if(fSigned  == kTRUE) sub2 = -sub2;
  Int_t sub = sub2 - sub1;     // order of subtraction is important
 
/*
  Mem().Init(max + 1,min);      //buffer: enough space for pre-comma, post-comma according to accuracy
  Mem().AssignFormatted(Max(sub,-sub)); 
  Mem().SetSign(sub);
  //Mem().FastInit(max+1,min,sub);
  return Mem();*/
  
  alu.Init(max + 1,min);
  alu.AssignFormatted(Max(sub,-sub)); 
  alu.SetSign(sub);

  return alu;

} 

//_____________________________________________________________________________  
AliTRDtrapAlu AliTRDtrapAlu::operator*(const AliTRDtrapAlu& binary){
  // * operator
  
    AliTRDtrapAlu alu;

  Int_t binPre   = binary.GetPre();
  Int_t binPost  = binary.GetPost();


  Int_t min      = Min(binPost,fPostCom);
  Int_t max      = Max(binPre,fPreCom);

  
  Int_t mult1 = binary.GetValue();
  Int_t mult2 = fValue;
  Int_t shift  = (Int_t)(fPostCom + binPost - min);
  Double_t fmult1 = (Double_t)mult1;
  Double_t fmult2 = (Double_t)mult2;
  (fmult1 > fmult2) ? fmult1 = fmult1/LUT(shift) : fmult2 = fmult2/LUT(shift);
  
    
  if (binary.GetSign() == kTRUE) fmult1 = -fmult1;
  if (fSigned  == kTRUE) fmult2 = -fmult2;
  Double_t fmult  = fmult1*fmult2;
  Int_t mult = (Int_t)fmult;
  Int_t sign = 1;
  if(mult<0) sign = -1;
  mult = Max(mult,-mult);
  //Int_t shift = fPostCom + binPost - min;
  //mult = mult>>(shift);

/*
  Mem().Init(2 * max + 1, min); // +1 to consider the borrow from the past-comma part; accuracy of past-comma part is determined by the minimum; therefore, for the result not more accuracy is guaranteed
  // be aware that this only works if 2*max+1+min <= 32!! adjusting the pre-comma place to the value would consume too much time

  Mem().AssignFormatted(mult);
  Mem().SetSign(sign);
  //mult = sign*mult;
  //Mem().FastInit(2*max+1,min,mult);
  return Mem();*/
  
  alu.Init(2 * max + 1, min); // +1 to consider the borrow from the past-comma part; accuracy of past-comma part is determined by the minimum; therefore, for the result not more accuracy is guaranteed
// be aware that this only works if 2*max+1+min <= 32!! adjusting the pre-comma place to the value would consume too much time

  alu.AssignFormatted(mult);
  alu.SetSign(sign);
  
  return alu;
}
  
//_____________________________________________________________________________  
AliTRDtrapAlu AliTRDtrapAlu::operator/(const AliTRDtrapAlu& binary){
  // / operator
 
    AliTRDtrapAlu alu;

  Int_t binPre  = binary.GetPre();
  Int_t binPost = binary.GetPost();
  Int_t min              = Min(binPost,fPostCom);
  Int_t max              = Max(binPre,fPreCom);
  
  Int_t div1             = binary.GetValue(); //value in integer format
  Int_t div2             = fValue;
  
  // this approach does not always work because it can exceed the range of integers
  //Int_t numerator     = div2 * LUT(min);
  Int_t numerator     = div2;
  if (fSigned == kTRUE) numerator = numerator*(-1);
  Int_t denominator   = div1;
  if (binary.GetSign() == kTRUE) denominator = denominator*(-1);
  Double_t fdiv       = 0.0;
  Double_t fLUT       = 0.0;
  Int_t div           = 0;
  

  if (div1 == 0){
      /*Mem().Init(max + 1,min);
    Mem().AssignFormatted(LUT(max+min+1)-1); // division by 0: set to max value
    //Mem().FastInit(max+1,min,div1);
    return Mem();*/
      alu.Init(max + 1,min);
      alu.AssignFormatted(LUT(max+min+1)-1); // division by 0: set to max value
      return alu;
  }
      
  fdiv = (Double_t)numerator/denominator;
  
  Int_t shift = fPostCom - binPost;
  
  if(shift>0){
      //denominator = denominator * LUT(shift);
      fLUT = (Double_t)LUT(min)/LUT(shift);
  }
  else {
      if(shift<0) {
	  shift = -shift;
	  //numerator =  numerator * LUT(shift);
	  fLUT = (Double_t)LUT(min)*LUT(shift);
      }
      else {
	  fLUT = (Double_t)LUT(min);
      }
  }

  fdiv = fdiv*fLUT;
  div = (Int_t)fdiv;
  
  Int_t sign = (div>=0) ? 1 : -1;
  div = Max(div,-div);
  
  // chose min as past-comma part because from a division of integers you can't get only an integer
  
  /*Mem().Init(max + 1,min); // max+1+min must <= 32!!
  Mem().SetSign(sign);
  Mem().AssignFormatted(div);
  
  return Mem();*/

  alu.Init(max + 1,min); // max+1+min must <= 32!!
  alu.SetSign(sign);
  alu.AssignFormatted(div);
  
  return alu;

  
}

//_____________________________________________________________________________  
Int_t AliTRDtrapAlu::MakePower(const Int_t& base,const  Int_t& exponent)const{
// calculate "base" to the power of "exponent"
  Int_t result = 1;
  
    for(Int_t i = 1; i <= exponent; i++){
    result = result * base;
  }
  return result;
}

//_____________________________________________________________________________  
Int_t AliTRDtrapAlu::LUT(const Int_t& index){   
  // simple look-up table for base=2
  
    AliTRDtrapAlu alu;

  static Bool_t fLUT = kFALSE;
  static Int_t gLUT[30];
  if (fLUT == kFALSE) {
    gLUT[0] = 1;
    for(Int_t i = 1; i<30; i++) {
      gLUT[i] = gLUT[i-1] * 2;
    }
  fLUT = kTRUE;
  } 
  if (index >=0 && index < 30){
    return gLUT[index];
  }
  else {
    
      //return Mem().MakePower(2,index);
      return alu.MakePower(2,index);
  }
}
