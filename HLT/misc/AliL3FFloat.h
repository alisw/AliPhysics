#ifndef ALIL3FFLOAT_H
#define ALIL3FFLOAT_H

#include "AliL3RootTypes.h"

#define DEFDIG 100
#define DEFMIN -1000000
#define DEFMAX  1000000

#define CALCSTATS
//#define FASTWITHROUNDINDERROS

class AliL3FFloat {
 public:
  AliL3FFloat(Double_t val=0)       {Set(val);}
  AliL3FFloat(const AliL3FFloat &f) {Set(f);}
  virtual ~AliL3FFloat();

  AliL3FFloat& operator = (const AliL3FFloat &f)  {Set(f); return *this;}
  AliL3FFloat& operator = (const Double_t f)      {Set(f); return *this;}

  operator const Double_t () const {return fVal;}
  operator const Float_t  () const {return (Float_t)fVal;}
  operator const Int_t  () const {return (Int_t)fVal;}

  friend ostream& operator<<(ostream &os, const AliL3FFloat &f);

  friend AliL3FFloat operator + (const AliL3FFloat &f1,const AliL3FFloat &f2);
  friend AliL3FFloat operator + (const AliL3FFloat &f1,const Double_t f2);
  friend AliL3FFloat operator + (const Double_t f1,    const AliL3FFloat &f2);
  friend AliL3FFloat operator + (const AliL3FFloat &f1);
  friend AliL3FFloat operator - (const AliL3FFloat &f1,const AliL3FFloat &f2);
  friend AliL3FFloat operator - (const AliL3FFloat &f1,const Double_t f2);
  friend AliL3FFloat operator - (const Double_t f1,    const AliL3FFloat &f2);
  friend AliL3FFloat operator - (const AliL3FFloat &f1);
  friend AliL3FFloat operator * (const AliL3FFloat &f1,const AliL3FFloat &f2);
  friend AliL3FFloat operator * (const AliL3FFloat &f1,const Double_t f2);
  friend AliL3FFloat operator * (const Double_t f1,    const AliL3FFloat &f2);
  friend AliL3FFloat operator / (const AliL3FFloat &f1,const AliL3FFloat &f2);
  friend AliL3FFloat operator / (const AliL3FFloat &f1,const Double_t f2);
  friend AliL3FFloat operator / (const Double_t f1,    const AliL3FFloat &f2);

  AliL3FFloat& operator += (const AliL3FFloat &f);
  AliL3FFloat& operator += (const Double_t f);
  AliL3FFloat& operator -= (const AliL3FFloat &f);
  AliL3FFloat& operator -= (const Double_t f);
  AliL3FFloat& operator *= (const AliL3FFloat &f);
  AliL3FFloat& operator *= (const Double_t f);
  AliL3FFloat& operator /= (const AliL3FFloat &f);
  AliL3FFloat& operator /= (const Double_t f);

  //==,!=,>=, ...
  //++,--

  static void PrintStat();
  static void SetParams(Int_t dig=DEFDIG,Int_t min=DEFMIN,Int_t max=DEFMAX);
  void Set(Double_t f=0);
  void Set(AliL3FFloat &f);
  inline Double_t GetVal()      const {return fVal;}
  inline Double_t GetExactVal() const {return fExactVal;}

 private:

  Double_t Round(Double_t f);
  Bool_t CheckUpperBound();
  Bool_t CheckLowerBound();
  Bool_t CheckBounds() {return (CheckUpperBound() && CheckLowerBound());}

  Double_t fVal;
  Double_t fExactVal;

  static Int_t fDigits;
  static Char_t fQuery[10];
  static Int_t fMax;
  static Int_t fMin;

#ifdef CALCSTATS
  static Int_t fN;
  static Int_t fNRounded;
  static Int_t fNOpAdds;
  static Int_t fNOpMults;
  static Int_t fNOpDivs;
  static Int_t fNOpSubs;
  static Int_t fNOverFlow;
  static Int_t fNUnderFlow;
  static Double_t fNDiff;
#endif

  ClassDef(AliL3FFloat,1)
};

#endif





