//$Id$

// Author: Constantin Loizides <mailto:loizides@fi.uib.no>
//*-- Copyright & copy CL

#include <stream.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "AliL3Logging.h"
#include "AliL3Logger.h"
#include "AliL3FFloat.h"

/** \class AliL3FFloat
//<pre>
//----------------------------------------------------
// AliL3FFloat
//
// Fixed Floating Point class for debugging purposes.
//
// The class behaves like a normal Double_t class but
// calculates everything in respect to fDigits (eg. 
// fDigits=100 -> 2 digits behind the comma). Further-
// more it keeps the exact value in floating precision
// and gathers some statistical information about 
// its usage.
//</pre>
 */


ClassImp(AliL3FFloat)

Int_t AliL3FFloat::fDigits =  DEFDIG;
Int_t AliL3FFloat::fMax = DEFMAX;
Int_t AliL3FFloat::fMin = DEFMIN;

Int_t AliL3FFloat::fN = 0;
Int_t AliL3FFloat::fNRounded = 0;
Int_t AliL3FFloat::fNOpAdds = 0;
Int_t AliL3FFloat::fNOpMults = 0;
Int_t AliL3FFloat::fNOpDivs = 0;
Int_t AliL3FFloat::fNOpSubs = 0;
Int_t AliL3FFloat::fNOverFlow = 0;
Int_t AliL3FFloat::fNUnderFlow = 0;

ostream& operator<<(ostream &os, const AliL3FFloat &f) 
{
  //  os << (Double_t)f << endl; 
  os << (Double_t)f << "(" << f.fExactVal << ")" << endl; 
  return os;
}

AliL3FFloat operator + (const AliL3FFloat &f1,const AliL3FFloat &f2)
{
  AliL3FFloat r(f1); 
  r+=f2; 
  return r;
}

AliL3FFloat operator + (const AliL3FFloat &f1,const Double_t f2)
{
  AliL3FFloat r(f1); 
  r+=f2; 
  return r;
}

AliL3FFloat operator + (const Double_t f1,const AliL3FFloat &f2)
{
  AliL3FFloat r(f1); 
  r+=f2; 
  return r;
}

AliL3FFloat operator + (const AliL3FFloat &f)
{
  AliL3FFloat r(f); 
  return r;
}


AliL3FFloat operator - (const AliL3FFloat &f1,const AliL3FFloat &f2)
{
  AliL3FFloat r(f1); 
  r+=f2; 
  return r;
}

AliL3FFloat operator - (const AliL3FFloat &f1,const Double_t f2)
{
  AliL3FFloat r(f1); 
  r-=f2; 
  return r;
}

AliL3FFloat operator - (const Double_t f1,const AliL3FFloat &f2)
{
  AliL3FFloat r(f1); 
  r-=f2; 
  return r;
}

AliL3FFloat operator - (const AliL3FFloat &f)
{
  AliL3FFloat r((-(Double_t)f)); 
  return r;
}

AliL3FFloat operator * (const AliL3FFloat &f1,const AliL3FFloat &f2)
{
  AliL3FFloat r(f1); 
  r*=f2; 
  return r;
}

AliL3FFloat operator / (const AliL3FFloat &f1,const AliL3FFloat &f2)
{
  AliL3FFloat r(f1); 
  r/=f2; 
  return r;
}

AliL3FFloat& AliL3FFloat::operator += (const AliL3FFloat &f)
{
  Double_t ev=fExactVal+f.GetExactVal();
  Set(fVal+(Double_t)f); 
  fExactVal=ev;
  fNOpAdds++; 
  return *this;
}

AliL3FFloat& AliL3FFloat::operator += (const Double_t f)     
{
  Double_t ev=fExactVal+f;
  Set(fVal+Round(f));   
  fExactVal=ev;
  fNOpAdds++; 
  return *this;
}

AliL3FFloat& AliL3FFloat::operator -= (const AliL3FFloat &f) 
{
  Double_t ev=fExactVal-f.GetExactVal();
  Set(fVal-(Double_t)f);
  fExactVal=ev;
  fNOpSubs++; 
  return *this;
}

AliL3FFloat& AliL3FFloat::operator -= (const Double_t f)
{
  Double_t ev=fExactVal-f;
  Set(fVal-Round(f)); 
  fExactVal=ev;
  fNOpSubs++; 
  return *this;
}

AliL3FFloat& AliL3FFloat::operator *= (const AliL3FFloat &f) 
{
  Double_t ev=fExactVal*f.GetExactVal();
  Set(fVal*(Double_t)f);
  fExactVal=ev;
  fNOpMults++;
  return *this;
}

AliL3FFloat& AliL3FFloat::operator *= (const Double_t f)     
{
  Double_t ev=fExactVal*f;
  Set(fVal*Round(f));   
  fExactVal=ev;
  fNOpMults++;
  return *this;
}

AliL3FFloat& AliL3FFloat::operator /= (const AliL3FFloat &f) 
{
  Double_t ev=fExactVal/f.GetExactVal();
  Set(fVal/(Double_t)f);
  fExactVal=ev;
  fNOpDivs++; 
  return *this;
}

AliL3FFloat& AliL3FFloat::operator /= (const Double_t f)     
{
  Double_t ev=fExactVal/f;
  Set(fVal/Round(f));   
  fExactVal=ev;
  fNOpDivs++; 
  return *this;
}

void AliL3FFloat::Set(Double_t val)
{
  fVal=Round(val);
  fExactVal=val;
  CheckBounds();
  fN++;
}

void AliL3FFloat::Set(AliL3FFloat &f)
{
  fVal=(Double_t)f;
  fExactVal=f.GetExactVal();
  CheckBounds();
  fN++;
}

Double_t AliL3FFloat::Round(Double_t val)
{
  Int_t dummy=Int_t(fDigits*val);
  Double_t ret=(Double_t)(dummy)/fDigits;
  if(ret!=val) fNRounded++;
  return ret;
}

Bool_t AliL3FFloat::CheckUpperBound()
{
  if(fVal>fMax){
    fVal=fMax;
    fNOverFlow++;
    return kFALSE;
  }

  return kTRUE;
}

Bool_t AliL3FFloat::CheckLowerBound()
{
  if(fVal<fMin){
    fVal=fMin;
    fNUnderFlow++;
    return kFALSE;
  }

  return kTRUE;
}

void AliL3FFloat::SetParams(Int_t dig,Int_t min,Int_t max)
{
  fDigits=dig;
  fMin=min;
  fMax=max;
}

void AliL3FFloat::PrintStat(){
  cout << "fN:            " << fN << endl;
  cout << "fNRounded:     " << fNRounded << endl;
  cout << "fNOpAdds:      " << fNOpAdds << endl;
  cout << "fNOpSubs:      " << fNOpSubs << endl;
  cout << "fNOpMults:     " << fNOpMults << endl;
  cout << "fNOpDivs:      " << fNOpDivs << endl;
  cout << "fNOpOverFlow:  " << fNOverFlow << endl;
  cout << "fNOpUnderFlow: " << fNUnderFlow << endl;
}
