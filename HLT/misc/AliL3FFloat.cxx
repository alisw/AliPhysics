//$Id$

// Author: Constantin Loizides <mailto:loizides@fi.uib.no>
//*-- Copyright & copy CL

#ifdef USEFFLOAT

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3Logger.h"

//calculate statistics according to usage and
//difference to floating point results
#define CALCSTATS

//use cast to int instead of sprintf and atoi
//#define FASTWITHROUNDINGERROS

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

Int_t AliL3FFloat::fDigits = DEFDIG;
Int_t AliL3FFloat::fMax    = DEFMAX;
Int_t AliL3FFloat::fMin    = DEFMIN;

#ifdef CALCSTATS
Int_t AliL3FFloat::fN          = 0;
Int_t AliL3FFloat::fNRounded   = 0;
Int_t AliL3FFloat::fNOpAdds    = 0;
Int_t AliL3FFloat::fNOpMults   = 0;
Int_t AliL3FFloat::fNOpDivs    = 0;
Int_t AliL3FFloat::fNOpSubs    = 0;
Int_t AliL3FFloat::fNOverFlow  = 0;
Int_t AliL3FFloat::fNUnderFlow = 0;
Double_t AliL3FFloat::fNDiff   = 0;
#endif

void AliL3FFloat::PrintStat(){
#ifdef CALCSTATS
  cout << "fN:            " << fN << endl;
  cout << "fNRounded:     " << fNRounded << endl;
  cout << "fNOpAdds:      " << fNOpAdds << endl;
  cout << "fNOpSubs:      " << fNOpSubs << endl;
  cout << "fNOpMults:     " << fNOpMults << endl;
  cout << "fNOpDivs:      " << fNOpDivs << endl;
  cout << "fNOpOverFlow:  " << fNOverFlow << endl;
  cout << "fNOpUnderFlow: " << fNUnderFlow << endl;
  if(fN) cout << "fNDiff:        " << fNDiff/fN << endl;
#else
  cerr << "Not compiled with #define CALCSTATS!" << endl;
#endif
}

AliL3FFloat::~AliL3FFloat()
{
#ifdef CALCSTATS
  Double_t diff=fabs(fVal-fExactVal);
  //  if(diff>10./fDigits) cout << diff << " Diff " << *this << endl;
  fNDiff+=diff;
#endif
}

ostream& operator<<(ostream &os, const AliL3FFloat &f) 
{
  os << (Double_t)f << "(" << f.fExactVal << ")"; 
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
  r-=f2; 
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

AliL3FFloat operator * (const AliL3FFloat &f1,const Double_t f2)
{
  AliL3FFloat r(f1); 
  r*=f2; 
  return r;
}

AliL3FFloat operator * (const Double_t f1,const AliL3FFloat &f2)
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

AliL3FFloat operator / (const AliL3FFloat &f1,const Double_t f2)
{
  AliL3FFloat r(f1); 
  r/=f2; 
  return r;
}

AliL3FFloat operator / (const Double_t f1,const AliL3FFloat &f2)
{
  AliL3FFloat r(f1); 
  r/=f2; 
  return r;
}

#ifdef USEINTS
void AliL3FFloat::SetParams(Int_t dig,Int_t min,Int_t max)
{
  fDigits=dig;
  fMin=min;
  fMax=max;
}

inline void AliL3FFloat::Set(Double_t val)
{
  Round(val);
  CheckBounds();
#ifdef CALCSTATS
  fN++;
#endif
}

inline void AliL3FFloat::Set(const AliL3FFloat &f)
{
  fVali=f.GetValInt();
  fVal=f.GetVal();
  fExactVal=f.GetExactVal();
  CheckBounds();
#ifdef CALCSTATS
  fN++;
#endif
}

inline void AliL3FFloat::Round(Double_t val)
{
  fExactVal=val;
  fVali=Fnt_t(val*fDigits);
  fVal=Double_t(fVali)/fDigits;
#ifdef CALCSTATS
  if(fVal!=fExactVal) fNRounded++;
#endif
}

inline Bool_t AliL3FFloat::CheckUpperBound()
{
  if(fVal>fMax){
    fVal=fMax;
    fVali=Fnt_t(fMax*fDigits);
#ifdef CALCSTATS
    fNOverFlow++;
#endif
    return kFALSE;
  }
  return kTRUE;
}

inline Bool_t AliL3FFloat::CheckLowerBound()
{
  if(fVal<fMin){
    fVal=fMin;
    fVali=Fnt_t(fMin*fDigits);
#ifdef CALCSTATS
    fNUnderFlow++;
#endif
    return kFALSE;
  }
  return kTRUE;
}

AliL3FFloat& AliL3FFloat::operator += (const AliL3FFloat &f)
{
  fExactVal+=f.GetExactVal();
  fVali+=f.GetValInt(); 
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpAdds++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator += (const Double_t f)     
{
  fExactVal+=f;
  fVali+=Fnt_t(f*fDigits);
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpAdds++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator -= (const AliL3FFloat &f) 
{
  fExactVal-=f.GetExactVal();
  fVali-=f.GetValInt(); 
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpSubs++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator -= (const Double_t f)
{
  fExactVal-=f;
  fVali-=Fnt_t(f*fDigits);
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpSubs++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator *= (const AliL3FFloat &f) 
{
  fExactVal*=f.GetExactVal();
  fVali=Fnt_t((fVali*f.GetValInt())/fDigits);
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpMults++;
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator *= (const Double_t f)     
{
  fExactVal*=f;
  fVali=Fnt_t(fVali*Fnt_t(f));
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpMults++;
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator /= (const AliL3FFloat &f) 
{
  fExactVal/=f.GetExactVal();
  fVali=Fnt_t(fVali*fDigits/f.GetValInt());
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpDivs++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator /= (const Double_t f)     
{
  fExactVal/=f;
  fVali=Fnt_t((fVali*fDigits)/(Int_t(f*fDigits)));
  fVal=Double_t(fVali)/fDigits;  
  CheckBounds();
#ifdef CALCSTATS
  fNOpDivs++; 
#endif
  return *this;
}

//--------------------------------------------------------
#else
//--------------------------------------------------------

Char_t AliL3FFloat::fQuery[10] = "%.2f";

inline void AliL3FFloat::Set(const Double_t val)
{
  fVal=Round(val);
  fExactVal=val;
  CheckBounds();
#ifdef CALCSTATS
  fN++;
#endif
}

inline void AliL3FFloat::Set(const AliL3FFloat &f)
{
  fVal=(Double_t)f;
  fExactVal=f.GetExactVal();
  CheckBounds();
#ifdef CALCSTATS
  fN++;
#endif
}

AliL3FFloat& AliL3FFloat::operator += (const AliL3FFloat &f)
{
  Double_t ev=fExactVal+f.GetExactVal();
  Set(fVal+(Double_t)f); 
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpAdds++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator += (const Double_t f)     
{
  Double_t ev=fExactVal+f;
  Set(fVal+Round(f));   
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpAdds++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator -= (const AliL3FFloat &f) 
{
  Double_t ev=fExactVal-f.GetExactVal();
  Set(fVal-(Double_t)f);
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpSubs++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator -= (const Double_t f)
{
  Double_t ev=fExactVal-f;
  Set(fVal-Round(f)); 
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpSubs++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator *= (const AliL3FFloat &f) 
{
  Double_t ev=fExactVal*f.GetExactVal();
  Set(fVal*(Double_t)f);
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpMults++;
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator *= (const Double_t f)     
{
  Double_t ev=fExactVal*f;
  Set(fVal*Round(f));   
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpMults++;
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator /= (const AliL3FFloat &f) 
{
  Double_t ev=fExactVal/f.GetExactVal();
  Set(fVal/(Double_t)f);
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpDivs++; 
#endif
  return *this;
}

AliL3FFloat& AliL3FFloat::operator /= (const Double_t f)     
{
  Double_t ev=fExactVal/f;
  Set(fVal/Round(f));   
  fExactVal=ev;
#ifdef CALCSTATS
  fNOpDivs++; 
#endif
  return *this;
}

inline Bool_t AliL3FFloat::CheckUpperBound()
{
  if(fVal>fMax){
    fVal=fMax;
#ifdef CALCSTATS
    fNOverFlow++;
#endif
    return kFALSE;
  }
  return kTRUE;
}

inline Bool_t AliL3FFloat::CheckLowerBound()
{
  if(fVal<fMin){
    fVal=fMin;
#ifdef CALCSTATS
    fNUnderFlow++;
#endif
    return kFALSE;
  }
  return kTRUE;
}

#ifdef FASTWITHROUNDINGERROS
inline Double_t AliL3FFloat::Round(Double_t val)
{
  Int_t dummy=Int_t(fDigits*val);
  Double_t ret=(Double_t)(dummy)/fDigits;
#ifdef CALCSTATS
  if(ret!=val) fNRounded++;
#endif
  return ret;
}
#else
inline Double_t AliL3FFloat::Round(Double_t val)
{
  static Char_t strnum[100];
  sprintf(strnum,fQuery,val);
  Double_t ret=atof(strnum);
#ifdef CALCSTATS
  if(ret!=val) fNRounded++;
#endif
  return ret;
}
#endif

void AliL3FFloat::SetParams(Int_t dig,Int_t min,Int_t max)
{
  fDigits=dig;
  Int_t prec=0;
  if(fDigits>0) prec=(Int_t)log10(fDigits);
  sprintf(fQuery,"%%.%df",prec);
  fMin=min;
  fMax=max;
}

#endif
#endif
