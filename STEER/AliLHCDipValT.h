#ifndef ALILHCDIPVALT_H
#define ALILHCDIPVALT_H

#include <typeinfo>
#include <TString.h>
#include <TObjString.h>
#include <TObject.h>
#include "AliLog.h"


/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
//  AliLHCDipValT: templated class to store the array from processed           //
//  LHC DIP data. Associated with the time stamp either of the value           //
//  acquisition or the last sample timestamp if many samples were processed    //
//                                                                             //
//  Author: ruben.shahoyan@cern.ch                                             //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

//_______________________________________________________________________________________
template<class Element> class AliLHCDipValT : public TObject
{
 public:
  //
  enum {
    kLastSpecial=BIT(14),         // last element of array is special (e.g. assigned error etc) 
    kProcessed1 =BIT(15),         // last element of array is special (e.g. assigned error etc) 
    kChar=BIT(22),
    kFloat=BIT(23)
  };
  //
  AliLHCDipValT(Int_t size=0, Double_t t=0);
  AliLHCDipValT(const AliLHCDipValT<Element> &src);
  virtual ~AliLHCDipValT() {delete[] fArray;}
  AliLHCDipValT& operator=(const AliLHCDipValT<Element> &src);
  Element&       operator[](Int_t i);
  Element        operator[](Int_t i)                          const;
  //
  void                SetSize(Int_t size);
  void                SetValue(Int_t i,Element v);
  void                SetValues(const Element *v, Int_t n);
  void                SetTimeStamp(Double_t v)                      {fTimeStamp = v;}
  // 
  Int_t               GetSizeTotal()                          const {return fSizeTot;}
  Element             GetValue(Int_t i=0)                     const;
  Element*            GetValues()                             const {return (Element*)fArray;}
  Double_t            GetTimeStamp()                          const {return fTimeStamp;}
  Char_t*             GetTimeAsString(Bool_t utc=kTRUE)       const {return TimeAsString(fTimeStamp,utc);}
  //
  void                SetProcessed1(Bool_t v=kTRUE)                 {SetBit(kProcessed1,v);}
  void                SetLastSpecial(Bool_t v=kTRUE)                {SetBit(kLastSpecial,v);}
  Bool_t              IsProcessed1()                          const {return TestBit(kProcessed1);}
  Bool_t              IsLastSpecial()                         const {return TestBit(kLastSpecial);}
  Int_t               GetSize()                               const {return IsLastSpecial() ? GetSizeTotal()-1:GetSizeTotal();}
  Bool_t              IsTypeC()                               const {return TestBit(kChar);}
  Bool_t              IsTypeF()                               const {return TestBit(kFloat);}
  Bool_t              IsTypeI()                               const {return !TestBit(kFloat|kChar);}
  //
  virtual void Clear(const Option_t *opt="");
  virtual void Print(const Option_t *opt="")                  const;
  //
  static Char_t*      TimeAsString(double t,Bool_t utc=kTRUE);
  //
 protected:
  //
  Double_t            fTimeStamp;        // timestamp of the entry
  Int_t               fSizeTot;          // vector total size (including special slots, like for errors)
  Element*            fArray;            //[fSizeTot] array of entries
  //
  ClassDef(AliLHCDipValT,1)
};


//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>::AliLHCDipValT(Int_t size,Double_t t) 
: fTimeStamp(t),fSizeTot(0),fArray(0)
{
  //def. constructor
  SetSize(size);
  if      (!strcmp(typeid(fArray).name(),typeid(Char_t*).name())) SetBit(kChar);
  else if (!strcmp(typeid(fArray).name(),typeid(Double_t*).name()) ||
	   !strcmp(typeid(fArray).name(),typeid(Float_t*).name() )) SetBit(kFloat);
  //
}

//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>::AliLHCDipValT(const AliLHCDipValT<Element> &src)
: TObject(src),fTimeStamp(src.fTimeStamp),fSizeTot(0),fArray(0)
{
  //copy constructor
  SetSize(src.GetSizeTotal());
  memcpy(fArray,src.fArray,GetSizeTotal()*sizeof(Element));
}

//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>& AliLHCDipValT<Element>::operator=(const AliLHCDipValT<Element> &src)
{
  //assingment
  if (this != &src) {
    ((TObject*)this)->operator=(src);
    if (GetSizeTotal()!=src.GetSizeTotal()) SetSize(src.GetSizeTotal());
    SetTimeStamp(src.GetTimeStamp());
    memcpy(fArray,src.fArray,GetSizeTotal()*sizeof(Element));    
  }
  return *this;
}

//__________________________________________________________________________
template<class Element>
Char_t* AliLHCDipValT<Element>::TimeAsString(double t, Bool_t utc)
{
  // time as string in UTC or local format
  static char buff[22];
  time_t tt = (time_t) t;
  struct tm *time = utc ? gmtime(&tt) : localtime(&tt);
  sprintf(buff,"%02d:%02d:%02d %02d/%02d/%04d",time->tm_hour,time->tm_min,time->tm_sec,
	  time->tm_mday,time->tm_mon+1,time->tm_year+1900);
  return (char*)buff;
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::SetValue(Int_t i,Element v)
{
  //assign value
  if (i>=GetSizeTotal() || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,GetSizeTotal()-1));
    return;
  }
  fArray[i] = v;
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::SetValues(const Element* v, Int_t n)
{
  //assign value
  if (n!=GetSizeTotal()) SetSize(n);
  memcpy(fArray,v,n*sizeof(Element));
}

//__________________________________________________________________________
template<class Element>
Element& AliLHCDipValT<Element>::operator[](Int_t i)
{
  //obtain value refeterence
  if (i>=GetSizeTotal() || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,GetSizeTotal()-1));
    return fArray[0];
  }
  return fArray[i];
}

//__________________________________________________________________________
template<class Element>
Element AliLHCDipValT<Element>::operator[](Int_t i) const
{
  //obtain value
  if (i>=GetSizeTotal() || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,GetSizeTotal()-1));
    return 0;
  }
  return fArray[i];
}

//__________________________________________________________________________
template<class Element>
Element AliLHCDipValT<Element>::GetValue(Int_t i) const
{
  //obtain value
  if (i>=GetSizeTotal() || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,GetSizeTotal()-1));
    return 0;
  }
  return fArray[i];
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::SetSize(Int_t sz)
{
  //resize
  Element* arr = 0;
  if (sz>0) {
    arr = new Element[sz];
    int nc = GetSizeTotal() > sz ? sz:GetSizeTotal(); // n elems to copy
    if (nc) memcpy(arr, fArray, nc*sizeof(Element));
    if (nc<sz) memset(arr+nc, 0, (sz-nc)*sizeof(Element));
    if (GetSizeTotal()) delete[] fArray;
    fArray = arr;
    fSizeTot = sz;
  }
  else {
    delete[] fArray;
    fArray = 0;
    fSizeTot = 0;
  }
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::Print(const Option_t *opt) const
{
  // print time and value
  TString str = opt; 
  str.ToLower();
  if (str.Contains("raw")) printf("%.1f ",GetTimeStamp());
  else printf("[%s] ",GetTimeAsString(!str.Contains("loc")));
  //
  TString tp = typeid(fArray).name();
  if ( tp==typeid(Char_t*).name() ) printf(": %s\n",(Char_t*)fArray);
  else {
    int sz = GetSize();
    if (sz>1) printf("\n");
    Bool_t eolOK = kFALSE;
    for (int i=0;i<sz;i++) {
      if      (tp == typeid(Int_t*).name()    || tp == typeid(UInt_t*).name() ) {
	if (!str.Contains("bit")) printf(" %6d |" ,(Int_t)fArray[i]);
	else {
	  printf(" ");
	  int val = (int)fArray[i];
	  for (int j=sizeof(int)*8;j--;) printf("%d",(val>>j)&0x1);
	  printf(" |");
	}

      }
      else if (tp == typeid(Double_t*).name() || tp == typeid(Float_t*).name()) printf(" %+.3e |",(Float_t)fArray[i]);
      else printf(" ");
      eolOK = kFALSE;
      if ( (i+1)%5 == 0) {printf("\n"); eolOK = kTRUE;}
    }
    if (IsLastSpecial()) {
      if (sz>1 && !eolOK) {printf("\n"); eolOK = kTRUE;}
      if (tp == typeid(Double_t*).name() || tp == typeid(Float_t*).name()) {
	printf(" Error: %+e\n",(Float_t)fArray[sz]);
	eolOK = kTRUE;
      }
    }
    if (!eolOK) printf("\n");
  }
  //
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::Clear(const Option_t *)
{
  // reset to 0 everything
  SetTimeStamp(0);
  memset(fArray,0,GetSizeTotal()*sizeof(Element));
}


typedef AliLHCDipValT<Double_t> AliLHCDipValD;
typedef AliLHCDipValT<Float_t>  AliLHCDipValF;
typedef AliLHCDipValT<Int_t>    AliLHCDipValI;
typedef AliLHCDipValT<Char_t>   AliLHCDipValC;



#endif
