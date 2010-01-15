#ifndef ALILHCDIPVALT_H
#define ALILHCDIPVALT_H

#include <typeinfo>
#include <TString.h>
#include <TObjString.h>
#include <TObject.h>
#include "AliDCSArray.h"
#include "AliLog.h"
class AliDCSArray;


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
  AliLHCDipValT(Int_t size=0, Double_t t=0,Int_t nv=1);
  AliLHCDipValT(AliLHCDipValT<Element> &src);
  virtual ~AliLHCDipValT() {delete[] fArray;}
  AliLHCDipValT& operator=(const AliLHCDipValT<Element> &src);
  AliLHCDipValT& operator+=(const AliLHCDipValT<Element> &src);
  AliLHCDipValT& operator+=(const AliDCSArray &src);
  Element&       operator[](Int_t i);
  Element        operator[](Int_t i)                          const;
  //
  void                SetSize(Int_t size);
  void                SetValue(Int_t i,Element v);
  void                SetValues(const Element *v, Int_t n);
  void                SetTimeStamp(Double_t v)                      {fTimeStamp = v;}
  void                SetNSamplesUsed(Int_t n)                      {SetUniqueID((UInt_t)n);}
  // 
  Int_t               GetSize()                               const {return fSize;}
  Element             GetValue(Int_t i=0)                     const;
  Element*            GetValues()                             const {return (Element*)fArray;}
  Double_t            GetTimeStamp()                          const {return fTimeStamp;}
  Int_t               GetNSamplesUsed()                       const {return GetUniqueID();}
  Char_t*             GetTimeAsString()                       const {return TimeAsString(fTimeStamp);}
  //
  virtual void Average(const Option_t *opt="");
  virtual void Clear(const Option_t *opt="");
  virtual void Print(const Option_t *opt="")                  const;
  //
  static Char_t*      TimeAsString(double t);
  //
 protected:
  //
  Double_t            fTimeStamp;        // timestamp of the entry
  Int_t               fSize;
  Element*            fArray;            //[fSize] array of entries
  //
  ClassDef(AliLHCDipValT,1)
};


//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>::AliLHCDipValT(Int_t size,Double_t t,Int_t nv) 
: fTimeStamp(t),fSize(0),fArray(0)
{
  //def. constructor
  SetNSamplesUsed(nv);
  SetSize(size);
}

//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>::AliLHCDipValT(AliLHCDipValT<Element> &src)
: TObject(src),fTimeStamp(src.fTimeStamp),fSize(0),fArray(0)
{
  //copy constructor
  SetSize(src.fSize);
  memcpy(fArray,src.fArray,fSize*sizeof(Element));
}

//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>& AliLHCDipValT<Element>::operator=(const AliLHCDipValT<Element> &src)
{
  //assingment
  if (this != &src) {
    ((TObject*)this)->operator=(src);
    SetTimeStamp(src.GetTimeStamp());
    if (fSize!=src.fSize) SetSize(src.fSize);
    memcpy(fArray,src.fArray,fSize*sizeof(Element));    
  }
  return *this;
}

//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>& AliLHCDipValT<Element>::operator+=(const AliLHCDipValT<Element> &src)
{
  //addition
  TString tp = typeid(fArray).name();
  if (tp == typeid(Char_t*).name()) {
    if (fSize<1 && src.fSize>0) *this=src; // for strings reassign only if this is empty
  }
  else if (fSize!=src.fSize) {
    AliError(Form("Sizes of arrays to add must be the same (%d vs %d)",fSize,src.fSize));
    return *this;    
  }
  else {
    double wtime0 = fTimeStamp*GetNSamplesUsed();
    double wtime1 = src.fTimeStamp*src.GetNSamplesUsed();
    int wtot = GetNSamplesUsed() + src.GetNSamplesUsed();
    SetNSamplesUsed(wtot);
    if (wtot>0) SetTimeStamp( (wtime0 + wtime1)/wtot );
    for (int i=fSize;i--;) fArray[i] += src.fArray[i];
  }
  return *this;
}

//__________________________________________________________________________
template<class Element>
AliLHCDipValT<Element>& AliLHCDipValT<Element>::operator+=(const AliDCSArray &src)
{
  //addition
  TString tp = typeid(fArray).name();
  int isstr = (tp==typeid(Char_t*).name()) + (src.GetType()==AliDCSArray::kString);
  if ( isstr == 1) {
    AliError("Only one of the sides is of the string type");
    return *this;
  }
  else if (isstr == 2) { // both are string type
    if (fSize<1 && src.GetNEntries()>0) {
      TString str;
      for (int i=0;i<src.GetNEntries();i++) {
	str += ((TObjString*)src.GetStringArray(i))->GetName();
	str += " ";
      }
      SetSize(str.Length()+1);
      memcpy(fArray,str.Data(),fSize*sizeof(char));
    }
  } else {
    if (fSize>src.GetNEntries()) {
      AliError(Form("Size of RHS (%d) is smaller than size of LHS (%d)",src.GetNEntries(),fSize));
      return *this;
    }
    for (int i=fSize;i--;) {
      Element &val = fArray[i];
      if      (src.GetType() == AliDCSArray::kDouble) val += src.GetDouble()[i];
      else if (src.GetType() == AliDCSArray::kFloat)  val += src.GetFloat()[i];
      else if (src.GetType() == AliDCSArray::kInt)    val += src.GetInt()[i];
      else if (src.GetType() == AliDCSArray::kUInt)   val += src.GetUInt()[i];
      else if (src.GetType() == AliDCSArray::kChar)   val += src.GetChar()[i];
      else if (src.GetType() == AliDCSArray::kBool)   val += src.GetBool()[i];
    }
  }
  //
  double wtime0 = fTimeStamp*GetNSamplesUsed();
  double wtime1 = src.GetTimeStamp();
  int wtot = GetNSamplesUsed() + 1;
  SetNSamplesUsed(wtot);
  if (wtot>0) SetTimeStamp( (wtime0 + wtime1)/wtot ); 
  //
  return *this;
}

//__________________________________________________________________________
template<class Element>
Char_t* AliLHCDipValT<Element>::TimeAsString(double t)
{
  time_t tt = (time_t) t;
  char* st = ctime(&tt);
  *(strchr(st,'\n')) = 0;
  return st;
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::SetValue(Int_t i,Element v)
{
  //assign value
  if (i>=fSize || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,fSize-1));
    return;
  }
  fArray[i] = v;
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::SetValues(const Element* v, Int_t n)
{
  //assign value
  if (n!=fSize) SetSize(n);
  memcpy(fArray,v,n*sizeof(Element));
}

//__________________________________________________________________________
template<class Element>
Element& AliLHCDipValT<Element>::operator[](Int_t i)
{
  //obtain value refeterence
  if (i>=fSize || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,fSize-1));
    return fArray[0];
  }
  return fArray[i];
}

//__________________________________________________________________________
template<class Element>
Element AliLHCDipValT<Element>::operator[](Int_t i) const
{
  //obtain value
  if (i>=fSize || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,fSize-1));
    return 0;
  }
  return fArray[i];
}

//__________________________________________________________________________
template<class Element>
Element AliLHCDipValT<Element>::GetValue(Int_t i) const
{
  //obtain value
  if (i>=fSize || i<0) {
    AliError(Form("Index %d is out of range 0:%d",i,fSize-1));
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
    int nc = fSize > sz ? sz:fSize; // n elems to copy
    if (nc) memcpy(arr, fArray, nc*sizeof(Element));
    if (nc<sz) memset(arr+nc, 0, (sz-nc)*sizeof(Element));
    if (fSize) delete[] fArray;
    fArray = arr;
    fSize = sz;
  }
  else {
    delete[] fArray;
    fSize = 0;
  }
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::Print(const Option_t *opt) const
{
  // print time and value
  TString str = opt; 
  str.ToLower();
  printf("Timestamp: ");
  if (str.Contains("raw")) printf("%f",GetTimeStamp());
  else printf("%s",GetTimeAsString());
  printf(" | samples %6d | Value: ",GetNSamplesUsed());
  //
  TString tp = typeid(fArray).name();
  if ( tp==typeid(Char_t*).name() ) printf(" => %s\n",(Char_t*)fArray);
  else {
    if (fSize>1) printf("\n");
    for (int i=0;i<fSize;i++) {
      if      (tp == typeid(Int_t*).name()    || tp == typeid(UInt_t*).name() ) printf("#%4d %+11d |",i,(Int_t)fArray[i]);
      else if (tp == typeid(Double_t*).name() || tp == typeid(Float_t*).name()) printf("#%4d %+e |",  i,(Double_t)fArray[i]);
      else printf(" ");
      if ( (i+1)%5 == 0) printf("\n");
    }
    if (fSize%5) printf("\n");
  }
  //
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::Clear(const Option_t *)
{
  // reset to 0 everything
  SetNSamplesUsed(0);
  SetTimeStamp(0);
  memset(fArray,0,fSize*sizeof(Element));
}

//__________________________________________________________________________
template<class Element>
void AliLHCDipValT<Element>::Average(const Option_t *)
{
  // average of samples
  int n = GetNSamplesUsed();
  if (n>0) for (int i=fSize;i--;) fArray[i] /= n;
  //
}


typedef AliLHCDipValT<Double_t> AliLHCDipValD;
typedef AliLHCDipValT<Float_t>  AliLHCDipValF;
typedef AliLHCDipValT<Int_t>    AliLHCDipValI;
typedef AliLHCDipValT<Char_t>   AliLHCDipValC;



#endif
