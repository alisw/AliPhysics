#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TDatime.h>
#include <TFile.h>
#include <TFile.h>
#include <TGrid.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TArrayF.h>
#include <Riostream.h>

class TInfo : public TObject {
 public:
  TInfo(Int_t rn=0) : fRunNo(rn), fMinT(159), fMaxT(159), fAvTime(0), fFirstTime(0), fLastTime(0) {;}
  virtual     ~TInfo() {;}
  UInt_t       AvTime()                   const { return fAvTime;}
  Float_t      Diff(Int_t ns )            const { return fMaxT.At(ns)-fMinT.At(ns); }
  UInt_t       FirstTime()                const { return fAvTime;}
  const char  *GetName()                  const { return Form("TempInfo_%d",fRunNo); }
  UInt_t       LastTime()                 const { return fAvTime;}
  Bool_t       IsValid(Int_t ns)          const { return ((fMinT.At(ns)!=0)&&(fMaxT.At(ns)!=0)); }
  Int_t        RunNo()                    const { return fRunNo; }
  Float_t      MinT(Int_t ns)             const { return fMinT.At(ns); } 
  Float_t      MaxT(Int_t ns)             const { return fMaxT.At(ns); } 
  TArrayF     &MinT()                           { return fMinT; }
  TArrayF     &MaxT()                           { return fMaxT; }
  Int_t        Nvalid()                   const { Int_t ret=0; for (Int_t i=0;i<159;++i) ret += IsValid(i); return ret;}
  void         Print(Option_t *option="") const { cout << "Runno: " << fRunNo << " with average time " << fAvTime << " and " << Nvalid() << " entries" << endl;
                                                  for (Int_t i=0;i<159;++i) {
						    if (IsValid(i)) 
						      cout << "  " << i << " minT=" << fMinT.At(i) << ", maxT=" << fMaxT.At(i) << ", diff=" << Diff(i) << endl;
						  } 
                                                 }
  void         Set(Int_t ns, Float_t min, Float_t max) { fMinT.SetAt(min,ns); fMaxT.SetAt(max,ns); }
  void         SetTime(UInt_t av, UInt_t f, UInt_t l)  { fAvTime=av; fFirstTime=f; fLastTime=l;}

  static       Int_t SM(Int_t ns)               { return ns / 8;}
 protected:
  Int_t        fRunNo;      // run number
  TArrayF      fMinT;       // min temperature per sensor
  TArrayF      fMaxT;       // max temperature per sensor 
  UInt_t       fAvTime;     // average start time
  UInt_t       fFirstTime;  // first time
  UInt_t       fLastTime;   // last time
  ClassDef(TInfo, 2); // Temperature info class
};
#endif
