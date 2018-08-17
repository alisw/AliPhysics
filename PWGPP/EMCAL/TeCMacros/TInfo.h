#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TArrayF.h>
#include <TString.h>

class TH2;

#ifndef _TINFO_
#define _TINFO_
class TInfo : public TObject {
 public:
  TInfo(Int_t rn=0) : fRunNo(rn), fAvgT(160), fRmsT(160), fMinT(160), fMaxT(160), fAvTime(0), fFirstTime(0), fLastTime(0), fNFaulty(25) {;}
  virtual     ~TInfo() {;}
  Float_t      AbsMinT(Int_t t=1)         const;
  Float_t      AbsMaxT(Int_t t=2)         const;
  TArrayF     &AvgT()                           { return fAvgT; }
  Float_t      AvgT(Int_t ns)             const { return fMinT.At(ns); }
  Float_t      AvgTempSM(Int_t sm, Int_t type=3) const;
  UInt_t       AvTime()                   const { return fAvTime; }
  Float_t      Diff(Int_t ns )            const { return fMaxT.At(ns)-fMinT.At(ns); }
  UInt_t       FirstTime()                const { return fFirstTime; }
  UInt_t       GetAverageTime()           const { return fAvTime;}
  UInt_t       GetFirstTime()             const { return fFirstTime;}
  TH2         *GetHist(Int_t type=1)      const;
  UInt_t       GetLastTime()              const { return fLastTime;}
  const char  *GetName()                  const { return Form("TempInfo_%d",fRunNo); }
  Int_t        GetRunNo()                 const { return fRunNo; }
  Int_t        GetNFaulty()               const { return fNFaulty; }
  Double_t     Fraction()                 const { return 1.*Nvalid()/(NSensors()-GetNFaulty()); }
  UInt_t       LastTime()                 const { return fAvTime; }
  Bool_t       IsValid(Int_t ns)          const { return ((fMinT.At(ns)>0)&&(fMaxT.At(ns)>0)); }
  Int_t        RunNo()                    const { return fRunNo; }
  Float_t      T(Int_t ns, Int_t t)       const;
  Float_t      MinT(Int_t ns)             const { return fMinT.At(ns); }
  Float_t      MaxT(Int_t ns)             const { return fMaxT.At(ns); }
  TArrayF     &MinT()                           { return fMinT; }
  TArrayF     &MaxT()                           { return fMaxT; }
  Int_t        Nvalid()                   const { Int_t ret=0; for (Int_t i=0;i<160;++i) ret += IsValid(i); return ret;}
  void         Print(Option_t *option="") const;
  TArrayF     &RmsT()                           { return fRmsT; }
  Float_t      RmsT(Int_t ns)             const { return fMaxT.At(ns); }
  void         Set(Int_t ns, Float_t avg, Float_t rms, Float_t min, Float_t max) { fAvgT.SetAt(avg,ns); fRmsT.SetAt(rms,ns); fMinT.SetAt(min,ns); fMaxT.SetAt(max,ns); }
  void         SetNFaulty(Int_t n)              { fNFaulty=n; }
  void         SetTime(UInt_t av, UInt_t f, UInt_t l)  { fAvTime=av; fFirstTime=f; fLastTime=l; }

  static const char *Type(Int_t t);
  static             Int_t GetBin(Int_t ns);
  static             Int_t SM(Int_t ns)               { return ns/8; }
  static             Int_t NSensors()                 { return 160; }
  static             Int_t SensId(Int_t sm, Int_t row, Int_t col);

 protected:
  Int_t        fRunNo;      // run number
  TArrayF      fAvgT;       // avg temperature per sensor
  TArrayF      fRmsT;       // rms temperature per sensor
  TArrayF      fMinT;       // min temperature per sensor
  TArrayF      fMaxT;       // max temperature per sensor
  UInt_t       fAvTime;     // average start time
  UInt_t       fFirstTime;  // first time
  UInt_t       fLastTime;   // last time
  Int_t        fNFaulty;    //!faulty sensors
  ClassDef(TInfo, 4); // Temperature info class
};
#endif
#endif
