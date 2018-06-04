#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TArrayF.h>
#include <TObject.h>
#include <TString.h>

class TH1;
class TH2;

#ifndef _LINFO_
#define _LINFO_

class LInfo : public TObject {
 public:
 LInfo(Int_t rn=0) : fRunNo(rn) { CreateHistograms(); }
  virtual     ~LInfo() {;}
  //  TH2         *GetHist(Int_t type=1)      const;
  const char  *GetName()                  const { return Form("LEDInfo_%d",fRunNo); }
  void         Print(Option_t *option="") const;
  void         FillStrip(Int_t mod,Int_t gain, Int_t strip, Double_t amp, Double_t rms);
  void         FillLed(Int_t mod,Int_t gain, Int_t col, Int_t row, Double_t amp, Double_t rms);
  void         Compute();
  static const Int_t kNSM = 20; 
  static Int_t NSM()     { return kNSM; }
  static Int_t NCol()    { return 48; }
  static Int_t NRow()    { return 24; }
  static Int_t NStrip()  { return NCol()/2; }

 protected:
  void         CreateHistograms();
  Int_t        fRunNo;                  // run number
  TH1         *fhStrip[kNSM][2];        // Ledmon info
  TH1         *fhStripCount[kNSM][2];   // Ledmon counts
  TH2         *fhLed[kNSM][2];          // Led info
  TH2         *fhLedCount[kNSM][2];     // Led counts
  TH2         *fhAmpOverMon[kNSM][2];   //! Led/Ledmon ratio

  ClassDef(LInfo, 1); // LED info class
};
#endif
#endif
