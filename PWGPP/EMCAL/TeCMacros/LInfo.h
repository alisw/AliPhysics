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
  void         Compute();
  //  TH2         *GetHist(Int_t type=1)      const;
  const char  *GetName()                  const { return Form("LEDInfo_%d",fRunNo); }
  void         Print(Option_t *option="") const;
  void         FillStrip(Int_t mod,Int_t gain, Int_t strip, Double_t amp, Double_t rms);
  void         FillLed(Int_t mod,Int_t gain, Int_t col, Int_t row, Double_t amp, Double_t rms);
  Double_t     FracStrips(Int_t sm, Int_t gain=0) const;

  static const Int_t kNSM = 20; 
  static Int_t NSM()     { return kNSM; }
  static Int_t NCol()    { return 48; } //eta direction
  static Int_t NRow()    { return 24; } //phi direction
  static Int_t NStrip()  { return NCol()/2; }

 protected:
  void         CreateHistograms();
  Int_t        fRunNo;                      // run number
  TH1         *fhStrip[kNSM][2];            // LedMon info
  TH1         *fhStripCount[kNSM][2];       // LedMon counts
  TH2         *fhLed[kNSM][2];              // Led info
  TH2         *fhLedCount[kNSM][2];         // Led counts
  TH2         *fhAmpOverMon[kNSM][2];       //! Led/LedMon ratio
  TH1         *fhStripRmsOverMean[kNSM][2]; //! RMS over Mean for LedMon
  TH2         *fhLedRmsOverMean[kNSM][2];   //! RMS over Mean for Led

  ClassDef(LInfo, 1); // LED info class
};
#endif
#endif
