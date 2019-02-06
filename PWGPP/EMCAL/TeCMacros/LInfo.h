#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TArrayF.h>
#include <TObject.h>
#include <TString.h>
#include <TProfile2D.h>
class TCanvas;

#ifndef _LINFO_
#define _LINFO_
class LInfo : public TObject {
 public:
  LInfo(Int_t rn=0) : fRunNo(rn), fIsComputed(0) { CreateHistograms(); }
  virtual     ~LInfo() {;}
  void         Compute();
  TCanvas     *DrawHist(Int_t which, Int_t gain=1, const char *opt=0) const;
  TProfile    *GetStripHist(Int_t sm, Int_t gain=1)         const { return fhStrip[sm][gain]; }
  TProfile    *GetStripRmsHist(Int_t sm, Int_t gain=1)      const { return fhStripCount[sm][gain]; }
  TProfile    *GetStripWeightedHist(Int_t sm, Int_t gain=1) const { return fhStripWeighted[sm][gain]; }
  TProfile    *GetLedMonHist(Int_t sm, Int_t gain=1)        const { return GetStripHist(sm,gain); }
  TProfile    *GetLedMonRmsHist(Int_t sm, Int_t gain=1)     const { return GetStripRmsHist(sm,gain); }
  TProfile    *GetLedMonWeightedHist(Int_t sm, Int_t gain=1)const { return GetStripWeightedHist(sm,gain); }
  TProfile2D  *GetLedHist(Int_t sm, Int_t gain=1)           const { return fhLed[sm][gain]; }
  TProfile2D  *GetLedRmsHist(Int_t sm, Int_t gain=1)        const { return fhLedCount[sm][gain]; }
  TProfile2D  *GetLedWeightedHist(Int_t sm, Int_t gain=1)   const { return fhLedWeighted[sm][gain]; }
  TH2         *GetLedOverMonHist(Int_t sm, Int_t gain=1)    const { return fhAmpOverMon[sm][gain]; }
  TH1         *GetLedMonDispHist(Int_t sm, Int_t gain=1)    const { return fhStripRmsOverMean[sm][gain]; }
  TH2         *GetLedDispHist(Int_t sm, Int_t gain=1)       const { return fhLedRmsOverMean[sm][gain]; }
  const char  *GetName()                                    const { return Form("LEDInfo_%d",fRunNo); }
  Int_t        GetRunNo()                                   const { return fRunNo; }
  void         FillLed(Int_t mod,Int_t gain, Int_t col, Int_t row, Double_t amp, Double_t rms);
  void         FillStrip(Int_t mod,Int_t gain, Int_t strip, Double_t amp, Double_t rms);
  Double_t     FracLeds(Int_t sm, Int_t gain=1) const;
  Double_t     FracStrips(Int_t sm, Int_t gain=1) const;
  void         Print(Option_t *option="") const;

  static const Int_t kNSM = 20; 
  static Int_t NSM()     { return kNSM; }
  static Int_t NCol()    { return 48; } //eta direction
  static Int_t NRow()    { return 24; } //phi direction
  static Int_t NStrip()  { return NCol()/2; }

 protected:
  void         CreateHistograms();
  Int_t        fRunNo;                      // run number
  TProfile    *fhStrip[kNSM][2];            // LedMon average
  TProfile    *fhStripCount[kNSM][2];       // LedMon rms
  TProfile    *fhStripWeighted[kNSM][2];    // LedMon weighted average
  TProfile2D  *fhLed[kNSM][2];              // Led average
  TProfile2D  *fhLedCount[kNSM][2];         // Led rms
  TProfile2D  *fhLedWeighted[kNSM][2];      // Led weighted average
  TH2         *fhAmpOverMon[kNSM][2];       //! Led/LedMon ratio
  TH1         *fhStripRmsOverMean[kNSM][2]; //! RMS over Mean for LedMon
  TH2         *fhLedRmsOverMean[kNSM][2];   //! RMS over Mean for Led
  Bool_t       fIsComputed;                 //! =1 if computed
  ClassDef(LInfo, 3); // LED info class
};
#endif
#endif
