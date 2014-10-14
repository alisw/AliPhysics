// This class is used in both AliRoot and libTRD. If you update it in
// one place you have to update also the other.

#ifndef ALITRDLTUPARAM_H
#define ALITRDLTUPARAM_H

#include "TObject.h"

class AliTRDltuParam : public TObject
{
 public:
  AliTRDltuParam();
  ~AliTRDltuParam();

  // configuration settings
  // called with special SCSN commands
  void SetPtMin(Int_t data)     { fPtMin = Float_t(data) / 1000.; }
  void SetMagField(Int_t data)  { fMagField = Float_t(data) / 1000.; }
  void SetOmegaTau(Int_t data)  { fOmegaTau = Float_t(data) / 1.e6; }
  void SetNtimebins(Int_t data) { fNtimebins = data; }
  void SetScaleQ0(Int_t data)   { fScaleQ0 = data; }
  void SetScaleQ1(Int_t data)   { fScaleQ1 = data; }
  void SetLengthCorrectionEnable(Int_t data) { fPidTracklengthCorr = Bool_t (data); }
  void SetTiltCorrectionEnable(Int_t data)   { fTiltCorr = Bool_t (data); }
  void SetPIDgainCorrectionEnable(Bool_t data)   { fPidGainCorr = data; }

  // set values directly
  void SetRawPtMin(Float_t data)     { fPtMin = data; }
  void SetRawMagField(Float_t data)  { fMagField = data; }
  void SetRawOmegaTau(Float_t data)  { fOmegaTau = data; }
  void SetRawNtimebins(Int_t data)   { fNtimebins = data; }
  void SetRawScaleQ0(Int_t data)     { fScaleQ0 = data; }
  void SetRawScaleQ1(Int_t data)     { fScaleQ1 = data; }
  void SetRawLengthCorrectionEnable(Bool_t data) { fPidTracklengthCorr = data; }
  void SetRawTiltCorrectionEnable(Bool_t data)   { fTiltCorr = data; }
  void SetRawPIDgainCorrectionEnable(Bool_t data)   { fPidGainCorr = data; }

  // retrieve the calculated information
  // which is written to the TRAPs
  Int_t GetDyCorrection(Int_t det, Int_t rob, Int_t mcm) const;
  void  GetDyRange(Int_t det, Int_t rob, Int_t mcm, Int_t ch, Int_t &dyMinInt, Int_t &dyMaxInt) const;
  void  GetCorrectionFactors(Int_t det, Int_t rob, Int_t mcm, Int_t ch,
			     UInt_t &cor0, UInt_t &cor1, Float_t gain = 1.) const;
  Int_t GetNtimebins() const;

  Float_t GetX(Int_t det, Int_t rob, Int_t mcm) const;
  Float_t GetLocalY(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const;
  Float_t GetLocalZ(Int_t det, Int_t rob, Int_t mcm) const;

  Float_t GetDist(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const;
  Float_t GetElongation(Int_t det, Int_t rob, Int_t mcm, Int_t ) const;
  Float_t GetPhi(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const;
  Float_t GetPerp(Int_t det, Int_t rob, Int_t mcm, Int_t ch) const;

 protected:
  // geometry constants
  static Float_t fgZrow[6][5];             // z-position of pad row edge
  static Float_t fgX[6];                   // x-position for all layers
  static Float_t fgTiltingAngle[6];	   // tilting angle for every layer
  static Float_t fgWidthPad[6];            // pad width for all layers
  static Float_t fgLengthInnerPadC0;       // inner pad length C0 chamber
  static Float_t fgLengthOuterPadC0;       // outer pad length C0 chamber
  static Float_t fgLengthInnerPadC1[6];    // inner pad length C1 chambers
  static Float_t fgLengthOuterPadC1[6];    // outer pad length C1 chambers
  static Float_t fgScalePad;		   // scaling factor for pad width
  static Float_t fgDriftLength;		   // length of the  parse gaintbl Krypton_2009-01 drift region
  static Float_t fgBinDy;		   // bin in dy (140 um)
  static Int_t   fgDyMax;		   // max dy for a tracklet (hard limit)
  static Int_t   fgDyMin;		   // min dy for a tracklet (hard limit)

  // settings
  Float_t fMagField;		// magnetic field
  Float_t fOmegaTau;		// omega tau, i.e. tan(Lorentz angle)
  Float_t fPtMin;		// min. pt for deflection cut
  Int_t   fNtimebins;	        // drift time in units of timebins << 5n
  UInt_t  fScaleQ0;	        // scale factor for accumulated charge Q0
  UInt_t  fScaleQ1;		// scale factor for accumulated charge Q1
  Bool_t  fPidTracklengthCorr;	// enable tracklet length correction
  Bool_t  fTiltCorr;		// enable tilt correction
  Bool_t  fPidGainCorr;         // enable MCM gain correction factor for PID

  ClassDef(AliTRDltuParam, 1);
};

#endif
