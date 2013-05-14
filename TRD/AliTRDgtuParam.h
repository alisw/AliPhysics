#ifndef ALITRDGTUPARAM_H
#define ALITRDGTUPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDgtuParam.h 27496 2008-07-22 08:35:45Z cblume $ */

// --------------------------------------------------------
//
// Singleton class to hold the parameters steering the GTU
// tracking
//
// --------------------------------------------------------

#include "TObject.h"

class AliTRDgeometry;

class AliTRDgtuParam : public TObject {
 public:
  virtual ~AliTRDgtuParam();

  static AliTRDgtuParam *Instance(); // Singleton
  static void Terminate();

  static Int_t GetNLinks() { return fgkNLinks; }
  static Int_t GetNLayers() { return fgkNLinks/2; }
  static Int_t GetNZChannels() { return fgkNZChannels; }
  static Int_t GetNRefLayers() { return fgkNRefLayers; }

  static Float_t GetChamberThickness() { return 3.0; }

  // ----- Bin widths (granularity) -----
  static Float_t GetBinWidthY() { return fgkBinWidthY; }
  static Float_t GetBinWidthdY() { return fgkBinWidthdY; }

  // ----- Bit Widths (used for internal representation) -----
  static Int_t GetBitWidthY() { return fgkBitWidthY; }
  static Int_t GetBitWidthdY() { return fgkBitWidthdY; }
  static Int_t GetBitWidthYProj() { return fgkBitWidthYProj; }
  static Int_t GetBitExcessY() { return fgkBitExcessY; }
  static Int_t GetBitExcessAlpha() { return fgkBitExcessAlpha; }
  static Int_t GetBitExcessYProj() { return fgkBitExcessYProj; }

  Float_t GetInnerPadLength(Int_t stack, Int_t layer) const {
      return (stack == 2) ? 9. : fgkInnerPadLength[layer];
  }
  Float_t GetOuterPadLength(Int_t stack, Int_t layer) const {
      return (stack == 2) ? 8. : fgkOuterPadLength[layer];
  }
  Float_t GetZrow(Int_t stack, Int_t layer, Int_t padrow) const {
    Float_t zRowCorrected = fgkRow0Pos[layer][stack] - GetOuterPadLength(stack, layer) + GetInnerPadLength(stack, layer);
    return zRowCorrected - (0.5 + padrow) * GetInnerPadLength(stack, layer);
  }

  AliTRDgeometry* GetGeo() const { return fGeo; }
  Float_t GetVertexSize() const { return fVertexSize; }
  Int_t GetCiAlpha(Int_t layer) const;
  Int_t GetCiYProj(Int_t layer) const;
  Int_t GetYt(Int_t stack, Int_t layer, Int_t zrow) const;
  Int_t GetDeltaY() const { return fgDeltaY; }
  Int_t GetDeltaAlpha() const { return fgDeltaAlpha; }
  Int_t GetZSubchannel(Int_t stack, Int_t layer, Int_t zchannel, Int_t zpos) const;
  static Int_t GetRefLayer(Int_t refLayerIdx);
//  Bool_t GetFitParams(TVectorD &rhs, Int_t k); // const
  Bool_t GetIntersectionPoints(Int_t k, Float_t &x1, Float_t &x2); // const
  static Int_t GetPt(Int_t layerMask, Int_t a, Float_t b, Float_t x1, Float_t x2, Float_t magField);

  Bool_t IsInZChannel(Int_t stack, Int_t layer, Int_t zchannel, Int_t zpos) const;

  void SetVertexSize(Float_t vertexsize) { fVertexSize = vertexsize; }

  static void SetDeltaY(Int_t dy) { fgDeltaY = dy; }
  static void SetDeltaAlpha(Int_t da) { fgDeltaAlpha = da; }

  static void SetUseGTUconst(Bool_t b) { fgUseGTUconst = b; }
  static Bool_t GetUseGTUconst() { return fgUseGTUconst; }

  static void SetUseGTUmerge(Bool_t b) { fgUseGTUmerge = b; }
  static Bool_t GetUseGTUmerge() { return fgUseGTUmerge; }

  static void SetLimitNoTracklets(Bool_t b) { fgLimitNoTracklets = b; }
  static Bool_t GetLimitNoTracklets() { return fgLimitNoTracklets; }

  // z-channel map
  Int_t GenerateZChannelMap(); // could have different modes (for beam-beam, cosmics, ...)
  Bool_t DisplayZChannelMap(Int_t zchannel = -1, Int_t subch = 0) const;

  // variables for pt-reconstruction (not used at the moment)
  Bool_t GenerateRecoCoefficients(Int_t trackletMask);
  Int_t   GetAki(Int_t k, Int_t i);
  Float_t GetBki(Int_t k, Int_t i);
  Float_t GetCki(Int_t k, Int_t i);
//  Float_t GetD(Int_t k) const;

  // B-field
  void SetMagField(Float_t field) { fMagField = field; }
  Float_t GetMagField() const { return fMagField; }

  static const Int_t fgkNZChannels = 3; // No. of z-channels
  static const Int_t fgkNLinks = 12;	// No. of links
  static const Int_t fgkFixLayer = 2;	// which layer is fixed for the generation of the z-channel map
  static const Int_t fgkNRefLayers = 3;	 // no. of reference layers

  static const Float_t fgkBinWidthY; // bin width for y-position
  static const Float_t fgkBinWidthdY; // bin width for deflection length

  static const Int_t fgkBitWidthY; // bit width for y-position
  static const Int_t fgkBitWidthdY; // bit width for deflection length
  static const Int_t fgkBitWidthYProj; // bit width for projected y-position
  static const Int_t fgkBitExcessY; // excess bits for y-position
  static const Int_t fgkBitExcessAlpha; // excess bits for alpha
  static const Int_t fgkBitExcessYProj; // excess bits for projected y-position

  static const Int_t fgkPtInfinity; // infinite pt as obtained when a == 0

 protected:
  static       Int_t fgDeltaY;    	// accepted deviation in y_proj, default: 9
  static       Int_t fgDeltaAlpha;      // accepted deviation in alpha, default: 11

  static       Int_t fgRefLayers[3];    // reference layers for track finding

  static       Bool_t fgUseGTUconst;    // use constants as in the GTU for the calculations
					       // instead of geometry derived quantities
  static       Bool_t fgUseGTUmerge;    // use merge algorithm exactly as in hardware
  static       Bool_t fgLimitNoTracklets; // limit the number of tracklets per layer
  static const Bool_t fgZChannelMap[5][16][6][16]; // z-channel tables as in GTU
  static const Float_t fgkRadius[6];    // layer radius as used in the GTU code
  static const Float_t fgkThickness;    // drift length as used in the GTU code
  static const Float_t fgkRow0Pos[6][5]; // geometry constant from GTU implementation
  static const Float_t fgkInnerPadLength[6]; // geometry constant from GTU implementation
  static const Float_t fgkOuterPadLength[6]; // geometry constant from GTU implementation
  static const Float_t fgkAcoeff[32][6]; // geometry constant from GTU implementation
  static const Int_t   fgkMaskID[64]; // geometry constant from GTU implementation

  Float_t fVertexSize;		// assumed vertex size (z-dir.) for the z-channel map

  Int_t fZChannelMap[5][16][6][16];		  // must be changed
  Int_t fZSubChannel[5][fgkNZChannels][6][16];    // must be changed

  Int_t fCurrTrackletMask; // current tracklet mask for which the coefficients have been calculated
  Float_t fAki[6]; // coefficients used for the fit, calculated for the current tracklet mask
  Float_t fBki[6]; // coefficients used for the fit, calculated for the current tracklet mask
  Float_t fCki[6]; // coefficients used for the fit, calculated for the current tracklet mask

  Float_t fMagField;            // magnetic field in T

  AliTRDgeometry *fGeo;		//! pointer to the TRD geometry

  static AliTRDgtuParam *fgInstance; // instance pointer

 private:
  AliTRDgtuParam();			     // instance only via Instance()
  AliTRDgtuParam(const AliTRDgtuParam &rhs); // not implemented
  AliTRDgtuParam& operator=(const AliTRDgtuParam &rhs); // not implemented

  ClassDef(AliTRDgtuParam, 1);
};

#endif
