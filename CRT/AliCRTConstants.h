#ifndef ALICRTCONSTANTS_H
#define ALICRTCONSTANTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliCRTConstants class
//
// This class serves to group constants needed by ACORDE detector in 1
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
//
// Author: Arturo Fernandez, Enrique Gamez
//         FCFM-UAP, Mexico.
//
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>

enum ECRMode {
  kSingleMuons,
  kMuonBundle,
  kMuonFlux
};

class AliCRTConstants {
public:
  virtual ~AliCRTConstants();

  static AliCRTConstants* Instance();

  const Float_t CageLenght() const;
  const Float_t CageWidth() const;
  const Float_t CageHeight() const;

  const Float_t SinglePaletteLenght() const;
  const Float_t SinglePaletteWidth() const;
  const Float_t SinglePaletteHeight() const;

  const Float_t ActiveAreaGap() const;
  const Float_t ActiveAreaLenght() const;
  const Float_t ActiveAreaWidth() const;
  const Float_t ActiveAreaHeight() const;

  const Float_t MagnetWidth() const;
  const Float_t MagnetLenght() const;
  const Float_t MagMinRadius() const;
  const Float_t MagMaxRadius() const;

  const Float_t Depth() const;

protected:
  AliCRTConstants() {}
  AliCRTConstants(const AliCRTConstants& ct) {}
  AliCRTConstants& operator=(const AliCRTConstants& ct) {return *this;}

  static AliCRTConstants* fgInstance;

  static const Float_t fgkCageLenght;
  static const Float_t fgkCageWidth;
  static const Float_t fgkCageHeight;

  static const Float_t fgkSinglePaletteLenght;
  static const Float_t fgkSinglePaletteWidth;
  static const Float_t fgkSinglePaletteHeight;

  static const Float_t fgkActiveAreaGap;

  static const Float_t fgkActiveAreaLenght;
  static const Float_t fgkActiveAreaWidth;
  static const Float_t fgkActiveAreaHeight;

  static const Float_t fgkMagnetWidth;
  static const Float_t fgkMagnetLenght;
  static const Float_t fgkMagMinRadius;
  static const Float_t fgkMagMaxRadius;

  static const Float_t fgkDepth;

 private:
  ClassDef(AliCRTConstants, 0)   // CRT(ACORDE) global constants
};
#endif // ALICRTCONSTANTS_H
