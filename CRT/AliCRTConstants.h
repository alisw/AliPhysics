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
//
// Author: Arturo Fernandez, Enrique Gamez
//         FCFM-UAP, Mexico.
//
/////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliCRTConstants {
 public:
  // Modules
  static const Float_t fgCageLenght; // Module lenght
  static const Float_t fgCageWidth;  // Module width
  static const Float_t fgCageHeight; // Module height

// The dimensions of the active area of a single palette are:
  static const Float_t fgSinglePaletteLenght; // 
  static const Float_t fgSinglePaletteWidth; // 
  static const Float_t fgSinglePaletteHeight; // 
  
  static const Float_t fgActiveAreaGap; // 
  
  // Aproximate dimensions of the active area of the module.
  static const Float_t fgActiveAreaLenght;
  static const Float_t fgActiveAreaWidth; // 
  static const Float_t fgActiveAreaHeight;
  
  // Magnet
  static const Float_t fgMagnetWidth; // 
  static const Float_t fgMagnetLenght; // 
  
  AliCRTConstants() {}
  virtual ~AliCRTConstants() {}

 private:
  ClassDef(AliCRTConstants, 0)   // CRT(ACORDE) global constants
};
    
#endif // ALICRTCONSTANTS_H
