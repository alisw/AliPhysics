 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
// AliCRTConstants class
//
// This class serves to group constants needed by CRT detector in 1
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
// Author: Arturo Fernandez, Enrique Gamez
//         FCFM-UAP, Mexico.
//
////////////////////////////////////////////////////////////////////////

#include "AliCRTConstants.h"

AliCRTConstants* AliCRTConstants::fgInstance = 0;

const Float_t AliCRTConstants::fgkCageLenght          = 477.6;
const Float_t AliCRTConstants::fgkCageWidth           = 166.7;
const Float_t AliCRTConstants::fgkCageHeight          =  10.7;
const Float_t AliCRTConstants::fgkSinglePaletteLenght = 363.0;
const Float_t AliCRTConstants::fgkSinglePaletteWidth  =  19.7;
const Float_t AliCRTConstants::fgkSinglePaletteHeight =   1;
const Float_t AliCRTConstants::fgkActiveAreaGap       = 0.7;
const Float_t AliCRTConstants::fgkActiveAreaLenght    = AliCRTConstants::fgkSinglePaletteLenght;
const Float_t AliCRTConstants::fgkActiveAreaWidth     = 156.7;
const Float_t AliCRTConstants::fgkActiveAreaHeight    = 2*AliCRTConstants::fgkSinglePaletteHeight + AliCRTConstants::fgkActiveAreaGap;
const Float_t AliCRTConstants::fgkMagnetWidth         = 654.4;
const Float_t AliCRTConstants::fgkMagnetLenght        = 1200;
const Float_t AliCRTConstants::fgkMagMinRadius        = 790;
const Float_t AliCRTConstants::fgkMagMaxRadius        = AliCRTConstants::fgkMagMinRadius + 20;
const Float_t AliCRTConstants::fgkDepth               =4420; // cm

ClassImp(AliCRTConstants)

//_____________________________________________________________________________
AliCRTConstants::AliCRTConstants()
  : TObject()
{
  // Default constructor
}

//_____________________________________________________________________________
AliCRTConstants::AliCRTConstants(const AliCRTConstants& ct)
  : TObject(ct)
{
  // Copy constructor
}

//_____________________________________________________________________________
AliCRTConstants& AliCRTConstants::operator=(const AliCRTConstants&)
{
  // Asingment operator
  return *this;
}

//_____________________________________________________________________________
AliCRTConstants* AliCRTConstants::Instance()
{
  if ( !fgInstance ) {
    fgInstance = new AliCRTConstants;
  }
  return fgInstance;
}

//_____________________________________________________________________________
AliCRTConstants::~AliCRTConstants()
{
  fgInstance = 0;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::CageLenght() const
{
  // Module lenght
  return fgkCageLenght;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::CageWidth() const
{
  // Module width
  return fgkCageWidth;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::CageHeight() const
{
  // Module height
  return fgkCageHeight;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::SinglePaletteLenght() const
{
  // Lenght of the scintillator active zone for a single counter
  return fgkSinglePaletteLenght;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::SinglePaletteWidth() const
{
  // Width of the scintillator active zone for a single counter
  return fgkSinglePaletteWidth;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::SinglePaletteHeight() const
{
  // Height of the scintillator active zone for a single counter
  return fgkSinglePaletteHeight;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::ActiveAreaGap() const
{ 
  // Gap betwen scintillators
  return fgkActiveAreaGap;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::ActiveAreaLenght() const
{
  // Lenght of the scintillator active zone
  return fgkActiveAreaLenght;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::ActiveAreaWidth() const
{
  // Width of the scintillator active zone
  return fgkActiveAreaWidth;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::ActiveAreaHeight() const
{
  // Height of the scintillator active zone
  return fgkActiveAreaHeight;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::MagnetWidth() const
{
  // Magnet  width
  return fgkMagnetWidth;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::MagnetLenght() const
{
  // Magnet lenght
  return fgkMagnetLenght;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::MagMinRadius() const
{
  // Magnet Inner radius
  return fgkMagMinRadius;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::MagMaxRadius() const
{
  // Magnet outer radius
  return fgkMagMaxRadius;
}

//_____________________________________________________________________________
Float_t AliCRTConstants::Depth() const
{
  // Alice IP depth
  return fgkDepth;
}
