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
// AliACORDEConstants class
//
// This class serves to group constants needed by ACORDE detector in 1
// easily accessible place. All constants are public const static data 
// members. The class is never instatiated.
// Author: Arturo Fernandez, Enrique Gamez
//         FCFM-UAP, Mexico.
//
////////////////////////////////////////////////////////////////////////

#include "AliACORDEConstants.h"

AliACORDEConstants* AliACORDEConstants::fgInstance = 0;

const Float_t AliACORDEConstants::fgkCageLenght          = 477.6;
const Float_t AliACORDEConstants::fgkCageWidth           = 166.7;
const Float_t AliACORDEConstants::fgkCageHeight          =  10.7;
const Float_t AliACORDEConstants::fgkSinglePaletteLenght = 363.0;
const Float_t AliACORDEConstants::fgkSinglePaletteWidth  =  19.7;
const Float_t AliACORDEConstants::fgkSinglePaletteHeight =   1;
const Float_t AliACORDEConstants::fgkActiveAreaGap       = 0.7;
const Float_t AliACORDEConstants::fgkActiveAreaLenght    = AliACORDEConstants::fgkSinglePaletteLenght;
const Float_t AliACORDEConstants::fgkActiveAreaWidth     = 156.7;
const Float_t AliACORDEConstants::fgkActiveAreaHeight    = 2*AliACORDEConstants::fgkSinglePaletteHeight + AliACORDEConstants::fgkActiveAreaGap;
const Float_t AliACORDEConstants::fgkMagnetWidth         = 654.4;
const Float_t AliACORDEConstants::fgkMagnetLenght        = 1200;
const Float_t AliACORDEConstants::fgkMagMinRadius        = 790;
const Float_t AliACORDEConstants::fgkMagMaxRadius        = AliACORDEConstants::fgkMagMinRadius + 20;
const Float_t AliACORDEConstants::fgkDepth               =4420; // cm

ClassImp(AliACORDEConstants)

//_____________________________________________________________________________
AliACORDEConstants::AliACORDEConstants()
  : TObject()
{
  // Default constructor
}

//_____________________________________________________________________________
AliACORDEConstants::AliACORDEConstants(const AliACORDEConstants& ct)
  : TObject(ct)
{
  // Copy constructor
}

//_____________________________________________________________________________
AliACORDEConstants& AliACORDEConstants::operator=(const AliACORDEConstants&)
{
  // Asingment operator
  return *this;
}

//_____________________________________________________________________________
AliACORDEConstants* AliACORDEConstants::Instance()
{
  if ( !fgInstance ) {
    fgInstance = new AliACORDEConstants;
  }
  return fgInstance;
}

//_____________________________________________________________________________
AliACORDEConstants::~AliACORDEConstants()
{
  fgInstance = 0;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::CageLenght() const
{
  // Module lenght
  return fgkCageLenght;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::CageWidth() const
{
  // Module width
  return fgkCageWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::CageHeight() const
{
  // Module height
  return fgkCageHeight;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::SinglePaletteLenght() const
{
  // Lenght of the scintillator active zone for a single counter
  return fgkSinglePaletteLenght;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::SinglePaletteWidth() const
{
  // Width of the scintillator active zone for a single counter
  return fgkSinglePaletteWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::SinglePaletteHeight() const
{
  // Height of the scintillator active zone for a single counter
  return fgkSinglePaletteHeight;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ActiveAreaGap() const
{ 
  // Gap betwen scintillators
  return fgkActiveAreaGap;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ActiveAreaLenght() const
{
  // Lenght of the scintillator active zone
  return fgkActiveAreaLenght;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ActiveAreaWidth() const
{
  // Width of the scintillator active zone
  return fgkActiveAreaWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::ActiveAreaHeight() const
{
  // Height of the scintillator active zone
  return fgkActiveAreaHeight;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::MagnetWidth() const
{
  // Magnet  width
  return fgkMagnetWidth;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::MagnetLenght() const
{
  // Magnet lenght
  return fgkMagnetLenght;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::MagMinRadius() const
{
  // Magnet Inner radius
  return fgkMagMinRadius;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::MagMaxRadius() const
{
  // Magnet outer radius
  return fgkMagMaxRadius;
}

//_____________________________________________________________________________
Float_t AliACORDEConstants::Depth() const
{
  // Alice IP depth
  return fgkDepth;
}
