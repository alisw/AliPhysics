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
//
// Author: Arturo Fernandez, Enrique Gamez
//         FCFM-UAP, Mexico.
//
////////////////////////////////////////////////////////////////////////

#include "AliCRTConstants.h"

const Float_t AliCRTConstants::fgCageLenght          = 477.6;
const Float_t AliCRTConstants::fgCageWidth           = 166.7;
const Float_t AliCRTConstants::fgCageHeight          =  10.7;
const Float_t AliCRTConstants::fgSinglePaletteLenght = 363.0;
const Float_t AliCRTConstants::fgSinglePaletteWidth  =  19.7;
const Float_t AliCRTConstants::fgSinglePaletteHeight =   1;
const Float_t AliCRTConstants::fgActiveAreaGap       = 0.7;
const Float_t AliCRTConstants::fgActiveAreaLenght    = AliCRTConstants::fgSinglePaletteLenght;
const Float_t AliCRTConstants::fgActiveAreaWidth     = 156.7;
const Float_t AliCRTConstants::fgActiveAreaHeight    = 2*AliCRTConstants::fgSinglePaletteHeight + AliCRTConstants::fgActiveAreaGap;
const Float_t AliCRTConstants::fgMagnetWidth         = 654.4;
const Float_t AliCRTConstants::fgMagnetLenght        = 1200;
const Float_t AliCRTConstants::fgMagMinRadius        = 790;
const Float_t AliCRTConstants::fgMagMaxRadius        = AliCRTConstants::fgMagMinRadius + 20;
const Float_t AliCRTConstants::fgDepth               =4420; // cm

ClassImp(AliCRTConstants)
