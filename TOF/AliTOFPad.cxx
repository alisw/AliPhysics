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

////////////////////////////////////////////////////////////////////////
//
// AliTOFPad class (class used in TOF Reconstruction)
//
//
// Authors: Bologna-ITEP-Salerno Group
// 
// Description: Physical description of the elementary TOF sensitive 
// volume (PAD) containing timing and charge induced data.
//
// Member variable summary description:
// - location of the pad according to the current pad numbering scheme
// - simulated time of flight and GEANT time of flight if the pad
//   has fired
// - matching flags with the tracks that have fired the pad
//
////////////////////////////////////////////////////////////////////////

#include "AliTOFPad.h"

ClassImp(AliTOFPad)
  
AliTOFPad::AliTOFPad() 
{
  //
  // Default ctor
  //
  fSector       = 0;
  fPlate        = 0;
  fStrip        = 0;
  fPixel        = 0;
  fTrack        =-1;
  fTrackMatched =-1;
  fState        = 0;
  fRealTime     = 0;
  fGeantTime    = 0;
  fCharge       = 1;
  fAverageTime  = 0;
  fHit          =-1;
}

//___________________________________________
AliTOFPad::AliTOFPad(Int_t sector, Int_t plate, Int_t strip ,Int_t pixel) 
{
  //
  // Par ctor
  //
  fSector       = sector;
  fPlate        = plate;
  fStrip        = strip;
  fPixel        = pixel;
  fTrack        =-1;
  fTrackMatched =-1;
  fState        = 0;
  fRealTime     = 0;
  fGeantTime    = 0;
  fCharge       = 1;
  fAverageTime  = 0;
  fHit          =-1;
}
//___________________________________________
void AliTOFPad::SetGeom(Int_t sector, Int_t plate, Int_t strip, Int_t pixel) 
{
  //
  // Set the pad location in TOF detector as
  // sector plate strip pixel
  //
  fSector = sector;
  fPlate  = plate;
  fStrip  = strip;
  fPixel  = pixel;   
}   

//___________________________________________
void AliTOFPad::SetTofChargeHit(Float_t realTime, Float_t charge, Float_t geantTime, Int_t hitnum)
{
  //
  // Set the realTime given by the PAD (if activated), the charge, 
  // the true time given by GEANT and the hit number
  
  fRealTime     = realTime;
  fGeantTime    = geantTime;
  fCharge       = charge;
  fHit          = hitnum;
}
