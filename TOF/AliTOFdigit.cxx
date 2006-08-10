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

//_________________________________________________________________________//
//                                                                         //
//  TOF digit: member variables                                            //
//  fSector  : TOF sector                                                  //
//  fPlate   : TOF plate                                                   //
//  fStrip   : strips number                                               //
//  fPadx    : pad number along x                                          //
//  fPadz    : pad number along z                                          //
//  fTdc     : TDC                                                         //
//  fAdc     : ADC                                                         //
//                                                                         //
//  Getters, setters and member functions  defined here                    //
//                                                                         //
// -- Authors: F. Pierella, A. Seganti, D. Vicinanza                       //
//_________________________________________________________________________//

 
#include "Riostream.h"

#include "AliTOFdigit.h"
#include "AliTOFGeometry.h"

ClassImp(AliTOFdigit)

//______________________________________________________________________________
AliTOFdigit::AliTOFdigit(Int_t *tracks, Int_t *vol,Float_t *digit)
:AliDigit(tracks)
{
//
// Constructor of digit object
//

  fSector = vol[0];
  fPlate  = vol[1];
  fStrip  = vol[2];
  fPadx   = vol[3];
  fPadz   = vol[4];
  fTdc    = digit[0];
  fTdcND  = digit[3];
  fAdc    = digit[1];
  fToT    = digit[2];

}

//____________________________________________________________________________
AliTOFdigit::AliTOFdigit(const AliTOFdigit & digit)
:AliDigit(digit)
{
  // 
  // copy ctor for AliTOFdigit object
  //

  Int_t i ;
  for ( i = 0; i < 3 ; i++)
    fTracks[i]  = digit.fTracks[i] ;
  fSector = digit.fSector;
  fPlate  = digit.fPlate;
  fStrip  = digit.fStrip;
  fPadx   = digit.fPadx;
  fPadz   = digit.fPadz;
  fTdc    = digit.fTdc;
  fTdcND    = digit.fTdcND;
  fAdc    = digit.fAdc;
  fToT = digit.fToT;

}

//______________________________________________________________________________
AliTOFdigit::AliTOFdigit(Int_t sector, Int_t plate, Int_t strip, Int_t padx,
Int_t padz, Float_t tdc, Float_t adc)
{
//
// Constructor for sdigit
//
  fSector = sector;
  fPlate  = plate;
  fStrip  = strip;
  fPadx   = padx;
  fPadz   = padz;  
  fTdc    = tdc;   
  fTdcND    = 0;   
  fAdc    = adc;     
  fToT = 0;
}
   
//______________________________________________________________________________
void AliTOFdigit::GetLocation(Int_t *Loc) const
{
//
// Get the cohordinates of the digit
// in terms of Sector - Plate - Strip - Pad
//

   Loc[0]=fSector;
   Loc[1]=fPlate;
   Loc[2]=fStrip;
   Loc[3]=fPadx;
   Loc[4]=fPadz;
}

//______________________________________________________________________________
Int_t AliTOFdigit::GetTotPad(AliTOFGeometry *tofGeom) const
{
//
// Get the "total" index of the pad inside a Sector
// starting from the digits data.
//

  Int_t before=0;

  switch(fPlate){ 
  case 0:
    //before = 0;
    break;
  case 1:
    before = tofGeom->NStripC();
    break;
  case 2:
    before = tofGeom->NStripC() +   AliTOFGeometry::NStripB();
    break;
  case 3:
    before = tofGeom->NStripC() +   AliTOFGeometry::NStripB() + AliTOFGeometry::NStripA();
    break;
  case 4:
    before = tofGeom->NStripC() + 2*AliTOFGeometry::NStripB() + AliTOFGeometry::NStripA();
    break;
  }
  
  Int_t pad = 2*fPadx + fPadz;
  //Int_t pad = fPadx+AliTOFGeometry::NpadX()*fPadz;
  Int_t strip  = fStrip + before;
  Int_t padTot = AliTOFGeometry::NpadXStrip()*strip + pad;

  return padTot;
}

//______________________________________________________________________________
void AliTOFdigit::AddTrack(Int_t track)
{
//
// Add a new and different track to the digit 
//
  if (track==fTracks[0] || track==fTracks[1] || track==fTracks[2]) return;
   if (fTracks[1]==0){
      fTracks[1] = track;
   }else if (fTracks[2]==0){
      fTracks[2] = track;
   }else{
   // printf("AliTOFdigit::AddTrack ERROR: Too many Tracks (>3) \n");
   }
}

// Overloading of Streaming, Sum and Comparison operators

//______________________________________________________________________________
Bool_t AliTOFdigit::operator==(AliTOFdigit const &digit) const
{
//
// Overloading of Comparison operator
//   
 if (fSector==digit.fSector &&
     fPlate==digit.fPlate &&
     fStrip==digit.fStrip &&
     fPadx==digit.fPadx &&
     fPadz==digit.fPadz &&
     fTdc==digit.fTdc &&
     fTdcND==digit.fTdcND &&
     fAdc==digit.fAdc &&
     fToT==digit.fToT ) return kTRUE;
     else return kFALSE;
}

//______________________________________________________________________________
AliTOFdigit& AliTOFdigit::operator+(AliTOFdigit const &digit)
{
//
// Overloading of Sum operator
// Note: Some convolution 
// between the two digit variables has to be inserted
//
if  (fSector==digit.fSector &&
     fPlate==digit.fPlate &&
     fStrip==digit.fStrip &&
     fPadx==digit.fPadx &&
     fPadz==digit.fPadz) {
                            // convolution to be inserted here
                             fTdc+=digit.fTdc;
                             fAdc+=digit.fAdc;
                           } else
                AliTOFdigit(fSector,fPlate,fStrip,fPadx,fPadz,fTdc,fAdc);
  return *this;
}

//______________________________________________________________________________
ostream& operator << (ostream& out, const AliTOFdigit &digit)
{
  //
  // Output streamer: output of the digit data
  //

  out << "Sector " << digit.fSector << ", Plate " << digit.fPlate << ", Strip " << digit.fStrip << endl;
  out << "Padx" << digit.fPadx << ", Padz " << digit.fPadz << endl;
  out << "TDC " << digit.fTdc << ", ADC "<< digit.fAdc << endl;

  return out;
}
