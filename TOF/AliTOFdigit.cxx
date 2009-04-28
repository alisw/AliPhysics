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
AliTOFdigit::AliTOFdigit()
  :AliDigit(),
   fSector(-1),
   fPlate(-1),
   fStrip(-1),
   fPadx(-1),
   fPadz(-1),
   fTdc(0),
   fTdcND(0),
   fAdc(0),
   fToT(0)
{
}
//______________________________________________________________________________
AliTOFdigit::AliTOFdigit(Int_t *tracks, Int_t *vol, Int_t *digit)
  :AliDigit(tracks),
   fSector(vol[0]),
   fPlate(vol[1]),
   fStrip(vol[2]),
   fPadx(vol[3]),
   fPadz(vol[4]),
   fTdc(digit[0]),
   fTdcND(digit[3]),
   fAdc(digit[1]),
   fToT(digit[2])
{
//
// Constructor of digit object
//
}

//____________________________________________________________________________
AliTOFdigit::AliTOFdigit(const AliTOFdigit & digit)
  :AliDigit(digit),
   fSector(digit.fSector),
   fPlate(digit.fPlate),
   fStrip(digit.fStrip),
   fPadx(digit.fPadx),
   fPadz(digit.fPadz),
   fTdc(digit.fTdc),
   fTdcND(digit.fTdcND),
   fAdc(digit.fAdc),
   fToT(digit.fToT)
{
  // 
  // copy ctor for AliTOFdigit object
  //

  Int_t i ;
  for ( i = 0; i < 3 ; i++)
    fTracks[i]  = digit.fTracks[i] ;

}

//______________________________________________________________________________
AliTOFdigit::AliTOFdigit(Int_t sector, Int_t plate, Int_t strip, Int_t padx,
			 Int_t padz, Int_t tdc, Int_t adc):
   fSector(sector),
   fPlate(plate),
   fStrip(strip),
   fPadx(padx),
   fPadz(padz),
   fTdc(tdc),
   fTdcND(0),
   fAdc(adc),
   fToT(0)
{
//
// Constructor for sdigit
//
}
   
//______________________________________________________________________________
void AliTOFdigit::GetLocation(Int_t *loc) const
{
//
// Get the cohordinates of the digit
// in terms of Sector - Plate - Strip - Pad
//

   loc[0] = fSector;
   loc[1] = fPlate;
   loc[2] = fStrip;
   loc[3] = fPadx;
   loc[4] = fPadz;
}

//______________________________________________________________________________
Int_t AliTOFdigit::GetTotPad() const
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
    before = AliTOFGeometry::NStripC();
    break;
  case 2:
    before = AliTOFGeometry::NStripC() +   AliTOFGeometry::NStripB();
    break;
  case 3:
    before = AliTOFGeometry::NStripC() +   AliTOFGeometry::NStripB() + AliTOFGeometry::NStripA();
    break;
  case 4:
    before = AliTOFGeometry::NStripC() + 2*AliTOFGeometry::NStripB() + AliTOFGeometry::NStripA();
    break;
  }
  
  Int_t pad = AliTOFGeometry::NpadZ()*fPadx + fPadz;
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
   if (fTracks[1]==-1)
      fTracks[1] = track;
   else if (fTracks[2]==-1)
      fTracks[2] = track;
   else
     printf("W-AliTOFdigit::AddTrack: Too many tracks (>3) that contribute to the same TOF digit\n");

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
AliTOFdigit AliTOFdigit::operator+(const AliTOFdigit &digit)
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
