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

//===================================================================
//  Segment element position in local coordinates of the detection element   
//        Gines MARTINEZ, SUBATECH July 04                
//  This class is one of the basic component of 
//  AliMUONSegmentationDetectionElement and contains al the 
//  info about a segment (pad or strip):
//          Id-indetectionelement,  x_local, y_local 
//  Detailed information in Alice Technical Note xxxxxxxx (2004)
//====================================================================

#include <TMath.h>

#include "AliMUONSegmentPosition.h"

//___________________________________________
ClassImp(AliMUONSegmentPosition)

Float_t AliMUONSegmentPosition::fUnit = 0.3 ; // 3mm unit for generation of segmenposition name

//
//___________________________________________
AliMUONSegmentPosition::AliMUONSegmentPosition() : TNamed()
{
//Default constructor
  fChannelId = 0;
  fX = 0.;
  fY = 0.;  
}
//___________________________________________
AliMUONSegmentPosition::AliMUONSegmentPosition(const Int_t channelId, const Float_t x, const Float_t y, const Int_t cathode) : TNamed()
{
  // Constructor to be used
  fName = Name(x,y,cathode);
  fTitle= Name(x,y,cathode);
  fChannelId = channelId;
  fX = x;
  fY = y;
  fCathode=cathode;
  fPadSizeX=0.;
  fPadSizeY=0.;  
}
//_______________________________________________
AliMUONSegmentPosition::~AliMUONSegmentPosition()
{
 // Destructor
}
//___________________________________________
Int_t AliMUONSegmentPosition::Compare(const TObject *obj) const
{
  // Comparison of two AliMUONSegmentPosition objects
  AliMUONSegmentPosition * myobj = ( AliMUONSegmentPosition *) obj;
  return (fChannelId > myobj->GetChannelId()) ? 1 : -1;
}

//___________________________________________
Float_t AliMUONSegmentPosition::Distance(Float_t x, Float_t y)
{
  return TMath::Sqrt( (fX-x)*(fX-x) + (fY-y)*(fY-y) ) ;
}
//___________________________________________
TString AliMUONSegmentPosition::Name(Float_t x, Float_t y, Int_t cathode) 
{
  // Definition of the name of AliMUONSegmentPosition
  // Our convention since the smaller pad pich is 5 mm is to choice a 3mm unit:
  // So for a position pair x,y and cathode plane icathode the name will be:
  // xp = TMath::Nint(x*10./3.);
  // yp = TMath::Nint(y+10./3.);
  // sprintf(name,"%d-%d-%d",xp,yp,cathode);
  Int_t xp = TMath::Nint(x/fUnit);
  Int_t yp = TMath::Nint(y/fUnit);
  char name[15];
  sprintf(name,"%d-%d-%d",xp,yp,cathode);
  return TString(name);
}
//___________________________________________
void AliMUONSegmentPosition::Print() const
{
  // Printing AliMUONSegmentManuIndex information
  Info("Print","Name=%s Id=%d X=%f Y=%f Cathode=%d\n",fName.Data(),fChannelId, fX, fY,fCathode);   
}
