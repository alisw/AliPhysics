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

#include "AliMUONSegmentPosition.h"

//___________________________________________
ClassImp(AliMUONSegmentPosition)

//
//___________________________________________
AliMUONSegmentPosition::AliMUONSegmentPosition() : TNamed()
{
  fChannelId = 0;
  fX = 0.;
  fY = 0.;  
}
//___________________________________________
AliMUONSegmentPosition::AliMUONSegmentPosition(const Int_t channelId, const Float_t x, const Float_t y, const Int_t cathode) : TNamed()
{
  char name[10];
  sprintf(name,"%5.2f-%5.2f",x,y);
  fName = name;
  fTitle= name;
  fChannelId = channelId;
  fX = x;
  fY = y;
  fCathode=cathode;
}
//_______________________________________________
AliMUONSegmentPosition::~AliMUONSegmentPosition()
{

}
//___________________________________________
Int_t AliMUONSegmentPosition::Compare(const TObject *obj) const
{
  AliMUONSegmentPosition * myobj = ( AliMUONSegmentPosition *) obj;
  return (fChannelId > myobj->GetChannelId()) ? 1 : -1;
}
//___________________________________________
void AliMUONSegmentPosition::Print() const
{
  printf("%s id=%d x=%f y=%f cathode=%d\n",fName.Data(),fChannelId, fX, fY,fCathode);   
}
