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
//  Segment element indexing in a detection element    
//        Gines MARTINEZ, SUBATECH July 04                
//  This class is the basic component of 
//  AliMUONSegmentationDetectionElement and contains al the 
//  info about a segment (pad or strip):
//          Id-indetectionelement,  ix ,iy 
//  Detailed information in Alice Technical Note xxxxxxxx (2004)
//====================================================================

#include "AliMUONSegmentIndex.h"

//___________________________________________
ClassImp(AliMUONSegmentIndex)

//
//___________________________________________
AliMUONSegmentIndex::AliMUONSegmentIndex() : TNamed()
{
  fChannelId = 0;
  fPadX = 0;
  fPadY = 0;
  fCathode=0;  
}
//___________________________________________
AliMUONSegmentIndex::AliMUONSegmentIndex(const Int_t channelId, const Int_t padX, const Int_t padY, const Int_t cathode) : TNamed()
{
  char name[10];
  sprintf(name,"%d-%d",padX,padY);
  fName = name;
  fTitle = name;
  fChannelId = channelId;
  fPadX = padX;
  fPadY = padY;
  fCathode=cathode;  
}
//_______________________________________________
AliMUONSegmentIndex::~AliMUONSegmentIndex()
{

}
//___________________________________________
Int_t AliMUONSegmentIndex::Compare(const TObject *obj) const
{
  AliMUONSegmentIndex * myobj = ( AliMUONSegmentIndex *) obj;
  return (fChannelId > myobj->GetChannelId()) ? 1 : -1;
}
//___________________________________________
void AliMUONSegmentIndex::Print() const
{
  printf("%s id=%d ix=%d iy=%d cathode=%d\n",fName.Data(),fChannelId,fPadX,fPadY,fCathode);   
}
