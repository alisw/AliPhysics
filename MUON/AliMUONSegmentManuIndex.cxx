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
//  Segment element indexing in a detection element for electronics   
//        Gines MARTINEZ, SUBATECH July 04                
//  This class is the basic component of 
//  AliMUONSegmentationDetectionElement and contains al the 
//  info about a segment (pad or strip):
//          Id-indetectionelement,  #manu, #manuchannel 
//  Detailed information in Alice Technical Note xxxxxxxx (2004)
//====================================================================

#include <TString.h>

#include "AliMUONSegmentManuIndex.h"
#include "AliLog.h"

//___________________________________________
ClassImp(AliMUONSegmentManuIndex)

//
//___________________________________________
AliMUONSegmentManuIndex::AliMUONSegmentManuIndex() 
{
  //Default constructor
  fChannelId= 0;; // Id of the channel within the detection element
  fManuId= 0;; // Manu id in the detection element
  fBusPatchId= 0;; // BusPatchId in the detection element up to 4 for slats
  fManuChannelId= 0;; 
}
//___________________________________________
AliMUONSegmentManuIndex::AliMUONSegmentManuIndex(Int_t channelId, Int_t manuId, Int_t busPatchId,  Int_t manuChannelId) : TNamed()
{  
  // Constructor to be used
  fName = Name(manuId, manuChannelId).Data();
  fTitle= Name(manuId, manuChannelId).Data();
  fChannelId     = channelId;
  fManuId        = manuId;
  fBusPatchId    = busPatchId;
  fManuChannelId = manuChannelId;  
}
//_______________________________________________
AliMUONSegmentManuIndex::~AliMUONSegmentManuIndex()
{
  // Destructor
}
//___________________________________________
Int_t AliMUONSegmentManuIndex::Compare(const TObject *obj) const
{
  // Comparison of two AliMUONSegmentManuIndex objects
 AliMUONSegmentManuIndex * myobj = ( AliMUONSegmentManuIndex *) obj;
  return (fChannelId > myobj->GetChannelId()) ? 1 : -1;
}
//___________________________________________
TString AliMUONSegmentManuIndex::Name(Int_t manuid, Int_t manuchannel) 
{
  // Definition of the name for TMap indexing
  char name[15];
  Int_t absid = manuid*64 + manuchannel;
  sprintf(name,"%d",absid);
  return TString(name);
}
//___________________________________________
void AliMUONSegmentManuIndex::Print(const char* /*opt*/) const
{
  // Printing AliMUONSegmentManuIndex information
  AliInfo(Form("Name=%s Id=%d BusPatch=%d ManuId=%d ManuChannelId=%d\n",fName.Data(),fChannelId,fBusPatchId,fManuId,fManuChannelId));   
}
