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

/*
$Log$
 authors: Roberto Preghenella, preghenella@bo.infn.it
          with contribution from Chiara Zampolli, zampolli@bo.infn.it 
*/

////////////////////////////////////////////////////////////////////////
//                                                                    //
//     This class provides the basic object to store just-decoded     //
//     raw data                                                       //
//                                                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "AliTOFHitData.h"

ClassImp(AliTOFHitData)

AliTOFHitData::AliTOFHitData():
  TObject(),fDDLID(-1),fSlotID(-1),fACQ(-1),fChain(-1),fPS(-1),fTDC(-1),fChan(-1),fTime(-1),fTimeBin(-1),fTOT(-1),fTOTBin(-1),fDeltaBunchID(-1),fDeltaEventCounter(-1)
{
  //ctor
}
//-----------------------------------------------------------------------------

AliTOFHitData::AliTOFHitData(const AliTOFHitData &source):
  TObject(),fDDLID(-1),fSlotID(-1),fACQ(-1),fChain(-1),fPS(-1),fTDC(-1),fChan(-1),fTime(-1),fTimeBin(-1),fTOT(-1),fTOTBin(-1),fDeltaBunchID(-1),fDeltaEventCounter(-1){ 
  // copy constructor 
  for (Int_t i = 0; i < 5; i++) this->fVolume[i]=source.fVolume[i];
  this->fDDLID=source.fDDLID;
  this->fSlotID=source.fSlotID;
  this->fACQ=source.fACQ;
  this->fChain=source.fChain;
  this->fPS=source.fPS;
  this->fTDC=source.fTDC;
  this->fTOT=source.fTOT;
  this->fTOTBin=source.fTOTBin;
  this->fChan=source.fChan;
  this->fTime=source.fTime;
  this->fTimeBin=source.fTimeBin;
  this->fDeltaBunchID=source.fDeltaBunchID;
  this->fDeltaEventCounter=source.fDeltaEventCounter;
}

//-----------------------------------------------------------------------------
AliTOFHitData& AliTOFHitData::operator=(const AliTOFHitData & source) { 
// assignment operator
  for (Int_t i = 0; i < 5; i++)
    this->fVolume[i]=source.fVolume[i];
  this->fDDLID=source.fDDLID;
  this->fSlotID= source.fSlotID;
  this->fACQ= source.fACQ;
  this->fChain= source.fChain;
  this->fPS= source.fPS;
  this->fTDC= source.fTDC;
  this->fChan= source.fChan;
  this->fTime= source.fTime;
  this->fTimeBin= source.fTimeBin;
  this->fTOT= source.fTOT;
  this->fTOTBin= source.fTOTBin;
  this->fDeltaBunchID=source.fDeltaBunchID;
  this->fDeltaEventCounter=source.fDeltaEventCounter;
  return *this;
}
//----------------------------------------------------------------------------
void AliTOFHitData::SetVolume(Int_t *Volume)
{
  // setting the TOF volume index 
  fVolume[0]=Volume[0];
  fVolume[1]=Volume[1];
  fVolume[2]=Volume[2];
  fVolume[3]=Volume[3];
  fVolume[4]=Volume[4];
}


  
