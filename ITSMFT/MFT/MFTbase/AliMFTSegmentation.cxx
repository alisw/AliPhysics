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

//=====================================================================================================================
//
//      Class for the virtual segmentation of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//=====================================================================================================================

#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "AliMFTPlane.h"
#include "AliMFTLadderSegmentation.h"
#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTSegmentation.h"
#include "TXMLEngine.h"

ClassImp(AliMFTSegmentation)

//=====================================================================================================================

AliMFTSegmentation::AliMFTSegmentation(): 
  TNamed(),
  fMFTHalves(NULL)
{ 

  // default constructor

}

//=====================================================================================================================

AliMFTSegmentation::AliMFTSegmentation(const Char_t *nameGeomFile): 
  TNamed(),
  fMFTHalves(NULL)
{ 
  AliDebug(2,"Entering");

  // constructor
  
  fMFTHalves = new TClonesArray("AliMFTHalfSegmentation", 2);
  fMFTHalves -> SetOwner(kTRUE);


  
  AliMFTHalfSegmentation *halfBottom = new AliMFTHalfSegmentation(nameGeomFile, kBottom);
  AliMFTHalfSegmentation *halfTop    = new AliMFTHalfSegmentation(nameGeomFile, kTop);

  new ((*fMFTHalves)[kBottom]) AliMFTHalfSegmentation(*halfBottom);
  new ((*fMFTHalves)[kTop])    AliMFTHalfSegmentation(*halfTop);

  delete halfBottom;
  delete halfTop;

  AliDebug(1,"MFT segmentation set!\n");

}

//=====================================================================================================================

AliMFTSegmentation::~AliMFTSegmentation() {

  if (fMFTHalves) fMFTHalves->Delete();
  delete fMFTHalves; 
  
}

//=====================================================================================================================

void AliMFTSegmentation::Clear(const Option_t* /*opt*/) {

  if (fMFTHalves) fMFTHalves->Delete();
  delete fMFTHalves; 
  fMFTHalves = NULL;
  
}

//=====================================================================================================================
Bool_t AliMFTSegmentation::Hit2PixelID(Double_t xHit, Double_t yHit, Double_t zHit, Int_t half, Int_t disk, Int_t ladder, Int_t sensor, Int_t &xPixel, Int_t &yPixel){

  Double_t master[3] = {xHit, yHit, zHit};
  Double_t local[3];
  AliMFTHalfSegmentation * halfSeg = ((AliMFTHalfSegmentation*)fMFTHalves->At(half));
  if(!halfSeg) return kFALSE;
  AliMFTHalfDiskSegmentation * diskSeg = halfSeg->GetHalfDisk(disk);
  if(!diskSeg) return kFALSE;
  AliMFTLadderSegmentation * ladderSeg = diskSeg->GetLadder(ladder);
  if(!ladderSeg) return kFALSE;
  AliMFTChipSegmentation * chipSeg = ladderSeg->GetSensor(sensor);
  if(!chipSeg) return kFALSE;

  AliDebug(1,Form(" ->  Global %f %f %f",master[0],master[1],master[2]));
  halfSeg->GetTransformation()->MasterToLocal(master, local);
  AliDebug(2,Form(" ->  Half %f %f %f",local[0],local[1],local[2]));
  for (int i=0; i<3; i++) master[i] = local[i];
  diskSeg->GetTransformation()->MasterToLocal(master, local);
  AliDebug(2,Form(" ->  Disk %f %f %f",local[0],local[1],local[2]));
  for (int i=0; i<3; i++) master[i] = local[i];
  ladderSeg->GetTransformation()->MasterToLocal(master, local);
  AliDebug(2,Form(" ->  Ladder %f %f %f",local[0],local[1],local[2]));
  for (int i=0; i<3; i++) master[i] = local[i];
  chipSeg->GetTransformation()->MasterToLocal(master, local);
  AliDebug(1,Form(" ->  Chip Pos %f %f %f",local[0],local[1],local[2]));
  
  
   return (chipSeg->Hit2PixelID(local[0], local[1], xPixel, yPixel));

}
//=====================================================================================================================

Int_t AliMFTSegmentation::GetDetElemLocalID(Int_t half, Int_t disk, Int_t ladder, Int_t sensor) const {

  Int_t localId =0;
  
  
  if (half==1) localId += GetHalf(0)->GetHalfDisk(disk)->GetNChips();
  
  for (int iLad=0; iLad<GetHalf(half)->GetHalfDisk(disk)->GetNLadders(); iLad++) {
    if (iLad<ladder) localId += GetHalf(half)->GetHalfDisk(disk)->GetLadder(iLad)->GetNSensors();
    else{
      for (int iSens=0; iSens<GetHalf(half)->GetHalfDisk(disk)->GetLadder(iLad)->GetNSensors(); iSens++) {
        if(iSens==sensor) return localId;
        localId++;
     }
    }
  }
  return -1;
}
