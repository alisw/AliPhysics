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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMFTHalfDisk
///
/// Class describing geometry of one half of a MFT disk
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"

#include "AliMFTHalfDisk.h"
#include "AliMFTGeometry.h"
#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTLadder.h"
#include "AliMFTHeatExchanger.h"
#include "AliMFTSupport.h"

/// \cond CLASSIMP
ClassImp(AliMFTHalfDisk);
/// \endcond

//=============================================================================================
/// \brief Default constructor

AliMFTHalfDisk::AliMFTHalfDisk():
TNamed(), fMFTSupport(NULL),fMFTHeatExchanger(NULL),fSegmentation(NULL){
  
}
//=============================================================================================
/// \brief Constructor
AliMFTHalfDisk::AliMFTHalfDisk(AliMFTHalfDiskSegmentation *segmentation):TNamed(segmentation->GetName(),segmentation->GetName()),
  fMFTSupport(NULL),
  fMFTHeatExchanger(NULL),
  fSegmentation(segmentation)
{
//  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  SetUniqueID(fSegmentation->GetUniqueID());
//  Int_t halfDiskID = mftGeom->GetHalfDiskID(GetUniqueID());
//  SetName(Form("D%d",halfDiskID));

  AliDebug(1,Form("Creating Half-Disk: %s Unique ID = %d ", GetName(), GetUniqueID()));

  fHalfDiskVolume = new TGeoVolumeAssembly(GetName());
  

  // Building Heat Exchanger Between faces
	TGeoVolumeAssembly * heatExchangerVol = CreateHeatExchanger();
	fHalfDiskVolume->AddNode(heatExchangerVol,1);
	
  //   Building Front Face of the Half Disk
  CreateLadders();

}


//=============================================================================================

AliMFTHalfDisk::~AliMFTHalfDisk() {
  delete fMFTSupport;
  delete fMFTHeatExchanger;
}

//=============================================================================================
/// \brief Build Heat exchanger
/// \return Pointer to the volume assembly holding the heat exchanger

TGeoVolumeAssembly * AliMFTHalfDisk::CreateHeatExchanger(){
  
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();

  AliMFTHeatExchanger * fMFTHeatExchanger = new AliMFTHeatExchanger();
  
  TGeoVolumeAssembly * vol = fMFTHeatExchanger->Create(mftGeom->GetHalfMFTID(GetUniqueID()), mftGeom->GetHalfDiskID(GetUniqueID()));
  
  return vol;
  
}

//=============================================================================================
/// \brief Build Ladders on the Half-disk
void AliMFTHalfDisk::CreateLadders(){
  AliDebug(1,"Start Building Ladders" );
	for (Int_t iLadder=0; iLadder<fSegmentation->GetNLadders(); iLadder++) {

    AliMFTLadderSegmentation * ladderSeg = fSegmentation->GetLadder(iLadder);
    if(!ladderSeg) AliFatal(Form("No Segmentation found for ladder %d ",iLadder));
    AliMFTLadder * ladder = new AliMFTLadder(ladderSeg);
    TGeoVolume * ladVol = ladder->CreateVolume();
    
    // Position of the center on the ladder volume in the ladder coordinate system
    TGeoBBox* shape = (TGeoBBox*)ladVol->GetShape();
    Double_t center[3];
    center[0] = shape->GetDX();
    center[1] = shape->GetDY();
    center[2] = shape->GetDZ();

    Double_t master[3];
    ladderSeg->GetTransformation()->LocalToMaster(center, master);
    Int_t ladderId = AliMFTGeometry::Instance()->GetLadderID(ladderSeg->GetUniqueID());
    
    fHalfDiskVolume->AddNode(ladVol,ladderId,new TGeoCombiTrans(master[0],master[1],master[2],ladderSeg->GetTransformation()->GetRotation()));
    
    delete ladder;
  }

}
