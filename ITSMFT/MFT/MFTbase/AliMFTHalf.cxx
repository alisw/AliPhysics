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
/// \class AliMFTHalf
///
/// Class Building the geometry of one half of the ALICE Muon Forward Tracker
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "TGeoMatrix.h"

#include "AliLog.h"

#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTHalfDisk.h"
#include "AliMFTGeometry.h"
#include "AliMFTHalf.h"


/// \cond CLASSIMP
ClassImp(AliMFTHalf);
/// \endcond

//=============================================================================================
/// \brief Default constructor

AliMFTHalf::AliMFTHalf():
TNamed(){
  
}

//=============================================================================================
/// \brief Constructor
AliMFTHalf::AliMFTHalf(AliMFTHalfSegmentation *seg):TNamed(),
fSegmentation(seg)
{
  
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  
  SetUniqueID(fSegmentation->GetUniqueID());
  
  SetName(Form("MFT_H_%d",mftGeom->GetHalfMFTID(GetUniqueID())));
  
  
  AliInfo(Form("Creating : %s ",GetName()));

  fHalfVolume = new TGeoVolumeAssembly(GetName());
  
  CreateHalfDisks();
}


//=============================================================================================

AliMFTHalf::~AliMFTHalf() {

  
}

//=============================================================================================
/// \brief Creates the Half-disks composing the Half-MFT 
void AliMFTHalf::CreateHalfDisks(){
  AliInfo(Form("Creating  %d Half-Disk ",fSegmentation->GetNHalfDisks()));

	for (Int_t iDisk = 0 ; iDisk < fSegmentation->GetNHalfDisks(); iDisk++) {
    AliMFTHalfDiskSegmentation * halfDiskSeg = fSegmentation->GetHalfDisk(iDisk);
    
    AliMFTHalfDisk * halfDisk = new AliMFTHalfDisk(halfDiskSeg);
		 Int_t halfDiskId = AliMFTGeometry::Instance()->GetHalfDiskID(halfDiskSeg->GetUniqueID());
    fHalfVolume->AddNode(halfDisk->GetVolume(),halfDiskId,halfDiskSeg->GetTransformation());
    delete halfDisk;
  }
  

}