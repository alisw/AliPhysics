
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
/// \class AliMFTGeometryBuilder
///
/// Class describing MFT Geometry Builder
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "TGeoVolume.h"
#include "TGeoManager.h"

#include "AliLog.h"

#include "AliMFTGeometryBuilder.h"
#include "AliMFTGeometry.h"
#include "AliMFTSegmentation.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTHalf.h"
#include "AliMFTHalfCone.h"


/// \cond CLASSIMP
ClassImp(AliMFTGeometryBuilder);
/// \endcond

//=============================================================================================
/// Default constructor

AliMFTGeometryBuilder::AliMFTGeometryBuilder():
TNamed(){
  
}


//=============================================================================================

AliMFTGeometryBuilder::~AliMFTGeometryBuilder() {
  
}


//=============================================================================================
/// \brief Build the MFT Geometry
void AliMFTGeometryBuilder::BuildGeometry(){
  
  AliMFTGeometry *mftGeo = AliMFTGeometry::Instance();
  
  AliMFTSegmentation * seg = mftGeo->GetSegmentation();
  
  TGeoVolume *volMFT = new TGeoVolumeAssembly("MFT");

  for (int iHalf=0; iHalf<2; iHalf++) {
    AliMFTHalfSegmentation * halfSeg = seg->GetHalf(iHalf);
    AliMFTHalf * halfMFT = new AliMFTHalf(halfSeg);
    volMFT->AddNode(halfMFT->GetVolume(),iHalf,halfSeg->GetTransformation());
    delete halfMFT;
  }
  
  
  /// \todo Add the service, Barrel, etc Those objects will probably be defined into the COMMON ITSMFT area.
  
  AliMFTHalfCone * halfCone = new AliMFTHalfCone();
  TGeoVolumeAssembly * fHalfCone1 = halfCone->CreateHalfCone(0);
  TGeoVolumeAssembly * fHalfCone2 = halfCone->CreateHalfCone(1);
  volMFT->AddNode(fHalfCone1,1);
  volMFT->AddNode(fHalfCone2,1);
  
  gGeoManager->GetVolume("ALIC")->AddNode(volMFT,0);
    
}
