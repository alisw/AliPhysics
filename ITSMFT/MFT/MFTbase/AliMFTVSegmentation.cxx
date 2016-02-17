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
/// \class AliMFTVSegmentation
///
/// Abstract base class for MFT Segmentation description
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "AliMFTVSegmentation.h"
#include "TGeoMatrix.h"

/// \cond CLASSIMP
ClassImp(AliMFTVSegmentation);
/// \endcond

//-----------------------------------------------------------------------------

AliMFTVSegmentation::AliMFTVSegmentation():TNamed(),
fTransformation(new TGeoCombiTrans())
{
  /// Default constructor
}
//-----------------------------------------------------------------------------

AliMFTVSegmentation::AliMFTVSegmentation(const AliMFTVSegmentation& input): TNamed(),
fTransformation(input.fTransformation)
{
  /// Copy constructor
  
  SetUniqueID(input.GetUniqueID());
  SetName(input.GetName());
  
}

//-----------------------------------------------------------------------------
void AliMFTVSegmentation::SetRotationAngles(const Double_t *ang){
  /// Set Rotation Angles
  if(!fTransformation) fTransformation = new TGeoCombiTrans();
  TGeoRotation *rot = new TGeoRotation();
  rot->SetAngles(ang[0], ang[1], ang[2]); // all angles in degrees
  fTransformation->SetRotation(rot);
  
};