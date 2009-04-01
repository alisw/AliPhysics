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
// $MpId: AliMpConstants.cxx,v 1.11 2006/05/24 13:58:29 ivana Exp $
// Category: basic

//-----------------------------------------------------------------------------
// Class AliMpConstants
// --------------------
// Class for globally used constants definition.
// Included in AliRoot: 2003/05/02
// Authors: David Guez, Ivana Hrivnacova; IPN Orsay
//-----------------------------------------------------------------------------

#include "AliMpConstants.h"

#include <TMath.h>

/// \cond CLASSIMP
ClassImp(AliMpConstants)
/// \endcond

// static data
const Double_t AliMpConstants::fgkLengthTolerance = 1e-05; // 0.1 mum
const Double_t AliMpConstants::fgkLengthStep = 0.1; // 1 mm
const Int_t    AliMpConstants::fgkStartPadIndex = 1;
const Int_t    AliMpConstants::fgkNofCathodes = 2;
const Int_t    AliMpConstants::fgkNofChambers = 14;
const Int_t    AliMpConstants::fgkNofTrackingChambers = 10;
const Int_t    AliMpConstants::fgkNofGeomModules = 20;
const Int_t    AliMpConstants::fgkNofLocalBoards = 234;
const Int_t    AliMpConstants::fgkTotalNofLocalBoards = 242;
const Int_t    AliMpConstants::fgkNonBendingManuMask(1<<10);
const Int_t    AliMpConstants::fgkManuNofChannels(64);
const Int_t    AliMpConstants::fgkLocalBoardNofChannels(16);

//_____________________________________________________________________________
AliMpConstants::AliMpConstants()
  : TObject() 
{
/// Default constructor  
}

//_____________________________________________________________________________
AliMpConstants::~AliMpConstants() 
{
///Destructor
}

//_____________________________________________________________________________
Bool_t  AliMpConstants::IsEqual(Double_t length1, Double_t length2)
{
/// Compare lengths within the length tolerance.

  return TMath::Abs(length1 - length2) < fgkLengthTolerance;
}  

//_____________________________________________________________________________
Bool_t  AliMpConstants::IsEqual(Double_t v1x, Double_t v1y, 
                                Double_t v2x, Double_t v2y)
{
/// Compare x, y vector coordinates within the length tolerance.

  return (  TMath::Abs(v1x - v2x) 
          + TMath::Abs(v1y - v2y)) < 2.*fgkLengthTolerance;
}

//_____________________________________________________________________________
Int_t AliMpConstants::ManuMask(AliMp::PlaneType planeType)
{
/// The manuIDs get an offset if they are in the non-bending plane

  return ( planeType == AliMp::kNonBendingPlane ) ? fgkNonBendingManuMask : 0;
}

//_____________________________________________________________________________
Int_t AliMpConstants::NofTriggerChambers() 
{ 
/// Return number of trigger chambers

  return fgkNofChambers - fgkNofTrackingChambers;
}
