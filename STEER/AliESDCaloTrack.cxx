/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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
 *                                                                        *
 **************************************************************************/

/* $Id$ */
#include "AliESDCaloTrack.h"
#include "AliPHOSRecParticle.h"
#include "AliEMCALRecParticle.h"

//-------------------------------------------------------------------------
//   Class AliESDCaloTrack
//   This is the class to deal with during the physical analysis of data
//   It converts calorimeter (PHOS or EMCAL) reconstructed particles   
//   into event summary data object
//-------------------------------------------------------------------------

ClassImp(AliESDCaloTrack)

AliESDCaloTrack::AliESDCaloTrack(AliPHOSRecParticle* recpart)
{
  // Convert AliPHOSRecParticle to AliESDCaloTrack
  fPx = recpart->Px();
  fPy = recpart->Py();
  fPz = recpart->Pz();
}

AliESDCaloTrack::AliESDCaloTrack(AliEMCALRecParticle* recpart)
{
  // Convert AliEMCALRecParticle to AliESDCaloTrack
  fPx = recpart->Px();
  fPy = recpart->Py();
  fPz = recpart->Pz();
}

