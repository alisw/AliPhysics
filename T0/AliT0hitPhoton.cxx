/**************************************************************************
 * Copyright(c) 2001, ALICE Experiment at CERN, All rights reserved.      *
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

#include "AliT0hitPhoton.h"

ClassImp(AliT0hitPhoton)
//___________________________________________
AliT0hitPhoton::AliT0hitPhoton(Int_t shunt, Int_t track, Int_t *vol, Float_t *hit)
  :AliHit(shunt, track)
{
// Constructor for object AliT0Cerenkov
    fArray  = vol[0];
    fPmt    = vol[1];
    fRadius = hit[0];
    fMomX   = hit[1];
    fMomY   = hit[2];
    fMomZ   = hit[3];
    fTime   = hit[4];
    fEtot   = hit[5];
}
