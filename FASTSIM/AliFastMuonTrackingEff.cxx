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
Revision 1.1.2.1  2003/04/15 15:57:31  hristov
Merging with v3-09-08

Revision 1.1  2003/01/06 10:13:09  morsch
First commit.

*/

#include "AliFastMuonTrackingEff.h"
#include "AliMUONFastTracking.h"

ClassImp(AliFastMuonTrackingEff)


AliFastMuonTrackingEff::AliFastMuonTrackingEff() :
    AliFastResponse("Efficiency", "Muon Tracking Efficiency")
{
    SetBackground();
}

void AliFastMuonTrackingEff::Init()
{
    fFastTracking = AliMUONFastTracking::Instance();
    fFastTracking->Init(fBackground);
}



Float_t AliFastMuonTrackingEff::Evaluate(Float_t pt, Float_t theta, Float_t phi)
{
    Float_t p = pt / TMath::Sin(theta*TMath::Pi()/180.);
    Float_t eff =  fFastTracking->Efficiency(p, theta, phi, Int_t(fCharge));
    return eff;
}
