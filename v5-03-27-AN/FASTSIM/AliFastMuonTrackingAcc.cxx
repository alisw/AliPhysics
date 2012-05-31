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

/* $Id$ */

// Realisation of AliFastResponse for the
// fast simulation of the muon spectrometer acceptance.
// The acceptance depends on the muon 3-vector which can be passed as (pt, theta, phi), 
// where pt is the transverse momentum, theta the polar angle and phi the azimuthal angle.
// Author: Andreas Morsch
// andreas.morsch@cern.ch 

#include <TMath.h>

#include "AliFastMuonTrackingAcc.h"
#include "AliMUONFastTracking.h"

ClassImp(AliFastMuonTrackingAcc)


AliFastMuonTrackingAcc::AliFastMuonTrackingAcc() :
    AliFastResponse("Acceptance", "Muon Tracking Acceptance"),
    fBackground(1.),
    fCharge(1.),
    fFastTracking(0)
{
    // Default Constructor
}

AliFastMuonTrackingAcc::AliFastMuonTrackingAcc(const AliFastMuonTrackingAcc & acc)
    :AliFastResponse(acc),
     fBackground(1.),
     fCharge(1.),
     fFastTracking(0)
{
// Copy constructor
    acc.Copy(*this);
}

void AliFastMuonTrackingAcc::Init()
{
    fFastTracking = AliMUONFastTracking::Instance();
    fFastTracking->Init(fBackground);
}



Float_t AliFastMuonTrackingAcc::Evaluate(Float_t /*charge*/, Float_t pt, Float_t theta, Float_t phi)
{
// Evaluate the tracking acceptance for 3-vector pt, theta, phi
    Float_t p = pt / TMath::Sin(theta*TMath::Pi()/180.);
    Float_t eff =  fFastTracking->Acceptance(p, theta, phi, Int_t(fCharge));
    return eff;
}

AliFastMuonTrackingAcc& AliFastMuonTrackingAcc::operator=(const  AliFastMuonTrackingAcc& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}
