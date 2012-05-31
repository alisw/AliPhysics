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
//
// Class for fast simulation of the ALICE Muon Spectrometer
// Tracking Efficiency.
// The efficiency depends on trasverse momentum pt, polar angle theta and azimuthal angle phi.
//
// Author: Alessandro de Falco 
// alessandro.de.falco@ca.infn.it
// 

#include <TMath.h>

#include "AliFastMuonTrackingEff.h"
#include "AliMUONFastTracking.h"

ClassImp(AliFastMuonTrackingEff)


AliFastMuonTrackingEff::AliFastMuonTrackingEff() :
    AliFastResponse("Efficiency", "Muon Tracking Efficiency"),
    fBackground(1.),
    fCharge(1.),
    fFastTracking(0)
{
//
// Constructor
}

AliFastMuonTrackingEff::AliFastMuonTrackingEff(const AliFastMuonTrackingEff& eff)
    :AliFastResponse(eff),
    fBackground(1.),
    fCharge(1.),
    fFastTracking(0)
{
// Copy constructor
    eff.Copy(*this);
}

void AliFastMuonTrackingEff::Init()
{
//
// Initialization
    fFastTracking = AliMUONFastTracking::Instance();
    fFastTracking->Init(fBackground);
}



Float_t AliFastMuonTrackingEff::Evaluate(Float_t /*charge*/, Float_t pt, Float_t theta, Float_t phi)
{
//
// Evaluate the efficience for muon with 3-vector (pt, theta, phi)
    Float_t p = pt / TMath::Sin(theta*TMath::Pi()/180.);
    Float_t eff =  fFastTracking->Efficiency(p, theta, phi, Int_t(fCharge));
    return eff;
}


AliFastMuonTrackingEff& AliFastMuonTrackingEff::operator=(const  AliFastMuonTrackingEff& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

