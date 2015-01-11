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

//-----------------------------------------------------------------------------
// Class AliMUONChamber
// -----------------------
// MUON tracking chamber class
// now only providing DisIntegration function
//-----------------------------------------------------------------------------

// --- ROOT includes ---
#include <TRandom.h>
#include <TMath.h>
#include "AliRun.h"


// --- MUON includes ---
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONHit.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONChamber)	
/// \endcond

//_______________________________________________________
AliMUONChamber::AliMUONChamber()
  : TObject(), 
    fId(0),
    fCurrentCorrel(1), // to avoid mistakes if ChargeCorrelInit is not called
    fResponse(0),
    fMUON(0)
{
/// Default constructor

  AliDebug(1, Form("default (empty) ctor this = %p", this));
}

//_______________________________________________________
AliMUONChamber::AliMUONChamber(Int_t id) 
  : TObject(), 
    fId(id),
    fCurrentCorrel(1), // to avoid mistakes if ChargeCorrelInit is not called
    fResponse(0),
    fMUON(0)
{
/// Standard constructor

    // muon
    fMUON = (AliMUON*)gAlice->GetModule("MUON");
    if (!fMUON) {
      AliFatal("MUON detector not defined.");
      return;
    }  

  AliDebug(1, Form("ctor this = %p", this) ); 
}

//_______________________________________________________
AliMUONChamber::~AliMUONChamber() 
{
/// Destructor

  AliDebug(1, Form("dtor this = %p", this));
  delete fResponse;
}

//_____________________________________________________
void AliMUONChamber::ChargeCorrelationInit() 
{
/// Initialisation of charge correlation for current hit
/// the value is stored, and then used by Disintegration

  // exponential is here to avoid eventual problems in 0 
  // factor 2 because chargecorrel is q1/q2 and not q1/qtrue
    fCurrentCorrel = TMath::Exp(gRandom->Gaus(0,fResponse->ChargeCorrel()/2));
}

//_____________________________________________________________________________
void
AliMUONChamber::SetResponseModel(const AliMUONResponse& thisResponse)
{
  delete fResponse;
  fResponse = static_cast<AliMUONResponse*>(thisResponse.Clone());
}

