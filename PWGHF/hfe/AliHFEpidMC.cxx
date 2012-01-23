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
//
// Class for TRD PID
// Implements the abstract base class AliHFEpidBase
// Make PID does the PID decision
// Class further contains TRD specific cuts and QA histograms
//
// Authors:
//   Markus Fasel <M.Fasel@gsi.de>
//
#include <TParticle.h>

#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
//#include "AliVParticle.h"

#include "AliHFEpidMC.h"
#include "AliHFEtools.h"

ClassImp(AliHFEpidMC)

//___________________________________________________________________
AliHFEpidMC::AliHFEpidMC():
  AliHFEpidBase()
{
  //
  // Default constructor
  //
}

//___________________________________________________________________
AliHFEpidMC::AliHFEpidMC(const Char_t *name):
  AliHFEpidBase(name)
{
  //
  // Standard constructor
  //
}

//___________________________________________________________________
Bool_t AliHFEpidMC::InitializePID(Int_t /*run*/){
  // 
  // Implementation of the framework function InitializePID
  // Not yet anything to implement in case of MC PID
  //
  return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidMC::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager * /*pidqa*/) const {
  //
  // returns MC PDG Code
  // Functionality implemented in the base class
  // (necessary for PID QA)
  //
  return AliHFEtools::GetPdg(track->GetRecTrack());
}
