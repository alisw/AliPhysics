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

//---------------------------------------------------------------------
// Jet finder base class
// manages the search for jets 
// Authors: jgcn@mda.cinvestav.mx
//          andreas.morsch@cern.ch
//          magali.estienne@subatech.in2p3.fr
//          alexandre.shabetai@cern.ch
//---------------------------------------------------------------------

#include <TFile.h>

#include "AliJetFinder.h"
#include "AliUA1JetHeaderV1.h"
#include "AliAODJetEventBackground.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"

ClassImp(AliJetFinder)

///////////////////////////////////////////////////////////////////////

AliJetFinder::AliJetFinder():
  fHeader(0x0),
  fAODjets(0x0),
  fNAODjets(0),
  fAODEvBkg(0),
  fDebug(0),
  fCalTrkEvent(0x0)
{
  // Constructor
}

//-----------------------------------------------------------------------
AliJetFinder::~AliJetFinder()
{
  // Destructor
}

//-----------------------------------------------------------------------
void AliJetFinder::WriteHeader()
{
  // Write the Headers
  TFile* f = new TFile("jets_local.root", "recreate");
  WriteHeaderToFile();
  f->Close();

}

//-----------------------------------------------------------------------
void AliJetFinder::WriteHeaderToFile()
{
  // write reader header
  AliJetHeader *rh = GetJetHeader();
  rh->Write();

}

//-----------------------------------------------------------------------
Bool_t AliJetFinder::ProcessEvent()
{
  // Process one event

  // Find jets
  FindJets();

  Reset();
  return kTRUE;

}

//-----------------------------------------------------------------------
void AliJetFinder::AddJet(AliAODJet p)
{
  // Add new jet to the list
  if (fAODjets) { new ((*fAODjets)[fNAODjets++]) AliAODJet(p);}
  else { Warning("AliJetFinder::AddJet(AliAODJet p)","fAODjets is null!");}

}

//-----------------------------------------------------------------------
void AliJetFinder::ConnectAOD(const AliAODEvent* aod)
{
  // Connect to the AOD
  fAODjets = aod->GetJets();
  fAODEvBkg = (AliAODJetEventBackground*)(aod->FindListObject(AliAODJetEventBackground::StdBranchName()));

}

//-----------------------------------------------------------------------
void AliJetFinder::ConnectAODNonStd(AliAODEvent* aod,const char *bname)
{
  // Connect non standard AOD jet and jet background branches 
  fAODjets = dynamic_cast<TClonesArray*>(aod->FindListObject(bname));
  fAODEvBkg = (AliAODJetEventBackground*)(aod->FindListObject(Form("%s_%s",AliAODJetEventBackground::StdBranchName(),bname)));
  // how is this is reset? Cleared? -> by the UserExec!!

}

