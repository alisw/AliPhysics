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
//---------------------------------------------------------------------

#include <Riostream.h>
#include <TFile.h>

#include "AliJetFinder.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"
#include "AliJetUnitArray.h"
#include "AliJetReaderHeader.h"
#include "AliJetHeader.h"
#include "AliJetReader.h"
#include "AliAODJetEventBackground.h"

ClassImp(AliJetFinder)

AliJetFinder::AliJetFinder():
    fReader(0x0),
    fHeader(0x0),
    fAODjets(0x0),
    fNAODjets(0),
    fAODEvBkg(0),
    fDebug(0)
{
  //
  // Constructor
  //
  fAODjets = 0;
}

////////////////////////////////////////////////////////////////////////
AliJetFinder::~AliJetFinder()
{
  //
  // Destructor
  //
}



////////////////////////////////////////////////////////////////////////
void AliJetFinder::WriteRHeaderToFile()
{
  // write reader header
    AliJetReaderHeader *rh = fReader->GetReaderHeader();
    rh->Write();
}


////////////////////////////////////////////////////////////////////////
void AliJetFinder::ConnectTree(TTree* tree, TObject* data)
{
    // Connect the input file
    fReader->ConnectTree(tree, data);
}

////////////////////////////////////////////////////////////////////////
void AliJetFinder::WriteHeaders()
{
    // Write the Headers
    TFile* f = new TFile("jets_local.root", "recreate");
    WriteRHeaderToFile();
    WriteJHeaderToFile();
    f->Close();
}

////////////////////////////////////////////////////////////////////////
Bool_t AliJetFinder::ProcessEvent()
{
  //
  // Process one event
  // Charged only jets
  //

  Bool_t ok = fReader->FillMomentumArray();
  if (!ok) return kFALSE;
  // Jets
  FindJets(); // V1
  Reset();  
  return kTRUE;
}

////////////////////////////////////////////////////////////////////////
Bool_t AliJetFinder::ProcessEvent2()
{
  //
  // Process one event
  // Charged only or charged+neutral jets
  //

  TRefArray* ref = new TRefArray();
  Bool_t procid = kFALSE;
  Bool_t ok = fReader->ExecTasks(procid,ref);

  // Delete reference pointer  
  if (!ok) {delete ref; return kFALSE;}
  // Jets
  FindJets();
  
  Int_t nEntRef = ref->GetEntries();

  for(Int_t i=0; i<nEntRef; i++)
    { 
      // Reset the UnitArray content which were referenced
      ((AliJetUnitArray*)ref->At(i))->SetUnitTrackID(0);
      ((AliJetUnitArray*)ref->At(i))->SetUnitEnergy(0.);
      ((AliJetUnitArray*)ref->At(i))->SetUnitCutFlag(kPtSmaller);
      ((AliJetUnitArray*)ref->At(i))->SetUnitCutFlag2(kPtSmaller);
      ((AliJetUnitArray*)ref->At(i))->SetUnitSignalFlag(kBad);
      ((AliJetUnitArray*)ref->At(i))->SetUnitSignalFlagC(kTRUE,kBad);
      ((AliJetUnitArray*)ref->At(i))->SetUnitDetectorFlag(kTpc);
      ((AliJetUnitArray*)ref->At(i))->SetUnitFlag(kOutJet);
      ((AliJetUnitArray*)ref->At(i))->ClearUnitTrackRef();

      // Reset process ID
      AliJetUnitArray* uA = (AliJetUnitArray*)ref->At(i);
      uA->ResetBit(kIsReferenced);
      uA->SetUniqueID(0);     
    }

  // Delete the reference pointer
  ref->Delete();
  delete ref;

  Reset();

  return kTRUE;
}


void AliJetFinder::AddJet(AliAODJet p)
{
// Add new jet to the list
  new ((*fAODjets)[fNAODjets++]) AliAODJet(p);
}

void AliJetFinder::ConnectAOD(const AliAODEvent* aod)
{
// Connect to the AOD
    fAODjets = aod->GetJets();
    fAODEvBkg = (AliAODJetEventBackground*)(aod->FindListObject(AliAODJetEventBackground::StdBranchName()));
}

////////////////////////////////////////////////////////////////////////
void AliJetFinder::ConnectAODNonStd(AliAODEvent* aod,const char *bname)
{

  fAODjets = dynamic_cast<TClonesArray*>(aod->FindListObject(bname));
  fAODEvBkg = (AliAODJetEventBackground*)(aod->FindListObject(Form("%s_%s",AliAODJetEventBackground::StdBranchName(),bname)));
  // how is this is reset? Cleared? -> by the UserExec!!
}

