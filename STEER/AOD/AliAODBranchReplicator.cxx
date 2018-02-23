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

//
// Base class of an object used for the replication 
// (and possibly filtering) of one (or several) AOD branches.
//
// Author: L. Aphecetche (Subatech)
//
// Intended usage is to be able to produce, besides the standard AOD (AliAOD.root)
// some light-weight AODs, by filtering (or skipping completely) the unwanted branches.
// 
// Exemple usage (pseudo-code) :
// 
// AliAODHandler* aodH = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
//
// AliAODExtension* ext = aodH->AddFilteredAOD("AliAOD.filtered.root","filtered AOD");
//
// ext->DisableReferences();
//  
// ext->FilterBranch("cascades",0x0); // AliAOD.filtered.root will *not* contain the cascades branch
//
//  
// AliAODBranchReplicator* murep = new AliAODMuonReplicator("MuonReplicator",
//                                                           "remove non muon tracks and non primary or pileup vertices",
//                                                           new AliAnalysisNonMuonTrackCuts,
//                                                           new AliAnalysisNonPrimaryVertices);
// ext->FilterBranch("tracks",murep);   // both the tracks and vertices branches 
// ext->FilterBranch("vertices",murep); // will be filtered by the MuonReplicator
//  

#include "AliAODBranchReplicator.h"

ClassImp(AliAODBranchReplicator)

//______________________________________________________________________________
AliAODBranchReplicator::AliAODBranchReplicator(const char* name, const char* title)
: TNamed(name,title)
{
  // default ctor (nop)
}

//______________________________________________________________________________
AliAODBranchReplicator::~AliAODBranchReplicator()
{
  // dtor (nop)
}
