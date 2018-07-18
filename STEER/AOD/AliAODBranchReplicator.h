#ifndef ALIAODBRANCHREPLICATOR_H
#define ALIAODBRANCHREPLICATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

/// \class AliAODBranchReplicator
/// \brief Base class of an object used for the replication
///
/// Base class of an object used for the replication
/// (and possibly filtering) of one (or several) AOD branches
///
/// Intended usage is to be able to produce, besides the standard AOD (AliAOD.root)
/// some light-weight AODs, by filtering (or skipping completely) the unwanted branches.
///
/// Exemple usage (pseudo-code) :
/// ~~~{.cxx}
/// AliAODHandler* aodH = AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler();
/// AliAODExtension* ext = aodH->AddFilteredAOD("AliAOD.filtered.root","filtered AOD");
/// ext->DisableReferences();
/// ext->FilterBranch("cascades",0x0); // AliAOD.filtered.root will *not* contain the cascades branch
/// AliAODBranchReplicator* murep = new AliAODMuonReplicator("MuonReplicator",
///                                                           "remove non muon tracks and non primary or pileup vertices",
///                                                           new AliAnalysisNonMuonTrackCuts,
///                                                           new AliAnalysisNonPrimaryVertices);
/// ext->FilterBranch("tracks",murep);   // both the tracks and vertices branches
/// ext->FilterBranch("vertices",murep); // will be filtered by the MuonReplicator
/// ~~~
/// \author L. Aphecetche (Subatech)

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class AliAODEvent;

class AliAODBranchReplicator : public TNamed
{
public:
  AliAODBranchReplicator(const char* name="", const char* title="");
  
  virtual ~AliAODBranchReplicator();

  /// Return the list of object we manage
  virtual TList* GetList() const = 0;
  
  /// Replicate (and optionally filter) the given aod event
  virtual void ReplicateAndFilter(const AliAODEvent& source) = 0;

  ClassDef(AliAODBranchReplicator,1) // AOD branch replicator base class
};

#endif

#ifndef ALIAODBRANCHREPLICATOR_H
#define ALIAODBRANCHREPLICATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

//
// Base class of an object used for the replication 
// (and possibly filtering) of one (or several) AOD branches.
//

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class AliAODEvent;

class AliAODBranchReplicator : public TNamed
{
public:
  AliAODBranchReplicator(const char* name="", const char* title="");
  
  virtual ~AliAODBranchReplicator();

  /// Return the list of object we manage
  virtual TList* GetList() const = 0;
  
  /// Replicate (and optionally filter) the given aod event
  virtual void ReplicateAndFilter(const AliAODEvent& source) = 0;

  ClassDef(AliAODBranchReplicator,1) // AOD branch replicator base class
};

#endif

#ifndef ALIAODBRANCHREPLICATOR_H
#define ALIAODBRANCHREPLICATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

//
// Base class of an object used for the replication 
// (and possibly filtering) of one (or several) AOD branches.
//

#ifndef ROOT_TNamed
#  include "TNamed.h"
#endif

class AliAODEvent;

class AliAODBranchReplicator : public TNamed
{
public:
  AliAODBranchReplicator(const char* name="", const char* title="");
  
  virtual ~AliAODBranchReplicator();

  /// Return the list of object we manage
  virtual TList* GetList() const = 0;
  
  /// Replicate (and optionally filter) the given aod event
  virtual void ReplicateAndFilter(const AliAODEvent& source) = 0;

  ClassDef(AliAODBranchReplicator,1) // AOD branch replicator base class
};

#endif

