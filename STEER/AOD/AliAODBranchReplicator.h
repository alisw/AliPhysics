#ifndef ALIAODBRANCHREPLICATOR_H
#define ALIAODBRANCHREPLICATOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

// Base class of an object used for the replication 
// (and possibly filtering) of one (or several) AOD branches.
//
// Author L. Aphecetche (Subatech)

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

