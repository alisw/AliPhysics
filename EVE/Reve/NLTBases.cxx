// $Header$

#include "NLTBases.h"
#include "Reve/NLTPolygonSet.h"

using namespace Reve;

//______________________________________________________________________________
// NLTProjectable
//

ClassImp(NLTProjectable)

//______________________________________________________________________________
NLTProjectable::NLTProjectable()
{
  // Default constructor.
}

//______________________________________________________________________________
NLTProjectable::~NLTProjectable()
{
  // Destructor.
  // Force projected replicas to unreference *this.

  while ( ! fProjectedList.empty())
  {
    fProjectedList.front()->UnRefProjectable(this);
  }
}


//______________________________________________________________________________
// NLTProjected
//

ClassImp(NLTProjected)

//______________________________________________________________________________
NLTProjected::NLTProjected() :
  fProjector   (0),
  fProjectable (0),
  fDepth       (0)
{
  // Default constructor.
}

//______________________________________________________________________________
NLTProjected::~NLTProjected()
{
  // Destructor.
  // If fProjectable is non-null, *this is removed from its list of
  // projected replicas.

  if (fProjectable) fProjectable->RemoveProjected(this);
}

//______________________________________________________________________________
void NLTProjected::SetProjection(NLTProjector* proj, NLTProjectable* model)
{
  fProjector   = proj;
  if (fProjectable) fProjectable->RemoveProjected(this);
  fProjectable = model;
  if (fProjectable) fProjectable->AddProjected(this);
}

//______________________________________________________________________________
void NLTProjected::UnRefProjectable(NLTProjectable* assumed_parent)
{
  static const Exc_t eH("NLTProjected::UnRefProjectable ");

  if (fProjectable != assumed_parent) {
    Warning(eH, "mismatch between assumed and real model. This is a bug.");
    assumed_parent->RemoveProjected(this);
    return;
  }

  if (fProjectable) {
    fProjectable->RemoveProjected(this);
    fProjectable = 0;
  }
}
