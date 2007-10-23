// $Header$

#include "NLTBases.h"
#include "Reve/NLTPolygonSet.h"

using namespace Reve;

//______________________________________________________________________
// NLTProjectable
//

ClassImp(NLTProjectable)

NLTProjectable::NLTProjectable()
{}

NLTProjectable::~NLTProjectable()
{
  while ( ! fProjectedList.empty())
  {
    fProjectedList.front()->UnRefProjectable(this);
  }
}

//______________________________________________________________________
// NLTProjected
//

ClassImp(NLTProjected)

NLTProjected::NLTProjected() :
  fProjector   (0),
  fProjectable (0),
  fDepth       (0)
{}

NLTProjected::~NLTProjected()
{
  if (fProjectable) fProjectable->RemoveProjected(this);
}

void NLTProjected::SetProjection(NLTProjector* proj, NLTProjectable* model)
{
  fProjector   = proj;
  if (fProjectable) fProjectable->RemoveProjected(this);
  fProjectable = model;
  fProjectable->AddProjected(this);
}

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
