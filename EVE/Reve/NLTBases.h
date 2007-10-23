// $Header$

#ifndef REVE_NLTBases_H
#define REVE_NLTBases_H

#include <Reve/Reve.h>

#include <list>

class TBuffer3D;

namespace Reve {

class NLTProjected;
class NLTProjector;

class NLTProjectable
{
private:
  NLTProjectable(const NLTProjectable&);            // Not implemented
  NLTProjectable& operator=(const NLTProjectable&); // Not implemented

protected:
  // Eventually, references to all projected instances.
  std::list<NLTProjected*> fProjectedList;

public:
  NLTProjectable();
  virtual ~NLTProjectable();

  virtual TClass* ProjectedClass() const = 0;

  virtual void AddProjected(NLTProjected* p)    { fProjectedList.push_back(p); }
  virtual void RemoveProjected(NLTProjected* p) { fProjectedList.remove(p); }

  ClassDef(NLTProjectable, 0);
}; // endclass NLTProjectable

/**************************************************************************/

class NLTProjected
{
private:
  NLTProjected(const NLTProjected&);            // Not implemented
  NLTProjected& operator=(const NLTProjected&); // Not implemented

protected:
  NLTProjector   *fProjector;
  NLTProjectable *fProjectable;

  Float_t         fDepth;

public:
  NLTProjected();
  virtual ~NLTProjected();

  virtual void SetProjection(NLTProjector* proj, NLTProjectable* model);
  virtual void UnRefProjectable(NLTProjectable* assumed_parent);

  virtual void SetDepth(Float_t d) { fDepth = d; }

  virtual void UpdateProjection() = 0;

  ClassDef(NLTProjected, 0);
}; // endclass NLTProjected

}

#endif
