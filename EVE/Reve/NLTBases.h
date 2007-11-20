// $Header$

#ifndef REVE_NLTBases_H
#define REVE_NLTBases_H

#include <Reve/Reve.h>

#include <list>

class TBuffer3D;

namespace Reve {

class NLTProjected;
class NLTProjector;

////////////////////////////////////////////////////////////////
//                                                            //
// NLTProjectable                                             //
//                                                            //
// Abstract base class for non-linear projectable objects.    //
//                                                            //
////////////////////////////////////////////////////////////////

class NLTProjectable
{
private:
  NLTProjectable(const NLTProjectable&);            // Not implemented
  NLTProjectable& operator=(const NLTProjectable&); // Not implemented

protected:
  std::list<NLTProjected*> fProjectedList; // references to projected instances.

public:
  NLTProjectable();
  virtual ~NLTProjectable();

  virtual TClass* ProjectedClass() const = 0;

  virtual void AddProjected(NLTProjected* p)    { fProjectedList.push_back(p); }
  virtual void RemoveProjected(NLTProjected* p) { fProjectedList.remove(p); }

  ClassDef(NLTProjectable, 0); // Abstract base class for non-linear projectable objects.  
}; // endclass NLTProjectable

////////////////////////////////////////////////////////////////
//                                                            //
// NLTProjected                                               //
//                                                            //
// Abstract base class for non-linear projected objects.      //
//                                                            //
////////////////////////////////////////////////////////////////

class NLTProjected
{
private:
  NLTProjected(const NLTProjected&);            // Not implemented
  NLTProjected& operator=(const NLTProjected&); // Not implemented

protected:
  NLTProjector   *fProjector;     // manager
  NLTProjectable *fProjectable;   // link to original object
  Float_t         fDepth;         // z coordinate

public:
  NLTProjected();
  virtual ~NLTProjected();

  virtual void SetProjection(NLTProjector* proj, NLTProjectable* model);
  virtual void UnRefProjectable(NLTProjectable* assumed_parent);

  virtual void SetDepth(Float_t d) { fDepth = d; }

  virtual void UpdateProjection() = 0;

  ClassDef(NLTProjected, 0); // Abstract base class for non-linear projected object. 
}; // endclass NLTProjected

}

#endif
