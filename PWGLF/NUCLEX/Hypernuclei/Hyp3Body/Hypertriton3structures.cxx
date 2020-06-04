#include "Hypertriton3structures.h"

ClassImp(RHyperTriton)
ClassImp(RHyperTriton3O2)
ClassImp(RHyperTriton3KF)
ClassImp(SHyperTriton<RHyperTriton3KF>)
ClassImp(SHyperTriton<RHyperTriton3O2>)

SHyperTriton3O2 __dummy_instanceO2__()
{
  SHyperTriton3O2 a;
  a.gPt = -12;
  return a;
}

SHyperTriton3KF __dummy_instanceKF__()
{
  SHyperTriton3KF a;
  a.gPt = -12;
  return a;
}
