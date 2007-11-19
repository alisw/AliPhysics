// $Header$

#include <TClass.h>

#include "PODs.h"

using namespace Reve;

void Reve::DisablePODTObjectStreamers()
{
  // Vector is not TObject

  // MCTrack derives from TParticle 
  TParticle::Class()->IgnoreTObjectStreamer(true);

  Hit::Class()->IgnoreTObjectStreamer(true);
  Cluster::Class()->IgnoreTObjectStreamer(true);

  RecTrack::Class()->IgnoreTObjectStreamer(true);
  // RecKink derives from RecTrack

  RecV0::Class()->IgnoreTObjectStreamer(true);

  GenInfo::Class()->IgnoreTObjectStreamer(true);
}

//______________________________________________________________________
// Vector
//

ClassImp(Reve::Vector)

Float_t Vector::Eta() const
{
  Float_t cosTheta = CosTheta();
  if (cosTheta*cosTheta < 1) return -0.5* TMath::Log( (1.0-cosTheta)/(1.0+cosTheta) );
  Warning("Eta","transverse momentum = 0! return +/- 10e10");
  return (z >= 0) ? 10e10 : -10e10;
}

Vector Vector::operator + (const Vector & b)
{
   return Vector(x + b.x, y + b.y, z + b.z);
}

Vector Vector::operator - (const Vector & b)
{
   return Vector(x - b.x, y - b.y, z - b.z);
}

Vector Vector::operator * (Float_t a)
{
   return Vector(a*x, a*y, a*z);
}
/**************************************************************************/
/**************************************************************************/

//______________________________________________________________________
// PathMark
//

ClassImp(Reve::PathMark)

const char* PathMark::type_name()
{
  switch (type)
  {
    case Daughter:  return "Daughter";
    case Reference: return "Reference";
    case Decay:     return "Decay";
    default:        return "Unknown";
  }
}

//ClassImp(Hit)
//ClassImp(RecTrack)
