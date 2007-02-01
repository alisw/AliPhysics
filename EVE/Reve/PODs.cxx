// $Header$

#include <TClass.h>

#include "PODs.h"

using namespace Reve;

void Reve::DisablePODTObjectStreamers()
{
  // Vector is not TObject

  // MCTrack derives from TParticle 
  TParticle::Class()->IgnoreTObjectStreamer(true);
  MCTrackRef::Class()->IgnoreTObjectStreamer(true);

  Hit::Class()->IgnoreTObjectStreamer(true);
  Cluster::Class()->IgnoreTObjectStreamer(true);

  RecTrack::Class()->IgnoreTObjectStreamer(true);
  // RecKink derives from RecTrack

  RecV0::Class()->IgnoreTObjectStreamer(true);

  GenInfo::Class()->IgnoreTObjectStreamer(true);
}

//______________________________________________________________________
// Point
//

ClassImp(Reve::Vector)

Float_t Vector::Eta() const
{
  Float_t cosTheta = CosTheta();
  if (cosTheta*cosTheta < 1) return -0.5* TMath::Log( (1.0-cosTheta)/(1.0+cosTheta) );
  Warning("Eta","transverse momentum = 0! return +/- 10e10");
  return (z >= 0) ? 10e10 : -10e10;
}

/**************************************************************************/
/**************************************************************************/

//ClassImp(Hit)
//ClassImp(RecTrack)
