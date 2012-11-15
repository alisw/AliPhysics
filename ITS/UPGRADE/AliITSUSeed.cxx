#include <TMath.h>
#include "AliITSUSeed.h"
#include "AliITSUAux.h"
using namespace AliITSUAux;
using namespace TMath;

ClassImp(AliITSUSeed)

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed() 
: fMass(kPionMass)
{
  // def c-tor
}

//_________________________________________________________________________
AliITSUSeed::~AliITSUSeed()
{
  // d-rot
}

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed(const AliITSUSeed& src) 
  : AliExternalTrackParam(src)
  , fMass(src.fMass)
{
  // def c-tor
}

//_________________________________________________________________________
AliITSUSeed &AliITSUSeed::operator=(const AliITSUSeed& src) 
{
  // def c-tor
  if (this == &src) return *this;
  fMass = src.fMass;
  AliExternalTrackParam::operator=(src);
  return *this;
}
