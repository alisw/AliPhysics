#include <TMath.h>
#include "AliITSUSeed.h"
#include "AliITSUAux.h"
using namespace AliITSUAux;
using namespace TMath;

ClassImp(AliITSUSeed)

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed() 
:  fClID(0)
  ,fParent(0)
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
  :AliExternalTrackParam(src)
  ,fClID(src.fClID)
  ,fParent(src.fParent) 
{
  // def c-tor
}

//_________________________________________________________________________
AliITSUSeed &AliITSUSeed::operator=(const AliITSUSeed& src) 
{
  // def c-tor
  if (this == &src) return *this;
  fClID = src.fClID;
  fParent = src.fParent;
  AliExternalTrackParam::operator=(src);
  return *this;
}
