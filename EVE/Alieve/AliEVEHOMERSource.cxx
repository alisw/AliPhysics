// $Header$

#include "AliEVEHOMERSource.h"

//______________________________________________________________________
// AliEVEHOMERSource
//

ClassImp(AliEVEHOMERSource)

AliEVEHOMERSource::AliEVEHOMERSource(const Text_t* n, const Text_t* t) :
  TEveElement(),
  TNamed(n, t),
  fSource(0)
{}

AliEVEHOMERSource::AliEVEHOMERSource(AliHLTHOMERSourceDesc* src, const Text_t* n, const Text_t* t) :
  TEveElement(),
  TNamed(n, t),
  fSource(src)
{}

