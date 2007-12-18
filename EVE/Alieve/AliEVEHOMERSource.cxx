// $Header$

#include "AliEVEHOMERSource.h"

using namespace Reve;
//using namespace Alieve;

//______________________________________________________________________
// AliEVEHOMERSource
//

ClassImp(AliEVEHOMERSource)

AliEVEHOMERSource::AliEVEHOMERSource(const Text_t* n, const Text_t* t) :
  Reve::RenderElement(),
  TNamed(n, t),
  fSource(0)
{}

AliEVEHOMERSource::AliEVEHOMERSource(AliHLTHOMERSourceDesc* src, const Text_t* n, const Text_t* t) :
  Reve::RenderElement(),
  TNamed(n, t),
  fSource(src)
{}

