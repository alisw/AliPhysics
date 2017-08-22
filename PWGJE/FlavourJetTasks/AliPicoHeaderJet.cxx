#include <TMap.h>
#include <TString.h>
#include <TParameter.h>

#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"

#include "AliPicoHeaderJet.h"

ClassImp(AliPicoHeaderJet)

//_____________________________________________________________________________
AliPicoHeaderJet::AliPicoHeaderJet(const TString s) :
AliPicoHeaderV0(s),
fBkgRho(-999.)
{
//
// AliPicoHeaderJet::AliPicoHeaderJet
//
}

//_____________________________________________________________________________
AliPicoHeaderJet::AliPicoHeaderJet(const AliPicoHeaderJet &src) :
AliPicoHeaderV0(src),
fBkgRho(src.fBkgRho)
{
//
// AliPicoHeaderJet::AliPicoHeaderJet
//
}

//_____________________________________________________________________________
AliPicoHeaderJet& AliPicoHeaderJet::operator=(const AliPicoHeaderJet &src)
{
//
// AliPicoHeaderJet::operator=
//

  if (&src==this) return *this;
  AliPicoHeaderV0::operator=(src);
//=============================================================================

  fBkgRho = src.fBkgRho;
//=============================================================================

  return *this;
}

//_____________________________________________________________________________
AliPicoHeaderJet::~AliPicoHeaderJet()
{
//
// AliPicoHeaderJet::~AliPicoHeaderJet
//
}

//_____________________________________________________________________________
void AliPicoHeaderJet::Reset()
{
//
//  AliPicoHeaderJet::Reset
//

  AliPicoHeaderV0::Reset();

  fBkgRho = -999.;

  return;
}
