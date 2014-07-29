#include <TLorentzVector.h>

#include "AliEmcalJet.h"

#include "AliPicoJet.h"

ClassImp(AliPicoJet)

//_____________________________________________________________________________
AliPicoJet::AliPicoJet() :
TObject(),
fKine(),
fArea(0.),
fLeadingPt()
{
//
// AliPicoJet::AliPicoJet
//
}

//_____________________________________________________________________________
AliPicoJet::AliPicoJet(AliEmcalJet const *pJet, Double_t dLeadingPt) :
TObject(),
fKine(pJet->Px(), pJet->Py(), pJet->Pz(), pJet->E()),
fArea(pJet->Area()),
fLeadingPt(dLeadingPt)
{
//
// AliPicoJet::AliPicoJet
//
}

//_____________________________________________________________________________
AliPicoJet::AliPicoJet(const AliPicoJet &src) :
TObject(src),
fKine(src.fKine),
fArea(src.fArea),
fLeadingPt(src.fLeadingPt)
{
//
// AliPicoJet::AliPicoJet
//
}

//_____________________________________________________________________________
AliPicoJet& AliPicoJet::operator=(const AliPicoJet &src)
{
//
// AliPicoJet::operator=
//

  if (&src==this) return *this;

  TObject::operator=(src);

  fKine = src.fKine;
  fArea = src.fArea;

  fLeadingPt = src.fLeadingPt;

  return *this;
}

//_____________________________________________________________________________
AliPicoJet::~AliPicoJet()
{
//
// AliPicoJet::~AliPicoJet
//
}
