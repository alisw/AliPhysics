#include "AliPicoDQ.h"

ClassImp(AliPicoDQ)

//_____________________________________________________________________________
AliPicoDQ::AliPicoDQ() :
TObject(),
fPt(0.),
fEta(0.),
fPhi(0.),
fM(0.)
{
//
//  AliPicoDQ::AliPicoDQ()
//
}

//_____________________________________________________________________________
AliPicoDQ::AliPicoDQ(const TLorentzVector &v) :
TObject(),
fPt(v.Pt()),
fEta(v.Eta()),
fPhi(v.Phi()),
fM(v.M())
{
//
// AliPicoDQ::AliPicoDQ(const TLorentzVector &v)
//
}

//_____________________________________________________________________________
AliPicoDQ::AliPicoDQ(const AliPicoDQ &a) :
TObject(a),
fPt(a.fPt),
fEta(a.fEta),
fPhi(a.fPhi),
fM(a.fM)
{
//
//  AliPicoDQ::AliPicoDQ(const AliPicoDQ &a)
//
}

//_____________________________________________________________________________
AliPicoDQ &AliPicoDQ::operator=(const AliPicoDQ &a)
{
//
//  AliPicoDQ &AliPicoDQ::operator=(const AliPicoDQ &a)
//

  if (&a==this) return *this;
//=============================================================================

  TObject::operator=(a);
//=============================================================================

  fPt  = a.fPt;
  fEta = a.fEta;
  fPhi = a.fPhi;
  fM   = a.fM;
//=============================================================================

  return *this;
}

//_____________________________________________________________________________
AliPicoDQ::~AliPicoDQ()
{
//
//  AliPicoDQ::~AliPicoDQ()
//
}
