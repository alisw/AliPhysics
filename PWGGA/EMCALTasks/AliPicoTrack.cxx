// $Id$
//
// Track class with minimal number of information 
// (targets at selection of primary tracks)
//

#include "AliPicoTrack.h"

//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack() :
  AliVTrack(),
  fPt(0), fEta(0), fPhi(0), fQ(0), fLabel(-1), fEtaEmc(0), fPhiEmc(0), fEmcal(0)
{
  // Default constructor.
}

//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack(Double_t pt, Double_t eta, Double_t phi, Byte_t q, Byte_t lab, 
                           Double_t etaemc, Double_t phiemc, Bool_t ise) :
  AliVTrack(),
  fPt(pt), fEta(eta), fPhi(phi), fQ(q), fLabel(lab), 
  fEtaEmc(etaemc), fPhiEmc(phiemc), fEmcal(ise)
{
  // Constructor.
}
  
//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack(const AliPicoTrack &pc) :
  AliVTrack(pc),
  fPt(pc.fPt), fEta(pc.fEta), fPhi(pc.fPhi), 
  fQ(pc.fQ), fLabel(pc.fLabel), 
  fEtaEmc(pc.fEtaEmc), fPhiEmc(pc.fPhiEmc), fEmcal(pc.fEmcal)
{
  // Constructor.
}

//_________________________________________________________________________________________________
AliPicoTrack &AliPicoTrack::operator=(const AliPicoTrack &pc)
{
  // Assignment operator.

  if (this!=&pc) {
    AliVTrack::operator=(pc);
    fPt     = pc.fPt;
    fEta    = pc.fEta;
    fPhi    = pc.fPhi;
    fQ      = pc.fQ;
    fLabel  = pc.fLabel;
    fEtaEmc = pc.fEtaEmc;
    fPhiEmc = pc.fPhiEmc;
    fEmcal  = pc.fEmcal;
  }

  return *this;
}

//_________________________________________________________________________________________________
Int_t AliPicoTrack::Compare(const TObject* obj) const
{
  // Compare this class with an other instance of this class used in a 
  // TCollection::Sort()/TClonesArray::Sort() which is descending.
  // Returns 0 when equal, 1 when this is smaller and -1 when bigger.

  const AliPicoTrack *t = dynamic_cast<const AliPicoTrack*>(obj);
  if (!t) 
    return -1;
  if (t->Pt()>Pt())
    return 1;
  if (t->Pt()<Pt())
    return -1;
  return 0;
}


