// $Id$
//
// Track class with minimal number of information 
// (targets at selection of primary tracks).
//
// Author: C.Loizides

#include "AliPicoTrack.h"
#include "AliExternalTrackParam.h"
#include "AliVCluster.h"

//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack() :
  AliVTrack(),
  fPt(0), fEta(0), fPhi(0), fQ(0), fLabel(-1), fTrackType(0), fEtaEmc(0), fPhiEmc(0), fEmcal(0), fClusId(-1)
{
  // Default constructor.
}

//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack(Double_t pt, Double_t eta, Double_t phi, Byte_t q, Int_t lab, Byte_t type,
                           Double_t etaemc, Double_t phiemc, Bool_t ise) :
  AliVTrack(),
  fPt(pt), fEta(eta), fPhi(phi), fQ(q), fLabel(lab), fTrackType(type), 
  fEtaEmc(etaemc), fPhiEmc(phiemc), fEmcal(ise), fClusId(-1)
{
  // Constructor.
}
  
//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack(const AliPicoTrack &pc) :
  AliVTrack(pc),
  fPt(pc.fPt), fEta(pc.fEta), fPhi(pc.fPhi), 
  fQ(pc.fQ), fLabel(pc.fLabel), fTrackType(pc.fTrackType),  
  fEtaEmc(pc.fEtaEmc), fPhiEmc(pc.fPhiEmc), fEmcal(pc.fEmcal),
  fClusId(pc.fClusId)
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
    fTrackType = pc.fTrackType;
    fEtaEmc = pc.fEtaEmc;
    fPhiEmc = pc.fPhiEmc;
    fEmcal  = pc.fEmcal;
    fClusId = pc.fClusId;
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

//_________________________________________________________________________________________________
void AliPicoTrack::GetEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
  // Calculate phi and eta difference between track and cluster.
 
  phidiff = 999;
  etadiff = 999;

  if (!t||!v)
    return;

  if (!t->IsEMCAL())
    return;

  Double_t veta = t->GetTrackEtaOnEMCal();
  Double_t vphi = t->GetTrackPhiOnEMCal();

  Float_t pos[3] = {0};
  v->GetPosition(pos);  
  TVector3 cpos(pos); 
  Double_t ceta     = cpos.Eta();
  Double_t cphi     = cpos.Phi();
  etadiff=veta-ceta;
  phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}
