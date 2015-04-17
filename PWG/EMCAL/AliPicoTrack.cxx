//
// Track class with minimal number of information 
// (targets at selection of primary tracks).
//
// Author: C.Loizides

#include "AliPicoTrack.h"
#include "AliExternalTrackParam.h"
#include "AliVCluster.h"
#include "AliAODTrack.h"

//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack() :
  AliVTrack(),
  fPt(0), fEta(0), fPhi(0), fM(0.13957), fQ(0), fLabel(-1), fTrackType(0), 
  fEtaEmc(0), fPhiEmc(0), fPtEmc(0), fEmcal(0), fFlag(0), fGeneratorIndex(-1), fClusId(-1), fOrig(0)
{
  // Default constructor.
}

//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack(Double_t pt, Double_t eta, Double_t phi, Byte_t q, Int_t lab, Byte_t type,
                           Double_t etaemc, Double_t phiemc, Double_t ptemc, Bool_t ise, Double_t mass) :
  AliVTrack(),
  fPt(pt), fEta(eta), fPhi(phi), fM(mass), fQ(q), fLabel(lab), fTrackType(type), 
  fEtaEmc(etaemc), fPhiEmc(phiemc), fPtEmc(ptemc), fEmcal(ise), fFlag(0), fGeneratorIndex(-1), fClusId(-1), fOrig(0)
{
  // Constructor.
}
  
//_________________________________________________________________________________________________
AliPicoTrack::AliPicoTrack(const AliPicoTrack &pc) :
  AliVTrack(pc),
  fPt(pc.fPt), fEta(pc.fEta), fPhi(pc.fPhi), fM(pc.fM),
  fQ(pc.fQ), fLabel(pc.fLabel), fTrackType(pc.fTrackType),  
  fEtaEmc(pc.fEtaEmc), fPhiEmc(pc.fPhiEmc), fPtEmc(pc.fPtEmc), fEmcal(pc.fEmcal), fFlag(pc.fFlag), fGeneratorIndex(pc.fGeneratorIndex),
  fClusId(pc.fClusId), fOrig(pc.fOrig)
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
    fM      = pc.fM;
    fQ      = pc.fQ;
    fLabel  = pc.fLabel;
    fTrackType = pc.fTrackType;
    fEtaEmc = pc.fEtaEmc;
    fPhiEmc = pc.fPhiEmc;
    fPtEmc  = pc.fPtEmc;
    fEmcal  = pc.fEmcal;
    fFlag  = pc.fFlag;
    fGeneratorIndex  = pc.fGeneratorIndex;
    fClusId = pc.fClusId;
    fOrig = pc.fOrig;
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
void AliPicoTrack::GetEtaPhiDiff(const AliVTrack *t, const AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
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

//_________________________________________________________________________________________________
Byte_t AliPicoTrack::GetTrackType(const AliVTrack *t)
{
  // Get track type encoded from bits 20 and 21.

  Byte_t ret = 0;
  if (t->TestBit(BIT(22)) && !t->TestBit(BIT(23)))
    ret = 1;
  else if (!t->TestBit(BIT(22)) && t->TestBit(BIT(23)))
    ret = 2;
  else if (t->TestBit(BIT(22)) && t->TestBit(BIT(23)))
    ret = 3;
  return ret;
}

//________________________________________________________________________
Byte_t AliPicoTrack::GetTrackType(const AliAODTrack *aodTrack, UInt_t filterBit1, UInt_t filterBit2)
{
  // Return track type: 0 = filterBit1, 1 = filterBit2 && ITS, 2 = filterBit2 && !ITS.
  // Returns 3 if filterBit1 and filterBit2 do not test.
  // WARNING: only works with AOD tracks and AOD filter bits must be provided. Otherwise will always return 0.

  Int_t res = 0;

  if (aodTrack->TestFilterBit(filterBit1)) {
    res = 0;
  }
  else if (aodTrack->TestFilterBit(filterBit2)) {
    if ((aodTrack->GetStatus()&AliVTrack::kITSrefit)!=0) {
      res = 1;
    }
    else {
      res = 2;
    }
  }
  else {
    res = 3;
  }

  return res;
}
