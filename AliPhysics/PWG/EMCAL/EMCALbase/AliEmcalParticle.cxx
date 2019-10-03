//
// Emcal particle class, which can contain either an AliVTrack or an AliVCluster
//
// Author: S.Aiola

#include "AliEmcalParticle.h"
#include "AliVCluster.h"
#include "AliLog.h"

//_________________________________________________________________________________________________
AliEmcalParticle::AliEmcalParticle() :
  AliVParticle(),
  fTrack(0), 
  fCluster(0),
  fNMatched(0),
  fId(-1),
  fPhi(0),
  fEta(0),
  fPt(0),
  fM(0),
  fMatchedPtr(0)
{
  // Default constructor.

  ResetMatchedObjects();
}

//_________________________________________________________________________________________________
AliEmcalParticle::AliEmcalParticle(AliVTrack* track, Int_t id) :
  AliVParticle(),
  fTrack(track), 
  fCluster(0),
  fNMatched(0),
  fId(id),
  fPhi(0),
  fEta(0),
  fPt(0),
  fM(0),
  fMatchedPtr(0)
{
  // Constructor.

  if (!track) {
    AliWarning("Null pointer passed as particle.");
    return;
  }

  fEta = fTrack->Eta();
  fPhi = fTrack->Phi();
  fPt  = fTrack->Pt();
  fM   = fTrack->M();

  ResetMatchedObjects();
}

//_________________________________________________________________________________________________
AliEmcalParticle::AliEmcalParticle(AliVCluster* cluster, Int_t id, Double_t vx, Double_t vy, Double_t vz, Int_t ude) :
  AliVParticle(),
  fTrack(0), 
  fCluster(cluster),
  fNMatched(0),
  fId(id),
  fPhi(0),
  fEta(0),
  fPt(0),
  fM(0),
  fMatchedPtr(0)
{
  // Constructor.

  if (!cluster) {
    AliWarning("Null pointer passed as particle.");
    return;
  }
  
  Double_t vtx[3]; vtx[0]=vx;vtx[1]=vy;vtx[2]=vz;
  TLorentzVector vect;
  if (ude >= 0 && ude <= AliVCluster::kLastUserDefEnergy) {
    fCluster->GetMomentum(vect, vtx, (AliVCluster::VCluUserDefEnergy_t)ude);
  }
  else {
    fCluster->GetMomentum(vect, vtx);
  }
  fEta = vect.Eta();
  fPhi = vect.Phi();
  fPt  = vect.Pt();
  fM   = 0.;

  ResetMatchedObjects();
}
  
//_________________________________________________________________________________________________
AliEmcalParticle::AliEmcalParticle(const AliEmcalParticle &p) :
  AliVParticle(p),
  fTrack(p.fTrack),
  fCluster(p.fCluster), 
  fNMatched(p.fNMatched),
  fId(p.fId),
  fPhi(p.fPhi),
  fEta(p.fEta),
  fPt(p.fPt),
  fM(p.fM),
  fMatchedPtr(p.fMatchedPtr)
{
  // Copy constructor.

  ResetMatchedObjects();

  memcpy(fMatchedIds, p.fMatchedIds, sizeof(UShort_t) * fSizeMatched);
  memcpy(fMatchedDist, p.fMatchedDist, sizeof(Double_t) * fSizeMatched);
}

//_________________________________________________________________________________________________
AliEmcalParticle::~AliEmcalParticle()
{
  // Destructor.
}

//_________________________________________________________________________________________________
AliEmcalParticle &AliEmcalParticle::operator=(const AliEmcalParticle &p)
{
  // Assignment operator.

  if (this != &p) {
    fTrack      = p.fTrack;
    fCluster    = p.fCluster;
    fNMatched   = p.fNMatched;
    fId         = p.fId;
    fPhi        = p.fPhi;
    fEta        = p.fEta;
    fPt         = p.fPt;
    fM          = p.fM;
    fMatchedPtr = p.fMatchedPtr;

    ResetMatchedObjects();
    memcpy(fMatchedIds,  p.fMatchedIds,  sizeof(UShort_t) * fSizeMatched);
    memcpy(fMatchedDist, p.fMatchedDist, sizeof(Double_t) * fSizeMatched);
  }

  return *this;
}

//_________________________________________________________________________________________________
void AliEmcalParticle::ResetMatchedObjects()
{
  // Reset matched objects.

  for (Int_t i = 0; i < fSizeMatched; i++) {
    fMatchedIds[i] = -1;
    fMatchedDist[i] = 999;
  }
}

//_________________________________________________________________________________________________
void AliEmcalParticle::AddMatchedObj(Int_t id, Double_t d)
{
  // Add a matched object.

  Int_t i = 0;
  while (i < fNMatched && d > fMatchedDist[i])
    ++i;
  
  if (i < fNMatched) {
    memmove(fMatchedIds  + i + 1, fMatchedIds  + i, sizeof(UShort_t) * (fNMatched - i));
    memmove(fMatchedDist + i + 1, fMatchedDist + i, sizeof(Double_t) * (fNMatched - i));
  }
  
  fMatchedIds[i]  = id;
  fMatchedDist[i] = d;
  ++fNMatched;

  if (fNMatched >= fSizeMatched)
    fNMatched = fSizeMatched - 1;
}

//_________________________________________________________________________________________________
TLorentzVector &AliEmcalParticle::GetLorentzVector(const Double_t *vertex) const
{
  // Make a TLorentzVector and return it.

  static TLorentzVector vect;

  if (fTrack) {
    vect.SetPtEtaPhiM(fPt, fEta, fPhi, fM);
  }
  else if (fCluster && vertex) {
    fCluster->GetMomentum(vect, const_cast<Double_t*>(vertex));
  }

  return vect;
}
