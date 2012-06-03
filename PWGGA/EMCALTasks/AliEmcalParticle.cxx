// $Id$
//
// Emcal particle class, which can contain either an AliVTrack or an AliVCluster
//
// Author: S.Aiola

#include "AliEmcalParticle.h"
#include "AliVCluster.h"
#include "AliLog.h"

//_________________________________________________________________________________________________
AliEmcalParticle::AliEmcalParticle() :
  TObject(),
  fTrack(0), 
  fCluster(0),
  fNMatched(0),
  fId(-1)
{
  // Default constructor.

  ResetMatchedObjects();
}

//_________________________________________________________________________________________________
AliEmcalParticle::AliEmcalParticle(TObject *particle, Int_t id) :
  TObject(),
  fTrack(0), 
  fCluster(0),
  fNMatched(0),
  fId(id)
{
  // Constructor.

  if (!particle)
    AliWarning("Null pointer passed as particle.");


  fTrack = dynamic_cast<AliVTrack*>(particle);
  if (!fTrack)
    fCluster = dynamic_cast<AliVCluster*>(particle);

  if (!fTrack && !fCluster)
    AliWarning("Particle type not recognized (not AliVTrack nor AliVCluster).");

  ResetMatchedObjects();
}
  
//_________________________________________________________________________________________________
AliEmcalParticle::AliEmcalParticle(const AliEmcalParticle &p) :
  TObject(),
  fTrack(p.fTrack),
  fCluster(p.fCluster), 
  fNMatched(p.fNMatched),
  fId(p.fId)
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
    fTrack = p.fTrack;
    fCluster = p.fCluster;
    fNMatched    = p.fNMatched;
    
    ResetMatchedObjects();
  
    memcpy(fMatchedIds, p.fMatchedIds, sizeof(UShort_t) * fSizeMatched);
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
    i++;
  
  if (i < fNMatched) {
    memmove(fMatchedIds  + i + 1, fMatchedIds  + i, sizeof(UShort_t) * (fNMatched - i));
    memmove(fMatchedDist + i + 1, fMatchedDist + i, sizeof(Double_t) * (fNMatched - i));
  }
  
  fMatchedIds[i]  = id;
  fMatchedDist[i] = d;
  fNMatched++;

  if (fNMatched >= fSizeMatched)
    fNMatched = fSizeMatched - 1;
}

//_________________________________________________________________________________________________
TLorentzVector &AliEmcalParticle::GetLorentzVector(const Double_t *vertex) const
{
  // Make a TLorentzVector and return it.

  static TLorentzVector vect;

  if (fTrack) {
    vect.SetPtEtaPhiM(fTrack->Pt(), fTrack->Eta(), fTrack->Phi(), M());
  }
  else if (fCluster && vertex) {
    fCluster->GetMomentum(vect, const_cast<Double_t*>(vertex));
  }

  return vect;
}
