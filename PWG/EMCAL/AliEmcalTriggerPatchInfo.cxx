// $Id$
//
// Emcal particle trigger class, which can contain either
//
// Author: J.Kral

#include "AliEmcalTriggerPatchInfo.h"
#include "AliLog.h"

//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo::AliEmcalTriggerPatchInfo() :
  TObject(),
  fCenterGeo(),
  fCenterMass(),
  fEdge1(),
  fEdge2(),
  fADCAmp(-1),
  fTriggerBits(0)
{
  // Default constructor.
  
}

  
//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo::AliEmcalTriggerPatchInfo(const AliEmcalTriggerPatchInfo &p) :
  TObject(p),
  fCenterGeo(p.fCenterGeo),
  fCenterMass(p.fCenterMass),
  fEdge1(p.fEdge1),
  fEdge2(p.fEdge2),
  fADCAmp(p.fADCAmp),
  fTriggerBits(p.fTriggerBits)
{
  // Copy constructor.
}

//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo::~AliEmcalTriggerPatchInfo()
{
  // Destructor.
}

//_________________________________________________________________________________________________
AliEmcalTriggerPatchInfo &AliEmcalTriggerPatchInfo::operator=(const AliEmcalTriggerPatchInfo &p)
{
  // Assignment operator.

  if (this != &p) {
    fCenterGeo = p.fCenterGeo;
    fCenterMass = p.fCenterMass;
    fEdge1 = p.fEdge1;
    fEdge2 = p.fEdge2;
    fADCAmp = p.fADCAmp;
    fTriggerBits = p.fTriggerBits;
  }

  return *this;
}

//_________________________________________________________________________________________________
void AliEmcalTriggerPatchInfo::SetLorentzVector( TLorentzVector &lv, TVector3 &v, Double_t e ){
  // sets the vector
  Double_t r = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2] ) ; 
  
  lv.SetPxPyPzE( e*v[0]/r,  e*v[1]/r,  e*v[2]/r,  e) ;   
}

