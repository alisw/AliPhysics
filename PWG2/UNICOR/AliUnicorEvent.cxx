/************************************************************************* 
* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. * 
*                                                                        * 
* Author: The ALICE Off-line Project.                                    * 
* Contributors are mentioned in the code where appropriate.              * 
*                                                                        * 
* Permission to use, copy, modify and distribute this software and its   * 
* documentation strictly for non-commercial purposes is hereby granted   * 
* without fee, provided that the above copyright notice appears in all   * 
* copies and that both the copyright notice and this permission notice   * 
* appear in the supporting documentation. The authors make no claims     * 
* about the suitability of this software for any purpose. It is          * 
* provided "as is" without express or implied warranty.                  * 
**************************************************************************/

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// parent class of all events; analyzers access data via this class
//=============================================================================

#include <TMath.h>
#include "AliUnicorEvent.h"

ClassImp(AliUnicorEvent)

//=============================================================================
void AliUnicorEvent::RP(Double_t &qx, Double_t &qy, Int_t harmonic) const 
{
  // simplest flow vector

  qx=0; 
  qy=0;
  for (int i=0; i<NParticles(); i++) {
    if (!ParticleGood(i,0)) continue;
    double pt = ParticlePt(i);
    if (pt>2.0) pt = 2.0; // from 2 GeV flow saturates anyway
    qx += pt*cos(harmonic*ParticlePhi(i));
    qy += pt*sin(harmonic*ParticlePhi(i));
  }
}
//=============================================================================
Double_t AliUnicorEvent::ParticleEta(Int_t i) const 
{
  // pseudorapidity

  double the = ParticleTheta(i); 
  if (the<0.0001) return 10; 
  else if (the>TMath::Pi()-0.0001) return -10;
  return -log(tan(the/2));
}
//=============================================================================
Double_t AliUnicorEvent::ParticleY(Int_t i, Double_t mass) const 
{
  // rapidity

  double pp = ParticleP(i);
  double ee = sqrt(fabs(mass*mass + pp*pp));
  double pz = ParticlePz(i);
  double yy = log((ee+pz)/(ee-pz))/2;
  return yy;
}
//=============================================================================

