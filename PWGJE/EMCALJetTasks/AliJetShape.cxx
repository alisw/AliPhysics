// $Id$
//
// Author: M. Verweij

#define AliJetShape_CXX
#include "AliJetShape.h"

/*
Double32_t AliJetShapeGRNum::result(const fastjet::PseudoJet &jet) const {

#ifdef FASTJET_VERSION
  if (!jet.has_constituents())
    throw Error("Angular structure can only be applied on jets for which the constituents are known.");

  Double_t A = 0.;
  vector<PseudoJet> constits = jet.constituents();
  for (vector<PseudoJet>::iterator ci = constits.begin(); ci!=constits.end(); ci++){
    for (vector<PseudoJet>::iterator cj = ci+1; cj!=constits.end(); cj++){
      Double_t dphi = ci.phi()-cj.phi();
      if(dphi<0.) dphi+=TMath::TwoPi();
      if(dphi>TMath::TwoPi()) dphi-=TMath::TwoPi();
      Double_t dr2 = (ci.eta()-cj.eta())*(ci.eta()-cj.eta()) + dphi*dphi;
      if(dr2>0.) {
	Double_t dr = TMath::Sqrt(dr2);
	Double_t x = fR-dr;
	//noisy function
	Double_t noise = TMath::Exp(-x*x/(2*fDRStep*fDRStep))/(TMath::Sqrt(2.*TMath::Pi())*fDRStep);
	//error function
	// Double_t erf = 0.5*(1.+TMath::Erf(x/(TMath::Sqrt(2.)*fDRStep)));
	A += ci.perp()*cj.perp()*dr2*noise;
      }
    }
  }
  return A;
#endif
  return 0;
}
*/

Double32_t AliJetShapeGRDen::result(const fastjet::PseudoJet &jet) const {

#ifdef FASTJET_VERSION
  if (!jet.has_constituents())
    throw Error("Angular structure can only be applied on jets for which the constituents are known.");

  Double_t A = 0.;
  vector<PseudoJet> constits = jet.constituents();
  for (vector<PseudoJet>::iterator ci = constits.begin(); ci!=constits.end(); ci++){
    for (vector<PseudoJet>::iterator cj = ci+1; cj!=constits.end(); cj++){
      Double_t dphi = ci.phi()-cj.phi();
      if(dphi<0.) dphi+=TMath::TwoPi();
      if(dphi>TMath::TwoPi()) dphi-=TMath::TwoPi();
      Double_t dr2 = (ci.eta()-cj.eta())*(ci.eta()-cj.eta()) + dphi*dphi;
      if(dr2>0.) {
	Double_t dr = TMath::Sqrt(dr2);
	Double_t x = fR-dr;
	//error function
	 Double_t erf = 0.5*(1.+TMath::Erf(x/(TMath::Sqrt(2.)*fDRStep)));
	A += ci.perp()*cj.perp()*dr2*erf;
      }
    }
  }
  return A;
#endif
  return 0;
}
