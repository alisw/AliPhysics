#ifndef AliJetShape_H
#define AliJetShape_H

// $Id$

#include <vector>
#include <TString.h>
#include "fastjet/PseudoJet.hh"
#ifdef FASTJET_VERSION
#include "fastjet/FunctionOfPseudoJet.hh"
#endif

#include "TMath.h"

using namespace std;

#ifdef FASTJET_VERSION
class AliJetShapeMass : public fastjet::FunctionOfPseudoJet<Double32_t>
{
 public:
  virtual std::string description() const{return "jet mass";}
  Double32_t result(const fastjet::PseudoJet &jet) const{ return jet.m();}
};

class AliJetShapeGRNum : public fastjet::FunctionOfPseudoJet<Double32_t>
{
 public:
  // default ctor
  AliJetShapeGRNum(Double_t r = 0.2, Double_t wr = 0.04) : fR(r),fDRStep(wr){}
  virtual std::string description() const{return "Numerator angular structure function";}
  //  static Int_t GetBin(Double_t x) {Int_t bin = TMath::FloorNint((x-kxmin)/mdx); return bin;}
  Double32_t result(const fastjet::PseudoJet &jet) const {
    if (!jet.has_constituents())
      return 0; //AliFatal("Angular structure can only be applied on jets for which the constituents are known.");

  Double_t A = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    /* Int_t uid = constits[ic].user_index(); */
    /* if (uid == -1) //skip ghost particle */
    /*   continue; */
    for(UInt_t jc = ic+1; jc < constits.size(); ++jc) {
      /* Int_t uid = constits[jc].user_index(); */
      /* if (uid == -1) //skip ghost particle */
      /* 	continue; */
      Double_t dphi = constits[ic].phi()-constits[jc].phi();
      if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
      if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
      Double_t dr2 = (constits[ic].eta()-constits[jc].eta())*(constits[ic].eta()-constits[jc].eta()) + dphi*dphi;
      if(dr2>0.) {
	Double_t dr = TMath::Sqrt(dr2);
	Double_t x = fR-dr;
	//noisy function
	Double_t noise = TMath::Exp(-x*x/(2*fDRStep*fDRStep))/(TMath::Sqrt(2.*TMath::Pi())*fDRStep);
	A += constits[ic].perp()*constits[jc].perp()*dr2*noise;
      }
    }
  }
  return A;
  }

 protected:
  Double_t fR;
  Double_t fDRStep;
};

class AliJetShapeGRDen : public fastjet::FunctionOfPseudoJet<Double32_t>
{
 public:
  // default ctor
  AliJetShapeGRDen(Double_t r = 0.2, Double_t wr = 0.04) : fR(r),fDRStep(wr){}
  virtual std::string description() const{return "Denominator angular structure function";}
  Double32_t result(const fastjet::PseudoJet &jet) const {
  if (!jet.has_constituents())
    return 0; //AliFatal("Angular structure can only be applied on jets for which the constituents are known.");

  Double_t A = 0.;
  std::vector<fastjet::PseudoJet> constits = jet.constituents();
  for(UInt_t ic = 0; ic < constits.size(); ++ic) {
    /* Int_t uid = constits[ic].user_index(); */
    /* if (uid == -1) //skip ghost particle */
    /*   continue; */
    for(UInt_t jc = ic+1; jc < constits.size(); ++jc) {
      /* Int_t uid = constits[jc].user_index(); */
      /* if (uid == -1) //skip ghost particle */
      /* 	continue; */
      Double_t dphi = constits[ic].phi()-constits[jc].phi();
      if(dphi<-1.*TMath::Pi()) dphi+=TMath::TwoPi();
      if(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
      Double_t dr2 = (constits[ic].eta()-constits[jc].eta())*(constits[ic].eta()-constits[jc].eta()) + dphi*dphi;
      if(dr2>0.) {
	Double_t dr = TMath::Sqrt(dr2);
	Double_t x = fR-dr;
	//error function
	Double_t erf = 0.5*(1.+TMath::Erf(x/(TMath::Sqrt(2.)*fDRStep)));
	A += constits[ic].perp()*constits[jc].perp()*dr2*erf;
      }
    }
  }
  return A;
  }

 protected:
  Double_t fR;
  Double_t fDRStep;
};

#endif
#endif
