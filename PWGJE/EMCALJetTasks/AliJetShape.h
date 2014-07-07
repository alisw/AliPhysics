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

/*
const Int_t knbins = 130;
const Double_t kxmin = -0.8;
const Double_t kxmax = 0.5;
static Int_t GetBin(Double_t x) {Double_t mdx = (kxmax-kxmin)/(double)knbins; Int_t bin = TMath::FloorNint((x-kxmin)/mdx); return bin;}
static Double_t gfcn[knbins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000005, 0.000020, 0.000069, 0.000226, 0.000696, 0.002014, 0.005473, 0.013977, 0.033528, 0.075556, 0.159953, 0.318105, 0.594298, 1.043025, 1.719657, 2.663457, 3.875307, 5.296916, 6.801375, 8.204024, 9.296377, 9.895942, 9.895942, 9.296377, 8.204024, 6.801375, 5.296916, 3.875307, 2.663457, 1.719657, 1.043025, 0.594298, 0.318105, 0.159953, 0.075556, 0.033528, 0.013977, 0.005473, 0.002014, 0.000696, 0.000226, 0.000069, 0.000020, 0.000005, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000};

static Double_t efcn[knbins] = {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001, 0.000002, 0.000006, 0.000019, 0.000053, 0.000144, 0.000369, 0.000889, 0.002020, 0.004332, 0.008774, 0.016793, 0.030396, 0.052081, 0.084566, 0.130295, 0.190787, 0.265986, 0.353830, 0.450262, 0.549738, 0.646170, 0.734014, 0.809213, 0.869705, 0.915434, 0.947919, 0.969604, 0.983207, 0.991226, 0.995668, 0.997980, 0.999111, 0.999631, 0.999856, 0.999947, 0.999981, 0.999994, 0.999998, 0.999999, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000};
*/

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
    if (!jet.has_constituents()) {
      Printf("Angular structure can only be applied on jets for which the constituents are known.");
      return 0; //AliFatal("Angular structure can only be applied on jets for which the constituents are known.");
    }
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
	//	Printf("Calculating numerator between %d-%d for r: %f dr: %f",ic,jc,fR,dr);
	//noisy function
	Double_t noise = TMath::Exp(-x*x/(2*fDRStep*fDRStep))/(TMath::Sqrt(2.*TMath::Pi())*fDRStep);
	/* Int_t bin = GetBin(x); */
	/* Double_t noise = 0; */
	/* if(bin<0) noise = gfcn[0]; */
	/* else if(bin>=knbins) noise = gfcn[knbins-1]; */
	/* else  noise = gfcn[bin]; */
	/* Printf("fR: %f dr: %f x=%f bin=%d noise: %f  %f",fR,dr,x,bin,noise,TMath::Exp(-x*x/(2*fDRStep*fDRStep))/(TMath::Sqrt(2.*TMath::Pi())*fDRStep)); */
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
	/* Int_t bin = GetBin(x); */
	/* Double_t erf = 0; */
	/* if(bin<0) erf = efcn[0]; */
	/* else if(bin>=knbins) erf = efcn[knbins-1]; */
	/* else  erf = efcn[bin]; */
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
