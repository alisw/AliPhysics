#ifndef AliJetShape_H
#define AliJetShape_H

// $Id$

#include <vector>
#include <TString.h>
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#ifdef FASTJET_VERSION
#include "fastjet/FunctionOfPseudoJet.hh"
#endif

#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVector3.h"
#include "TVector2.h"
using namespace std;

#ifdef FASTJET_VERSION
//________________________________________________________________________
class AliJetShapeMass : public fastjet::FunctionOfPseudoJet<Double32_t>
{
 public:
  virtual std::string description() const{return "jet mass";}
  Double32_t result(const fastjet::PseudoJet &jet) const{ return jet.m();}
};

//________________________________________________________________________
class AliJetShapeGRNum : public fastjet::FunctionOfPseudoJet<Double32_t>
{
 public:
  // default ctor
  AliJetShapeGRNum(Double_t r = 0.2, Double_t wr = 0.04) : fR(r),fDRStep(wr){}
  virtual std::string description() const{return "Numerator angular structure function";}
  //  static Int_t GetBin(Double_t x) {Int_t bin = TMath::FloorNint((x-kxmin)/mdx); return bin;}
  Double32_t result(const fastjet::PseudoJet &jet) const;

 protected:
  Double_t fR;
  Double_t fDRStep;
};

//________________________________________________________________________
class AliJetShapeGRDen : public fastjet::FunctionOfPseudoJet<Double32_t>
{
 public:
  // default ctor
  AliJetShapeGRDen(Double_t r = 0.2, Double_t wr = 0.04) : fR(r),fDRStep(wr){}
  virtual std::string description() const{return "Denominator angular structure function";}
  Double32_t result(const fastjet::PseudoJet &jet) const;

 protected:
  Double_t fR;
  Double_t fDRStep;
};

//________________________________________________________________________
class AliJetShapeAngularity : public fastjet::FunctionOfPseudoJet<Double32_t>{
public:
  virtual std::string description() const{return "Angularity:radial moment";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

//________________________________________________________________________
class AliJetShapepTD : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "pTD";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShapeConstituent : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "constituents";}
  Double_t result(const fastjet::PseudoJet &jet) const {
    if (!jet.has_constituents())
      return 0; 
    Double_t num = 0.;
    std::vector<fastjet::PseudoJet> constits = jet.constituents();
    num=1.*constits.size();  
    return num;
  }
};

//________________________________________________________________________
class AliJetShapeCircularity : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "circularity denominator";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

//________________________________________________________________________
class AliJetShapeSigma2 : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "cms sigma2";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

//________________________________________________________________________
class AliJetShapeLeSub : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "leading mins subleading";}
  Double32_t result(const fastjet::PseudoJet &jet) const {
    if (!jet.has_constituents())
      return 0;
    std::vector<fastjet::PseudoJet> constits = jet.constituents();
    std::vector<fastjet::PseudoJet> sortedconstits=sorted_by_pt(constits); 
    if(sortedconstits.size()<2) return 0;
    Double_t num=TMath::Abs(sortedconstits[0].perp()-sortedconstits[1].perp());
    return num;
  }
};

//__________________________________________________________________________
class AliJetShape1subjettiness_kt : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "1subJettiness kt exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape2subjettiness_kt : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "2subJettiness kt exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape3subjettiness_kt : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "3subJettiness kt exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShapeOpeningAngle_kt : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "Opening Angle of Subjet Axes kt exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape1subjettiness_ca : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "1subJettiness ca exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape2subjettiness_ca : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "2subJettiness ca exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShapeOpeningAngle_ca : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "Opening Angle of Subjet Axes ca exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape1subjettiness_akt02 : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "1subJettiness akt02 exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape2subjettiness_akt02 : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "2subJettiness akt02 exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShapeOpeningAngle_akt02 : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "Opening Angle of Subjet Axes akt02 exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape1subjettiness_onepassca : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "1subJettiness ca sd exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShape2subjettiness_onepassca : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "2subJettiness ca sd exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};

class AliJetShapeOpeningAngle_onepassca : public fastjet::FunctionOfPseudoJet<Double32_t>{
 public:
  virtual std::string description() const{return "Opening Angle of Subjet Axes ca sd exclusive";}
  Double32_t result(const fastjet::PseudoJet &jet) const;
};






#endif
#endif

