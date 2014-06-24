#ifndef AliJetShape_H
#define AliJetShape_H

// $Id$

#include <vector>
#include <TString.h>
#include "fastjet/PseudoJet.hh"
#ifdef FASTJET_VERSION
#include "fastjet/FunctionOfPseudoJet.hh"
#endif

#ifdef FASTJET_VERSION
class AliJetShapeMass : public fastjet::FunctionOfPseudoJet<Double32_t>
{
 public:
  virtual std::string description() const{return "jet mass";}
  Double32_t result(const fastjet::PseudoJet &jet) const{ return jet.m();}
};
#endif
#endif
