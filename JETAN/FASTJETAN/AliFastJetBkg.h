#ifndef ALIFASTJETBKG_H
#define ALIFASTJETBKG_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//---------------------------------------------------------------------
// Class to calculate the background per unit area
// manages the search for jets 
// Authors: Elena Bruna elena.bruna@yale.edu
//          Sevil Salur ssalur@lbl.gov
//  
// 2011 :
// renamed from AliJetBkg to AliFastJetBkg as this class uses only FASTJET based algos        
//---------------------------------------------------------------------
#include "fastjet/PseudoJet.hh"
#include <vector>
using std::vector;
#ifdef __CINT__
namespace fastjet {
  class PsuedoJet;
}
#endif
class TString;
class TClonesArray;
class AliFastJetInput;

class AliFastJetBkg : public TObject
{
 public:
  AliFastJetBkg();
  AliFastJetBkg(const AliFastJetBkg &input);
  AliFastJetBkg& operator=(const AliFastJetBkg& source);
  virtual          ~AliFastJetBkg() {;}
  void             SetHeader(AliJetHeader *header)  {fHeader=header;}
  void             SetFastJetInput(AliFastJetInput *fjinput)  {fInputFJ=fjinput;}
  void             BkgFastJetb(Double_t& x,Double_t& y, Double_t& z);
  void             BkgFastJetWoHardest(Double_t& x,Double_t& y, Double_t& z);
  Float_t          BkgFastJet();
  Float_t          BkgChargedFastJet();
  Float_t          BkgStat();
  Float_t          BkgFastJetCone(TClonesArray* fAODJets);

  Bool_t           EmcalAcceptance(Float_t eta, Float_t phi, Float_t radius) const;
  static Double_t  BkgFunction(Double_t *x,Double_t *par);
    
 private:
  Double_t         CalcRho(vector<fastjet::PseudoJet> input_particles,Double_t RparamBkg,TString method);
  void             CalcRhob(Double_t& median, Double_t& sigma, Double_t& meanarea,
			    vector<fastjet::PseudoJet> input_particles,Double_t RparamBkg,TString method);
  void             CalcRhoWoHardest(Double_t& median, Double_t& sigma, Double_t& meanarea,
				    vector<fastjet::PseudoJet> input_particles,Double_t RparamBkg,TString method);

  AliJetHeader*    fHeader;  //! header
  AliFastJetInput* fInputFJ; //! input particles

  ClassDef(AliFastJetBkg, 2)   //  Fastjet backgroud analysis
 
};
 
#endif
