#ifndef ALIGENPERFORMANCE_H
#define ALIGENPERFORMANCE_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// Perormance generator for Jet  according generic functions
//  
//  TF1 *   fF1Momentum;           // momentum distribution function inGeV
//  TF1 *   fFPhi;                // phi distribution function in rad
//  TF1 *   fFTheta;              // theta distribution function in rad
//  TF3 *   fFPosition;           // position distribution function in cm
//  TF1 *   fFPdg;                // pdg distribution function  
//  We assume that the moment, postion and PDG code of particles are independent  
//  Only tracks/particle crossing the reference radius at given z range
//
// Origin: marian.ivanov@cern.ch


#include "AliGenerator.h"
class TF3;
class TTreeSRedirector;
class TChain;

class AliGenPerformance : public AliGenerator{
public:

  AliGenPerformance();
  AliGenPerformance(const AliGenPerformance& perf);
  AliGenPerformance &operator=(const AliGenPerformance& perf);
  virtual ~AliGenPerformance() {}
  virtual void Generate();
  virtual void Init();
  void SetNJets(Double_t nJets){fNJets=nJets;}
  Double_t GetNJets()const {return fNJets;}
  void SetFunctions(TF1 * momentum, TF1 *fphi=0, TF1 *ftheta=0, TF3 * position=0, TF1* pdg=0);
  void SetStreamer(TTreeSRedirector *pcstream){fTestStream=pcstream;}
  TTreeSRedirector * GetStreamer(){return fTestStream;}
  static void TestAliGenPerformance(Int_t nEvents, TF1 *f1pt, TF1 *fpdg);
  static TChain *  MakeKineChain();
private:
  Float_t fNJets;                 // mean number of jets
  TF1 *   fF1Momentum;            // momentum distribution function
  TF1 *   fFPhi;                  // phi distribution function
  TF1 *   fFTheta;                // theta distribution function
  TF3 *   fFPosition;             // position distribution function 
  TF1 *   fFPdg;                  // pdg distribution function  
  TTreeSRedirector *fTestStream;  // test stream
  //
  ClassDef(AliGenPerformance,1) // performance generator
};

#endif
