#ifndef ALIJETCONTROLPLOTS_H
#define ALIJETCONTROLPLOTS_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
 
//---------------------------------------------------------------------
// Jet Control Plots class 
// manages histograms with control plots of jet searching
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------


#include <TObject.h>

class TFile;
class TClonesArray;
class TH1I;
class TH1D;
class AliJetReader;
class AliJet;

class AliJetControlPlots : public TObject
{
 public:
  AliJetControlPlots();
  ~AliJetControlPlots();

  // setter
  // getters
  TH1I *GetNJetsH() {return fNJetsH;}
  TH1I *GetMultH() {return fMultH;}
  TH1D *GetEneH() {return fEneH;}
  TH1D *GetPtH() {return fPtH;}
  TH1D *GetEtaH() {return fEtaH;}
  TH1D *GetFragH() {return fFragH;}
  TH1D *GetFragLnH() {return fFragLnH;}
  TH1D *GetPhiH() {return fPhiH;}
  TH1D *GetFractionInJetH() {return fInJetH;}
  
  // others
  void FillHistos(AliJet *j,AliJetReader *r);
  void PlotHistos();
  void SetProperties(TH1* h,const char* x, const char* y) const;
  void Normalize();

 protected:
  TH1I *fNJetsH;    // distribution of number of jets
  TH1I *fMultH;    // jet multiplicity
  TH1D *fPtH;       // pt spectra
  TH1D *fEtaH;      // eta distribution
  TH1D *fEneH;      // energy distribution
  TH1D *fFragH;      // jet fragmentation
  TH1D *fFragLnH;      // jet fragmentation in ln scale
  TH1D *fPhiH;      // phi distribution
  TH1D *fInJetH;    // percentage of input particles in a jet
  Int_t fNJetT;     // total number of jets for normalization

  ClassDef(AliJetControlPlots,1)
};
 
#endif

