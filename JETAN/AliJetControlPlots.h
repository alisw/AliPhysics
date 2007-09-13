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
class TH1F;
class TH1;

class AliJetReader;
class AliJet;

class AliJetControlPlots : public TObject
{
 public:
  AliJetControlPlots();
  ~AliJetControlPlots();

  // setter
  // getters
  TH1I *GetNJetsH()         const {return fNJetsH;}
  TH1I *GetMultH()          const {return fMultH;}
  TH1D *GetPhiH()           const {return fPhiH;}
  TH1D *GetFractionInJetH() const {return fInJetH;}
  TH1D *GetEneH()           const {return fEneH;}
  TH1D *GetPtH()            const {return fPtH;}
  TH1D *GetEtaH()           const {return fEtaH;}
  TH1D *GetFragH()          const {return fFragH;}
  TH1D *GetFragLnH()        const {return fFragLnH;}
  TH1D *GetFragrH()         const {return fFragrH;}
  TH1D *GetFragLnrH()       const {return fFragLnrH;}
  TH1D *GetShapeH()         const {return fShapeH;}
  TH1D *GetShaperH()        const {return fShaperH;}  
  
  // others
  void FillHistos(AliJet *j);
  void PlotHistos();
  void SetProperties(TH1* h,const char* x, const char* y) const;
  void Normalize();

 protected:
  AliJetControlPlots(const AliJetControlPlots& rControlPlots);
  AliJetControlPlots& operator = (const AliJetControlPlots& rcp);

  TH1I *fNJetsH;   // distribution of number of jets
  TH1I *fMultH;    // jet multiplicity
  TH1D *fPtH;      // pt spectra
  TH1D *fEtaH;     // eta distribution
  TH1D *fEneH;     // energy distribution
  TH1D *fFragH;    // leading jet fragmentation (selected part)
  TH1D *fFragLnH;  // leading jet fragmentation in ln scale
  TH1D *fFragrH;   // leading jet fragmentation (rejected part)
  TH1D *fFragLnrH; // leading jet fragmentation in ln scale
  TH1D *fShapeH;   // leading jet shape (selected part)
  TH1D *fShaperH;  // leading jet shape (rejected part)  
  TH1D *fPhiH;     // phi distribution
  TH1D *fInJetH;   // percentage of input particles in a jet
  Int_t fNJetT;    // total number of jets for normalization

  ClassDef(AliJetControlPlots,1)
};
 
#endif

