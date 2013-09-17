#ifndef ALIJETHISTOS_H
#define ALIJETHISTOS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------
// Jet Histos class
// Creates and fills a few cummon histograms for jet analysis
//
//---------------------------------------------------------------------

class TList;
class TClonesArray;
class TH1I;
class TH1F;

class AliJetHistos : public TObject {
 public:
  AliJetHistos();
  ~AliJetHistos();

  void CreateHistos();
  void AddHistosToList(TList *list) const;
  void FillHistos(TClonesArray *jets);
  
 private:
  void SetProperties(TH1* h,const char* x, const char* y) const;
  AliJetHistos(const AliJetHistos &det);
  AliJetHistos &operator=(const AliJetHistos &det);

 private:

  TH1I *fNJetsH;           // distribution of number of jets
  TH1F *fPtH;              // pt spectra
  TH1F *fEtaH;             // eta distribution
  TH1F *fEneH;             // energy distribution
  TH1F *fPhiH;             // phi distribution

  ClassDef(AliJetHistos,2) // some jet histos

};

#endif

    
