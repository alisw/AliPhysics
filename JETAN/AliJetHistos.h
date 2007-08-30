#ifndef ALIJETHISTOS_H
#define ALIJETHISTOS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class Tlist;
class TClonesArray;
class TH1I;
class TH1F;

class AliJetHistos : public TObject {
 public:
  AliJetHistos();
  ~AliJetHistos();

  void AddHistosToList(TList *list);
  void FillHistos(TClonesArray *jets);

 protected:
  void SetProperties(TH1* h,const char* x, const char* y) const;
 private:

  TH1I *fNJetsH;   // distribution of number of jets
  TH1F *fPtH;      // pt spectra
  TH1F *fEtaH;     // eta distribution
  TH1F *fEneH;     // energy distribution
  TH1F *fPhiH;     // phi distribution


  ClassDef(AliJetHistos,1)
};

#endif

