#ifndef ALIJETFINDER_H
#define ALIJETFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet finder base class
// manages the search for jets 
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <TObject.h>

class TFile;
class AliJet;
class AliJetReader;
class AliJetControlPlots;
class AliLeading;

class AliJetFinder : public TObject 
{
 public:

  AliJetFinder();
  virtual ~AliJetFinder();

  // getters
  virtual AliJet *GetJets() {return fJets;}
  virtual Bool_t GetPlotMode() const {return fPlotMode;}
  virtual TFile* GetOutputFile() {return fOut;}
  // setters
  virtual void SetPlotMode(Bool_t b);
  virtual void SetOutputFile(const char *name="jets.root");
  virtual void SetJetReader(AliJetReader* r) {fReader=r;}

  // others
  virtual void PrintJets();
  virtual void Run();
  virtual void WriteJetsToFile(Int_t i);
  virtual void WriteRHeaderToFile();
  // the following have to be implemented for each specific finder
  virtual void Init() { }
  virtual void Reset() { }
  virtual void FindJets() { }
  virtual void WriteJHeaderToFile() { }
  virtual void GetGenJets();

 protected:
  AliJetFinder(const AliJetFinder& rJetFinder);
  AliJetFinder& operator = (const AliJetFinder& rhsf);

  Bool_t fPlotMode;              // do you want control plots?
  AliJet* fJets;                 // pointer to jet class
  AliJet* fGenJets;              // pointer to generated jets
  AliLeading* fLeading;          // pointer to leading particle data 
  AliJetReader* fReader;         // pointer to reader
  AliJetControlPlots* fPlots;    // pointer to control plots
  TFile* fOut;                   // output file

  ClassDef(AliJetFinder,1)
};

#endif
