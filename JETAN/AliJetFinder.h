#ifndef ALIJETFINDER_H
#define ALIJETFINDER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet finder base class
// manages the search for jets 
// Authors: jgcn@mda.cinvestav.mx
//          andreas.morsch@cern.ch
//---------------------------------------------------------------------

#include <TObject.h>
#include "AliAODJet.h"
#include "AliJetReader.h"

class TFile;
class TTree;
class AliJet;
class AliJetHeader;
class AliJetControlPlots;
class AliLeading;
class AliAODJet;
class AliAODEvent;

class AliJetFinder : public TObject 
{
 public:
  AliJetFinder();
  virtual ~AliJetFinder();

  // getters
  virtual AliJetReader *GetReader() {return fReader;}
  virtual AliJetHeader *GetHeader() {return fHeader;}
  virtual AliJet *GetJets() {return fJets;}
  virtual Bool_t GetPlotMode() const {return fPlotMode;}
  virtual TFile* GetOutputFile() {return fOut;}
  // setters
  virtual void SetPlotMode(Bool_t b);
  virtual void SetOutputFile(const char *name="jets.root");
  virtual void SetJetReader(AliJetReader* r) {fReader=r;}
  virtual void SetJetHeader(AliJetHeader* h) {fHeader=h;}
  // others
  virtual void   AddJet(AliAODJet jet);
  virtual void   PrintJets();
  virtual void   Run();
  virtual void   WriteRHeaderToFile();  
  // the following have to be implemented for each specific finder
  virtual void Init() {}
  virtual void Reset() {fNAODjets = 0;}
  virtual void FindJets() {}
  virtual void FindJetsTPC(){}
  virtual void WriteJHeaderToFile() { }
  // some methods to allow steering from the outside
  virtual Bool_t ProcessEvent(Long64_t entry);
  virtual void   FinishRun();
  virtual void   ConnectTree(TTree* tree, TObject* data);
  virtual void   ConnectAOD(AliAODEvent* aod);
  virtual TTree* MakeTreeJ(char* name);
  virtual void   WriteHeaders();
  virtual void   WriteJetsToFile() {;}
  virtual void   TestJet() {;}
 protected:
  AliJetFinder(const AliJetFinder& rJetFinder);
  AliJetFinder& operator = (const AliJetFinder& rhsf);
  TTree* fTreeJ;                 //! pointer to jet tree
  Bool_t fPlotMode;              //! do you want control plots?
  AliJet* fJets;                 //! pointer to jet class
  AliJet* fGenJets;              //! pointer to generated jets
  AliLeading*   fLeading;        //! pointer to leading particle data 
  AliJetReader* fReader;         //  pointer to reader
  AliJetHeader* fHeader;         //  pointer to header
  TClonesArray* fAODjets;        //! reconstructed jets
  Int_t         fNAODjets;       //! number of reconstructed jets
  AliJetControlPlots* fPlots;    //! pointer to control plots
  TFile* fOut;                   //! output file
  ClassDef(AliJetFinder,2)
};

#endif
