#ifndef ALIJETSELECTOR_H
#define ALIJETSELECTOR_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//---------------------------------------------------------------------
// Jet finder base class
// manages the search for jets 
// Author: jgcn@mda.cinvestav.mx
//---------------------------------------------------------------------

#include <TSelector.h>

class TTree;
class AliJetFinder;

class AliJetSelector : public TSelector 
{
 public:
  AliJetSelector(TTree* tree = 0);
  virtual ~AliJetSelector();
  virtual void   Config();
  virtual Int_t  Version() const {return 1;}
  virtual void   Begin(TTree* tree) ;
  virtual void   SlaveBegin(TTree* tree);
  virtual void   Init(TTree* tree);
  virtual Bool_t Notify();
  virtual Bool_t Process(Long64_t entry);
  virtual void   SetOption(const char *option) { fOption = option; }
  virtual void   SetObject(TObject *obj) { fObject = obj; }
  virtual void   SetInputList(TList *input) {fInput = input;}
  virtual TList* GetOutputList() const { return fOutput; }
  virtual void   SlaveTerminate();
  virtual void   Terminate();
  
 protected:
  AliJetFinder* fJetFinder;
  
  ClassDef(AliJetSelector, 1)
};

#endif
