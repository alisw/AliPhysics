#ifndef ALICDFJETHEADER_H
#define ALICDFJETHEADER_H

/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 *
*/

// settings for jet finder process

#include "AliJetHeader.h"

class AliCdfJetHeader : public AliJetHeader
{
 public:
  AliCdfJetHeader();
  virtual ~AliCdfJetHeader() { }

  // Getters
  Double_t GetJetPtCut() const      { return fJetPtCut; }
  Int_t    GetMinPartJet() const    { return fMinPartJet; }
  Bool_t   GetAnalyseJets() const   { return fkAnalyseJets; }

  // Setters
  void     SetJetPtCut(Double_t jetptcut)          { fJetPtCut = jetptcut; }
  void     SetAODwrite(Bool_t aodwrite)            { fAODwrite = aodwrite; }
  void     SetAODtracksWrite(Bool_t aodtrackswrite){ fAODtracksWrite = aodtrackswrite; }
  void     SetMinPartJet(Int_t npart)              { fMinPartJet = npart; }
  void     SetAnalyseJets(Bool_t flag = kTRUE)     { fkAnalyseJets = flag; }

  Bool_t   IsAODwrite() const       { return fAODwrite; }
  Bool_t   IsAODtracksWrite() const { return fAODtracksWrite; }

 protected:
  // Parameters of algorithm
  Int_t    fMinPartJet;           // minimum number of particles in jet
  Double_t fJetPtCut;             // pt cut of jets

  Bool_t   fAODwrite;             // flag for writing to AOD
  Bool_t   fAODtracksWrite;       // flag for writing tracks to AOD

  Bool_t   fkAnalyseJets;         // analyse jets 

  ClassDef ( AliCdfJetHeader, 2 ) // CDF jet header class

};

#endif
