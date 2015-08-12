#ifndef ALIGENEXTEXEC_H
#define ALIGENEXTEXEC_H

/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Event generator that can runs an external executable
// to produce events which are read as HepMC.
// 
// Author: Jochen Klein <Jochen.Klein@cern.ch>

#include "AliGenExtFile.h"

class AliGenReader;
class TTree;

class AliGenExtExec : public AliGenExtFile
{
public:
  AliGenExtExec(const TString &scriptpath = "./gen.sh");
  AliGenExtExec(Int_t npart, const TString &scriptpath = "./gen.sh");

  virtual ~AliGenExtExec();

  virtual void SetPathScript(const TString &path = "./gen.sh");
  virtual void SetPathFIFO(const TString &path = "gen.hepmc");

  virtual void Init();

protected:
  Bool_t StartGen();
  Bool_t StopGen();

  TString fPathScript;    // path to executable script
  TString fPathFIFO;      // path used for FIFO

  Int_t fPID;             // PID of running generator

private:
  AliGenExtExec(const AliGenExtExec &ext);              // not implemented
  AliGenExtExec & operator=(const AliGenExtExec & rhs); // not implemented

  ClassDef(AliGenExtExec, 1)
};

#endif
