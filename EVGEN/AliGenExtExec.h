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

  enum GenExtMode_t {
    kFIFO = 0,
    kAlternatingFiles
  };

  enum GenExtInput_t {
    kHepMC = 0,
    kEPOSroot
  };

  virtual void SetPathScript(const TString &path = "./gen.sh");
  virtual void SetPathFIFO(const TString &path = "gen.hepmc");
  virtual void SetPathFile1(const TString &path = "gen1.root");
  virtual void SetPathFile2(const TString &path = "gen2.root");
  virtual void SetMode(GenExtMode_t mode) { fMode = mode; }
  virtual void SetInput(GenExtInput_t input) { fInput = input; }

  virtual void Init();
  virtual void Generate();

protected:
  Bool_t StartGen();
  Bool_t StopGen();

  // configuration settings
  TString fPathScript;           // path to executable script
  TString fPathFIFO;             // path used for FIFO
  TString fPathFile1;            // path used for file 1
  TString fPathFile2;            // path used for file 2
  Int_t fEventNumberInFileMax;   // max number of events in file
  GenExtMode_t fMode;            // mode for external generator
  GenExtInput_t fInput;          // input type to choose reader

  // transient variables
  Int_t fPID;                    //! PID of running generator
  Int_t fFileConnected;          //! which file is connected (for alternating files)
  Int_t fEventNumberInFile;      //! event number in file

private:
  AliGenExtExec(const AliGenExtExec &ext);              // not implemented
  AliGenExtExec & operator=(const AliGenExtExec & rhs); // not implemented

  ClassDef(AliGenExtExec, 2)
};

#endif
