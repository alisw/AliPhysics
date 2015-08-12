/**************************************************************************
 * Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Event generator that can runs an external executable
// to produce events which are read as HepMC.
// 
// Author: Jochen Klein <Jochen.Klein@cern.ch>

#include <sys/types.h>
#include <sys/stat.h>

#include "TSystem.h"

#include "AliLog.h"
#include "AliGenReaderHepMC.h"

#include "AliGenExtExec.h"

AliGenExtExec::AliGenExtExec(const TString &scriptpath) :
  AliGenExtFile(),
  fPathScript(scriptpath),
  fPathFIFO("gen.hepmc"),
  fPID(0)
{
  // default constructor

  gSystem->ExpandPathName(fPathScript);
}

AliGenExtExec::AliGenExtExec(Int_t npart, const TString &scriptpath) :
  AliGenExtFile(npart),
  fPathScript(scriptpath),
  fPID(0)
{
  // constructor

  gSystem->ExpandPathName(fPathScript);
}

AliGenExtExec::~AliGenExtExec()
{
  // destructor

  StopGen();
}

void AliGenExtExec::SetPathScript(const TString &path)
{
  // set and expand path to executable script

  fPathScript = path;
  gSystem->ExpandPathName(fPathScript);
}

void AliGenExtExec::SetPathFIFO(const TString &path)
{
  // set and expand path used for FIFO

  fPathFIFO = path;
  gSystem->ExpandPathName(fPathFIFO);
}


void AliGenExtExec::Init()
{
  // initialization:
  // start generator, setup reader, and initialize AliGenExtFile

  // setup and run the generator
  StartGen();

  // create and set HepMC reader
  AliGenReaderHepMC *reader = new AliGenReaderHepMC();
  reader->SetFileName(fPathFIFO);
  SetReader(reader);

  // proceed with init using FIFO
  AliGenExtFile::Init();
}

Bool_t AliGenExtExec::StartGen()
{
  // prepare FIFO and start the generator (if not yet running)
  // by forking the script provided at fPathScript

  if (fPID != 0) {
    AliError(Form("generator already running with PID %i", fPID));
    return kFALSE;
  }

  // create FIFO
  mkfifo(fPathFIFO, 0666);

  // fork generator
  fPID = fork();
  if (fPID < 0) {
    AliError("forking generator failed");
    return kFALSE;
  } else if (fPID == 0) {
    execl("/bin/bash", "bash", "-c", (fPathScript + " " + fPathFIFO + " > gen.log 2>&1").Data(), (char *) 0);
  } else {
    AliInfo(Form("running generator with PID %i", fPID));
  }

  return kTRUE;
}

Bool_t AliGenExtExec::StopGen()
{
  // kill generator if running

  if (fPID == 0) {
    AliError("generator not running, nothing killed");
    return kFALSE;
  }

  gSystem->Exec(Form("kill %i", fPID));
  gSystem->Unlink(fPathFIFO);

  fPID = 0;

  return kTRUE;
}
