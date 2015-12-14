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
// to produce events which are read in HepMC or ROOT format.
//
// Author: Jochen Klein <Jochen.Klein@cern.ch>

#include <sys/types.h>
#include <sys/stat.h>

#include "TSystem.h"

#include "AliLog.h"
#include "AliGenReaderHepMC.h"
#include "AliGenEposReader.h"

#include "AliGenExtExec.h"

AliGenExtExec::AliGenExtExec(const TString &scriptpath) :
  AliGenExtFile(),
  fPathScript(scriptpath),
  fPathFIFO("gen.hepmc"),
  fPathFile1("gen1.root"),
  fPathFile2("gen2.root"),
  fEventNumberInFileMax(100),
  fMode(kFIFO),
  fInput(kHepMC),
  fPID(0),
  fFileConnected(0),
  fEventNumberInFile(0)
{
  // default constructor

  gSystem->ExpandPathName(fPathScript);
}

AliGenExtExec::AliGenExtExec(Int_t npart, const TString &scriptpath) :
  AliGenExtFile(npart),
  fPathScript(scriptpath),
  fPathFIFO("gen.hepmc"),
  fPathFile1("gen1.root"),
  fPathFile2("gen2.root"),
  fEventNumberInFileMax(100),
  fMode(kFIFO),
  fInput(kHepMC),
  fPID(0),
  fFileConnected(0),
  fEventNumberInFile(0)
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

void AliGenExtExec::SetPathFile1(const TString &path)
{
  // set and expand path used for file 1

  fPathFile1 = path;
  gSystem->ExpandPathName(fPathFile1);
}

void AliGenExtExec::SetPathFile2(const TString &path)
{
  // set and expand path used for file 2

  fPathFile2 = path;
  gSystem->ExpandPathName(fPathFile2);
}

void AliGenExtExec::Init()
{
  // initialization:
  // start generator, setup reader, and initialize AliGenExtFile

  // setup and run the generator
  StartGen();

  // create and set HepMC reader
  if (fInput == kHepMC) {
    AliGenReaderHepMC *reader = new AliGenReaderHepMC();
    SetReader(reader);
  } else if (fInput == kEPOSroot) {
    AliGenEposReader *reader = new AliGenEposReader();
    SetReader(reader);
  } else {
    AliFatal(Form("Invalid input type: %i", fInput));
  }

  if (fMode == kFIFO) {
    Reader()->SetFileName(fPathFIFO);
    AliGenExtFile::Init();
  }
}

void AliGenExtExec::Generate()
{
  // generate next event
  // depending on mode read from FIFO or files

  switch (fMode) {
  case kFIFO:
    AliGenExtFile::Generate();
    break;

  case kAlternatingFiles:
    // connect a new file if needed
    // (i.e. no file yet or reached the end)
    while ((fFileConnected == 0) ||
	   (fEventNumberInFile == fEventNumberInFileMax)) {

      if ((fFileConnected != 1) &&
	  !gSystem->AccessPathName(fPathFile1)) {
	if (fInput == kEPOSroot)
	  ((AliGenEposReader*) Reader())->ChangeFile(fPathFile1);
	else
	  AliFatal(Form("Alternating files not supported for input type: %i", fInput));

	if (fFileConnected == 2)
	  gSystem->Unlink(fPathFile2);
	fFileConnected = 1;
	fEventNumberInFile = 0;
      } else if ((fFileConnected != 2) &&
		 !gSystem->AccessPathName(fPathFile2)) {
	if (fInput == kEPOSroot)
	  ((AliGenEposReader*) Reader())->ChangeFile(fPathFile2);
	else
	  AliFatal(Form("Alternating files not supported for input type: %i", fInput));

	if (fFileConnected == 1)
	  gSystem->Unlink(fPathFile1);
	fFileConnected = 2;
	fEventNumberInFile = 0;
      } else {
	// wait for further input to become available
	sleep(1);
      }
    }

    fEventNumberInFile++;

    AliGenExtFile::Generate();

    break;

  default:
    AliFatal(Form("Unknown mode for external generator: %i", fMode));
  }
}

Bool_t AliGenExtExec::StartGen()
{
  // prepare FIFO and start the generator (if not yet running)
  // by forking the script provided at fPathScript

  if (fPID != 0) {
    AliError(Form("generator already running with PID %i", fPID));
    return kFALSE;
  }

  switch (fMode) {
  case kFIFO:
    // create FIFO and attach reader to it
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
    break;

  case kAlternatingFiles:
    // fork generator
    fPID = fork();
    if (fPID < 0) {
      AliError("forking generator failed");
      return kFALSE;
    } else if (fPID == 0) {
      execl("/bin/bash", "bash", "-c", (fPathScript + " " + fPathFile1 + " " + fPathFile2 + " > gen.log 2>&1").Data(), (char *) 0);
    } else {
      AliInfo(Form("running generator with PID %i", fPID));
    }
    break;

  default:
    AliFatal(Form("Unknown mode for external generator: %i", fMode));
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
