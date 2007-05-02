/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


#if !defined(__CINT__) || defined(__MAKECINT__)
//#include "iostream.h"

#include <TClassTable.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TTree.h>

#include "AliHeader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliMagF.h"

#include "AliMUON.h"
#include "AliMUONDisplay.h"
#endif



void MUONdisplay (Int_t nevent=0, 
                  TString fileName="galice.root",
                  TString fileNameSim="galice_sim.root") {

  // set off mag field 
  AliMagF::SetReadField(kFALSE);
 
  // Getting runloader 
  AliRunLoader * RunLoaderSim = AliRunLoader::Open(fileNameSim.Data(),"MUONFolderSim","READ");
  if (RunLoaderSim == 0x0) {
    Error("MUONdisplay","Inut file %s error!",fileName.Data());
    return;   
  }
  RunLoaderSim->LoadHeader();
  RunLoaderSim->LoadKinematics("READ");

  //  if (RunLoader->GetAliRun() == 0x0) 
  RunLoaderSim->LoadgAlice();
  gAlice = RunLoaderSim->GetAliRun();

  // Getting MUONloader 
  AliLoader * MUONLoaderSim  = RunLoaderSim->GetLoader("MUONLoader");
  MUONLoaderSim->LoadHits("READ");
  MUONLoaderSim->LoadDigits("READ");
  
  // Getting runloader 
  AliRunLoader * RunLoaderRec = AliRunLoader::Open(fileName.Data(),"MUONFolder","READ");
  if (RunLoaderRec == 0x0) {
    Error("MUONdisplay","Inut file %s error!",fileName.Data());
    return;   
  }
  AliLoader * MUONLoaderRec  = RunLoaderRec->GetLoader("MUONLoader");
  MUONLoaderRec->LoadRecPoints("READ");
  MUONLoaderRec->LoadTracks("READ");


  // Create Event Display object
  AliMUONDisplay *muondisplay = new AliMUONDisplay(750, MUONLoaderSim, MUONLoaderRec);

  // Display first event
  RunLoaderSim->GetEvent(nevent);
  RunLoaderRec->GetEvent(nevent);
  muondisplay->ShowNextEvent(0);
}
