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



void MUONdisplay (Int_t nevent=0, TString fileName="galice.root") {

  // set off mag field 
  AliMagF::SetReadField(kFALSE);
 
  // Getting runloader 
  AliRunLoader * RunLoader = AliRunLoader::Open(fileName.Data(),"MUONFolder","READ");
  if (RunLoader == 0x0) {
    Error("MUONdisplay","Inut file %s error!",fileName.Data());
    return;   
  }
  RunLoader->LoadHeader();
  RunLoader->LoadKinematics("READ");

  //  if (RunLoader->GetAliRun() == 0x0) 
  RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();

  // Getting MUONloader 
  AliLoader * MUONLoader  = RunLoader->GetLoader("MUONLoader");
  MUONLoader->LoadHits("READ");
  MUONLoader->LoadDigits("READ");
  MUONLoader->LoadRecPoints("READ");
  MUONLoader->LoadTracks("READ");


  // Create Event Display object
  AliMUONDisplay *muondisplay = new AliMUONDisplay(750, MUONLoader);

  // Display first event
  RunLoader->GetEvent(nevent);
  muondisplay->ShowNextEvent(0);
}
