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

/* $Id$ */

// Macro MUONTracker.C (TO BE COMPILED)
// for testing the C++ reconstruction code
// Output is using aliroot standard output MUON.Tracks.root
// The output is a TClonesArray of AliMUONTracks.
// Gines MARTINEZ, Subatech, sep 2003

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TClonesArray.h>

#include "AliRun.h"
#include "AliESD.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONEventReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackHit.h"
#include "AliMUONTrackParam.h"
#endif

void MUONRecoTest(Text_t *FileName = "galice.root")
{
  //
  cout << "MUONRecoTest" << endl;
  cout << "FileName ``" << FileName << "''" << endl;
  
  // Creating Run Loader and openning file containing Hits, Digits and RecPoints
  AliRunLoader * RunLoader = AliRunLoader::Open(FileName,"Event","UPDATE");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",FileName);
    return;
  }
  // Loading AliRun master
  RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();
 
  // Loading MUON subsystem
  AliMUON* MUON = (AliMUON *) gAlice->GetDetector("MUON");

  //AliESD *event = new AliESD(); 
  MUON->Reconstruct();
  //  MUON->FillESD(event);
}
