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
#include <TClonesArray.h>

#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#endif

void MUONHit2SDigit(Int_t iEveMin, Int_t iEveMax, Text_t *FileName = "galice.root")
{
  //
  cout << "MUONSDigitTest" << endl;
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
  AliLoader    * gime; 
  gime        = RunLoader->GetLoader("MUONLoader");
 
  gime->LoadHits("READ");
  gime->LoadSDigits("RECREATE");

  // Loading MUON subsystem
  AliMUON* MUON = (AliMUON *) gAlice->GetDetector("MUON");

  Int_t nEvents = RunLoader->GetNumberOfEvents();
  Int_t nEveMax = TMath::Min(iEveMax,nEvents-1);
  for (Int_t iEvent = iEveMin; iEvent <= nEveMax; iEvent++){
    RunLoader->GetEvent(iEvent);
    MUON->Hits2SDigits();
  }
}
