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

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliMUON.h"
#include "AliMUONData.h"
#endif

AliRun * gAlice;

//get trigger decision and write it in TreeR of MUON.RecPoints.root

void MUONtrigger (char* filename="galice.root", 
		  Int_t evNumber1=0, Int_t evNumber2=9999)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Creating Run Loader and openning file containing Hits
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","UPDATE");
  if (RunLoader ==0x0) {
    printf(">>> Error : Error Opening %s file \n",filename);
    return;
  }

// Loading AliRun master
  RunLoader->UnloadgAlice();
  RunLoader->LoadgAlice();
  gAlice = RunLoader->GetAliRun();

// Loading MUON subsystem
  AliMUON * MUON = (AliMUON *) gAlice->GetDetector("MUON");
  AliLoader * MUONLoader = RunLoader->GetLoader("MUONLoader");
  AliMUONData * muondata = MUON->GetMUONData();

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  MUONLoader->LoadDigits("READ");
  MUONLoader->LoadRecPoints("UPDATE");

// Testing if trigger has already been done
  RunLoader->GetEvent(0);
  if (MUONLoader->TreeR()) {
      if (muondata->IsTriggerBranchesInTree()) {
	  MUONLoader->UnloadRecPoints();
	  MUONLoader->LoadRecPoints("RECREATE");
	  printf("Recreating RecPoints files\n");
      }
  }
  
//   Loop over events
  if (evNumber2>nevents) evNumber2=nevents;
  for (Int_t ievent=evNumber1; ievent<evNumber2; ievent++) { // event loop
      printf("event %d\n",ievent);
      RunLoader->GetEvent(ievent);       
      if (MUONLoader->TreeR() == 0x0)  MUONLoader->MakeRecPointsContainer();
      
      muondata->MakeBranch("GLT");
      muondata->SetTreeAddress("D,GLT");
      MUON->Trigger(ievent); 
      muondata->ResetDigits();
      muondata->ResetTrigger();
  } 
  MUONLoader->UnloadDigits();
  MUONLoader->UnloadRecPoints();
}














