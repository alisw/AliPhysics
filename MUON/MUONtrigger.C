#include "iostream.h"
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
  muondata->SetLoader(MUONLoader);

  Int_t ievent, nevents;
  nevents = RunLoader->GetNumberOfEvents();

  MUONLoader->LoadDigits("READ");
  MUONLoader->LoadRecPoints("UPDATE");
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (evNumber2>nevents) evNumber2=nevents;
   for (Int_t ievent=evNumber1; ievent<evNumber2; ievent++) { // event loop
       printf("event %d\n",ievent);
       RunLoader->GetEvent(ievent);       
       if (MUONLoader->TreeR() == 0x0) { 
	 MUONLoader->MakeTree("R");
       }
       muondata->MakeBranch("GLT");
       muondata->SetTreeAddress("D,GLT");
       MUON->Trigger(ievent); 
       muondata->ResetDigits();
       muondata->ResetTrigger();
   } // event loop 
   MUONLoader->UnloadDigits();
   MUONLoader->UnloadRecPoints();
}














