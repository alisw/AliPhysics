// 0 = all
// 1 = not pion
// 2 = not kaon
// 3 = not proton
// 4 = not muon
// 5 = not electron
// 6 = not neutron


Int_t particle_type=0;

#include "iostream.h"

void RICHdigit (Int_t nEvents = 1,Int_t type = 0)
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    }else {
      delete gAlice;
      gAlice = 0;
   }

    if(type)
      {
	AliRunDigitizer * manager = new AliRunDigitizer(2,1);
	manager->SetInputStream(0,"galice.root");
	manager->SetInputStream(1,"bgr.root");
	manager->SetOutputFile("galice.root");
	AliRICHDigitizer *dRICH  = new AliRICHDigitizer(manager);
	manager->Exec("deb");
      }
    else
      {
	AliRunDigitizer * manager = new AliRunDigitizer(1,1);
	manager->SetInputStream(0,"galice.root");
	manager->SetNrOfEventsToWrite(nEvents);
	AliRICHDigitizer *dRICH  = new AliRICHDigitizer(manager);
	manager->Exec("deb");
      }

   //delete gAlice;
  printf("\nEnd of Macro  *************************************\n");
}



