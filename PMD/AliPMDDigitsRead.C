// ----------------------------------------------------//
//                                                     //
//    This is a macro  to read PMD.Digits.root         //
//                                                     //
// ----------------------------------------------------//

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include <stdlib.h>

void AliPMDDigitsRead(Int_t nevt = 1) 
{
  TStopwatch timer;
  timer.Start();
  
  //  FILE *fp = fopen("junk.dat","w");

  AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root");
				  
  if (!fRunLoader)
    {
      Error("Open","Can not open session for file %s.",file);
    }
  
  fRunLoader->LoadgAlice();
  fRunLoader->LoadHeader();
  gAlice = fRunLoader->GetAliRun();
  
  if (gAlice)
    {
      printf("AliRun object found on file.\n");
    }
  else
    {
      printf("Could not find AliRun object.\n");
    }
  fPMD  = (AliPMD*)gAlice->GetDetector("PMD");
  fPMDLoader = fRunLoader->GetLoader("PMDLoader");
  if (fPMDLoader == 0x0)
    {
      cerr<<"OpengAlice : Can not find PMD or PMDLoader\n";
    }
  
  fPMDLoader->LoadDigits("READ");
  TClonesArray *fDigits; 

  // -------------------------------------------------------------- //

  Int_t    det = 0,smn = 0;
  Int_t    xpos,ypos;
  Float_t  adc;
  Int_t    isup;
  Int_t    idet;
  Float_t  clusdata[7];
  
  for (Int_t ievt = 2; ievt <nevt; ievt++)
    {
      fRunLoader->GetEvent(ievt);
  
      fTreeD = fPMDLoader->TreeD();
      if (fTreeD == 0x0)
	{
	  cout << " Can not get TreeD" << endl;
	}
      AliPMDdigit  *pmddigit;
      TBranch *branch = fTreeD->GetBranch("PMDDigit");
      branch->SetAddress(&fDigits);
  
      Int_t nmodules = (Int_t) fTreeD->GetEntries();

      cout << " Total number of modules in an event = " << nmodules << endl;

      for (Int_t imodule = 0; imodule < nmodules; imodule++)
	{
	  fTreeD->GetEntry(imodule); 
	  Int_t nentries = fDigits->GetLast();
	  cout << "nentries = " << nentries << endl;
	  for (Int_t ient = 0; ient < nentries+1; ient++)
	    {
	      pmddigit = (AliPMDdigit*)fDigits->UncheckedAt(ient);
	      
	      det    = pmddigit->GetDetector();
	      smn    = pmddigit->GetSMNumber();
	      xpos   = pmddigit->GetRow();
	      ypos   = pmddigit->GetColumn();
	      adc    = pmddigit->GetADC();
	      Int_t trno   = pmddigit->GetTrackNumber();
	      
	      //fprintf(fp,"%d %d %d %d %f \n ",ievt,smn,xpos,ypos,adc);
	      
	    }
	} // modules
    }

  timer.Stop();
  timer.Print();
}

