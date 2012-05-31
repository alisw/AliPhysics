// ----------------------------------------------------//
//                                                     //
//    This is a macro  to read PMD.SDigits.root         //
//                                                     //
// ----------------------------------------------------//
#include <Riostream.h>

void AliPMDSDigitsRead(Int_t nevt = 1) 
{
  TStopwatch timer;
  timer.Start();

  //FILE *fpw = fopen("sdigits.dat","w");

  TH2F *h2 = new TH2F("h2","Y vs. X",200,-100.,100.,200,-100.,100.);

  AliRunLoader *fRunLoader = AliRunLoader::Open("galice.root");

  if (!fRunLoader)
    { 
      Error("Open","Can not open session for file %s.",file);
    }
  

  Int_t nEvent = fRunLoader->GetNumberOfEvents();

  fPMDLoader = fRunLoader->GetLoader("PMDLoader");

  if (fPMDLoader == 0x0)
    {
      cerr<<"OpengAlice : Can not find PMD or PMDLoader\n";
    }
  
  fPMDLoader->LoadSDigits("READ");
  TClonesArray *fSDigits = 0x0; 
  
  // -------------------------------------------------------------- //
  
  Int_t    det = 0,smn = 0;
  Int_t    xpos, ypos;
  Int_t    xpad, ypad;
  Float_t  edep;
  Float_t  xx,yy;

  AliPMDUtility cc;  


  AliPMDsdigit  *pmdsdigit = 0x0;
  TBranch *branch;

  if (nevt > nEvent)
    {
      nevt = nEvent;
      printf("Input number is more than the number of events, set to number of events = %d\n",nevt);
    }
  
  for (Int_t ievt = 0; ievt <nevt; ievt++)
    {
      fRunLoader->GetEvent(ievt);
      fTreeS = fPMDLoader->TreeS();
      if (fTreeS == 0x0)
	{
	  cout << " Can not get TreeD" << endl;
	}

      branch = fTreeS->GetBranch("PMDSDigit");
      branch->SetAddress(&fSDigits);
      
      Int_t nmodules = (Int_t) fTreeS->GetEntries();
      
      for (Int_t imodule = 0; imodule < nmodules; imodule++)
	{
	  fTreeS->GetEntry(imodule); 
	  Int_t nentries = fSDigits->GetLast();
	  for (Int_t ient = 0; ient < nentries+1; ient++)
	    {
	      pmdsdigit = (AliPMDsdigit*)fSDigits->UncheckedAt(ient);
	      
	      det    = pmdsdigit->GetDetector();
	      smn    = pmdsdigit->GetSMNumber();
	      xpos   = pmdsdigit->GetRow();
	      ypos   = pmdsdigit->GetColumn();
	      edep   = pmdsdigit->GetCellEdep();
	      Int_t trno   = pmdsdigit->GetTrackNumber();

	      //if (ievt == 0) fprintf(fpw,"Track # = %d\n",trno);
	      
	      if(smn <12)
		{
		  xpad = ypos;
		  ypad = xpos;
		}
	      else if(smn >=12 && smn < 24)
		{
		  xpad = xpos;
		  ypad = ypos;
		}
	      
	      if(det == 1)
		{
		  cc.RectGeomCellPos(smn,xpad,ypad,xx,yy);
		  h2->Fill(xx,yy);  
		}
	    }
	} // modules

    }   
  h2->Draw();

  timer.Stop();
  timer.Print();
}

