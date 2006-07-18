// ----------------------------------------------------//
//                                                     //
//    This is a macro  to read PMD.Digits.root         //
//                                                     //
// ----------------------------------------------------//

void AliPMDDigitsRead(Int_t nevt = 1) 
{
  TStopwatch timer;
  timer.Start();
  
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
  AliPMDUtility cc;

  TH2F *h2 = new TH2F("h2","Y vs. X",200,-100.,100.,200,-100.,100.);

  // -------------------------------------------------------------- //

  Int_t    det = 0,smn = 0;
  Int_t    xpos, ypos, xpad, ypad;
  Float_t  adc;
  Float_t  xx,yy;

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
	      
	      if(det==1)
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

