Int_t AliTOFanalyzeSDigitsV2(Int_t ndump=15, Int_t iEvNum=0)
{
  //
  // Analyzes the TOF sdigits and fills QA-histograms 
  // report problems to pierella@bo.infn.it
  // iEvNum=0 means all events in the file
  //
  // Updated to the new I/O by: A. De Caro, C. Zampolli

  Int_t rc = 0;

  // adc and tdc
  TH1F *htdc  = new TH1F("htdc","TDC [bin]",5000,0.,150000.);
  TH1F *hadc   = new TH1F("hadc","ADC [bin]",100,0., 3000.);
  // TOF sdigit volumes
  TH1F *hsector  = new TH1F("hsector","Sector",20,0.,20.);
  TH1F *hplate   = new TH1F("hplate","Plate ", 6,0., 6.);
  TH1F *hstrip   = new TH1F("hstrip","Strip ",25,0.,25.);
  TH1F *hpadz    = new TH1F("hpadz","Pad along z ",3,0.,3.);
  TH1F *hpadx    = new TH1F("hpadx","Pad along x",50,0.,50.);
  // ADC-TDC correlation
  TH2F *h2tdcVSadc = new TH2F("h2tdcVSadc","TDC [bin] VS ADC [bin]",500,0.,150000.,100,0.,3000.);
  
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  if (gAlice) 
    {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
   }
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::fgkDefaultEventFolderName,"read");
  if (!rl)
    {
      cerr<<"Can't load RunLoader from file"<<"!\n";
      rc = 1;
      return rc;
    }

  rl->LoadgAlice();
  gAlice = rl->GetAliRun();

  if (!gAlice)
    {
      cerr << "<AliTOFanalyzeSDigits> AliRun object not found on file\n ";
      rc = 2;
      return rc;
    }

  AliLoader *tofl = rl->GetLoader("TOFLoader");
  AliTOF *tof = (AliTOF *) rl->GetAliRun()->GetDetector("TOF");

  if (tof==0x0 || tofl==0x0)
    {
      cerr << "<AliTOFanalyzeSDigits> no TOF detector found" << endl;
      rc = 3;
      return rc;
    }
  
  cout << "First " << ndump << " SDigits found in TOF TreeS branch have:" << endl;

  if (iEvNum == 0)
    {
      rl->LoadHeader();
      TTree *TE = rl->TreeE();
      iEvNum = (Int_t)TE->GetEntries();
    }

  AliTOFSDigit *tofsdigit;
  
  for (Int_t ievent = 0; ievent < iEvNum; ievent++) {
    printf ("Processing event %d \n", ievent);
    rl->GetEvent(ievent);
    
    // Get the pointer SDigit tree
    tofl->LoadSDigits();
    TTree *TS=tofl->TreeS();
    tof->SetTreeAddress();
    
    if(!TS)
      {
	cout << "<AliTOFanalyzeSDigits> No TreeS found" << endl;
	rc = 4;
	return rc;
      }
    
    TClonesArray * TOFsdigits = new TClonesArray("AliTOFSDigit",1000);
    TOFsdigits = tof->SDigits();
    TOFsdigits = TS->GetBranch("TOF")->SetAddress(&TOFsdigits); 

    Int_t nEntries = TS->GetEntries(); 

    for (Int_t iEntry = 0; iEntry < nEntries; iEntry ++) 
      {
	tof->ResetDigits();
	TS->GetEvent(iEntry);
	Int_t ndig = TOFsdigits->GetEntriesFast();
	cout << "<AliTOFanalyzeSDigits> found " << ndig
	     << " TOF sdigits for event " << ievent << endl;
	
	for (Int_t k=0; k<ndig; k++) {
	  tofsdigit = (AliTOFSDigit*) TOFsdigits->UncheckedAt(k);
	  Float_t firstTDC = tofsdigit->GetTdc(0);
	  Float_t firstADC = tofsdigit->GetAdc(0);
	  htdc->Fill(firstTDC);
	  hadc->Fill(firstADC);
	  // TOF sdigit volumes
	  Int_t sector = tofsdigit->GetSector(); // range [1-18]
	  Int_t plate  = tofsdigit->GetPlate();  // range [1- 5]
	  Int_t strip  = tofsdigit->GetStrip();  // range [1-20]
	  Int_t padz   = tofsdigit->GetPadz();   // range [1- 2]
	  Int_t padx   = tofsdigit->GetPadx();   // range [1-48]
	  // it is QA, then I perform QA!
	  Bool_t isSDigitBad = (sector<1 || sector>18 || plate<1 || plate >5 || padz<1 || padz>2 || padx<1 || padx>48);
	  
	  if (isSDigitBad) {
	    cout << "<AliTOFanalyzeSDigits>  strange sdigit found" << endl;
	    rc = 4;
	    return rc;
	  }
	  
	  if(k<ndump){
	    cout << k << "-th | " << "Sector " << sector << " | Plate " << plate << " | Strip " << strip << " | PadZ " << padz << " | PadX " << padx << endl;
	    cout << k << "-th | ADC " << firstADC << " [bin] | TDC " << firstTDC << " [bin]" << endl;
	    cout << "----------------------------------------------------"<< endl;
	  }
	  
	  // filling sdigit volume histos
	  hsector->Fill(sector);
	  hplate->Fill(plate);
	  hstrip->Fill(strip);
	  hpadx->Fill(padx);
	  hpadz->Fill(padz);
	  h2tdcVSadc->Fill(firstTDC,firstADC);
	  
	}
      }
    
    tofl->UnloadSDigits();
  
  } // end loop on events
  
  rl->UnloadHeader();
  rl->UnloadgAlice();

  TFile *fout = new TFile("TOF_sdigitsQA.root","RECREATE");
  htdc->Write();
  hadc->Write();
  h2tdcVSadc->Write();
  hsector->Write();
  hplate->Write();
  hstrip->Write();
  hpadz->Write();
  hpadx->Write();
  fout->Close();

  delete htdc;
  delete hadc;
  delete h2tdcVSadc;
  delete hsector;
  delete hplate;
  delete hstrip;   
  delete hpadz;
  delete hpadx;

  if (gAlice) 
    {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  return rc;
  
}

