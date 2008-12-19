Int_t AliTOFanalyzeSDigitsV2(Int_t ndump=0, Int_t iEvNum=-1, Int_t nEvent=0)
{

  /////////////////////////////////////////////////////////////////////////
  //
  // Analyzes the TOF sdigits and fills QA-histograms 
  // iEvNum=-1 and nEvent=0 means all events in the file
  // 
  // root[0] .L AliTOFanalyzeSDigitsV2.C
  // root[1] AliTOFanalyzeSDigitsV2()
  //
  // If you want analyze only the sdigits in the 2th event
  // (existing int the header file), see in the following:
  //
  // root[0] .L AliTOFanalyzeSDigitsV2.C
  // root[1] AliTOFanalyzeSDigitsV2(0,2,1)
  //
  // Updated to the new I/O by: A. De Caro, C. Zampolli
  //
  // Report problems to: decaro@sa.infn.it
  //
  /////////////////////////////////////////////////////////////////////////

  Int_t rc = 0;

  // adc and tdc
  TH1F *htdc     = new TH1F("htdc","TDC [bin]",5000,0.,150000.);
  TH1F *hadc     = new TH1F("hadc","ADC [bin]",100,0., 3000.);
  // TOF sdigit volumes
  TH1F *hsector  = new TH1F("hsector","Sector",18,0.,18.);
  TH1F *hplate   = new TH1F("hplate","Plate ", 5,0., 5.);
  TH1F *hstrip   = new TH1F("hstrip","Strip ",20,0.,20.);
  TH1F *hpadz    = new TH1F("hpadz","Pad along z ",2,0.,2.);
  TH1F *hpadx    = new TH1F("hpadx","Pad along x",48,0.,48.);
  // ADC-TDC correlation
  TH2F *h2tdcVSadc = new TH2F("h2tdcVSadc","TDC [bin] VS ADC [bin]",500,0.,150000.,100,0.,3000.);
  
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  if (gAlice) 
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
   }
  
  AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),"read");
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

  rl->LoadHeader();

  AliLoader *tofl = rl->GetLoader("TOFLoader");
  AliTOF *tof = (AliTOF *) rl->GetAliRun()->GetDetector("TOF");

  if (tof==0x0 || tofl==0x0)
    {
      cerr << "<AliTOFanalyzeSDigits> no TOF detector found" << endl;
      rc = 3;
      return rc;
    }
  
  cout << "First " << ndump << " SDigits found in TOF TreeS branch have: \n";

  Int_t upperLimit;
  Int_t bottomLimit;

  if (iEvNum<0) bottomLimit=0;
  else bottomLimit = iEvNum;

  if (nEvent == 0) upperLimit = (Int_t)(rl->GetNumberOfEvents());
  else upperLimit = nEvent+bottomLimit;

  AliTOFSDigit *tofsdigit;
  
  for (Int_t ievent = bottomLimit; ievent < upperLimit; ievent++) {

    rl->GetEvent(ievent);
    printf ("Processing event %d \n", ievent);
    
    // Get the pointer SDigit tree
    tofl->LoadSDigits("read");
    TTree *TS=tofl->TreeS();
    tof->SetTreeAddress();
    
    if(!TS)
      {
	cout << "<AliTOFanalyzeSDigits> No TreeS found \n";
	rc = 4;
	return rc;
      }
    
    TClonesArray * TOFsdigits = new TClonesArray("AliTOFSDigit",1000);
    //TOFsdigits = tof->SDigits();
    tof->SDigits();
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
	  Int_t sector = tofsdigit->GetSector(); // range [0-17]
	  Int_t plate  = tofsdigit->GetPlate();  // range [0- 4]
	  Int_t strip  = tofsdigit->GetStrip();  // range [0-19]
	  Int_t padz   = tofsdigit->GetPadz();   // range [0- 1]
	  Int_t padx   = tofsdigit->GetPadx();   // range [0-47]
	  // it is QA, then I perform QA!
	  Bool_t isSDigitBad = (sector<0 || sector>17 || plate<0 || plate >4 || padz<0 || padz>1 || padx<0 || padx>47);
	  

	  if (isSDigitBad) {
	    cout << "<AliTOFanalyzeSDigits>  strange sdigit found \n";
	    rc = 4;
	    return rc;
	  }
	  
	  if(k<ndump){
	    cout << k << "-th | Sector " << sector << " | Plate " << plate << " | Strip " << strip << " | PadZ " << padz << " | PadX " << padx << endl;
	    cout << k << "-th | ADC " << firstADC << " [bin] | TDC " << firstTDC << " [bin] \n";
	    cout << "---------------------------------------------------- \n";
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
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  return rc;
  
}

