Int_t AliTOFanalyzeDigits(TString headersFile, Int_t ndump=15, Int_t iEvNum=0)
{
  //
  // Analyzes the TOF digits and fills QA-histograms 
  // report problems to pierella@bo.infn.it
  // iEvNum=0 means all events in the file
  // Author: F. Pierella (Bologna University)

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice = 0;
  }
  
  Int_t rc=0;

  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(headersFile.Data());
  if(file){
    cout<<"headerFile already open \n";
  }
  else {
    if(!file)file=TFile::Open(headersFile.Data());
  }
  
  // Get AliRun object from file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
  }
  

  if (iEvNum == 0) iEvNum = (Int_t) gAlice->TreeE()->GetEntries();


  AliTOFdigit *tofdigit;

  AliTOF * tof = (AliTOF *) gAlice->GetDetector("TOF") ;

  if (!tof) {
    cout << "<AliTOFanalyzeDigits> No TOF detector found" << endl;
    rc = 2;
    return rc;
  }

  // adc and tdc
  TH1F *htdc  = new TH1F("htdc","TDC [bin]",500,0.,15000.);
  TH1F *hadc   = new TH1F("hadc","ADC [bin]",100,0., 3000.);

  // TOF digit volumes
  TH1F *hsector  = new TH1F("hsector","Sector",20,0.,20.);
  TH1F *hplate   = new TH1F("hplate","Plate ", 6,0., 6.);
  TH1F *hstrip   = new TH1F("hstrip","Strip ",25,0.,25.);
  TH1F *hpadz    = new TH1F("hpadz","Pad along z ",3,0.,3.);
  TH1F *hpadx    = new TH1F("hpadx","Pad along x",50,0.,50.);
  // ADC-TDC correlation
  TH2F *h2tdcVSadc = new TH2F("h2tdcVSadc","TDC [bin] VS ADC [bin]",100,0.,3000.,500,0.,15000.);

  cout << "First " << ndump << " Digits found in TOF TreeD branch have:" << endl;

  for (Int_t ievent = 0; ievent < iEvNum; ievent++) {

    gAlice->GetEvent(ievent) ;
    if(gAlice->TreeD()==0) {
      cout << "<AliTOFanalyzeDigits> No  TreeD found" << endl;
      rc = 4;
      return rc;
    }
    
    
    
    Int_t ndig, k;
    gAlice->ResetDigits();
    gAlice->TreeD()->GetEvent(ievent);
    TClonesArray * TOFdigits   = tof->Digits();
    
    ndig=TOFdigits->GetEntries();
    
    cout << "<AliTOFanalyzeDigits> found " << ndig
	 << " TOF digits for event " << ievent << endl;
    
    for (k=0; k<ndig; k++) {
      tofdigit= (AliTOFdigit*) TOFdigits->UncheckedAt(k);
      Float_t tdc=tofdigit->GetTdc();
      Float_t adc=tofdigit->GetAdc();
      htdc->Fill(tdc);
      hadc->Fill(adc);
      // TOF digit volumes
      Int_t sector    = tofdigit->GetSector(); // range [1-18]
      Int_t plate     = tofdigit->GetPlate();  // range [1- 5]
      Int_t strip     = tofdigit->GetStrip();  // range [1-20]
      Int_t padz      = tofdigit->GetPadz();   // range [1- 2]
      Int_t padx      = tofdigit->GetPadx();   // range [1-48]
      // it is QA, then I perform QA!
      Bool_t isDigitBad = (sector<1 || sector>18 || plate<1 || plate >5 || padz<1 || padz>2 || padx<1 || padx>48);

      if (isDigitBad) {
	cout << "<AliTOFanalyzeDigits>  strange digit found" << endl;
	rc = 3;
	return rc;
      }
      
      if(k<ndump){
	cout << k << "-th | " << "Sector " << sector << " | Plate " << plate << " | Strip " << strip << " | PadZ " << padz << " | PadX " << padx << endl;
	cout << k << "-th | ADC " << adc << " [bin] | TDC " << tdc << " [bin]" << endl;
	cout << "----------------------------------------------------"<< endl;
      }

      // filling digit volume histos
      hsector->Fill(sector);
      hplate->Fill(plate);
      hstrip->Fill(strip);
      hpadx->Fill(padx);
      hpadz->Fill(padz);
      h2tdcVSadc->Fill(adc,tdc);

    }
  
  } // end loop on events

  TFile *fout = new TFile("TOF_digitsQA.root","RECREATE");
  htdc->Write();
  hadc->Write();
  h2tdcVSadc->Write();
  hsector->Write();
  hplate->Write();
  hstrip->Write();
  hpadz->Write();
  hpadx->Write();
  fout->Close(); 


  return rc;
  
}

