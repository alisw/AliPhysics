Int_t AliTOFanalyzeSDigits(TString headersFile, Int_t iEvNum=0)
{
  //
  // Analyzes the TOF sdigits and fills QA-histograms 
  //
  // iEvNum=0 means all events in the file

  Int_t rc=0;

  TFile * file = (TFile*) gROOT->GetFile(headersFile.Data() ) ;

  //File was not opened yet

  if(file == 0){
    if(headersFile.Contains("rfio"))
      file =	TFile::Open(headersFile,"update") ;
    else
      file = new TFile(headersFile.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }


  if (iEvNum == 0) iEvNum = (Int_t) gAlice->TreeE()->GetEntries();


  AliTOFSDigit *tofsdigit;

  AliTOF * tof = (AliTOF *) gAlice->GetDetector("TOF") ;

  if (!tof) {
    cout << "<AliTOFanalyzeSDigits> No TOF detector found" << endl;
    rc = 2;
    return rc;
  }

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


  for (Int_t ievent = 0; ievent < iEvNum; ievent++) {

    gAlice->GetEvent(ievent) ;
    if(gAlice->TreeS()==0) {
      cout << "<AliTOFanalyzeSDigits> No  TreeS found" << endl;
      rc = 4;
      return rc;
    }
    
    
    
    Int_t ndig, k;
    gAlice->ResetDigits();
    gAlice->TreeS()->GetEvent(ievent);
    TClonesArray * TOFdigits   = tof->SDigits();
    
    ndig=TOFdigits->GetEntries();
    
    cout << "<AliTOFanalyzeSDigits> found " << ndig
	 << " TOF sdigits for event " << ievent << endl;
    
    for (k=0; k<ndig; k++) {
      tofsdigit= (AliTOFSDigit*) TOFdigits->UncheckedAt(k);
      Float_t firstTDC=tofsdigit->GetTdc(0);
      Float_t firstADC=tofsdigit->GetAdc(0);
      htdc->Fill(firstTDC);
      hadc->Fill(firstADC);
      // TOF sdigit volumes
      Int_t sector    = tofsdigit->GetSector(); // range [1-18]
      Int_t plate     = tofsdigit->GetPlate();  // range [1- 5]
      Int_t strip     = tofsdigit->GetStrip();  // range [1-20]
      Int_t padz      = tofsdigit->GetPadz();   // range [1- 2]
      Int_t padx      = tofsdigit->GetPadx();   // range [1-48]
      // it is QA, then I perform QA!
      Bool_t isSDigitBad = (sector<1 || sector>18 || plate<1 || plate >5 || padz<1 || padz>2 || padx<1 || padx>48);

      if (isSDigitBad) {
	cout << "<AliTOFanalyzeHits>  strange hit found" << endl;
	rc = 3;
	return rc;
      }

      // filling sdigit volume histos
      hsector->Fill(sector);
      hplate->Fill(plate);
      hstrip->Fill(strip);
      hpadx->Fill(padx);
      hpadz->Fill(padz);
      h2tdcVSadc->Fill(firstTDC,firstADC);

      //cout << "firstTDC " << firstTDC << " firstADC " << firstADC << endl;
    }
  
  }

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


  return rc;
  
}

