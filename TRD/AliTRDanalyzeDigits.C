Int_t AliTRDanalyzeDigits()
{
  //
  // Analyzes the digits
  //

  Int_t rc = 0;

  if (!gAlice) {
    cout << "<AliTRDanalyzeDigits> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  Int_t nPart = gAlice->GetEvent(0);

  // Get the pointer to the TRD detector 
  AliTRD *TRD = (AliTRD *) gAlice->GetDetector("TRD");
  if (!TRD) {
    cout << "<AliTRDanalyzeDigits> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Get the digitizer object
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  AliTRDdigitizer *Digitizer = (AliTRDdigitizer *) file->Get("digitizer");
  if (!Digitizer) {
    cout << "<AliTRDanalyzeDigits> No digitizer object found" << endl;
    rc = 3;
    return rc;
  }

  // Define the histograms
  Int_t adcRange = ((Int_t) Digitizer->GetADCoutRange());
  TH1F *hAmpAll   = new TH1F("hAmpAll"  ,"Amplitude of the digits (all)"
                            ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);
  TH1F *hAmpEl    = new TH1F("hAmpEl"   ,"Amplitude of the digits (electrons)"
                            ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);
  TH1F *hAmpPi    = new TH1F("hAmpPi"   ,"Amplitude of the digits (pions)"
                            ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);
  TH1F *hAmpNoise = new TH1F("hAmpNoise","Amplitude of the digits (noise)"
                            ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);

  // Get the pointer to the geometry object
  AliTRDgeometry *TRDgeometry;
  if (TRD) {
    TRDgeometry = TRD->GetGeometry();
  }
  else {
    cout << "<AliTRDanalyzeDigits> No TRD geometry found" << endl;
    rc = 4;
    return rc;
  }

  // Create the digits manager
  AliTRDdigitsManager *DigitsManager = new AliTRDdigitsManager();

  // Read the digits from the file
  if (!(DigitsManager->ReadDigits())) {
    cout << "<AliTRDanalyzeDigits> Cannot read the digits" << endl;
    rc = 5;
    return rc;
  }
  
  // Define the detector matrix for one chamber
  Int_t iSec = TRD->GetSensSector();
  Int_t iCha = TRD->GetSensChamber();
  Int_t iPla = 0;
  Int_t  rowMax = TRDgeometry->GetRowMax(iPla,iCha,iSec);
  Int_t  colMax = TRDgeometry->GetColMax(iPla);
  Int_t timeMax = TRDgeometry->GetTimeMax();
  cout << "<AliTRDanalyzeDigits> Geometry: rowMax = "  <<  rowMax
                                      << " colMax = "  <<  colMax
                                      << " timeMax = " << timeMax << endl;
  AliTRDmatrix *TRDmatrix = new AliTRDmatrix(rowMax,colMax,timeMax,iSec,iCha,iPla);
  // Get the detector number
  Int_t iDet = TRDgeometry->GetDetector(iPla,iCha,iSec);
  cout << "<AliTRDanalyzeDigits> iSec = " << iSec
                            << " iCha = " << iCha 
                            << " iPla = " << iPla 
                            << " iDet = " << iDet << endl;

  // Loop through the detector pixel
  Int_t countDigits = 0;
  for (Int_t time = 0; time < timeMax; time++) {
    for (Int_t  col = 0;  col <  colMax;  col++) {
      for (Int_t  row = 0;  row <  rowMax;  row++) {

        AliTRDdigit *Digit = DigitsManager->GetDigit(row,col,time,iDet);
        Int_t amp = Digit->GetAmp();

        Int_t track0 = DigitsManager->GetTrack(0,row,col,time,iDet);
        Int_t track1 = DigitsManager->GetTrack(1,row,col,time,iDet);
        TParticle *Part = 0;
        if (track0 > -1) {
          Part = gAlice->Particle(track0);
	}

        if (amp > 0) {

          countDigits++;
          TRDmatrix->SetSignal(row,col,time,amp);

	}

	// Total spectrum
        hAmpAll->Fill(amp);

	// Noise spectrum
        if (track0 < 0) {
          hAmpNoise->Fill(amp);
	}          

	// Electron digit
        if ((Part) && (Part->GetPdgCode() ==   11) && (track1 < 0)) {
          hAmpEl->Fill(amp);
	}

        // Pion digit
        if ((Part) && (Part->GetPdgCode() == -211) && (track1 < 0)) {
          hAmpPi->Fill(amp);
	}

        delete Digit;

      }
    }
  }

  cout << "<AliTRDanalyzeDigits> Found " << countDigits << " digits in total" << endl;

  // Display the detector matrix
  TRDmatrix->Draw();
  TRDmatrix->ProjRow();
  TRDmatrix->ProjCol();
  TRDmatrix->ProjTime();

  TCanvas *cDigits = new TCanvas("cDigits","AliTRDanalyzeDigits",50,50,600,600);
  cDigits->Divide(2,2);
  cDigits->cd(1);
  gPad->SetLogy();
  hAmpAll->Draw();
  cDigits->cd(2);
  gPad->SetLogy();
  gPad->SetLogx();
  hAmpNoise->Draw();
  cDigits->cd(3);
  gPad->SetLogy();
  hAmpEl->Draw();
  cDigits->cd(4);
  gPad->SetLogy();
  hAmpPi->Draw();

  return rc;

}
