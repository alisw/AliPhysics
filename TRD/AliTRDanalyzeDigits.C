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
  gAlice->GetEvent(0);

  // Get the pointer to the TRD detector 
  AliTRD *TRD = (AliTRD *) gAlice->GetDetector("TRD");
  if (!TRD) {
    cout << "<AliTRDanalyzeDigits> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Define the histograms
  TH1F *hAmp = new TH1F("hAmp","Amplitude of the digits",256,-0.5,255.5);

  // Get the pointer to the geometry object
  AliTRDgeometry *TRDgeometry;
  if (TRD) {
    TRDgeometry = TRD->GetGeometry();
  }
  else {
    cout << "<AliTRDanalyzeDigits> No TRD geometry found" << endl;
    rc = 3;
    return rc;
  }

  // Create the digits manager
  AliTRDdigitsManager *DigitsManager = new AliTRDdigitsManager();

  // Read the digits from the file
  if (!(DigitsManager->ReadDigits())) {
    cout << "<AliTRDanalyzeDigits> Cannot read the digits" << endl;
    rc = 4;
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

        if (amp > 0) {
          countDigits++;
          hAmp->Fill(amp);
          TRDmatrix->SetSignal(row,col,time,amp);
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
  gPad->SetLogy();
  hAmp->Draw();

  return rc;

}
