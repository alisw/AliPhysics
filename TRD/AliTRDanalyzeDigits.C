Int_t AliTRDanalyzeDigits()
{
  //
  // Analyzes the digits
  //

  Int_t rc = 0;

  const Int_t kNpla = AliTRDgeometry::Nplan();

  if (!gAlice) {
    cout << "<AliTRDanalyzeDigits> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  Int_t nPart = gAlice->GetEvent(0);

  // Get the pointer to the TRD detector 
  AliTRD *trd = (AliTRD *) gAlice->GetDetector("TRD");
  if (!trd) {
    cout << "<AliTRDanalyzeDigits> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Get the digitizer object
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  AliTRDdigitizer *digitizer = (AliTRDdigitizer *) file->Get("digitizer");
  if (!digitizer) {
    cout << "<AliTRDanalyzeDigits> No digitizer object found" << endl;
    rc = 3;
    return rc;
  }

  // Define the histograms
  Int_t adcRange = ((Int_t) digitizer->GetADCoutRange());
  TH1F *hAmpAll   = new TH1F("hAmpAll"  ,"Amplitude of the digits (all)"
                            ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);
  TH1F *hAmpEl    = new TH1F("hAmpEl"   ,"Amplitude of the digits (electrons)"
                            ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);
  TH1F *hAmpPi    = new TH1F("hAmpPi"   ,"Amplitude of the digits (pions)"
                            ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);
  TH1F *hAmpNoise = new TH1F("hAmpNoise","Amplitude of the digits (noise)"
                            ,5,-0.5,4.5);

  // Get the pointer to the geometry object
  AliTRDgeometry *geo;
  if (trd) {
    geo = trd->GetGeometry();
  }
  else {
    cout << "<AliTRDanalyzeDigits> No TRD geometry found" << endl;
    rc = 4;
    return rc;
  }

  // Create the digits manager
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();

  digitsManager->Open("TRD_test.root");

  // Read the digits from the file
  if (!(digitsManager->ReadDigits())) {
    cout << "<AliTRDanalyzeDigits> Cannot read the digits" << endl;
    rc = 5;
    return rc;
  }

  AliTRDmatrix *matrix;

  Int_t countDigits = 0;
  Int_t iSec        = trd->GetSensSector();
  Int_t iCha        = trd->GetSensChamber();
  Int_t timeMax     = geo->GetTimeTotal();

  TProfile *hAmpTimeEl = new TProfile("hAmpTimeEl","Amplitude of the digits (electrons)"
				      ,timeMax,-0.5,((Double_t) timeMax)-0.5);
  TProfile *hAmpTimePi = new TProfile("hAmpTimePi","Amplitude of the digits (pions)"
				      ,timeMax,-0.5,((Double_t) timeMax)-0.5);

  // Loop over all planes
  for (Int_t iPla = 0; iPla < kNpla; iPla++) {

    Int_t iDet   = geo->GetDetector(iPla,iCha,iSec);
    Int_t rowMax = geo->GetRowMax(iPla,iCha,iSec);
    Int_t colMax = geo->GetColMax(iPla);
  
    if (iPla == 0) {
      matrix = new AliTRDmatrix(rowMax,colMax,timeMax,iSec,iCha,iPla);
    }

    // Loop through the detector pixel
    for (Int_t time = 0; time < timeMax; time++) {
      for (Int_t  col = 0;  col <  colMax;  col++) {
        for (Int_t  row = 0;  row <  rowMax;  row++) {

          AliTRDdigit *digit    = digitsManager->GetDigit(row,col,time,iDet);
          Int_t        amp      = digit->GetAmp();
          Int_t        track0   = digitsManager->GetTrack(0,row,col,time,iDet);
          Int_t        track1   = digitsManager->GetTrack(1,row,col,time,iDet);
          TParticle   *particle = 0;
          if (track0 > -1) {
            particle = gAlice->Particle(track0);
	  }

          if (amp > 0) {
            countDigits++;
            if (iPla == 0) {
              matrix->SetSignal(row,col,time,amp);
	    }
	  }

	  // Total spectrum
          hAmpAll->Fill(amp);

	  // Noise spectrum
          if (track0 < 0) {
            hAmpNoise->Fill(amp);
	  }          

	  // Electron digit
          if ((particle) && (particle->GetPdgCode() ==   11) && (track1 < 0)) {
            hAmpEl->Fill(amp);
            hAmpTimeEl->Fill(time,amp);
	  }

          // Pion digit
          if ((particle) && (particle->GetPdgCode() == -211) && (track1 < 0)) {
            hAmpPi->Fill(amp);
            hAmpTimePi->Fill(time,amp);
	  }

          delete digit;

        }
      }
    }

  }

  cout << "<AliTRDanalyzeDigits> Found " << countDigits << " digits in total" << endl;

  // Display the detector matrix
  matrix->Draw();
  matrix->ProjRow();
  matrix->ProjCol();
  matrix->ProjTime();

  TCanvas *cDigits = new TCanvas("cDigits","AliTRDanalyzeDigits",50,50,600,800);
  cDigits->Divide(2,3);
  cDigits->cd(1);
  gPad->SetLogy();
  hAmpAll->SetXTitle("Amplitude (ADC-channels)");
  hAmpAll->SetYTitle("Entries");
  hAmpAll->Draw();
  cDigits->cd(2);
  gPad->SetLogy();
  hAmpNoise->SetXTitle("Amplitude (ADC-channels)");
  hAmpNoise->SetYTitle("Entries");
  hAmpNoise->Draw();
  cDigits->cd(3);
  gPad->SetLogy();
  hAmpEl->SetXTitle("Amplitude (ADC-channels)");
  hAmpEl->SetYTitle("Entries");
  hAmpEl->Draw();
  cDigits->cd(4);
  gPad->SetLogy();
  hAmpPi->SetXTitle("Amplitude (ADC-channels)");
  hAmpPi->SetYTitle("Entries");
  hAmpPi->Draw();
  cDigits->cd(5);
  hAmpTimeEl->SetXTitle("Timebin number");
  hAmpTimeEl->SetYTitle("Mean amplitude");
  hAmpTimeEl->Draw("HIST");
  cDigits->cd(6);
  hAmpTimePi->SetXTitle("Timebin number");
  hAmpTimePi->SetYTitle("Mean amplitude");
  hAmpTimePi->Draw("HIST");

  TFile *fileOut = new TFile("digits_test.root","RECREATE");
  hAmpAll->Write();
  hAmpNoise->Write();
  hAmpEl->Write();
  hAmpPi->Write();
  hAmpTimeEl->Write();
  hAmpTimePi->Write();
  fileOut->Close();

  return rc;

}
