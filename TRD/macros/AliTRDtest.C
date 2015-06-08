
Int_t AliTRDtest() 
{
  //
  // Test macro for the TRD code
  //

  Int_t rc = 0;

  AliSimulation sim;
  sim.SetConfigFile("$(ALICE_ROOT)/TRD/Macros/AliTRDconfig.C");
  sim.SetLoadAlignFromCDB(0);
  sim.Run();

  // Analyze the TRD hits
  if (rc = AliTRDanalyzeHits()) return rc;

  // Analyze the digits
  //if (rc = AliTRDanalyzeDigits()) return rc;

  return rc;

}

//_____________________________________________________________________________
Int_t AliTRDanalyzeHits()
{
  //
  // Analyzes the hits and fills QA-histograms 
  //

  Int_t rc = 0;

  //AliRunLoader *rl = gAlice->GetRunLoader();
  AliRunLoader *rl = AliRunLoader::Open("TRD_test.root"
				       ,AliConfig::GetDefaultEventFolderName());
  if (!rl) {
    cout << "<AliTRDanalyzeHits> No RunLoader found" << endl;
    rc = 1;
    return rc;
  }
  
  AliLoader* loader = rl->GetLoader("TRDLoader");
  if (!loader) {
    cout << "<AliTRDanalyzeHits> No TRDLoader found" << endl;
    rc = 2;
    return rc;
  }

  rl->LoadgAlice();
  rl->LoadHeader();
  rl->LoadKinematics();
  rl->LoadHits();

  // Get the pointer to the TRD detector 
  gAlice = rl->GetAliRun();
  AliTRD *trd = (AliTRD *) gAlice->GetDetector("TRD");
  if (!trd) {
    cout << "<AliTRDanalyzeHits> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Get the pointer to the geometry object
  AliTRDgeometry *geo;
  if (trd) {
    geo = (AliTRDgeometry *) trd->GetGeometry();
  }
  else {
    cout << "<AliTRDanalyzeHits> No TRD geometry found" << endl;
    rc = 3;
    return rc;
  }

  AliTRDCommonParam *par = AliTRDCommonParam::Instance();

  // Define the histograms
  TH1F *hQdedx  = new TH1F("hQdedx","Charge dedx-hits",100,0.0,1000.0);
  TH1F *hQtr    = new TH1F("hQtr"  ,"Charge TR-hits"  ,100,0.0,1000.0);

  Float_t rmin   = geo->Rmin();
  Float_t rmax   = geo->Rmax();
  Float_t length = geo->GetChamberLength(0,2);
  Float_t width  = geo->GetChamberWidth(0);
  Int_t   ncol   = par->GetColMax(0);
  Int_t   nrow   = par->GetRowMax(0,2,13);
  Int_t   ntime  = ((Int_t) (rmax - rmin) / 22.0);

  TH2F *hZY     = new TH2F("hZY"   ,"Y vs Z (chamber 0)", nrow,-length/2.,length/2.
                                                        ,ntime,      rmin,     rmax);
  TH2F *hXZ     = new TH2F("hXZ"   ,"Z vs X (plane 0)"  , ncol, -width/2., width/2.
                                                        , nrow,-length/2.,length/2.);

  // Get the pointer hit tree
  TTree *hitTree = loader->TreeH();  
  if (!hitTree) {
    cout << "<AliTRDanalyzeHits> No hit tree found" << endl;
    rc = 4;
    return rc;
  }

  Int_t countHits = 0;
  Int_t nBytes    = 0;

  // Get the number of entries in the hit tree
  // (Number of primary particles creating a hit somewhere)
  Int_t nTrack = (Int_t) hitTree->GetEntries();
  cout << "<AliTRDanalyzeHits> Found " << nTrack 
       << " primary particles with hits" << endl;

  // Loop through all entries in the tree
  for (Int_t iTrack = 0; iTrack < nTrack; iTrack++) {

    gAlice->ResetHits();
    nBytes += hitTree->GetEvent(iTrack);

    // Loop through the TRD hits  
    Int_t iHit = 0;
    AliTRDhit *hit = (AliTRDhit *) trd->FirstHit(-1);
    while (hit) {

      countHits++;
      iHit++;

      Float_t x     = hit->X();
      Float_t y     = hit->Y();
      Float_t z     = hit->Z();
      Float_t q     = hit->GetCharge();
      Int_t   track = hit->Track();
      Int_t   det   = hit->GetDetector();
      Int_t   plane = geo->GetPlane(det);

      if      (q > 0) {
        hQdedx->Fill(q);
        hZY->Fill(z,y);
        if (plane == 0) {
          hXZ->Fill(x,z);
	}
      }
      else if (q < 0) {
        hQtr->Fill(TMath::Abs(q));
      }

      hit = (AliTRDhit *) trd->NextHit();         

    }

  }

  cout << "<AliTRDanalyzeHits> Found " << countHits << " hits in total" << endl;

  TCanvas *cHits = new TCanvas("cHits","AliTRDanalyzeHits",50,50,600,600);
  cHits->Divide(2,2);
  cHits->cd(1);
  hXZ->Draw("COL");
  cHits->cd(2);
  hZY->Draw("COL");
  cHits->cd(3);
  gPad->SetLogy();
  hQdedx->Draw();
  cHits->cd(4);
  gPad->SetLogy();
  hQtr->Draw();

  return rc;

}

//_____________________________________________________________________________
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

  AliRunLoader *rl = gAlice->GetRunLoader();
  if (!rl) {
    cout << "<AliTRDanalyzeHits> No RunLoader found" << endl;
    rc = 2;
    return rc;
  }

  // Import the Trees for the event nEvent in the file
  rl->LoadDigits();
  
  AliLoader* loader = rl->GetLoader("TRDLoader");
  if (!loader) {
    cout << "<AliTRDanalyzeHits> No TRDLoader found" << endl;
    rc = 3;
    return rc;
  }

  // Get the pointer to the TRD detector 
  AliTRD *trd = (AliTRD *) gAlice->GetDetector("TRD");
  if (!trd) {
    cout << "<AliTRDanalyzeDigits> No TRD detector found" << endl;
    rc = 4;
    return rc;
  }

  // Get the parameter object
  AliTRDSimParam    *parSim = AliTRDSimParam::Instance();
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  AliTRDcalibDB     *cal    = AliTRDcalibDB::Instance();

  // Define the histograms
  Int_t adcRange = ((Int_t) parSim->GetADCoutRange());
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
    rc = 5;
    return rc;
  }

  // Create the digits manager
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  digitsManager->SetDebug(1);

  // Read the digits from the file
  if (!(digitsManager->ReadDigits(loader->TreeD()))) {
    cout << "<AliTRDanalyzeDigits> Cannot read the digits" << endl;
    rc = 6;
    return rc;
  }

  // Get the particle stack
  AliStack *kineStack = rl->Stack();
  if (!kineStack) {
    cout << "<AliTRDanalyzeDigits> Cannot find the KINE stack" << endl;
    rc = 7;
    return rc;
  }

  Int_t countDigits = 0;
  Int_t iSec = 0;
  Int_t iCha = 2;
  Int_t timeMax     = cal->GetNumberOfTimeBins();

  TProfile *hAmpTimeEl = new TProfile("hAmpTimeEl","Amplitude of the digits (electrons)"
				      ,timeMax,-0.5,((Double_t) timeMax)-0.5);
  TProfile *hAmpTimePi = new TProfile("hAmpTimePi","Amplitude of the digits (pions)"
				      ,timeMax,-0.5,((Double_t) timeMax)-0.5);

  // Loop over all planes
  for (Int_t iPla = 0; iPla < kNpla; iPla++) {

    Int_t iDet   = geo->GetDetector(iPla,iCha,iSec);
    Int_t rowMax = parCom->GetRowMax(iPla,iCha,iSec);
    Int_t colMax = parCom->GetColMax(iPla);
  
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
            particle = (TParticle *) kineStack->Particle(track0);
	  }

          if (amp > 0) {
            countDigits++;
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

  TCanvas *cDigits = new TCanvas("cDigits","AliTRDanalyzeDigits",100,100,600,800);
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

  return rc;

}

