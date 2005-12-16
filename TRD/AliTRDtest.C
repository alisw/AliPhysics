
Int_t AliTRDtest() 
{
  //
  // Test macro for the TRD code
  //

  Int_t rc = 0;

  TFile *file;

  // Initialize the test setup 
  gAlice->Init("$(ALICE_ROOT)/TRD/AliTRDconfig.C");

  // Run one event and create the hits
  gAlice->Run(1);

  //if (gAlice) delete gAlice;
  //file   = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  //gAlice = (AliRun *) file->Get("gAlice");

  // Analyze the TRD hits
  if (rc = AliTRDanalyzeHits()) return rc;

  // Run the digitization
  if (rc = AliTRDcreateDigits()) return rc;

//    if (gAlice) delete gAlice;
//    file   = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
//    gAlice = (AliRun *) file->Get("gAlice");

  // Analyze the digits
  if (rc = AliTRDanalyzeDigits()) return rc;

//    // Create the cluster
//    if (rc = AliTRDcreateCluster()) return rc;

//    if (gAlice) delete gAlice;
//    file   = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
//    gAlice = (AliRun *) file->Get("gAlice");

//    // Analyze the cluster
//    if (rc = AliTRDanalyzeCluster()) return rc;

  //file   = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  //file->Close();

  // Find the tracks
  //if (rc = AliTRDcreateTracks()) return rc;

  return rc;

}

//_____________________________________________________________________________
Int_t AliTRDanalyzeHits()
{
  //
  // Analyzes the hits and fills QA-histograms 
  //

  Int_t rc = 0;

  AliRunLoader *rl = gAlice->GetRunLoader();
  if (!rl) {
    cout << "<AliTRDanalyzeHits> No RunLoader found" << endl;
    rc = 1;
    return rc;
  }

  // Import the Trees for the event nEvent in the file
  rl->GetEvent(0);
  
  AliLoader* loader = rl->GetLoader("TRDLoader");
  if (!loader) {
    cout << "<AliTRDanalyzeHits> No TRDLoader found" << endl;
    rc = 2;
    return rc;
  }

  // Get the pointer to the TRD detector 
  AliTRD *trd = (AliTRD *) gAlice->GetDetector("TRD");
  if (!trd) {
    cout << "<AliTRDanalyzeHits> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Get the pointer to the geometry object
  AliTRDgeometryFull *geo;
  if (trd) {
    geo = (AliTRDgeometryFull *) trd->GetGeometry();
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
Int_t AliTRDcreateDigits()
{
  //
  // Creates the digits from the hits of the slow simulator
  //

  Int_t rc = 0;

  if (!gAlice) {
    cout << "<AliTRDcreateDigits> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  gAlice->GetEvent(0);

  // Create the TRD digitzer 
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("TRDdigitizer"
                                                  ,"TRD digitizer class");
  digitizer->Open("TRD_test.root");
  digitizer->InitDetector();
  digitizer->InitOutput(0);

  // Set the parameter
  digitizer->SetDebug(1);

  // Define the parameter object
  // If no external parameter object is defined, 
  // default parameter will be used
  AliTRDparameter *parameter = new AliTRDparameter("TRDparameter"
						  ,"TRD parameter class");
  digitizer->SetParameter(parameter);

  // Create the digits
  if (!(digitizer->MakeDigits())) {
    rc = 2;
    return rc;
  }

  // Write the digits into the input file
  //if (!(digitizer->MakeBranch())) {
  //  rc = 3;
  //  return rc;
  //}

  // Write the digits into the input file
  if (!(digitizer->WriteDigits())) {
    rc = 4;
    return rc;
  }

  // Save the parameter object in the AliROOT file
  if (!(parameter->Write())) {
    rc = 4;
    return rc;
  }

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
  rl->LoadKinematics();
  rl->GetEvent(0);
  rl->GetHeader();
  
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

  AliTRDmatrix *matrix;

  Int_t countDigits = 0;
  Int_t iSec        = trd->GetSensSector();
  Int_t iCha        = trd->GetSensChamber();
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
            particle = (TParticle *) kineStack->Particle(track0);
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

  return rc;

}

//_____________________________________________________________________________
Int_t AliTRDcreateCluster()
{
  //
  // Creates the cluster from the digits
  //

  Int_t rc = 0;

  // Create the clusterizer
  AliTRDclusterizerV1 *clusterizer =  new AliTRDclusterizerV1("TRDclusterizer"
                                                             ,"TRD clusterizer class");
  clusterizer->SetVerbose(1);

  // Define the parameter object
  // If no external parameter object is defined, 
  // default parameter will be used
  AliTRDparameter *parameter = new AliTRDparameter("TRDparameter"
						  ,"TRD parameter class");
  parameter->SetClusMaxThresh(0);
  parameter->SetClusSigThresh(0);
  clusterizer->SetParameter(parameter);
 
  // Open the file
  if (!(clusterizer->Open("TRD_test.root",0))) {
    rc = 1;
    return rc;
  }    

  // Load the digits
  if (!(clusterizer->ReadDigits())) {
    rc = 2;
    return rc;
  }    

  // Find the cluster
  if (!(clusterizer->MakeClusters())) {
    rc = 3;
    return rc;
  }

  // Write the cluster tree into the file 
  if (!(clusterizer->WriteClusters(-1))) {
    rc = 4;
    return rc;
  }

  // Save the clusterizer class in the file
  if (!(clusterizer->Write())) {
    rc = 5;
    return rc;
  }

  return rc;

}

//_____________________________________________________________________________
Int_t AliTRDanalyzeCluster()
{
  //
  // Analyzes the cluster
  //

  Int_t rc = 0;

  if (!gAlice) {
    cout << "<AliTRDanalyzeCluster> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  gAlice->GetEvent(0);

  // Get the pointer to the TRD detector 
  AliTRD *trd = (AliTRD *) gAlice->GetDetector("TRD");
  if (!trd) {
    cout << "<AliTRDanalyzeCluster> No TRD detector found" << endl;
    rc = 2;
    return rc;
  }

  // Define the histograms
  TH1F *hClusAll   = new TH1F("hClusAll"  ,"Amplitude of the cluster (all)"     
                                          ,501,-0.5,500.5);
  TH1F *hClusNoise = new TH1F("hClusNoise","Amplitude of the cluster (noise)"   
                                          , 11,-0.5, 10.5);
  TH1F *hClusEl    = new TH1F("hClusEl"   ,"Amplitude of the cluster (electron)"
                                          ,501,-0.5,500.5);
  TH1F *hClusPi    = new TH1F("hClusPi"   ,"Amplitude of the cluster (pion)"    
                                          ,501,-0.5,500.5);

  // Get the pointer to the geometry object
  AliTRDgeometry *geo;
  if (trd) {
    geo = trd->GetGeometry();
  }
  else {
    cout << "<AliTRDanalyzeCluster> No TRD geometry found" << endl;
    rc = 3;
    return rc;
  }

  // Get the pointer to the hit-tree
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  TTree *clusterTree = (TTree *) file->Get("TreeR0_TRD");
  if (!(clusterTree)) {
    cout << "<AliTRDanalyzeCluster> No tree with clusters found" << endl;
    rc = 4;
    return rc;
  }

  // Get the pointer to the hit container
  TObjArray *clusterArray = trd->RecPoints();
  if (!(clusterArray)) {
    cout << "<AliTRDanalyzeCluster> No clusterArray found" << endl;
    rc = 5;
    return rc;
  }

  // Set the branch address
  clusterTree->GetBranch("TRDcluster")->SetAddress(&clusterArray);
  Int_t nEntries = clusterTree->GetEntries();
  cout << "<AliTRDanalyzeCluster> Number of entries in the cluster tree = " 
       << nEntries 
       << endl;

  Int_t countCluster = 0;
  Int_t countOverlap = 0;

  // Loop through all entries in the tree
  Int_t nbytes;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

    // Import the tree
    nbytes += clusterTree->GetEvent(iEntry);

    // Get the number of points in the detector 
    Int_t nCluster = clusterArray->GetEntriesFast();

    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) {

      // Get the information for this digit
      AliTRDcluster *cluster = (AliTRDcluster *) clusterArray->UncheckedAt(iCluster);
      Int_t    detector = cluster->GetDetector();      
      Int_t    sector   = geo->GetSector(detector);
      Int_t    plane    = geo->GetPlane(detector);
      Int_t    chamber  = geo->GetChamber(detector);
      Float_t  energy   = cluster->GetQ();
      Int_t    track0   = cluster->GetLabel(0);
      Int_t    track1   = cluster->GetLabel(1);
      Int_t    track2   = cluster->GetLabel(2);
      TParticle *particle = 0;
      if (track0 > -1) {
        particle = gAlice->Particle(track0);
      }

      countCluster++;
      if (!cluster->Isolated()) countOverlap++;

      // Total spectrum
      hClusAll->Fill(energy);

      if (cluster->Isolated()) {

        // Noise spectrum
        if (track0 < 0) {
          hClusNoise->Fill(energy);
        }          

        // Electron cluster
        if ((particle) && (particle->GetPdgCode() ==   11) && (track1 < 0)) {
          hClusEl->Fill(energy);
        }

        // Pion cluster
        if ((particle) && (particle->GetPdgCode() == -211) && (track1 < 0)) {
          hClusPi->Fill(energy);
        }

      }

    }

  }

  cout << "<AliTRDanalyzeCluster> Found " << countCluster << " cluster in total"    << endl;
  cout << "<AliTRDanalyzeCluster> Found " << countOverlap << " overlapping cluster" << endl;
  cout << endl;

  TCanvas *cCluster = new TCanvas("cCluster","AliTRDanalyzeCluster",50,50,600,600);
  cCluster->Divide(2,2);

  TF1 *fun;
  cCluster->cd(1);
  gPad->SetLogy();
  hClusAll->Fit("landau","0");
  fun = (TF1 *) hClusAll->GetListOfFunctions()->First();
  Float_t meanAll = fun->GetParameter(1);
  hClusAll->Draw();
  fun->SetLineColor(2);
  fun->Draw("SAME");
  
  cCluster->cd(2);
  gPad->SetLogy();
  Float_t meanNoise = hClusNoise->GetMean();
  hClusNoise->Draw();

  cCluster->cd(3);
  gPad->SetLogy();
  hClusEl->Fit("landau","0");
  fun = (TF1 *) hClusEl->GetListOfFunctions()->First();
  fun->SetLineColor(2);
  Float_t meanEl = fun->GetParameter(1);
  hClusEl->Draw();
  fun->Draw("SAME");

  cCluster->cd(4);
  gPad->SetLogy();
  hClusPi->Fit("landau","0");
  fun = (TF1 *) hClusPi->GetListOfFunctions()->First();
  fun->SetLineColor(2);
  Float_t meanPi = fun->GetParameter(1);
  hClusPi->Draw();
  fun->Draw("SAME");

  cout << endl;
  cout << "##################################################################" << endl;
  cout << "    Mean all       = " << meanAll   << endl;
  cout << "    Mean noise     = " << meanNoise << endl;
  cout << "    Mean electrons = " << meanEl    << endl;
  cout << "    Mean pions     = " << meanPi    << endl;
  cout << "##################################################################" << endl;
  cout << endl;

  return rc;

}

//_____________________________________________________________________________
Int_t AliTRDcreateTracks()
{
  //
  // Creates the tracks
  //

  Int_t rc = 0;

  // Create the tracker
  AliTRDtracker *tracker = new AliTRDtracker("TRDtracker","TRD tracker");

  // Read in the kine tree and the cluster
  tracker->GetEvent("TRD_test.root","TRD_test.root");
 
  // Find the tracks
  TH1F *hs = new TH1F("hs","hs",100,0.0,1.0);
  TH1F *hd = new TH1F("hd","hd",100,0.0,1.0);
  tracker->Clusters2Tracks(hs,hd);
 
  // Store the tracks
  tracker->WriteTracks("TRD_test.root");

  return rc;                                                                    

}
