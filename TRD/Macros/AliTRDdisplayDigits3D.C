//_____________________________________________________________________________
Int_t AliTRDdisplayDigits3D(Int_t event = 0, Int_t thresh = 4
                          , Bool_t sdigits = kFALSE) 
{
  //  
  //  TRD digits display
  //
  //  Input parameter:
  //    <event>   : Event number 
  //    <thresh>  : Threshold to suppress the noise
  //    <sdigits> : If kTRUE it will display summable digits, normal digits otherwise.
  //                The signal event is displayed in yellow.
  //

  Char_t *inputFile = "galice.root";

  // Define the objects
  AliTRDv1       *trd;
  AliTRDgeometry *geo;

  TString       evfoldname = AliConfig::GetDefaultEventFolderName();
  AliRunLoader *runLoader  = AliRunLoader::GetRunLoader(evfoldname);
  if (!runLoader) {
    runLoader = AliRunLoader::Open(inputFile
                                  ,AliConfig::GetDefaultEventFolderName()
                                  ,"UPDATE");
  }
  if (!runLoader) {
    printf("Can not open session for file %s.",inputFile);
    return kFALSE;
  }
   
  if (!runLoader->GetAliRun()) {
    runLoader->LoadgAlice();
  }
  gAlice = runLoader->GetAliRun();  
  if (!gAlice) {
    printf("Could not find AliRun object.\n");
    return kFALSE;
  }

  runLoader->GetEvent(event);
  
  AliLoader *loader = runLoader->GetLoader("TRDLoader");
  if (!loader) {
    printf("Can not get TRD loader from Run Loader");
  }
  loader->LoadDigits();
  
  // Get the pointer to the detector object
  trd = (AliTRDv1*) gAlice->GetDetector("TRD");

  // Get the pointer to the geometry object
  if (trd) {
    geo = trd->GetGeometry();
  }
  else {
    printf("Cannot find the geometry\n");
    return 1;
  }

  AliCDBManager *cdbManager  = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliTRDcalibDB *calibration = AliTRDcalibDB::Instance();
  calibration->SetRun(0);

  TCanvas *c1 = new TCanvas("digits","TRD digits display",0,0,700,730);
  TView   *v  = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  v->SetParallel();
  c1->Clear();
  c1->SetFillColor(1);
  c1->SetTheta(90.0);
  c1->SetPhi(0.0);

  Int_t markerColorSignal = 2;
  Int_t markerColorBgnd   = 7;
  Int_t markerColorMerged = 5;
  Int_t mask              = 10000000;

  // Create the digits manager
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  digitsManager->SetSDigits(sdigits);

  // Read the digits from the file
  if (sdigits) {
    digitsManager->ReadDigits(loader->TreeS());
  }
  else {
    if (!loader->TreeD()) {
      printf("mist\n");
      return kFALSE;
    }
    digitsManager->ReadDigits(loader->TreeD());
  }

  Int_t totalsignal = 0;
  Int_t totalbgnd   = 0;
  Int_t totalmerged = 0;
  Int_t timeMax     = calibration->GetNumberOfTimeBins();

  // Loop through all detectors
  for (Int_t idet = 0; idet < geo->Ndet(); idet++) {

    printf("<AliTRDdisplayDigits3D> Loading detector %d\n",idet);
    AliTRDdataArrayI *digits  = digitsManager->GetDigits(idet);
    digits->Expand();

    Int_t isec    = geo->GetSector(idet);
    Int_t icha    = geo->GetChamber(idet);
    Int_t ipla    = geo->GetPlane(idet);
    AliTRDpadPlane *padPlane = new AliTRDpadPlane(ipla,icha);
    Int_t  rowMax = padPlane->GetNrows();
    Int_t  colMax = padPlane->GetNcols();

    Int_t ndigits = digits->GetOverThreshold(thresh);

    if (ndigits > 0) {

      TPolyMarker3D *pmSignal = new TPolyMarker3D(ndigits);
      Int_t isignal = 0;

      for (Int_t time = 0; time < timeMax; time++) {
        for (Int_t  col = 0;  col <  colMax;  col++) {
          for (Int_t  row = 0;  row <  rowMax;  row++) {

            Int_t amp = digits->GetDataUnchecked(row,col,time);

            if (amp > thresh) {
          
              Double_t glb[3];
              Double_t loc[3];

              loc[0] = row;
              loc[1] = col;
              loc[2] = time;
              geo->Local2Global(idet,loc,glb);
              Double_t x = glb[0];
              Double_t y = glb[1];
              Double_t z = glb[2];

              pmSignal->SetPoint(isignal,x,y,z);
              isignal++;
              totalsignal++;

	    }

	  }
        }
      }

      digits->Compress(1,0);

      pmSignal->SetMarkerSize(1); 
      pmSignal->SetMarkerColor(markerColorSignal);
      pmSignal->SetMarkerStyle(1);
      pmSignal->Draw();
   
    }

  }

  delete padPlane;

  TGeometry *geoAlice = gAlice->GetGeometry();
  TNode     *main     = (TNode *) ((geoAlice->GetListOfNodes())->First());
  TIter      next(main->GetListOfNodes());
  TNode     *module   = 0;
  while ((module = (TNode *) next())) {
    Char_t ch[100];
    sprintf(ch,"%s\n",module->GetTitle());
    if ((ch[0] == 'T') && ((ch[1] == 'R') || (ch[1] == 'P'))) {
      module->SetVisibility( 3);
    }
    else {
      module->SetVisibility(-1);
    }
  }
      
  geoAlice->Draw("same");

  c1->Modified(); 
  c1->Update(); 

  return 0;

}

