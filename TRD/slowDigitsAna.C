{

/////////////////////////////////////////////////////////////////////////
//
// Example macro for the analysis of the TRD digits and the use
// of the AliTRDmatrix class.
//
/////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }

  // Input file name
//Char_t *alifile = "galice_d_v1.root"; 
  Char_t *alifile = "galice_c_v1.root"; 

  // Event number
  Int_t   nEvent  = 0;

  // Define the objects
  AliTRDv1     *TRD;
  TClonesArray *TRDDigits;
  AliTRDdigit  *OneTRDDigit;

  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits
  TFile *gafl = (TFile*) gROOT->GetListOfFiles()->FindObject(alifile);
  if (!gafl) {
    cout << "Open the ALIROOT-file " << alifile << endl;
    gafl = new TFile(alifile);
  }
  else {
    cout << alifile << " is already open" << endl;
  }

  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*) gafl->Get("gAlice");
    if (gAlice)  
      cout << "AliRun object found on file" << endl;
    else
      gAlice = new AliRun("gAlice","Alice test program");
  }

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) break;
  
  // Get the pointer to the hit-tree
  TTree *DigitsTree = gAlice->TreeD();

  // Get the pointer to the detector classes
  TRD = (AliTRDv1*) gAlice->GetDetector("TRD");
  // Get the pointer to the hit container
  if (TRD) TRDDigits = TRD->Digits();

  DigitsTree->GetBranch("TRD")->SetAddress(&TRDDigits);

  // Define the detector matrix for one chamber (Sector 6, Chamber 3, Plane 1)
  const Int_t iSec = 6;
  const Int_t iCha = 3;
  const Int_t iPla = 1;
  Int_t  rowMax = TRD->GetRowMax(iPla,iCha,iSec);
  Int_t  colMax = TRD->GetColMax(iPla);
  Int_t timeMax = TRD->GetTimeMax();
  cout << " rowMax = "  << rowMax
       << " colMax = "  << colMax
       << " timeMax = " << timeMax << endl;
  AliTRDmatrix *TRDMatrix = new AliTRDmatrix(rowMax,colMax,timeMax,iSec,iCha,iPla);

  Int_t nEntries = DigitsTree->GetEntries();
  cout << "Number of entries in digits tree = " << nEntries << endl; 

  // Loop through all entries in the tree
  Int_t nbytes;
  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {

    cout << "iEntry = " << iEntry << endl;

    // Import the tree
    gAlice->ResetDigits();
    nbytes += DigitsTree->GetEvent(iEntry);

    // Get the number of digits in the detector 
    Int_t nTRDDigits = TRDDigits->GetEntriesFast();
    cout << " nTRDDigits = " << nTRDDigits << endl;    

    // Loop through all TRD digits
    for (Int_t iTRDDigits = 0; iTRDDigits < nTRDDigits; iTRDDigits++) {

      // Get the information for this digit
      OneTRDDigit = (AliTRDdigit*) TRDDigits->UncheckedAt(iTRDDigits);
      Int_t signal    = OneTRDDigit->fSignal;
      Int_t   sector  = OneTRDDigit->fSector;
      Int_t   chamber = OneTRDDigit->fChamber;
      Int_t   plane   = OneTRDDigit->fPlane;
      Int_t   row     = OneTRDDigit->fRow;
      Int_t   col     = OneTRDDigit->fCol;
      Int_t   time    = OneTRDDigit->fTime;

      // Fill the detector matrix
      if (signal > 1) {
        TRDMatrix->SetSignal(row,col,time,signal);
      }

    }

  }

  // Display the detector matrix
  TRDMatrix->Draw();
  TRDMatrix->DrawRow(18);
  TRDMatrix->DrawCol(58);
  TRDMatrix->DrawTime(20);

}
