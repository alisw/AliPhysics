void slowDigitsAna() {

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
  Char_t *alifile = "galice_d_v1.root"; 

  // Event number
  Int_t   nEvent  = 0;

  // Define the objects
  AliTRDv1       *TRD;
  AliTRDgeometry *TRDgeometry;

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
  gAlice = (AliRun*) gafl->Get("gAlice");
  if (gAlice)  
    cout << "AliRun object found on file" << endl;
  else
    gAlice = new AliRun("gAlice","Alice test program");

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(nEvent);
  if (nparticles <= 0) break;
  
  // Get the pointer to the detector classes
  TRD = (AliTRDv1*) gAlice->GetDetector("TRD");
  // Get the pointer to the digits container and the geometry
  if (TRD) {
    TRDgeometry = TRD->GetGeometry();
  }
  else {
    cout << "Cannot find the geometry" << endl;
    break;
  }

  // Define the segment array for the digits
  AliTRDsegmentArray *DigitsArray = new AliTRDsegmentArray(540);

  // Load the digits from the tree
  DigitsArray->LoadArray("TRDdigits");

  // Define the detector matrix for one chamber
  const Int_t iSec = 13;
  const Int_t iCha = 3;
  const Int_t iPla = 3;
  Int_t  rowMax = TRDgeometry->GetRowMax(iPla,iCha,iSec);
  Int_t  colMax = TRDgeometry->GetColMax(iPla);
  Int_t timeMax = TRDgeometry->GetTimeMax();
  cout << "Geometry: rowMax = "  <<  rowMax
                << " colMax = "  <<  colMax
                << " timeMax = " << timeMax << endl;
  AliTRDmatrix *TRDmatrix = new AliTRDmatrix(rowMax,colMax,timeMax,iSec,iCha,iPla);

  // Get the digits for this detector
  Int_t iDet = TRDgeometry->GetDetector(iPla,iCha,iSec); 
  AliTRDdataArray *Digits = (AliTRDdataArray *) DigitsArray->At(iDet);
  Digits->Dump();
  // Expand the digits array
  //Digits->Expand();

  //Float_t signal = Digits->GetData(0,0,29);

  // Loop through the detector pixel
  for (Int_t time = 0; time < timeMax; time++) {
    for (Int_t  col = 0;  col <  colMax;  col++) {
      for (Int_t  row = 0;  row <  rowMax;  row++) {

        Float_t signal = Digits->GetData(row,col,time);
        if (signal != 0) {
	}
        TRDmatrix->SetSignal(row,col,time,signal);

      }
    }
  }

  // Display the detector matrix
  TRDmatrix->Draw();
  TRDmatrix->DrawRow(18);
  TRDmatrix->DrawCol(58);
  TRDmatrix->DrawTime(20);

}
