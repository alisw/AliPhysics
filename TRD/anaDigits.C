void anaDigits() 
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
  Char_t *alifile = "galice.root"; 

  // Event number
  Int_t   nEvent  = 0;

  // Define the objects
  AliTRDv1       *trd;
  AliTRDgeometry *geo;
  AliTRDdigit    *digit;

  Int_t           track;

  // Connect the AliRoot file containing Geometry, Kine, Hits, and digits
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
  
  // Get the pointer to the detector object
  trd = (AliTRDv1*) gAlice->GetDetector("TRD");

  // Get the pointer to the geometry object
  if (trd) {
    geo = trd->GetGeometry();
  }
  else {
    cout << "Cannot find the geometry" << endl;
    break;
  }

  // Create the digits manager
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  digitsManager->SetVerbose(1);

  // Read the digits from the file
  digitsManager->Open(alifile);
  digitsManager->ReadDigits();

  // Define the detector matrix for one chamber
  const Int_t iSec = 11;
  const Int_t iCha = 2;
  const Int_t iPla = 0;
  Int_t  rowMax = geo->GetRowMax(iPla,iCha,iSec);
  Int_t  colMax = geo->GetColMax(iPla);
  Int_t timeMax = geo->GetTimeMax();
  cout << "Geometry: rowMax = "  <<  rowMax
                << " colMax = "  <<  colMax
                << " timeMax = " << timeMax << endl;
  AliTRDmatrix *matrix = new AliTRDmatrix(rowMax,colMax,timeMax,iSec,iCha,iPla);

  // Get the detector number
  Int_t iDet = geo->GetDetector(iPla,iCha,iSec); 
  cout << " iDet = " << iDet << endl;

  // Loop through the detector pixel
  for (Int_t time = 0; time < timeMax; time++) {
    for (Int_t  col = 0;  col <  colMax;  col++) {
      for (Int_t  row = 0;  row <  rowMax;  row++) {

        digit = digitsManager->GetDigit(row,col,time,iDet);
        track = digitsManager->GetTrack(0,row,col,time,iDet);
        
        matrix->SetSignal(row,col,time,digit->GetAmp());

        delete digit;

      }
    }
  }

  // Display the detector matrix
  matrix->Draw();
  //matrix->DrawRow(18);
  //matrix->DrawCol(58);
  //matrix->DrawTime(20);
  matrix->ProjRow();
  matrix->ProjCol();
  matrix->ProjTime();

}
