void AliTRDanaDigits() 
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
  AliTRDv1        *trd;
  AliTRDgeometry  *geo;
  AliTRDdigit     *digit;
  AliTRDparameter *par;

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

  // The parameter object
  par = new AliTRDparameter("TRDparameter","TRD parameter class");

  // Create the digits manager
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  digitsManager->SetDebug(1);

  // Read the digits from the file
  digitsManager->Open(alifile);
  digitsManager->ReadDigits();

  // Get the detector number
  Int_t iDet = 514;
  cout << " iDet = " << iDet << endl;

  // Define the detector matrix for one chamber
  const Int_t iSec = geo->GetSector(iDet);
  const Int_t iCha = geo->GetChamber(iDet);
  const Int_t iPla = geo->GetPlane(iDet);
  Int_t  rowMax = par->GetRowMax(iPla,iCha,iSec);
  Int_t  colMax = par->GetColMax(iPla);
  Int_t timeMax = par->GetTimeMax();
  cout << "Geometry: rowMax = "  <<  rowMax
                << " colMax = "  <<  colMax
                << " timeMax = " << timeMax << endl;
  AliTRDmatrix *matrix = new AliTRDmatrix(rowMax,colMax,timeMax,iSec,iCha,iPla);

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
  matrix->ProjRow();
  matrix->ProjCol();
  matrix->ProjTime();

}
