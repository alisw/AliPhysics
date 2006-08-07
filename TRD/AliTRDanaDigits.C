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

  AliRunLoader* rl = AliRunLoader::Open(alifile);
  AliLoader* loader = rl->GetLoader("TRDLoader");
  rl->LoadDigits();
  
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  
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
  digitsManager->ReadDigits(loader->TreeD());

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

  // Loop through the detector pixel
  for (Int_t time = 0; time < timeMax; time++) {
    for (Int_t  col = 0;  col <  colMax;  col++) {
      for (Int_t  row = 0;  row <  rowMax;  row++) {

        digit = digitsManager->GetDigit(row,col,time,iDet);
        track = digitsManager->GetTrack(0,row,col,time,iDet);
        
        delete digit;

      }
    }
  }

  delete rl;

}
