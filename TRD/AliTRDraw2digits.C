void AliTRDraw2digits()
{

  AliTRDrawData *raw = new AliTRDrawData();
  raw->SetDebug(3);
  raw->Raw2Digit();

  AliTRDdigitsManager *digitsManager = raw->GetDigitsManager();

  // The geometry object
  AliTRDgeometryFull *geo = new AliTRDgeometryFull();

  // The parameter object
  AliTRDparameter    *par = new AliTRDparameter("TRDparameter"
                                               ,"TRD parameter class");

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


