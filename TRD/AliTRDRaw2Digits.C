#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1.h>
#include <TCanvas.h>
#include "AliTRDrawData.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdigit.h"
#include "AliTRDgeometryFull.h"
#include "AliTRDparameter.h"
#include "AliTRDmatrix.h"
#include "AliRawReaderFile.h"
#include "AliRunLoader.h"
#endif

void AliTRDRaw2Digits(Int_t iEvent = 0, Int_t iDet = 0)
{

  AliTRDrawData *raw = new AliTRDrawData();
  raw->SetDebug(1);
  AliRawReaderFile rawReader(iEvent);
  AliTRDdigitsManager *digitsManagerRaw = raw->Raw2Digits(&rawReader);

  // The geometry object
  AliTRDgeometryFull *geo = new AliTRDgeometryFull();

  // The parameter object
  AliTRDparameter    *par = new AliTRDparameter("TRDparameter"
                                               ,"TRD parameter class");

  // Print the event and detector number
  cout << " iEvent = " << iEvent << endl;
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

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  AliLoader* loader = rl->GetLoader("TRDLoader");
  rl->GetEvent(iEvent);
  loader->LoadDigits();
  
  // Create the digits manager
  AliTRDdigitsManager *digitsManager = new AliTRDdigitsManager();
  digitsManager->SetDebug(1);

  // Read the digits from the file
  digitsManager->ReadDigits(loader->TreeD());

  TH1F *DigDiff = new TH1F("DigDiff","DigDiff",100,-10,+10);
  Int_t DigAmpRaw, DigAmp;

  // Loop through the detector pixel
  for (Int_t time = 0; time < timeMax; time++) {
    for (Int_t  col = 0;  col <  colMax;  col++) {
      for (Int_t  row = 0;  row <  rowMax;  row++) {

        AliTRDdigit* digit = digitsManagerRaw->GetDigit(row,col,time,iDet);
        
        matrix->SetSignal(row,col,time,digit->GetAmp());

	DigAmpRaw = digit->GetAmp();

        digit = digitsManager->GetDigit(row,col,time,iDet);

	DigAmp = digit->GetAmp();

	DigDiff->Fill((Float_t)DigAmp-(Float_t)DigAmpRaw);

        delete digit;

      }
    }
  }

  // Display the detector matrix
  matrix->Draw();
  matrix->ProjRow();
  matrix->ProjCol();
  matrix->ProjTime();

  TCanvas *c1 = new TCanvas("c1","Canvas 1",10,10,600,500);
  DigDiff->Draw();

}


