#if !defined( __CINT__) || defined(__MAKECINT__)

#include <TFile.h>
#include <TProfile2D.h>
#include <TTree.h>
#include <Riostream.h>
#include <TSystem.h>
#include "AliReconstruction.h"
#include "../TRD/AliTRDCalibra.h"
#include "AliCDBManager.h"
#include "TStopwatch.h"

#endif



void AliTRDFindTheCoef() 
{
  //
  // This macro takes a 2D histo or vector in the file TRD.calibration.root
  // tries to find the coeffficients
  // writes the result in the form of a tree in the file coeftest.root
 
  TStopwatch timer;
  timer.Start();

 
  
  ////Set the CDBManager(You have to use the same as during the reconstruction)*************************
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT"); 
  man->SetRun(0);

  ////Set the parameters of AliTRDCalibra***************
  AliTRDCalibra *calibra = AliTRDCalibra::Instance();
  
  ////Take the Histo2d or tree in the TRD.calibration.root file 
  TFile *file = new TFile("TRD.calibration.root","READ");
  //TProfile2D *h = (TProfile2D *) file->Get("PRF2d");
  TTree *h = (TTree *) file->Get("treePRF2d");
  //h->SetDirectory(0);
 
  
  //TProfile2D *h = (TProfile2D *) file->Get("PRF2d");
  TTree *h1 = (TTree *) file->Get("treePH2d");
  //h->SetDirectory(0);
 
  
  //TProfile2D *h = (TProfile2D *) file->Get("PRF2d");
  TTree *h2 = (TTree *) file->Get("treeCH2d");
  //h->SetDirectory(0);
 
 
  ////How many bins did you have?
  //calibra->SetNumberBinCharge(100);
  //calibra->SetNumberBinPRF(20);


  ////Which method do you want to use (It is always the default method that will be put in the database)
  //calibra->SetMeanChargeOn();
  //calibra->SetFitChargeBisOn();
  //calibra->SetFitPHOn();
  //calibra->SetPeriodeFitPH(10);

  //Some details?
  //calibra->SetRangeFitPRF(0.5);//fit from -0.5 and 0.5 with a gaussian the PRF
  //calibra->SetT0Shift(0.1433);//will always abstract 0.1433 mus to the result of the method for time 0

  ////What do you want to see?
  calibra->SetDebug(1);//0 (nothing to see), 1, 2, 3, or 4
  //calibra->SetDet(0,1,14);//in case of fDebug = 3 and 4
  //calibra->SetFitVoir(2);//in case of fDebug = 2

  ////How many statistics do you want to accept?
  calibra->SetMinEntries(10);// 1 entry at least to fit

  ////Do you want to write the result
  //calibra->SetWriteCoef(1);

  ////Do you want to change the name of the file (TRD.coefficient.root)
  //calibra->SetWriteNameCoef("coeftest.root");

  //Set the mode on the z and rphi direction for each calibration paramaters
  calibra->SetModeCalibrationFromTObject((TObject *)h2,0);
  calibra->SetModeCalibrationFromTObject((TObject *)h1,1);
  calibra->SetModeCalibrationFromTObject((TObject *)h,2);
  
  //Fit the pad response function 
  calibra->FitPRFOnline(h);  
   
  //Fit the avreage pulse height
  calibra->FitPHOnline(h1);  
  
  //Fit the deposited charge
  calibra->FitCHOnline(h2);  

  file->Close();
  
  timer.Stop();
  timer.Print();

}
