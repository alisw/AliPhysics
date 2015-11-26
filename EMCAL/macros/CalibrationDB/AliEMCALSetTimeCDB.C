/// \file PrintEMCALCalibTime.C
/// \brief Print time parameters in OCDB
///
/// Script to create calibration parameters and store them into CDB
/// Thre sets of calibration parameters can be created:
/// 1) equal parameters
/// 2) randomly distributed parameters for decalibrated detector silumations
/// 3) gausian  distributed parameters for decalibrated detector silumations
/// 
/// Execute like this:
///.x $ALICE_ROOT/EMCAL/macros/CalibrationDB/AliEMCALSetTimeCDB.C
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__)
#include <TControlBar.h>
#include <TString.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>

#include "AliRun.h"
#include "AliEMCALCalibTime.h"
#include "AliEMCALGeoParams.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif

///
/// Main method
/// When execution, menu appears
//------------------------------------------------------------------------
void AliEMCALSetTimeCDB()
{
  TControlBar *menu = new TControlBar("vertical","EMCAL CDB");
  menu->AddButton("Help to run EMCAL CDB","Help()",
		  "Explains how to use EMCAL CDS menus");
  menu->AddButton("Equal CC","SetCC(0)",
		  "Set equal calibration coefficients");
  menu->AddButton("Random De-calibration","SetCC(1)",
		  "Set random decalibration calibration coefficients");
  menu->AddButton("Gaussian De-calibration","SetCC(2)",
		  "Set gausian decalibration calibration coefficients");
  //  menu->AddButton("Read equal CC","GetCC(0)",
  //		  "Read initial equal calibration coefficients");
  //  menu->AddButton("Read random CC","GetCC(1)",
  //		  "Read random decalibration calibration coefficients");
  //  menu->AddButton("Read gaussian CC","GetCC(2)",
  //		  "Read gausian decalibration calibration coefficients");
  menu->Show();
}

///
//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\nSet calibration parameters and write them into ALICE CDB. Press button \"Equal CC\" to create equal pedestals and gain factors. Press button \"Decalibrate\" to create random pedestals and gain factors to imitate decalibrated detector\n";
  printf(string);
}

///
/// Writing calibration coefficients into the Calibration DB
/// Arguments:
///   flag=0: all calibration coefficients are equal
///   flag=1: all calibration coefficients random (decalibration)
///   flag=2: all calibration coefficients have Gaussian random distribution (decalibration)
//------------------------------------------------------------------------
void SetCC(Int_t flag=0)
{
  TString dbFolder;
  Int_t firstRun   =  0; // What is this
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  char* objFormat  = "";
  
  if      (flag == 0)
  {
    dbFolder  ="local://InitCalibDB";
    firstRun  =  0;
    lastRun   =  AliCDBRunRange::Infinity();
    objFormat = "EMCAL constant time shift of 600 ns";
  }
  else if (flag == 1)
  {
    dbFolder  ="local://DeCalibDB";
    firstRun  =  0;
    lastRun   = AliCDBRunRange::Infinity();
    objFormat = "EMCAL random time shift";
  }
  else if (flag == 2)
  {
    dbFolder  ="local://DeCalibGausDB"; // create directory DeCalibDB in current directory
    firstRun  =  0;
    lastRun   = AliCDBRunRange::Infinity();
    objFormat = "EMCAL random gaussian time shift";
  }
  
  AliEMCALCalibTime *calibti=new AliEMCALCalibTime("EMCAL");
  
  Float_t timeShift0[] = {600.,600.,600.,600.} ;  
  Float_t timeShift [] = {600.,600.,600.,600.} ;  
  
  TRandom rn;
  Int_t nSMod = AliEMCALCalibTime::fgkECALDCALModules;
  Int_t nCol  = AliEMCALGeoParams::fgkEMCALCols;
  Int_t nRow  = AliEMCALGeoParams::fgkEMCALRows;
 
  for(Int_t bc = 0; bc < 4; bc ++)
  {
    printf("--->>> BC %d\n",bc);
    for(Int_t supermodule = 0; supermodule < nSMod; supermodule++)
    {
      printf("\t --->>> Fill SM %d\n",supermodule);
      // Set all the channels even the known to not exist in 1/3 sm and DCAL
      for(Int_t column = 0; column < nCol; column++)
      {
        for(Int_t row = 0; row < nRow; row++)
        {
          if (flag == 1)
          {
            // Spread calibration coefficients uniformly 
            timeShift[bc]  = rn.Uniform(500,700);
          }
          else if (flag == 2)
          { // Gaussian
            timeShift[bc] = rn.Gaus(timeShift0[bc], 0.1 );
          }
          
          calibti->SetTimeChannel(supermodule,column,row,bc,timeShift[bc]);
          
          //cout<<"Set col "<<column<<" row "<<row<<" shift "<< timeShift[bc]<<endl;
        } // row
      } // col
    } // SM
  } // BC
  
  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Gustavo Conesa");
  
  AliCDBId id("EMCAL/Calib/Time",firstRun,lastRun); // create in EMCAL/Calib/Time dbFolder
  
  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBStorage* loc = man->GetStorage(dbFolder.Data());
  loc->Put(calibti, id, &md);
  
}

