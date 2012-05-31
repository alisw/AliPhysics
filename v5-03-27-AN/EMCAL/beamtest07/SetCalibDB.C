// Script to create calibration parameters and store them into CDB 
// from gain values determined for 2007 beam test
//
// JLK 2007-12-04
//
//
#if !defined(__CINT__)
#include "AliRun.h"
#include "AliEMCALCalibData.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif

void SetCalibDB() {

  Int_t firstRun   =  0; // What is this                                                
  Int_t lastRun    = 9999999;
  Int_t beamPeriod =  1;
  char* objFormat  = "";

  TString DBFolder  ="local://$ALICE_ROOT/OCDB/EMCAL/beamtest07";
  firstRun  =  0;
  lastRun   =  999999999;
  objFormat = "EMCAL beam test 2007 gain factors and pedestals";

  AliEMCALCalibData *calibda=new AliEMCALCalibData("EMCAL");

  Float_t fADCpedestal = 0.009;
  Float_t fADCchannel  = 0.0153;  // 250 GeV / (16*1024)                               
  Float_t ped = 0.;
  Float_t cc = fADCchannel;

  Int_t nSMod  = 12;
  Int_t nCol   = 48;
  Int_t nRow   = 24;
  Int_t nRow2  = 12; //Modules 11 and 12 are half modules                          

  Int_t colOffset = 40;
  Int_t rowOffset = 8;

  Double_t gain_ratios[8][8] =
    {
      15.9289, 16.2141, 16.1204, 15.9118,
      15.9363, 15.9402, 16.2257, 16.0097,
      16.058, 16.1116, 16.039, 16.4167,
      16.2148, 16.1399, 16.1515, 16.2194,
      15.9082, 16.0776, 16.0496, 16.2353,
      15.8054, 16.2158, 16.2344, 16.1023,
      15.8903, 16.2387, 16.13, 16.157,
      16.0685, 16.172, 16.3495, 16.3887,
      16.2842, 16.049, 16.4328, 16.3954,
      16.4226, 15.7254, 16.1634, 16.3182,
      16.4216, 16.1201, 16.0000, 16.2305,
      16.0266, 16.3573, 16.1382, 16.237,
      16.2981, 16.1796, 15.854, 16.4189,
      15.6425, 16.287, 16.3293, 16.6308,
      16.2469, 16.0412, 16.252, 16.3367,
      16.1412, 16.0646, 16.3996, 16.3479
    };

  
  Float_t gains[8][8] =
    {
      4.43274, 6.7283, 8.23733, 3.59882,
      4.2717, 2.85658, 4.86389, 2.71961,
      3.05523, 3.02552, 3.50615, 3.26494,
      6.69024, 2.51058, 8.42275, 2.83824,
      8.05074, 5.36051, 4.36794, 4.73468,
      9.9684, 5.5, 6.42999, 5.6,
      7.37306, 5.28314, 5.27662, 5.26982,
      3.29468, 5.23107, 6.40948, 4.06855,
      4.09685, 5.37323, 5.32816, 5.89487,
      9.2395, 5.3, 4.77239, 5.0,
      4.85923, 3.44063, 4.74517, 5.28772,
      3.80171, 4.84878, 5.12039, 4.59205,
      2.34745, 3.16971, 3.61231, 3.65195,
      3.43496, 3.4, 3.65678, 2.9,
      2.71648, 3.39577, 3.40896, 3.31741,
      3.24286, 3.51346, 2.61503, 3.44246
    };

  for(Int_t supermodule=0; supermodule < nSMod; supermodule++) {
    for(Int_t column=0; column< nCol; column++) {
      if(supermodule >= 10)
        nRow = nRow2;
      for(Int_t row=0; row< nRow; row++) {
	if(supermodule < 2 && column > 39 && row > 7 && row < 16) {
	  cc = 1./gain_ratios[column-colOffset][row-rowOffset]/gains[column-colOffset][row-rowOffset];
	  cout << "column = " << column << " column - colOffset = " << column-colOffset << " row = " << " row Offset = " << row-rowOffset << endl;
	}
        calibda->SetADCchannel(supermodule,column,row,cc);
        calibda->SetADCpedestal(supermodule,column,row,ped);
      }
    }
  }

  //Store calibration data into database
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("David Silvermyr");
  AliCDBId id("EMCAL/Calib/Data",firstRun,lastRun); // create in
						    // EMCAL/Calib/Data DBFolder                     

  AliCDBManager* man = AliCDBManager::Instance();
  AliCDBStorage* loc = man->GetStorage(DBFolder.Data());
  loc->Put(calibda, id, &md);

}







