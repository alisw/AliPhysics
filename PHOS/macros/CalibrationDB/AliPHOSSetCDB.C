/* $Id$ */

// Script to create calibration parameters and store them into CDB
// Two sets of calibration parameters can be created:
// 1) equal parameters
// 2) randomly distributed parameters for decalibrated detector silumations

void AliPHOSSetCDB()
{
   menu = new TControlBar("vertical","PHOS CDB");
   menu->AddButton("Help to run PHOS CDB","Help()",
		   "Explains how to use PHOS CDS menus");
   menu->AddButton("Equal CC","SetCC(0)",
		   "Set equal calibration coefficients");
   menu->AddButton("Decalibrate","SetCC(1)",
		   "Set random calibration coefficients");
   menu->Show();
}

//------------------------------------------------------------------------
Help()
{
  char *string =
    "\nSet calibration parameters and write them into ALICE CDB.
Press button \"Equal CC\" to create equal pedestals and gain factors.
Press button \"Decalibrate\" to create random pedestals and gain factors to imitate decalibrated detector\n";
  printf(string);
}
//------------------------------------------------------------------------
SetCC(Int_t flag=0)
{
  // Writing calibration coefficients into the Calibration DB
  // Arguments:
  //   flag=0: all calibration coefficients are equal
  //   flag=1: all calibration coefficients random (decalibration)
  // Author: Boris Polishchuk (Boris.Polichtchouk@cern.ch)

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  char* objFormat;

  if      (flag == 0) {
    DBFolder  ="InitCalibDB";
    firstRun  =  0;
    lastRun   =  0;
    objFormat = "PHOS initial gain factors and pedestals";
  }
  else if (flag == 1) {
    DBFolder  ="DeCalibDB";
    firstRun  =  0;
    lastRun   = 10;
    objFormat = "PHOS random pedestals and ADC gain factors (5x64x56)";
  }

  // create DB directory
  if(!gSystem->OpenDirectory(DBFolder)){
    printf("Warning: folder %s does not exist, I will create it!",
	   DBFolder.Data());
    TString command = "mkdir "+ DBFolder;
    gSystem->Exec(command.Data());
  }

  AliPHOSCalibData *calibda=new AliPHOSCalibData("PHOS");
  
  Float_t fADCpedestalEmc = 0.005;
  Float_t fADCchanelEmc   = 0.0015;

  TRandom rn;

  for(Int_t module=1; module<6; module++) {
    for(Int_t column=1; column<57; column++) {
      for(Int_t row=1; row<65; row++) {
	if (flag == 0) {
	  // Decalibration:
	  // Spread calibration coefficients uniformly with
	  // Cmax/Cmin = 5, (Cmax-Cmin)/2 = 0.0015
	  // and pedestals 0.005 +-10%
	  fADCchanelEmc  =rn.Uniform(0.00075,0.00375);
	  fADCpedestalEmc=rn.Uniform(0.0045,0.0055);
	}
	calibda->SetADCchannelEmc (module,column,row,fADCchanelEmc);
	calibda->SetADCpedestalEmc(module,column,row,fADCpedestalEmc);
      }
    }
  }

  //Store calibration data into database

  AliCDBMetaData md("PHOS/Calib/GainFactors_and_Pedestals",
		    firstRun,lastRun,beamPeriod,
		    objFormat,
		    "B. Polishchuk", 
		    "PHOS calibration");

  AliCDBLocal *loc = new AliCDBLocal(DBFolder.Data());
  AliCDBStorage::Instance()->Put(calibda, md);
  AliCDBStorage::Instance()->Delete();
}
