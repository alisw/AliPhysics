/*
Implemented modes:
0 - write the same value for all towers
1 - use individual V30 settings
*/

int 
//____________________________________________________________________
void WriteBiasFix(const int mode = 0, const int biasSetting = 390) 
{ 

  if (mode == 0) { // fixed values
    SetAll(390);
  }
  else if (mode == 1) {
  }

}

  gSystem->Load("AliEMCALCalibAPD_cxx");
  AliEMCALCalibAPD *calibAPD = new AliEMCALCalibAPD();

  calibAPD->ReadCalibAPDInfo(10000, "dilan-APD-database.csv");
  //calibAPD->ReadCalibAPDInfo(10000, "paola-APD-database.csv");
  calibAPD->WriteCalibAPDInfo("dummy.txt");

  int fNCalibAPD = calibAPD->GetNCalibAPD();
  AliEMCALCalibAPD::AliEMCALCalibAPDData * fData = calibAPD->GetCalibAPDData();
  for (int i=0; i<fNCalibAPD; i++) {
    cout << " i " << i
	 << " fAPDNum " << fData[i].fAPDNum 
	 << " fSerialNum " << fData[i].fSerialNum 
	 << " fDarkCurrent " << fData[i].fDarkCurrent << endl;
  }

  /*
  // 1: create a dummy file
  calibAPD->GenerateDummyAPDInfo(nAPD);
  */

  /*
  // 2: test I/O
  calibAPD->ReadCalibAPDInfo(nAPD, "dummy.txt");
  calibAPD->WriteCalibAPDInfo("dummy2.txt");
  */

  /*
  // 3: see if it works ok if we genarate values first, and then try to read others/overwriting
  calibAPD->GenerateDummyAPDInfo(2*nAPD); // some extra APDs
  calibAPD->ReadCalibAPDInfo(nAPD, "dummy.txt");
  calibAPD->WriteCalibAPDInfo("dummy3.txt");
  */

  /*
  // 4: other way around from #3
  calibAPD->ReadCalibAPDInfo(nAPD, "dummy.txt");
  calibAPD->GenerateDummyAPDInfo(2*nAPD);
  calibAPD->WriteCalibAPDInfo("dummy4.txt");
  */
  
}

