void ConvertCalib(const char * name = "gains.dat.final.gain0"){
  //macro to convert ASCII file with calibration parametrs 
  //to AliRoot structure AliPHOSCalibrationData.

  AliPHOSCalibrator * c = new AliPHOSCalibrator() ;
  c->SetConTableDB("ConeTableDB.root") ;
  c->ReadFromASCII(name) ;
  cout << "Write Ped " << endl ;
  c->WritePedestals("v1") ;
  cout << "Write Ped OK" << endl ;
  c->WriteGains("v1") ;
  cout << "Write Gains OK" << endl ;
  delete c ;
}
