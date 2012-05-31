void makeCalDB(void){
  //Read Connection Table DB
  TFile f("ConTableDB.root") ;
  AliPHOSConTableDB * ctdb = (AliPHOSConTableDB *)f.Get("AliPHOSConTableDB") ;
  f.Close() ;

  AliPHOSCalibrationDB * c = new AliPHOSCalibrationDB("beamtest.root") ;
  c->SetConTableDB(ctdb) ;
  c->ReadCalibrationParameters("gains.dat","gains") ;
  c->ReadCalibrationParameters("pedestals.dat","pedestals") ;
  TFile g("beamtest.root","recreate") ;
  g.cd() ;
  c->Write(0) ;
  g.Close() ;
}
