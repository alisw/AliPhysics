void makecontDB(void){
  AliPHOSConTableDB * c = new AliPHOSConTableDB("Beamtest2002") ;
  c->SetNRaws(8) ;
  c->SetNCols(8) ;
  c->BuildDB() ;
  TFile f("ConTableDB.root","recreate") ;
  f.cd() ;
  c->Write(0) ;
  f.Close() ;

}
