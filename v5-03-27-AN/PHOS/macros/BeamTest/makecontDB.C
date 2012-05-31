void makecontDB(void){
  AliPHOSConTableDB * c = new AliPHOSConTableDB("Beamtest2002") ;
  c->SetNRaws(16) ;
  c->SetNCols(16) ;
  c->BuildDB() ;
  TFile f("ConTableDB.root","recreate") ;
  f.cd() ;
  c->Write(0) ;
  f.Close() ;

}
