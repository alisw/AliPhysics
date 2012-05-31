void raw2digits(const char * inputfile = "RUN_4400.fz"){
  //Read Connection table DB
  TFile f("ConTableDB.root") ;
  AliPHOSConTableDB * cdb = f.Get("AliPHOSConTableDB") ;
  f.Close() ;

  AliPHOSRaw2Digits * r = new AliPHOSRaw2Digits(inputfile) ;
  //Set position of the target in the given run.
  //The reference system is following
  //Z axis along beam direction, from target to prototype (0-surface of prototype)
  //X axis along columns of prototype (0-center of prototype)
  //Y axis along raws of prototype    (0-center of prototype)
  Double_t pos[3]={0,0,-452.} ;
  r->SetTargetPosition(pos) ;
  r->SetConTableDB(cdb) ;
//   r->SetDebugLevel(1) ;
  r->ExecuteTask() ;
  r->Print("deb") ;
	

}
