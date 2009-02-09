//
// Macro to be invoked before Calibration analysis 
// Setup TPC OCDB entries
// 
// This is just example macro  - some path are hardwired
//  TO BE MODIFIED BY USERS 



void ConfigOCDB(Float_t bfield){
  // 
  //
  // import geometry
  //

  printf("SETUP OCBD for PROOF\n");
  TGeoManager::Import("/u/miranov/proof/geometry.root");
  AliGeomManager::LoadGeometry("/u/miranov/proof/geometry.root");
  //
  //
  // Setup magnetic field
  //
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", 2, 1., 1., 10., AliMagF::k5kG));
  //
  //
  //
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Parameters","local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/ClusterParam","local://$ALICE_ROOT/OCDB");
  //  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/PadTime0","local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/PadTime0","local:///u/miranov/OCDB0");

  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data","local:///lustre_alpha/alice/alien/alice/data/2008/LHC08d/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Temperature","local:///lustre_alpha/alice/alien/alice/data/2008/LHC08d/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Goofie","local:///lustre_alpha/alice/alien/alice/data/2008/LHC08d/OCDB/");


  AliCDBManager::Instance()->SetRun(1);

  AliTPCClusterParam * paramCl = AliTPCcalibDB::Instance()->GetClusterParam(); 
  AliTPCParam   * paramTPC = AliTPCcalibDB::Instance()->GetParameters();
  paramCl->SetInstance(paramCl);
  //paramTPC->Dump();
  printf("\n\nSET EXB FIELD\t%f\n\n", bfield);
  AliTPCcalibDB::Instance()->SetExBField(bfield);
  //
  //
  //
  printf("END of SETUP OCBD for PROOF\n");
}


void ConfigAlien(){
  //
  // Setup-activate alien
  //

  //myvar=342
  //while [ $myvar -ne 360 ] ; do  echo enable alien on lxb$myvar; lsrun -m lxb$myvar  /u/miranov/.aliensetup;  myvar=$(( $myvar + 1 )) ; echo $myvar ; done 
  gSystem->Exec("/u/miranov/.aliensetup >setup.log"); 
  //ifstream in;
  //in.open("path.txt");
  
  TString envString;
  
  gSystem->Setenv("LD_LIBRARY_PATH",envString.Data());
  gSystem->Setenv("GBBOX_ENVFILE","/tmp/xxxxxxx");
  printf("LOAD LIBRARIES start\n\n\n");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libSTAT.so");
  gSystem->Load("libTPCcalib.so");
  //
  gSystem->Load("libXrdClient.so");
  gSystem->Load("libNetx.so");
  printf("LOAD LIBRARIES end\n\n\n");
  TGrid * alien = TGrid::Connect("alien://",0,0,"t");
  if (alien) {
    printf("Alien activated\n");
  }else{
    printf("Alien not activated\n");
  }
}
