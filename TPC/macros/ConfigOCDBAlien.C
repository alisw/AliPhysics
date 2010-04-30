//
// Macro to Setup OCDB  
// This is just example macro - using the alien OCDB
// Responsible: marian.ivanov@cern.ch
// To be used:
// 1. Before invocation of the calibration - in the calibration trains
// 2. To setup calibration viewer.
//  
// ConfigOCDB  - setup default and specific data storage
// SetupCustom - user sepcific configuration 
//             - Values in local cache of OCDB are overwritten



void SetupCustom(Int_t run);

void ConfigOCDBAlien(Int_t crun=-1){
  // 
  TGrid * alien =     TGrid::Connect("alien://",0,0,"t"); 
  printf("SETUP OCBD for TPC\n");
  //
  AliCDBManager::Instance()->SetDefaultStorage("raw://");
   
  Int_t run =crun;
  if (run<0) run =0;
  AliCDBManager::Instance()->SetRun(run);
  AliTPCcalibDB::Instance()->SetRun(run);
  SetupCustom(run);
}


void SetupCustom(Int_t run){
  //
  //
  // Custom part - to be empty once we are happy with the calibration
  //
  //
  // Setup magnetic field - In future this should be part of GRP functionality
  //
  AliGRPObject *grp = AliTPCcalibDB::GetGRP(run);
  Float_t current = 0;
  Float_t bz      = 0;
  if (grp){
    current = grp->GetL3Current((AliGRPObject::Stats)0);
    bz = 5*current/30000.;
    printf("Run%d\tL3 current%f\tBz\t%f\n",run,current,bz);
  }
  else{
    printf("Run%d\tL3 current%f\tBz\t%f\n",run,current,bz);
  }
  AliMagF::BMap_t smag = AliMagF::k5kG;
  Double_t bzfac = bz/5;
  Double_t bzfacOrig=bzfac;
  if (TMath::Abs(bzfac)<0.01) {  // force default magnetic field if 0 field used
    bzfac=1;
    bz=5;
  }
  AliMagF * magF = new AliMagF("Maps","Maps", 2, bzfac, 1., smag);
  TGeoGlobalMagField::Instance()->SetField(magF);
  printf("\n\nSET EXB FIELD\t\n\n");
  AliTPCcalibDB::Instance()->SetExBField(magF);
  
  AliTPCClusterParam * paramCl = AliTPCcalibDB::Instance()->GetClusterParam(); 
  AliTPCParam   * paramTPC = AliTPCcalibDB::Instance()->GetParameters();
  paramCl->SetInstance(paramCl);

   if (TMath::Abs(bzfacOrig)<0.05){
    tpcRecoParam->SetUseExBCorrection(kFALSE);
  }
  //
  //
  //
   printf("END of SETUP OCBD for TPC\n");
}


