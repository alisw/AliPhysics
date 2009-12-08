//
// Macro to Setup OCDB  
// This is just example macro
// Responsible: marian.ivanov@cern.ch
// To be used:
// 1. Before invocation of the calibration - in the calibration trains
// 2. To setup calibration viewer.
//  
// ConfigOCDB  - setup default and specific data storage
// SetupCustom - user sepcific configuration 
//             - Values in local cache of OCDB are overwritten



void SetupCustom(Int_t run);

void ConfigOCDB(Int_t crun=-1){
  // 
  printf("SETUP OCBD for TPC\n");
  //
  Int_t run =crun;
  if (run<0) run =0;
  AliCDBManager::Instance()->SetDefaultStorage("local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetRun(run);
  //
  //
  // custom calibration to test before committing
  //
  TString ocdbStorage="local://";
  ocdbStorage+=gSystem->Getenv("workdir");
  ocdbStorage+="/calibNoDrift/OCDB";
  printf("OCDB storage\t%s\n",ocdbStorage.Data());
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/TimeGain",ocdbStorage.Data());
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/TimeDrift",ocdbStorage.Data());

  AliCDBManager::Instance()->Get("TPC/Calib/TimeDrift",run)->Dump();
  SetupCustom(run);
  AliTPCcalibDB::Instance()->SetRun(run);
}


void SetupCustom(Int_t run){
  //
  //
  // Custom part - to be empty once we are happy with the calibration
  //
  //
  // Setup magnetic field
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
  AliMagF* magF= new AliMagF("Maps","Maps", bzfac, 1., smag);
  TGeoGlobalMagField::Instance()->SetField(magF);  
  printf("\n\nSET EXB FIELD\t\n\n");
  AliTPCcalibDB::Instance()->SetExBField(magF);
  //
  //
  // import geometry
  //
  //
  TGeoManager::Import("/u/miranov/proof/geometry.root");
  AliGeomManager::LoadGeometry("/u/miranov/proof/geometry.root");

  AliTPCClusterParam * paramCl = AliTPCcalibDB::Instance()->GetClusterParam(); 
  AliTPCParam   * paramTPC = AliTPCcalibDB::Instance()->GetParameters();
  paramCl->SetInstance(paramCl);

  //
  // Setup reco param
  //
  AliTPCTransform *transform     = AliTPCcalibDB::Instance()->GetTransform() ;
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  transform->SetCurrentRecoParam(tpcRecoParam);
  tpcRecoParam->SetUseRPHICorrection(kTRUE); 
  tpcRecoParam->SetUseTOFCorrection(kFALSE);
  //
  tpcRecoParam->SetUseGainCorrectionTime(1);
  tpcRecoParam->SetUseDriftCorrectionTime(1);
  tpcRecoParam->SetUseDriftCorrectionGY(1);
  //
  tpcRecoParam->SetUseRadialCorrection(kFALSE);
  tpcRecoParam->SetUseQuadrantAlignment(kTRUE);
  //
  tpcRecoParam->SetUseSectorAlignment(kFALSE);
  tpcRecoParam->SetUseFieldCorrection(kFALSE);
  tpcRecoParam->SetUseExBCorrection(kTRUE);
  if (TMath::Abs(bzfacOrig)<0.05){
    tpcRecoParam->SetUseExBCorrection(kFALSE);
  }
  //
  //
  //
  TFile fposcor("~/OCDB/calibUnlin.root");
  AliTPCPointCorrection *pcorr = fposcor.Get("correction");
  if (pcorr) pcorr->SetInstance(pcorr); 
  //
  //
  //
  printf("END of SETUP OCBD for TPC\n");
}


