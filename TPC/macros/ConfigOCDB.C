/// \file ConfigOCDB.C
///
/// Macro to Setup OCDB
/// This is just example macro
/// \author marian.ivanov@cern.ch
/// To be used:
/// 1. Before invocation of the calibration - in the calibration trains
/// 2. To setup calibration viewer.
///
/// ConfigOCDB  - setup default and specific data storage
/// SetupCustom - user sepcific configuration
///             - Values in local cache of OCDB are overwritten



void SetupCustom(Int_t run);

void ConfigOCDB(Int_t crun=-1){
  ///

  printf("SETUP OCBD for TPC\n");
  //
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Parameters","local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/ClusterParam","local:///u/miranov/OCDB/TPCcosmic2/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/PadTime0","local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Temperature","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/Goofie","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/HighVoltage","local:///lustre/alice/alien/alice/data/2009/OCDB/");
  Int_t run =crun;
  if (run<0) run =0;
  AliCDBManager::Instance()->SetRun(run);
  SetupCustom(run);
}


void SetupCustom(Int_t run){
  /// Custom part - to be empty once we are happy with the calibration
  ///
  /// Setup magnetic field

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
  if (bzfac==0) {  // force default magnetic field if 0 field used
    bzfac=1;
    bz=5;
  }
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", bzfac, 1., smag));

  printf("\n\nSET EXB FIELD\t%f\n\n", -bz);
  AliTPCcalibDB::Instance()->SetExBField(-bz);
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
  //
  tpcRecoParam->SetUseRadialCorrection(kFALSE);
  tpcRecoParam->SetUseQuadrantAlignment(kTRUE);
  //
  tpcRecoParam->SetUseSectorAlignment(kFALSE);
  tpcRecoParam->SetUseDriftCorrectionTime(kFALSE);
  tpcRecoParam->SetUseDriftCorrectionGY(kTRUE);
  tpcRecoParam->SetUseGainCorrectionTime(kFALSE);
  tpcRecoParam->SetUseFieldCorrection(kFALSE);
  tpcRecoParam->SetUseExBCorrection(kTRUE);
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


