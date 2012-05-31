/*  
  Example macro to build AliTPCClusterParam
  postprocessing the output of the calibration using tracks
  In the future this macro will be part of the Preprocesor
  ..

*/
void StoreObject(AliTPCClusterParam *param);

void MakeClusterParam(const char *fnresolc="Output.root", const char *fnresolg="Output.root"){
  gSystem->Load("libTPCcalib.so");
  TFile fresolc(fnresolc);
  TFile fresolg(fnresolg);
  AliTPCcalibTracks *calibtracks = (AliTPCcalibTracks*)fresolc.Get("calibTracks");
  AliTPCcalibTracksGain *calibtracksGain = (AliTPCcalibTracksGain*)fresolg.Get("calibTracksGain");

  AliTPCClusterParam clParam;
  //
  // Make a resolution tree
  //
  calibtracks->MakeResPlotsQTree(200,"plots");  
  TFile fres("plots/resol.root");
  TTree *treeres = (TTree*)fres.Get("Resol");
  // Fit the resolution parameterization
  clParam.FitResol(treeres);
  clParam.FitRMS(treeres);
  clParam.SetInstance(&clParam);
  TF1 f1z_z("f1z_z","AliTPCClusterParam::SGetError0Par(1,0,x,0)",0,250);
  //
  // angular effect calibration - usable only with the 
  // cosmic tracks 
  //
  calibtracksGain->UpdateClusterParam(&clParam);
  //
  //
  //
  TFile fclparam("TPCClusterParam.root","recreate");
  clParam->Write("Param");
  fclparam.Close();
  StoreObject(clParam);
  AliTPCClusterParam::SetInstance(&clParam);
}
 
void StoreObject(AliTPCClusterParam *clParam)
{ 
  //
  //
  //
  Int_t gkDummyRun = 0;
  char *gCDBpath   = "local://$ALICE_ROOT/OCDB";
  AliCDBMetaData *md1= new AliCDBMetaData(); 
  AliCDBId id1("TPC/Calib/ClusterParam", gkDummyRun, gkDummyRun); 
  AliCDBStorage* gStorLoc = 0;
  AliCDBManager *man = AliCDBManager::Instance();
  gStorLoc = man->GetStorage(gCDBpath);

  md1->SetObjectClassName("AliTPCClusterParam");
  md1->SetResponsible("Marian Ivanov");
  md1->SetBeamPeriod(1);
  md1->SetAliRootVersion("v5-08-Release"); //root version
  md1->SetComment("Calibration data using the MC cosmic");
  gStorLoc->Put(&clParam, id1, md1); 
}


