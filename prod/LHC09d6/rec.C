void rec() {

  AliReconstruction reco;
  
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetRunReconstruction("ITS VZERO");

  // ITS Efficiency
  reco.SetRunPlaneEff(kTRUE);
  AliITSRecoParam *itsRecoParam = AliITSRecoParam::GetPlaneEffParam(-1);
  itsRecoParam->SetOptTrackletsPlaneEff(kTRUE);


 //****** FIRST PHYSICS 2009 (same as COSMICS 2009) *********************
  //AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetLowFluxParam();
   itsRecoParam->SetClusterErrorsParam(2);
   // find independently ITS SA tracks
   itsRecoParam->SetSAUseAllClusters();
   itsRecoParam->SetOuterStartLayerSA(AliITSgeomTGeo::GetNLayers()-2);
   // to maximize efficiency
   itsRecoParam->SetAllowProlongationWithEmptyRoad();     // larger seach windows for SA (in case of large misalignments)
   itsRecoParam->SetNLoopsSA(33);
   itsRecoParam->SetFactorSAWindowSizes(20);
   // additional error due to misal (B off)
   itsRecoParam->SetClusterMisalErrorY(1.0,1.0,1.0,1.0,1.0,1.0); // [cm]
   itsRecoParam->SetClusterMisalErrorZ(1.0,1.0,1.0,1.0,1.0,1.0); // [cm]
   // additional error due to misal (B on)
   itsRecoParam->SetClusterMisalErrorYBOn(0.0,0.0,0.1,0.1,0.1,0.1); // [cm]
   itsRecoParam->SetClusterMisalErrorZBOn(0.1,0.1,0.1,0.1,0.1,0.1); // [cm]
   //***********************************************************************
   itsRecoParam->SetEventSpecie(AliRecoParam::kLowMult);

   reco.SetRecoParam("ITS",itsRecoParam);



  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
//   reco.SetSpecificStorage("GRP/GRP/Data",
// 			  Form("local://%s",gSystem->pwd()));
  // We store the object in AliEn during the simulation
  reco.SetSpecificStorage("GRP/GRP/Data",
			  "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");

			  
  reco.SetSpecificStorage("GRP/Calib/MeanVertexSPD","alien://folder=/alice/cern.ch/user/r/rgrosso/ShiftedVertex");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
