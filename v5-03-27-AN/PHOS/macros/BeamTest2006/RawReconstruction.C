void Rec(TString file="raw.root")
{
  // Reconstruction of RAW data from the input file raw.root
  // Boris Polichtchouk, 13 Mar 2008

  //AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local://BadMap");

  AliReconstruction rec ;
  rec.SetRunTracking("") ;
  rec.SetRunVertexFinder(kFALSE) ; 
  rec.SetRunLocalReconstruction("PHOS") ;
  rec.SetFillESD("PHOS") ;

  //Set rec. parameters different from the default ones.
  AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  recEmc->SetSubtractPedestals(kTRUE);
  recEmc->SetMinE(0.01);                     //Minimal Digit energy
  recEmc->SetClusteringThreshold(0.02);      //Minimal cluster seed energy
  recEmc->SetDecoderVersion("v1");           //Comment out to choose Max-Ped version
  recEmc->SetOldRCUFormat(kTRUE);

  AliPHOSReconstructor::SetRecoParamEmc(recEmc);

  rec.SetInput(file.Data());  // read RAW data
  rec.SetNumberOfEventsPerFile(5000);	
  rec.Run();


}
