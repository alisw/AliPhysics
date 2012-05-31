void RecRaw(char * file)
{
  // Reconstruction of RAW data from the input (raw) root file 
  // D.Peressounko after Boris Polichtchouk, 31 Aug 2007

   AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  // AliCDBManager::Instance()->SetDefaultStorage("local://./");
  // Provide here address of misalignment parametrs, Calibration or bad modules maps
  // AliCDBManager::Instance()->SetSpecificStorage("PHOS/*","local:///data/prsnko/");


  AliReconstruction rec ;
  rec.SetOption("PHOS","OldRCUFormat");
  rec.SetRunTracking("") ;
  rec.SetRunVertexFinder(kFALSE) ; 
  //Digits are produced as a by-product of local reconstruction...
  rec.SetRunLocalReconstruction("PHOS") ;
  //Here we do not want to produce ESD
  rec.SetFillESD("") ;

  //Uncomment following lines if you want to set rec. 
  //parameters other than default ones.
  // AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  // recEmc->SetSubtractPedestals(kTRUE);
  // //Minimal energy of digits used in clusteriztion
  // recEmc->SetMinE(0.01);
  // //Minimal energy of cluster seed
  // recEmc->SetClusteringThreshold(0.02);
  // //Choose here method of energy/time extraction:
  // //fitting of samples - "v1"
  // //maximal value extraction - ""
  // recEmc->SetDecoderVersion("v1");
  // AliPHOSReconstructor::SetRecoParamEmc(recEmc);

  rec.SetInput(file);  // read RAW data
  //If necessary, set first event and number of events to reconstruct
  //  rec.SetEventRange(0,20);
  rec.Run();

}
