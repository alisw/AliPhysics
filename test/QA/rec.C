void rec() {
  const char * kYear = "08" ; 
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08c/OCDB/");
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetSpecificStorage("EMCAL/*","local://DB");
  
  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetRecoParam("TPC",AliTPCRecoParam::GetLowFluxParam());
  reco.SetRecoParam("TRD",AliTRDrecoParam::GetLowFluxParam());
  reco.SetRecoParam("PHOS",AliPHOSRecoParam::GetDefaultParameters());
  reco.SetRecoParam("MUON",AliMUONRecoParam::GetLowFluxParam());
	
 	AliTPCReconstructor::SetStreamLevel(1);
//  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON T0 VZERO FMD PMD ZDC ACORDE");
  reco.SetRunReconstruction("ITS TPC TRD HMPID PHOS EMCAL MUON T0 VZERO FMD PMD ZDC ACORDE");
  reco.SetRunQA("ALL:ALL") ;
	 reco.SetInLoopQA() ; 
	  
	// AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
  AliQA::SetQARefStorage("local://$ALICE_ROOT") ;
  
// **** The field map settings must be the same as in Config.C !
  AliMagFMaps *field=new AliMagFMaps("Maps","Maps",2,1.,10.,AliMagFMaps::k5kG);
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}

