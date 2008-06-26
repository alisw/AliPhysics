void rec() {
  const char * kYear = "08" ; 
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("alien://Folder=/alice/simulation/2007/PDC07_v4-09-Rev-00/Ideal/CDB/");
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetSpecificStorage("EMCAL/*","local://DB");
  
  AliReconstruction reco;

  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetLowFluxParam();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(0);
  reco.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS EMCAL MUON T0 VZERO FMD PMD ZDC");
  //AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  //	recEmc->SetSubtractPedestals(kFALSE);
  //	AliPHOSReconstructor::SetRecoParamEmc(recEmc);  
	if (! reco.SetRunQA("ALL:ALL") ) {
		printf("SetRunQA ERROR\n") ; 
		exit ; 
	}
	reco.SetInLoopQA() ; 
	  
  AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
  AliQA::SetQARefDataDirName("Sim") ; //Data, Pedestals, BlackEvent, .....
  
// **** The field map settings must be the same as in Config.C !
  AliMagFMaps *field=new AliMagFMaps("Maps","Maps",2,1.,10.,AliMagFMaps::k5kG);
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);

  TStopwatch timer;
  timer.Start();
  gSystem->Exec("mv -f galice.root galice0.root");
  reco.Run();
  gSystem->Exec("mv -f galice0.root galice.root")
  timer.Stop();
  timer.Print();
}

