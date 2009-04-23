void rec() {
 
	AliReconstruction reco;
	reco.SetWriteESDfriend();
	reco.SetWriteAlignmentData();

	reco.SetRecoParam("ITS",AliITSRecoParam::GetLowFluxParam());
	reco.SetRecoParam("TPC",AliTPCRecoParam::GetLowFluxParam());
	reco.SetRecoParam("TRD",AliTRDrecoParam::GetLowFluxParam());
	reco.SetRecoParam("PHOS",AliPHOSRecoParam::GetDefaultParameters());
	reco.SetRecoParam("MUON",AliMUONRecoParam::GetLowFluxParam());
	reco.SetRecoParam("EMCAL",AliEMCALRecParam::GetLowFluxParam());

 	reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");


	/*
	  reco.SetSpecificStorage("GRP/GRP/Data",
	  "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
	*/
	// No write access to the OCDB => local specific storage
	reco.SetSpecificStorage("GRP/GRP/Data",
				Form("local://%s/../",gSystem->pwd()));


	// add TRD standalone tracks
	reco.SetOption("TRD", "sl_tr_1");

 	reco.SetRunQA("ALL:ALL");

 
	TStopwatch timer;
	timer.Start();
	reco.Run();
	timer.Stop();
	timer.Print();
}
