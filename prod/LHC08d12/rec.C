void rec() {

  AliReconstruction reco;
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();

  reco.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");

  reco.SetRecoParam("ITS",AliITSRecoParam::GetHighFluxParam());
  reco.SetRecoParam("TPC",AliTPCRecoParam::GetHighFluxParam());
  reco.SetRecoParam("TRD",AliTRDrecoParam::GetHighFluxParam());
  reco.SetRecoParam("PHOS",AliPHOSRecoParam::GetDefaultParameters());
  reco.SetRecoParam("EMCAL",AliEMCALRecParam::GetHighFluxParam());
  reco.SetRecoParam("MUON",AliMUONRecoParam::GetHighFluxParam());

  // No write access to the OCDB => specific storage
  reco.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));

  reco.SetRunQA("ALL:ALL");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
