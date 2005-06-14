void
Reconstruction (char *fileName = "galice.root")
{

  TPluginManager *pluginManager = gROOT->GetPluginManager ();
  pluginManager->AddHandler ("AliReconstructor", "MUON",
			     "AliMUONReconstructor", "MUON",
			     "AliMUONReconstructor()");

  AliReconstruction MuonRec (fileName);
  MuonRec.SetRunTracking ("");
  MuonRec.SetRunVertexFinder (kFALSE);
  MuonRec.SetRunLocalReconstruction ("MUON");
  MuonRec.SetFillESD ("MUON");
  MuonRec.Run ();
}
