Int_t AliITSStoreFindableTracks
(Int_t nMinClusters = 5, const Text_t *evname = "galice", Int_t evnum = 0)
{
	gSystem->SetIncludePath("-I- -I$ALICE_ROOT/ITS -I$ALICE_ROOT/STEER");

	gROOT->LoadMacro("$ALICE_ROOT/ITS/AliITSStoreFindableTracksCompiled.C+");
	AliITSStoreFindableTracksCompiled(nMinClusters, evname, evnum);
}
