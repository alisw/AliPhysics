void init_trd()
{
	TString macdir("$(REVESYS)/alice-macros");
	gSystem->ExpandPathName(macdir);
	gROOT->GetListOfBrowsables()->Add
	  (new TSystemDirectory(macdir.Data(), macdir.Data()));
	Reve::AssertMacro("region_marker.C");

	Alieve::TRDLoaderManager *trd=new Alieve::TRDLoaderManager("TRD manager", "Loader manager for TRD data monitoring");
	gReve->AddRenderElement(trd);
	gReve->AddToListTree(trd, kTRUE);
}
