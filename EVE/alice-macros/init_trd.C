void init_trd()
{
	TString macdir("$(REVESYS)/alice-macros");
  gSystem->ExpandPathName(macdir);
  gROOT->GetListOfBrowsables()->Add
    (new TSystemDirectory(macdir.Data(), macdir.Data()));
	Reve::AssertMacro("region_marker.C");
	
	Alieve::TRDLoader *trdl=new Alieve::TRDLoader();
	gReve->AddRenderElement(trdl);
	gReve->DrawRenderElement(trdl);
}