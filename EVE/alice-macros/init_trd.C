void init_trd()
{
	TString macdir("$(REVESYS)/alice-macros");
  gSystem->ExpandPathName(macdir);
  gROOT->GetListOfBrowsables()->Add
    (new TSystemDirectory(macdir.Data(), macdir.Data()));
	Reve::AssertMacro("region_marker.C");
	
	Alieve::TRDLoader *trdl=new Alieve::TRDLoader("Simulations", "TRDLoader for simulations");
	gReve->AddRenderElement(trdl);
	gReve->DrawRenderElement(trdl);
	
	Alieve::TRDLoader *trdld=new Alieve::TRDLoaderSingle("Digits", "TRDLoaderSingle for Digits");
	gReve->AddRenderElement(trdld);
	gReve->DrawRenderElement(trdld);
	
	Alieve::TRDLoader *trdlc=new Alieve::TRDLoaderSingle("Clusters", "TRDLoaderSingle for Clusters");
	gReve->AddRenderElement(trdlc);
	gReve->DrawRenderElement(trdlc);
}