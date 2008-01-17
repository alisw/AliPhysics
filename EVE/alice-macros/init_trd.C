void init_trd()
{
  TString macdir("$(REVESYS)/alice-macros");
  gSystem->ExpandPathName(macdir);
  gROOT->GetListOfBrowsables()->Add
    (new TSystemDirectory(macdir.Data(), macdir.Data()));
  AssertMacro("region_marker.C");

  Alieve::TRDLoaderManager *trd=new Alieve::TRDLoaderManager("TRD manager", "Loader manager for TRD data monitoring");
  gEve->AddElement(trd);
  gEve->AddToListTree(trd, kTRUE);
}
