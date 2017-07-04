// -*- C++ -*-


void ProcVdM(Int_t fillNumber)
{
  gROOT->LoadMacro("AliLuminousRegionFit.cxx+");

  gROOT->LoadMacro(Form("MakeLumiRegion_%d.C+", fillNumber));
  gROOT->LoadMacro(Form("PlotLumiRegion_%d.C", fillNumber));

  gSystem->Exec(Form("mkdir -p {root,pdf}/%d", fillNumber));

  MakeLumiRegion();
  PlotLumiRegion();

}
