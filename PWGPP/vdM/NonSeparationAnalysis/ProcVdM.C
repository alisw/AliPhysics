// -*- C++ -*-


void ProcVdM(Int_t fillNumber, Int_t bc=-1)
{
  gROOT->LoadMacro("AliLuminousRegionFit.cxx+");
  gROOT->LoadMacro("PlotBkgd.C+");

  gROOT->LoadMacro(Form("MakeLumiRegion_%d.C+", fillNumber));
  gROOT->LoadMacro(Form("PlotLumiRegion_%d.C", fillNumber));
  gROOT->LoadMacro(Form("PlotBkgdLumiRegion_%d.C+", fillNumber));

  gSystem->Exec(Form("mkdir -p {root,pdf}/%d", fillNumber));

//  MakeLumiRegion();
  if (bc == -1)
    PlotLumiRegion();
  else
    PlotLumiRegion(bc);
//  PlotBkgdLumiRegion();
}
