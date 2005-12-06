//
//  $Id$
// 
// Script to draw detail of the FMD
//
void DrawDetail()
{
  // gAlice->Init("FMD/scripts/ConfigInner.C");
  AliLog::SetModuleDebugLevel("FMD", 6);
  gAlice->Init("$(ALICE)/FMD/Config.C");
  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  gROOT->LoadMacro("$(ALICE)/FMD/ViewFMD.C");
  gInterpreter->ProcessLine("ViewFMD()");
  // gROOT->LoadMacro("VZERO/ViewVZERO.C");
  // gInterpreter->ProcessLine("ViewVZERO()");
  // gROOT->LoadMacro("START/ViewSTART.C");
  // gInterpreter->ProcessLine("ViewSTART()");
  // gROOT->LoadMacro("macros/ViewITS.C");
  // gInterpreter->ProcessLine("ViewITS()");
  // gROOT->LoadMacro("FMD/scripts/ViewPIPE.C");
  // gInterpreter->ProcessLine("ViewPIPE()");
  // gMC->Gsatt("ITSV", "seen", 1);
  // gMC->Gsatt("0STR", "seen", 1);
  // gMC->Gsatt("0STL", "seen", 1);
  // gMC->Gsatt("0SUP", "seen", 1);
  // gMC->Gsatt("FMD1", "seen", 1);
  // gMC->Gsatt("FMD2", "seen", 1);
  // gMC->Gsatt("FOAC", "seen", 1);
  // gMC->Gsatt("FOVF", "seen", 1);
  // gMC->Gsatt("FOVB", "seen", 1);
  // gMC->Gsatt("FIRG", "seen", 1);
  // gMC->Gsatt("FIRG", "colo", 1);
  // gMC->Gsatt("F3IH", "seen", -1);
  // gMC->Gsatt("F3II", "seen", -1);
  // gMC->Gsatt("F3IJ", "seen", -1);
  // gMC->Gsatt("F3IK", "seen", -1);
  gMC->Gsatt("FMD3", "seen", 1);
  gMC->Gsatt("FMD3", "colo", 3);
  gMC->Gsatt("FMD3", "lsty", 1);
  gMC->Gsatt("F3SL", "colo", 2);
  gMC->Gsatt("F3SN", "colo", 2);
  gMC->Gsatt("F3SB", "colo", 2);
  gMC->Gsatt("F3SF", "colo", 2);
  // gMC->Gsatt("FIAC", "lwid", 2);
  // gMC->Gsatt("FIAC", "colo", 3);
  // gMC->Gsatt("FOAC", "lwid", 2);
  // gMC->Gsatt("FOAC", "colo", 3);
  gMC->Gdopt("hide", "off");
  gMC->Gdopt("shad", "off");
#if 1
  gMC->Gdopt("hide", "on");
  // gMC->Gdopt("shad", "on");
  // gMC->Gsatt("*", "fill", 7);
#endif 
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 10000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 90, 0, 0, -5, 10, .15, .15);

  gPad->Modified();
  gPad->cd();
  gPad->Print("FMD3_detail.png");
}
//____________________________________________________________________
//
// EOF
//
