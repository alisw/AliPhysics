//
// Script to draw detail of the FMD
//
void DrawInner()
{
  gAlice->Init("FMD/scripts/ConfigInner.C");
  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  gROOT->LoadMacro("FMD/ViewFMD.C");
  gInterpreter->ProcessLine("ViewFMD()");
  gROOT->LoadMacro("VZERO/ViewVZERO.C");
  gInterpreter->ProcessLine("ViewVZERO()");
  gROOT->LoadMacro("START/ViewSTART.C");
  gInterpreter->ProcessLine("ViewSTART()");
  gROOT->LoadMacro("macros/ViewITS.C");
  gInterpreter->ProcessLine("ViewITS()");
  // gROOT->LoadMacro("FMD/scripts/ViewPIPE.C");
  // gInterpreter->ProcessLine("ViewPIPE()");
  // gMC->Gsatt("ITSV", "seen", 1);
  gMC->Gsatt("0STR", "seen", 1);
  gMC->Gsatt("0STL", "seen", 1);
  gMC->Gsatt("0SUP", "seen", 1);
  // gMC->Gsatt("FMD1", "seen", 1);
  // gMC->Gsatt("FMD2", "seen", 1);
  gMC->Gsatt("FMD3", "seen", 0);
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 60, 0, 0, 10, 10, .10, .10);
  // gMC->Gdhead(1111, "FMD3 detail");
  // gMC->Gdman(16, 10, "MAN");
  
  gPad->Modified();
  gPad->cd();
  gPad->Print("Inner.png");
}
