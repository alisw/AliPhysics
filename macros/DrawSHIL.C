void DrawSHIL()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewSHIL.C");
   gInterpreter->ProcessLine("ViewSHIL()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 30, 30, 0, 25, 18, .025, .025);
   gMC->Gdhead(1111, "Muon Shield");
   gMC->Gdman(16, 6, "MAN");
}
