void DrawABSO()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewABSO.C");
   gInterpreter->ProcessLine("ViewABSO()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 21.5, 15, .04, .04);
   gMC->Gdhead(1111, "Muon Absorber");
   gMC->Gdman(16, 6, "MAN");
}
