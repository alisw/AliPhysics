void DrawPIPE()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewPIPE.C");
   gInterpreter->ProcessLine("ViewPIPE()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 3, 5, .04, .04);
   gMC->Gdhead(1111, "Beam Pipe");
   gMC->Gdman(16, 6, "MAN");
}
