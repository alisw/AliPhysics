void DrawTOF0()
{
   geant3->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewTOF0.C");
   gInterpreter->ProcessLine("ViewTOF0()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
   gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .025, .025);
   gMC->Gdhead(1111, "Time Projection Chamber");
   gMC->Gdman(18, 4, "MAN");
}
