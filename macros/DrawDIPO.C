void DrawDIPO()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewDIPO.C");
   gInterpreter->ProcessLine("ViewDIPO()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
//   gMC->SetClipBox("*", 0, 1000, -1000, 1000, -3000, 3000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 30, 30, 0, 17, 13.5, .019, .019);
   gMC->Gdhead(1111, "Magnetic Dipole");
   gMC->Gdman(16, 4, "MAN");
}
