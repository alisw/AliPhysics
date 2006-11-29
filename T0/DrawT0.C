void DrawT0()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewT0.C");
   gInterpreter->ProcessLine("ViewT0()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   //  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
   gMC->DefaultRange();
   //gMC->Gdraw("alic", 40, 30, 0, -8, 2, 1, 1);
   gMC->Gdraw("alic", 40, 30, 0, 10, 9, 1., 1.);
   gMC->Gdhead(1111, "T0 Detector");
   //  gMC->Gdman(13, 9, "MAN");
}
