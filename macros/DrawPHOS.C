void DrawPHOS()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewPHOS.C");
   gInterpreter->ProcessLine("ViewPHOS()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 500, -800, 800, -200, 200);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 13, 20, .035, .035);
   gMC->Gdhead(1111, "Photon Detector");
   gMC->Gdman(18, 4, "MAN");
}
