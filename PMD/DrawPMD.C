void DrawPMD()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewPMD.C");
   gInterpreter->ProcessLine("ViewPMD()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 22, 15.5, .04, .04);
   gMC->Gdhead(1111, "Photon Multiplicity Detector");
   gMC->Gdman(17, 5, "MAN");
}
