//void DrawCPV()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->Macro("ViewCPV.C");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 13, 20, .025, .025);
   gMC->Gdhead(1111, "CPV Detector");
   gMC->Gdman(18, 4, "MAN");
}
