//void DrawMUON()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->Macro("ViewMUON.C");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 2000, -2000, 2000, 0, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 22, 16, .017, .017);
   gMC->Gdhead(1111, "Muon Arm");
   gMC->Gdman(16, 6, "MAN");
}
