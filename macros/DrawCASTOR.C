//void DrawCASTOR()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->Macro("ViewCASTOR.C");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 20, -20, 20, -1900, -1700);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, -191.5, -78, .19, .19);
   gMC->Gdhead(1111, "Castor");
   gMC->Gdman(15,-2, "MAN");
}
