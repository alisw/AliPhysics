void DrawRICH()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewRICH.C");
   gInterpreter->ProcessLine("ViewRICH()");
   gMC->Gdopt("hide", "off");
   gMC->Gdopt("shad", "off");
   gMC->Gsatt("*", "fill", 1);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 2000, 0, 2000, 0, 2000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 30, 40, 0, -2, 4, .03, .03);
   gMC->Gdhead(1111, "Ring Imaging Cherenkov");
   gMC->Gdman(16, 6, "MAN");
}
