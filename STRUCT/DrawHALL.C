void DrawHALL()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewHALL.C");
   gInterpreter->ProcessLine("ViewHALL()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 12, 7.5, .005, .005);
   gMC->Gdhead(1111, "Experimental Hall");
   gMC->Gdman(18, 2, "MAN");
}
