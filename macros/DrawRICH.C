void DrawHMPID()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewHMPID.C");
   gInterpreter->ProcessLine("ViewHMPID()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   //   gMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 120, 0, 10, 0, .03, .03);
   gMC->Gdhead(1111, "Ring Imaging Cherenkov");
   gMC->Gdman(16, 6, "MAN");
}
