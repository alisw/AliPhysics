void DrawCRT()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewCRT.C");
   gInterpreter->ProcessLine("ViewCRT()");
   //gMC->Gdopt("proj", "pers");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("ALIC", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 0, 90, 0, 10, 9.5, .009, .009);
   gMC->Gdhead(1111, "View of CRT (ACORDE)");
   gMC->Gdman(18, 4, "MAN");
}
