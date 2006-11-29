void DrawACORDE()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewACORDE.C");
   gInterpreter->ProcessLine("ViewACORDE()");
   //gMC->Gdopt("proj", "pers");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("ALIC", 0, 3000, -3000, 3000, -6000, 6000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 30, 40, 0, 10, 9.5, .0095, .0095);
   //gMC->Gdraw("alic", 30, 40, 0, -30, -60, .09, .09); //Zoom right side
   //gMC->Gdraw("alic", 30, 30, 0, 10, -10, .03, .03); // Zoom right side
   gMC->Gdhead(1111, "View of ACORDE (ACORDE)");
   gMC->Gdman(18, 4, "MAN");
}
