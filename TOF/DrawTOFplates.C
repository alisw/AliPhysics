void DrawTOFplates()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewTOFplates.C");
   gInterpreter->ProcessLine("ViewTOFplates()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 5);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 1000, 0, 1000, 0, 1000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 45, 40, 0, 10, 10, .015, .015);
   gMC->Gdhead(1111, "TOF Detector (Plates)");
   gMC->Gdman(18, 4, "MAN");             
}
