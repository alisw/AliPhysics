void DrawEMCAL(){
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewEMCAL.C");
   gInterpreter->ProcessLine("ViewEMCAL()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 600, -90, 2000, 0, 1000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 8, 7.5, .02, .02);
   gMC->Gdhead(1111, "Large E&M Calorimeter Detector");
   gMC->Gdman(16, 4, "MAN");
   gMC->Gdman(2, 4, "WM2");
}
