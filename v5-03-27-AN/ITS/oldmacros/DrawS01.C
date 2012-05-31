void DrawS01()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewS01.C");
   gInterpreter->ProcessLine("ViewS01()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 100, -1000, 1000, -1000, 1000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 10, 10, .1, .1);
   gMC->Gdhead(1111, "Inner Tracking System");
   gMC->Gdman(16, 6, "MAN");
}
