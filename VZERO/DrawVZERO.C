void DrawVZERO()

{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewVZERO.C");
   gInterpreter->ProcessLine("ViewVZERO()");

   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");

   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 13, 11, .05, .05);
   gMC->Gdhead(1111, "VZERO detector");
   gMC->Gdman(16, 6, "MAN");

}
