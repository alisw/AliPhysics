//void DrawSTART()
{
   geant3->Gsatt("*", "seen", -1);
   geant3->Gsatt("alic", "seen", 0);
   gROOT->Macro("ViewSTART.C");
   geant3->Gdopt("hide", "on");
   geant3->Gdopt("shad", "on");
   geant3->Gsatt("*", "fill", 7);
   geant3->SetClipBox(".");
   geant3->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
   geant3->DefaultRange();
   geant3->Gdraw("alic", 40, 30, 0, 6, 9, .08, .08);
   geant3->Gdhead(1111, "START Detector");
   geant3->Gdman(13, 9, "MAN");
}
