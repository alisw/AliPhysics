void DrawSTART()
{

 gMC->Gsatt("*", "seen", -1);
 gMC->Gsatt("alic", "seen", 0);
 
 gROOT->LoadMacro("ViewSTART.C");
 gInterpreter->ProcessLine("ViewSTART()");
 gMC->Gdopt("hide", "on");
 gMC->Gdopt("shad", "on");
 gMC->Gsatt("*", "fill", 7);
 //   gMC->SetClipBox(".");
 // gMC->SetClipBox("*", 0, 300, -300, 300, -300, 300);
 //   gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
 gMC->DefaultRange();
 gMC->Gdraw("ALIC", 40, 30, 0, -10, 2, 0.5, 0.5);
 gMC->Gdhead(1111, "START Detector");
 //  gMC->Gdman(13, 9, "MAN");
}
