void DrawSTART()
{

 gMC->Gsatt("*", "seen", -1);
 gMC->Gsatt("alic", "seen", 0);
 
 gROOT->LoadMacro("ViewSTART.C");
 gInterpreter->ProcessLine("ViewSTART()");
 gMC->Gdopt("hide", "on");
 gMC->Gdopt("shad", "on");
 gMC->Gsatt("*", "fill", 7);
 gMC->SetClipBox(".");
 gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
 gMC->DefaultRange();
 gMC->Gdraw("ALIC", 40, 30, 0, 6, 9, 0.05, 0.05);
 gMC->Gdhead(1111, "START Detector");
 gMC->Gdman(13, 9, "MAN");
}
