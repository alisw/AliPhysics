void DrawTPC()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewTPC.C");
   gInterpreter->ProcessLine("ViewTPC()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("TPMW",-300,300,-300,300,254.,270.);
   gMC->SetClipBox("TESR",-300,300,-300,300,254.,270.);
   gMC->SetClipBox("TSSW",-300,300,-300,300,283.,284.);
   gMC->SetClipBox("TSWC",-300,300,-300,300,283.,284.);  

   gMC->SetClipBox("*", 0, 300, -300, 300, -290, 290);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .025, .025);
   gMC->Gdhead(1111, "Time Projection Chamber");
   gMC->Gdman(18, 4, "MAN");
}
