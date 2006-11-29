void DrawALICE()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   gROOT->LoadMacro("ViewABSO.C");gInterpreter->ProcessLine("ViewABSO()");
   gROOT->LoadMacro("ViewCASTOR.C");gInterpreter->ProcessLine("ViewCASTOR()");
   gROOT->LoadMacro("ViewDIPO.C");gInterpreter->ProcessLine("ViewDIPO()");
   gROOT->LoadMacro("ViewFMD.C");gInterpreter->ProcessLine("ViewFMD()");
   gROOT->LoadMacro("ViewHALL.C");gInterpreter->ProcessLine("ViewHALL()");
   gROOT->LoadMacro("ViewITS.C");gInterpreter->ProcessLine("ViewITS()");
   gROOT->LoadMacro("ViewMAG.C");gInterpreter->ProcessLine("ViewMAG()");
   gROOT->LoadMacro("ViewMUON.C");gInterpreter->ProcessLine("ViewMUON()");
   gROOT->LoadMacro("ViewPHOS.C");gInterpreter->ProcessLine("ViewPHOS()");
   gROOT->LoadMacro("ViewPIPE.C");gInterpreter->ProcessLine("ViewPIPE()");
   gROOT->LoadMacro("ViewPMD.C");gInterpreter->ProcessLine("ViewPMD()");
   gROOT->LoadMacro("ViewHMPID.C");gInterpreter->ProcessLine("ViewHMPID()");
   gROOT->LoadMacro("ViewSHIL.C");gInterpreter->ProcessLine("ViewSHIL()");
   gROOT->LoadMacro("ViewSTART.C");gInterpreter->ProcessLine("ViewSTART()");
   gROOT->LoadMacro("ViewTOF.C");gInterpreter->ProcessLine("ViewTOF()");
   gROOT->LoadMacro("ViewTPC.C");gInterpreter->ProcessLine("ViewTPC()");
   gROOT->LoadMacro("ViewTRD.C");gInterpreter->ProcessLine("ViewTRD()");
   gROOT->LoadMacro("ViewFRAME.C");gInterpreter->ProcessLine("ViewFRAME()");
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 3000, -3000, 3000, -3500, 3500);
   gMC->DefaultRange();
   //   gMC->Gdraw("alic", 40, 30, 0, 12, 7.5, .005, .005);
   gMC->Gdraw("alic", 40, 30, 0, 13, 8, .012, .012);
   gMC->Gdhead(1111, "ALICE Detector");
   gMC->Gdman(18, 2, "MAN");
}
