//void DrawTRD()
{
   geant3->Gsatt("*", "seen", -1);
   geant3->Gsatt("alic", "seen", 0);
   AliTRD *TRD = gAlice->GetModule("TRD");
   if      (TRD->IsVersion() == 0)
     gROOT->Macro("ViewTRD0.C");
   else if (TRD->IsVersion() == 1)
     gROOT->Macro("ViewTRD1.C");
   else if (TRD->IsVersion() == 2)
     gROOT->Macro("ViewTRD2.C");
   geant3->Gdopt("hide", "on");
   geant3->Gdopt("shad", "on");
   geant3->Gsatt("*", "fill", 7);
   geant3->SetClipBox(".");
   geant3->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
   geant3->DefaultRange();
   geant3->Gdraw("alic", 40, 30, 0, 12, 9.4, .020, .020);
   geant3->Gdhead(1111, "Transition Radiation Detector");
   geant3->Gdman(18, 4, "MAN");
}
