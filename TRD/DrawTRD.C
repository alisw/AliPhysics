void DrawTRD()
{
   gMC->Gsatt("*", "seen", -1);
   gMC->Gsatt("alic", "seen", 0);
   AliTRD *TRD = (AliTRD *) gAlice->GetModule("TRD");
   AliTRDgeometry *Geo = TRD->GetGeometry(); 
   gROOT->LoadMacro("ViewTRD.C");
   if (Geo->IsVersion() == 0) {
     gInterpreter->ProcessLine("ViewTRDhole()");
   }
   else {
     if (Geo->GetPHOShole()) {
       if (Geo->GetRICHhole()) {
         gInterpreter->ProcessLine("ViewTRDfull3()");
       }
       else {
         gInterpreter->ProcessLine("ViewTRDfull1()");
       }
     }
     else {
       if (Geo->GetRICHhole()) {
         gInterpreter->ProcessLine("ViewTRDfull2()");
       }
       else {
         gInterpreter->ProcessLine("ViewTRDfull0()");
       }
     }
   }
   gMC->Gdopt("hide", "on");
   gMC->Gdopt("shad", "on");
   gMC->Gsatt("*", "fill", 7);
   gMC->SetClipBox(".");
   gMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
   gMC->DefaultRange();
   gMC->Gdraw("alic", 40, 30, 0, 12, 9.4, .020, .020);
   gMC->Gdhead(1111, "Transition Radiation Detector");
   gMC->Gdman(18, 4, "MAN");
}
