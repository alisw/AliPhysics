void ViewTOF()
{
   AliDetector *TOF  = gAlice->GetDetector("TOF");
   gMC->Gsatt("*","seen",0);
   gMC->Gsatt("BTO1","seen",0);
   gMC->Gsatt("BTO2","seen",0);
   gMC->Gsatt("BTO3","seen",0);
   gMC->Gsatt("FLT1","seen",1);
   gMC->Gsatt("FLT2","seen",1);
   gMC->Gsatt("FLT3","seen",1);
//   gMC->Gsatt("BTO1","fill",7);
//   gMC->Gsatt("BTO2","fill",7);
//   gMC->Gsatt("BTO3","fill",7);
//   gMC->Gsatt("FLT1","fill",7);
//   gMC->Gsatt("FLT2","fill",7);
//   gMC->Gsatt("FLT3","fill",7);
   gMC->Gsatt("FSTR","seen",1);
//   gMC->Gsatt("FPAD","seen",0);
}
