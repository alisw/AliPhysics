void ViewFMD()
{
   gMC->Gsatt("FMD1","seen",0);
   gMC->Gsatt("FMD2","seen",0);
   gMC->Gsatt("FMD3","seen",0);
   gMC->Gsatt("FRGI","seen",0);
   gMC->Gsatt("FRGO","seen",0);
   gMC->Gsatt("FVFI","seen",0);
   gMC->Gsatt("FVFO","seen",0);
   gMC->Gsatt("FVBI","seen",0);
   gMC->Gsatt("FVBO","seen",0);
   gMC->Gsatt("FACI","seen",1);
   gMC->Gsatt("FACO","seen",1);
   gMC->Gsatt("FPTI","seen",1);
   gMC->Gsatt("FPTO","seen",1);
   gMC->Gsatt("FPBI","seen",1);
   gMC->Gsatt("FPBO","seen",1);
}
