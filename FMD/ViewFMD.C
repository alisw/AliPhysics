void ViewFMD()
{
   gMC->Gsatt("FMD1","seen",0);
   gMC->Gsatt("FMD2","seen",0);
   gMC->Gsatt("FMD3","seen",0);
   gMC->Gsatt("RNGI","seen",0);
   gMC->Gsatt("RNGO","seen",0);
   gMC->Gsatt("VFI","seen",0);
   gMC->Gsatt("VFO","seen",0);
   gMC->Gsatt("VBI","seen",0);
   gMC->Gsatt("VBO","seen",0);
   gMC->Gsatt("ACTI","seen",1);
   gMC->Gsatt("ACTO","seen",1);
   gMC->Gsatt("PBTI","seen",1);
   gMC->Gsatt("PBTO","seen",1);
   gMC->Gsatt("PBBI","seen",1);
   gMC->Gsatt("PBBO","seen",1);
}
