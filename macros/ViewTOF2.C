//void ViewTOF2()
{
   gMC->Gsatt("FBAR","seen",0);
   gMC->Gsatt("FTOF","seen",0);
   gMC->Gsatt("FSTP","seen",0);
   gMC->Gsatt("FPET","seen",0);
   gMC->Gsatt("FST2","seen",1);
   gMC->Gsatt("FLT2","seen",0);
   gMC->Gsatt("FST3","seen",1);
   gMC->Gsatt("FLT3","seen",0);
   gMC->Gsatt("FST4","seen",1);
   gMC->Gsatt("FLT4","seen",1);
   gMC->Gsatt("FPG1","seen",0);
   gMC->Gsatt("FPG2","seen",0);
   gMC->Gsatt("FGP1","seen",0);
   gMC->Gsatt("FGP2","seen",0);
}
