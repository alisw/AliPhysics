void ViewTOF1()
{
   gMC->Gsatt("FBAR","seen",0);
   gMC->Gsatt("FGP1","seen",1);
   gMC->Gsatt("FGP2","seen",1);
   gMC->Gsatt("FPL1","seen",1);
   gMC->Gsatt("FPL2","seen",1);
   gMC->Gsatt("FPL3","seen",1);
   gMC->Gsatt("FPL4","seen",1);
}
