void ViewCASTOR()
{
   gMC->Gsatt("OCTA","seen",0);
   gMC->Gsatt("EM  ","seen",0);
   gMC->Gsatt("HAD ","seen",0);
   gMC->Gsatt("CAL ","seen",0);
   gMC->Gsatt("CALT","seen",1);
   gMC->Gsatt("OCT ","seen",0);
   gMC->Gsatt("SLEM","seen",1);
   gMC->Gsatt("SLHA","seen",1);
   gMC->Gsatt("SAHA","seen",1);
}
