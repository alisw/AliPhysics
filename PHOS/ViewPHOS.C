//void ViewPHOS()
{
  // Set drawing attributes for PHOS version #2

  // PHOS
  // 1st level
   gMC->Gsatt("PHOS","seen",1);

   // 2nd level
   gMC->Gsatt("PTXW","seen",1);

   // 3rd level
   gMC->Gsatt("PUFP","seen",1);
   gMC->Gsatt("PAIR","seen",0);

   // 4th level
   gMC->Gsatt("PTCB","seen",1);
   gMC->Gsatt("PUCP","seen",1);
   gMC->Gsatt("PASP","seen",1);
   gMC->Gsatt("PTIP","seen",1);
   gMC->Gsatt("PTXP","seen",1);

   // 5th level
   gMC->Gsatt("PCBL","seen",-2);

   // 6th level
   gMC->Gsatt("PROW","seen",0);

   // 7th level
   gMC->Gsatt("PCEL","seen",0);

   // 8th level
   gMC->Gsatt("PSTC","seen",0);

   // 9th level
   gMC->Gsatt("PPAP","seen",0);

   // 10th level
   gMC->Gsatt("PXTL","seen",0);
   gMC->Gsatt("PSUP","seen",0);

   // 11th level
   gMC->Gsatt("PPIN","seen",0);

   // PHOS CPV
   // 1st level
   gMC->Gsatt("PCPV","seen",-2);

   // 2nd level
   gMC->Gsatt("PCRO","seen",0);

   // 3rd level
   gMC->Gsatt("PCCE","seen",0);

}
