void ViewPHOS()
{
  // Set drawing attributes for PHOS version #2

  gMC->Gsatt("PHOS","seen",0);
  gMC->Gsatt("EMCA","seen",0);
  gMC->Gsatt("PTXW","seen",0);
  gMC->Gsatt("PUFP","seen",1);
  gMC->Gsatt("PAIR","seen",1);
  gMC->Gsatt("PTCB","seen",1);
  gMC->Gsatt("PCBL","seen",0);
  gMC->Gsatt("PROW","seen",0);
  gMC->Gsatt("PCEL","seen",0);
  gMC->Gsatt("PSTC","seen",0);
  gMC->Gsatt("PPAP","seen",0);
  gMC->Gsatt("PXTL","seen",0);
  gMC->Gsatt("PSUP","seen",0);
  gMC->Gsatt("PPIN","seen",0);
  gMC->Gsatt("PUCP","seen",1);
  gMC->Gsatt("PASP","seen",1);
  gMC->Gsatt("PTIP","seen",1);
  gMC->Gsatt("PTXP","seen",1);
  gMC->Gsatt("PPSD","seen",0);
  gMC->Gsatt("MPPS","seen",1);
  gMC->Gsatt("TLPS","seen",1);
  gMC->Gsatt("UPPS","seen",1);
  gMC->Gsatt("ANPS","seen",0);
  gMC->Gsatt("GGPS","seen",0);
  gMC->Gsatt("GROW","seen",0);
  gMC->Gsatt("GCEL","seen",0);
  gMC->Gsatt("CAPS","seen",1);
  gMC->Gsatt("PCPS","seen",1);
  gMC->Gsatt("LPPS","seen",1);
  gMC->Gsatt("UAPP","seen",0);
  gMC->Gsatt("LCPP","seen",0);
  gMC->Gsatt("LAPP","seen",1);
}
