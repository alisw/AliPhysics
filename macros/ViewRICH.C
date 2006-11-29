void ViewHMPID()
{
   gMC->Gsatt("HMPID","seen",0);
   gMC->Gsatt("SRIC","seen",0);
   gMC->Gsatt("HONE","seen",1);
   gMC->Gsatt("ALUM","seen",1);
   gMC->Gsatt("QUAR","seen",1);
   gMC->Gsatt("SPAC","seen",1);
   gMC->Gsatt("OQUA","seen",1);
   gMC->Gsatt("BARR","seen",1);
   gMC->Gsatt("META","seen",1);
   gMC->Gsatt("GAP ","seen",1);
   gMC->Gsatt("CSI ","seen",1);
   gMC->Gsatt("GRID","seen",1);
}
