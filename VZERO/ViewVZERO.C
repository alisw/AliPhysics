void ViewVZERO()
{
     gMC->Gsatt("v0RI","seen",0);
     
     gMC->Gsatt("v0CA","seen",0);
     gMC->Gsatt("v0IR","seen",1);
     gMC->Gsatt("v0ER","seen",1);
     
     gMC->Gsatt("v0R0","seen",0);     
     gMC->Gsatt("v0R1","seen",1);
     gMC->Gsatt("v0R2","seen",1);
     gMC->Gsatt("v0R3","seen",1);
     gMC->Gsatt("v0R4","seen",1); 
     gMC->Gsatt("V0R5","seen",1); 
     gMC->Gsatt("V0R6","seen",1);
     
     gMC->Gsatt("v0LE","seen",0);
     
     gMC->Gsatt("v0L0","seen",0);
     gMC->Gsatt("v0L1","seen",1);
     gMC->Gsatt("v0L2","seen",1);
     gMC->Gsatt("v0L3","seen",1);
     gMC->Gsatt("v0L4","seen",1);
}

