//void ViewTPC()
{
 gMC->Gsatt("TPC ","seen",0); 
 gMC->Gsatt("TGAS","seen",0); 
 gMC->Gsatt("TPSG","seen",0); 
 gMC->Gsatt("TPHV","seen",1); 
 gMC->Gsatt("TRCS","seen",1); 
 gMC->Gsatt("TSGA","seen",-1); 
 gMC->Gsatt("TRCL","seen",1); 
 gMC->Gsatt("TLGA","seen",-1); 
 gMC->Gsatt("TSWS","seen",1); 
 gMC->Gsatt("TPW1","seen",1); 
 gMC->Gsatt("TPW2","seen",1); 
 gMC->Gsatt("TPW3","seen",1); 
 gMC->Gsatt("TPW4","seen",1); 
 gMC->Gsatt("TSPI","seen",1); 
 gMC->Gsatt("TSP1","seen",0); 
 gMC->Gsatt("TSPO","seen",1); 
 gMC->Gsatt("TSP2","seen",0); 
 gMC->Gsatt("TSWH","seen",1); 
 gMC->Gsatt("TSW1","seen",0); 
 gMC->Gsatt("TCOV","seen",0); 
 gMC->Gsatt("TPOI","seen",1); 
 gMC->Gsatt("TPIV","seen",1); 
 gMC->Gsatt("TPVD","seen",1); 

}
