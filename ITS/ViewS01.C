void ViewS01(){
    gMC->Gsatt("ITSV","seen",0); // Air "Gen Air"
      gMC->Gsatt("ITSD","seen",0); // Air "Gen Air"
        gMC->Gsatt("IS01","seen",0); //  Air "Air"
	  gMC->Gsatt("ICY2","seen",1); // Carbon Fiber "SDD C (M55J)
	  gMC->Gsatt("I212","seen",1); // Carbon "ITS SandW C"
	  gMC->Gsatt("I211","seen",1); // PC board "G10FR4"
	  gMC->Gsatt("I217","seen",1); // Carbon Fiber "SDD/SSD rings"
	  gMC->Gsatt("I219","seen",1); // Carbon Fiber "SDD/SSD rings"
	  gMC->Gsatt("I214","seen",1); // PC board "G10FR4"
	  gMC->Gsatt("I213","seen",1); // PC board "G10FR4"
	  gMC->Gsatt("I215","seen",2); // Air "Air"
	  gMC->Gsatt("I216","seen",1); // Carbon "ITS Sandw C"
}
