void ViewITSSPD(Int_t version1,Int_t version2){
    gMC->Gsatt("ITSV","seen",0);// Air
      gMC->Gsatt("ITSD","seen",0);// Air
        gMC->Gsatt("IT12","seen",0);// Air
	if(version1==10){gMC->Gsatt("I651","seen",1);gMC->Gsatt("I651","colo",7);} // Services
          gMC->Gsatt("I12B","seen",0);// Air
            gMC->Gsatt("I10B","seen",0);// Air
              gMC->Gsatt("I105","seen",1);gMC->Gsatt("I105","colo",7);// SPD End ladder
              gMC->Gsatt("I108","seen",1);gMC->Gsatt("I108","colo",7);// SPD bus
              gMC->Gsatt("I109","seen",1);gMC->Gsatt("I109","colo",7);// SPD bus
              gMC->Gsatt("I107","seen",0);//Air
                gMC->Gsatt("I106","seen",1);gMC->Gsatt("I106","colo",6);// Silicon
                gMC->Gsatt("I101","seen",1);gMC->Gsatt("I101","colo",6);// Silicon
                  gMC->Gsatt("ITS1","seen",1);gMC->Gsatt("ITS1","colo",6);// Silicon
            gMC->Gsatt("I20B","seen",0);// Air
              gMC->Gsatt("I105","seen",1);gMC->Gsatt("I105","colo",7);// SPD End ladder
              gMC->Gsatt("I109","seen",1);gMC->Gsatt("I109","colo",7);// SPD bus
              gMC->Gsatt("I108","seen",1);gMC->Gsatt("I108","colo",7);// SPD bus
              gMC->Gsatt("I107","seen",0);// Air
                gMC->Gsatt("I1D6","seen",1);gMC->Gsatt("I1D6","colo",6);// Silicon
                gMC->Gsatt("I1D1","seen",1);gMC->Gsatt("I1D1","colo",6);// Silicon
                  gMC->Gsatt("ITS2","seen",1);gMC->Gsatt("ITS2","colo",6);// Silicon
            gMC->Gsatt("I123","seen",1);gMC->Gsatt("I123","colo",1);// carbon
            gMC->Gsatt("I121","seen",1);gMC->Gsatt("I121","colo",1);// carbon
            gMC->Gsatt("I122","seen",1);gMC->Gsatt("I122","colo",1);// carbon
            gMC->Gsatt("I120","seen",1);gMC->Gsatt("I120","colo",1);// carbon
            gMC->Gsatt("I144","seen",1);gMC->Gsatt("I144","colo",1);// carbon
            gMC->Gsatt("I143","seen",1);gMC->Gsatt("I143","colo",1);// carbon
            gMC->Gsatt("I142","seen",1);gMC->Gsatt("I142","colo",1);// carbon
            gMC->Gsatt("I141","seen",1);gMC->Gsatt("I141","colo",1);// carbon
            gMC->Gsatt("I140","seen",1);gMC->Gsatt("I140","colo",1);// carbon
            gMC->Gsatt("I139","seen",1);gMC->Gsatt("I139","colo",1);// carbon
            gMC->Gsatt("I138","seen",1);gMC->Gsatt("I138","colo",1);// carbon
            gMC->Gsatt("I137","seen",1);gMC->Gsatt("I137","colo",1);// carbon
            gMC->Gsatt("I136","seen",1);gMC->Gsatt("I136","colo",1);// carbon
            gMC->Gsatt("I135","seen",1);gMC->Gsatt("I135","colo",1);// carbon
            gMC->Gsatt("I134","seen",1);gMC->Gsatt("I134","colo",1);// carbon
            gMC->Gsatt("I133","seen",1);gMC->Gsatt("I133","colo",1);// carbon
            gMC->Gsatt("I132","seen",1);gMC->Gsatt("I132","colo",1);// carbon
            gMC->Gsatt("I131","seen",1);gMC->Gsatt("I131","colo",1);// carbon
            gMC->Gsatt("I130","seen",1);gMC->Gsatt("I130","colo",1);// carbon
            gMC->Gsatt("I129","seen",1);gMC->Gsatt("I129","colo",1);// carbon
            gMC->Gsatt("I128","seen",1);gMC->Gsatt("I128","colo",1);// carbon
            gMC->Gsatt("I126","seen",1);gMC->Gsatt("I126","colo",1);// carbon
            gMC->Gsatt("I125","seen",1);gMC->Gsatt("I125","colo",1);// carbon
            gMC->Gsatt("I124","seen",1);gMC->Gsatt("I124","colo",1);// carbon
            gMC->Gsatt("I113","seen",0);//Air
              gMC->Gsatt("I112","seen",1);gMC->Gsatt("I112","colo",1);// carbon
              gMC->Gsatt("I111","seen",1);gMC->Gsatt("I111","colo",1);// carbon
              gMC->Gsatt("I118","seen",1);gMC->Gsatt("I118","colo",2);// glue
              gMC->Gsatt("I110","seen",1);gMC->Gsatt("I110","colo",1);// carbon
              gMC->Gsatt("I114","seen",1);gMC->Gsatt("I114","colo",3);// INOX
              gMC->Gsatt("I115","seen",1);gMC->Gsatt("I115","colo",4);// water
              gMC->Gsatt("I116","seen",1);gMC->Gsatt("I116","colo",3);// INOX
                gMC->Gsatt("I117","seen",1);gMC->Gsatt("I117","colo",4);// Water
          gMC->Gsatt("I650","seen",0); //air
	    gMC->Gsatt("I666","seen",1);gMC->Gsatt("I666","colo",5);//Aluminum
            gMC->Gsatt("I667","seen",1);gMC->Gsatt("I667","colo",5);//Aluminum
              gMC->Gsatt("I668","seen",1);gMC->Gsatt("I668","colo",4);// water
            gMC->Gsatt("I669","seen",1);gMC->Gsatt("I669","colo",3);// INOX
              gMC->Gsatt("I670","seen",1);gMC->Gsatt("I670","colo",4);// water
            gMC->Gsatt("I671","seen",1);gMC->Gsatt("I671","colo",6);// Silicon
              gMC->Gsatt("I672","seen",1);gMC->Gsatt("I672","colo",4);// water
            gMC->Gsatt("I673","seen",1);gMC->Gsatt("I673","colo",6);// Silicon
              gMC->Gsatt("I674","seen",1);gMC->Gsatt("I674","colo",4);// water
              gMC->Gsatt("I675","seen",1);gMC->Gsatt("I675","colo",5);//Aluminum
            gMC->Gsatt("I676","seen",1);gMC->Gsatt("I676","colo",6);// Silicon
              gMC->Gsatt("I677","seen",1);gMC->Gsatt("I677","colo",4);// water
              gMC->Gsatt("I678","seen",1);gMC->Gsatt("I678","colo",5);//Aluminum
}
