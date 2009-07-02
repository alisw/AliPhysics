//#define VERBOSEARGS

{
	// set job and simulation variables as :
	// root run.C  --run <x> --event <y> --process <proc> --minhard <min> --maxhard <max> --minpt <minpt> ...
	
	int nrun = 0;
	int nevent = 0;
	int seed = 0;
	
	char sseed[1024];
	char srun[1024];
	char sevent[1024];
	char secms[1024];
	char sprocess[1024];
	char sminpthard[1024];
	char smaxpthard[1024];
	char sminptgammapi0[1024];
	char spi0gammafrag[1024];
	char squench[1024];
	char sqhat[1024];
	char smedlength[1024];
	char strig[1024];
	char syear[1024];
	
	sprintf(srun,"");
	sprintf(sevent,"");
	sprintf(secms,"");
	sprintf(sprocess,"");
	sprintf(sminpthard,"");
	sprintf(smaxpthard,"");
	sprintf(sminptgammapi0,"");
	sprintf(spi0gammafrag,"");
	sprintf(squench,"");
	sprintf(sqhat,"");
	//sprintf(smedlength,"");
	sprintf(strig,"");
	sprintf(syear,"");
	
	for (int i=0; i< gApplication->Argc();i++){
#ifdef VERBOSEARGS
		printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
#endif
		if (!(strcmp(gApplication->Argv(i),"--run")))
			nrun = atoi(gApplication->Argv(i+1));
		sprintf(srun,"%d",nrun);
	
		if (!(strcmp(gApplication->Argv(i),"--event")))
			nevent = atoi(gApplication->Argv(i+1));
		sprintf(sevent,"%d",nevent);
				
		if (!(strcmp(gApplication->Argv(i),"--process")))
			sprintf(sprocess, gApplication->Argv(i+1));
				
		if (!(strcmp(gApplication->Argv(i),"--year")))
			sprintf(syear, gApplication->Argv(i+1));

		if (!(strcmp(gApplication->Argv(i),"--ecms")))
			sprintf(secms,gApplication->Argv(i+1));
						
		if (!(strcmp(gApplication->Argv(i),"--minhard")))
			sprintf(sminpthard,gApplication->Argv(i+1));
							
		if (!(strcmp(gApplication->Argv(i),"--maxhard")))
			sprintf(smaxpthard,gApplication->Argv(i+1));
	
		if (!(strcmp(gApplication->Argv(i),"--pi0gammafrag")))
			sprintf(spi0gammafrag,gApplication->Argv(i+1));
																						
		if (!(strcmp(gApplication->Argv(i),"--minpt")))
			sprintf(sminptgammapi0,gApplication->Argv(i+1));
																	
		if (!(strcmp(gApplication->Argv(i),"--quench")))
			sprintf(squench,gApplication->Argv(i+1));
										
		if (!(strcmp(gApplication->Argv(i),"--qhat")))
			sprintf(sqhat,gApplication->Argv(i+1));
											
		//if (!(strcmp(gApplication->Argv(i),"--medlength")))
		//	sprintf(smedlength,gApplication->Argv(i+1));
											
		if (!(strcmp(gApplication->Argv(i),"--trigger")))
			sprintf(strig,gApplication->Argv(i+1));
													
	}
	
	//rec params for reconstruction
	//only for release < 16
	//char cmd[200] ; 
	//sprintf(cmd, ".! tar zxvf DBpp.tgz") ; 
	//gROOT->ProcessLine(cmd) ; 
  	
	seed = nrun * 100000 + nevent;
	sprintf(sseed,"%d",seed);
	
	if (seed==0) {
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stderr,"!!!!  WARNING! Seeding variable for MC is 0          !!!!\n");
		fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	} else {
		fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		fprintf(stdout,"!!!  MC Seed is %d \n",seed);
		fprintf(stdout,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	
	// set the seed environment variable
	gSystem->Setenv("CONFIG_SEED",sseed);
	gSystem->Setenv("RUN",srun);
	gSystem->Setenv("EVENT",sevent);
	gSystem->Setenv("DC_RUN_TYPE",sprocess);//"kPyGammaJetPHOS");
	gSystem->Setenv("YEAR",syear);//"14000");	
	gSystem->Setenv("ECMS",secms);//"14000");
	gSystem->Setenv("PTHARDMIN",sminpthard);//"20");
	gSystem->Setenv("PTHARDMAX",smaxpthard);//"30");
	gSystem->Setenv("PI0GAMMAINDET",spi0gammafrag);//"1");
	gSystem->Setenv("PTGAMMAPI0MIN",sminptgammapi0);//"1");
	gSystem->Setenv("QUENCHING",squench);
	gSystem->Setenv("QHAT",sqhat);
	//gSystem->Setenv("MEDLENGTH",smedlength);
	gSystem->Setenv("TRIGGER",strig);
	
	cout<< "SIMRUN:: Run = " << gSystem->Getenv("RUN") << "; Event = " << gSystem->Getenv("EVENT") 
	<< "; Process = "    << gSystem->Getenv("DC_RUN_TYPE")
	<< "; E cms = " << gSystem->Getenv("ECMS") 
	<< "; Trigger configuration = "<< gSystem->Getenv("TRIGGER")
	<< "; Year = "<<gSystem->Getenv("YEAR")	<< endl;
	
	cout<< "         minpthard = " << gSystem->Getenv("PTHARDMIN") 
	<< "; maxpthard =" << gSystem->Getenv("PTHARDMAX") 
	<< "; force pi0/gamma frag in detector = "     << gSystem->Getenv("PI0GAMMAINDET") 
	<< "; minpt = "     << gSystem->Getenv("PTGAMMAPI0MIN") 
	<< "; quenching option ="<<gSystem->Getenv("QUENCHING") 
	<< "; qhat = "<<gSystem->Getenv("QHAT") 
	//<< "; medium length = "<<gSystem->Getenv("MEDIUMLENGTH")
	<< endl;
	
	cout<<">>>>> SIMULATION <<<<<"<<endl;
	gSystem->Exec("aliroot -b -q sim.C > sim.log 2>&1");
	cout<<">>>>> RECONSTRUCTION <<<<<"<<endl;
	gSystem->Exec("aliroot -b -q rec.C > rec.log 2>&1");
	cout<<">>>>> TAG <<<<<"<<endl;
	gSystem->Exec("aliroot -b -q tag.C > tag.log 2>&1");
	cout<<">>>>> CHECK ESD <<<<<"<<endl;
	gSystem->Exec("aliroot -b -q CheckESD.C > check.log 2>&1");
	
}
