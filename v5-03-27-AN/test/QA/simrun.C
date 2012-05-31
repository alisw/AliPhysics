//#define VERBOSEARGS

{
  // set job and simulation variables as :
  // root run.C  --run <x> --event <y> --process <proc> --minhard <min> --maxhard <max> --minpt <minpt>
  
  int nrun = 0;
  int nevent = 0;
  int seed = 0;
  
  char sseed[1024];
  char srun[1024];
  char sevent[1024];
  char sprocess[1024];
  char sminpthard[1024];
  char smaxpthard[1024];
  char sminptgammapi0[1024];
  char squench[1024];
  char sqhat[1024];

  sprintf(srun,"");
  sprintf(sevent,"");
  sprintf(sprocess,"");
  sprintf(sminpthard,"");
  sprintf(smaxpthard,"");
  sprintf(sminptgammapi0,"");
  sprintf(squench,"");
  sprintf(sqhat,"");

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
    
    if (!(strcmp(gApplication->Argv(i),"--minhard")))
      sprintf(sminpthard,gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--maxhard")))
      sprintf(smaxpthard,gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--minpt")))
      sprintf(sminptgammapi0,gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--quench")))
      sprintf(squench,gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--qhat")))
      sprintf(sqhat,gApplication->Argv(i+1));

  }
  
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
  gSystem->Setenv("DC_RUN",srun);
  gSystem->Setenv("DC_EVENT",sevent);
  gSystem->Setenv("DC_RUN_TYPE",sprocess);//"kPyGammaJetPHOS");
  gSystem->Setenv("PTHARDMIN",sminpthard);//"20");
  gSystem->Setenv("PTHARDMAX",smaxpthard);//"30");

  gSystem->Setenv("PTGAMMAPI0MIN",sminptgammapi0);//"1");
  gSystem->Setenv("QUENCHING",squench);
  gSystem->Setenv("QHAT",sqhat);


  gSystem->Setenv("ECMS","14000");
  gSystem->Setenv("ALIMDC_RAWDB1","./mdc1");
  gSystem->Setenv("ALIMDC_RAWDB2","./mdc2");
  gSystem->Setenv("ALIMDC_TAGDB","./mdc1/tag");
  gSystem->Setenv("ALIMDC_RUNDB","./mdc1/meta");
  cout<< "SIMRUN:: Run " << gSystem->Getenv("DC_RUN") << " Event " << gSystem->Getenv("DC_EVENT") 
	  << " Process "    << gSystem->Getenv("DC_RUN_TYPE") 
	  << " minpthard " << gSystem->Getenv("PTHARDMIN") 
	  << " maxpthard " << gSystem->Getenv("PTHARDMAX") 
	  << " minpt "     << gSystem->Getenv("PTGAMMAPI0MIN") 
	  << endl;
  //gSystem->Exec("cp $ROOTSYS/etc/system.rootrc .rootrc");
  cout<<">>>>> SIMULATION <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q sim.C > sim.log 2>&1");
  cout<<">>>>> SIMULATION QA <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q simqa.C > simqa.log 2>&1");
  cout<<">>>>> RECONSTRUCTION <<<<<"<<endl;
  gSystem->Exec("rm galice.root");
  gSystem->Exec("Aliroot -b -q rec.C > rec.log 2>&1");
  cout<<">>>>> RECONSTRUCTION QA <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q recqa.C > recqa.log 2>&1");
  cout<<">>>>> TAG <<<<<"<<endl;
  if( gSystem->Getenv("ALIEN_JDL_OUTPUTDIR"))
    gSystem->Exec("aliroot -b -q tag.C > tag.log 2>&1");
  cout<<">>>>> CHECK ESD <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q CheckESD.C > check.log 2>&1");

}
