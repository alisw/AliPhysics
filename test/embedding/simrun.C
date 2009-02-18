// #define VERBOSEARGS
// simrun.C
{
// set job and simulation variables as :
// root.exe -b -q simrun.C  --run <x> --event <y> --process <kPythia6/kPhojet> --field <kNoField/k5kG> --energy <900/10000>

  int nrun = 0;
  int nevent = 0;
  int seed = 0;

  char sseed[1024];
  char srun[1024];
  char sevent[1024];
  char sprocess[1024];
  char sfield[1024];
  char senergy[1024];

  sprintf(srun,"");
  sprintf(sevent,"");
  sprintf(sprocess,"");
  sprintf(sfield,"");
  sprintf(senergy,"");

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

    if (!(strcmp(gApplication->Argv(i),"--field")))
      sprintf(sfield,gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--energy")))
      sprintf(senergy,gApplication->Argv(i+1));

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
  gSystem->Setenv("CONFIG_RUN_TYPE",sprocess); // kPythia6 or kPhojet
  gSystem->Setenv("CONFIG_FIELD",sfield);      // kNoField or k5kG
  gSystem->Setenv("CONFIG_ENERGY",senergy);    // 900 or 10000 (GeV)
  gSystem->Setenv("DC_RUN",srun); // Not used in Config.C
  gSystem->Setenv("DC_EVENT",sevent); // Not used in Config.C
  
// Needed to produce simulated RAW data
  gSystem->Setenv("ALIMDC_RAWDB1","./mdc1");
  gSystem->Setenv("ALIMDC_RAWDB2","./mdc2");
  gSystem->Setenv("ALIMDC_TAGDB","./mdc1/tag");
  gSystem->Setenv("ALIMDC_RUNDB","./mdc1/meta");
  cout    << "SIMRUN:: Run " << gSystem->Getenv("DC_RUN") << " Event " << gSystem->Getenv("DC_EVENT")
	  << endl;


  // Background simulation
  gSystem->Setenv("CONFIG_EMBEDDING","kBackground");

  cout<<">>>>> BACKGROUND SIMULATION <<<<<"<<endl;
  gSystem->Exec("mkdir BackgroundFull");
  gSystem->Exec("cp Config.C BackgroundFull/");
  gSystem->Exec("cp sim.C BackgroundFull/");
  gSystem->Exec("cp rec.C BackgroundFull/");
  gSystem->ChangeDirectory("BackgroundFull/");
  gSystem->Exec("aliroot -b -q 'sim.C(0)' > sim.log 2>&1");
  cout<<">>>>> BACKGROUND RECONSTRUCTION <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q 'rec.C(0)' > rec.log 2>&1");
  gSystem->ChangeDirectory("../");

  // Convert Raw to SDigits
  cout << ">>>>> CONVERTING RAW 2 SDIGITS <<<<<" << endl;
  gSystem->Exec("mkdir Background");
  gSystem->Exec("cp BackgroundFull/raw.root Background/");
  gSystem->Exec("cp BackgroundFull/AliESDs.root Background/");
  gSystem->Exec("cp -a BackgroundFull/GRP Background/");
  gSystem->Exec("cp sim.C Background/");
  gSystem->ChangeDirectory("Background/");
  gSystem->Exec("aliroot -b -q 'sim.C(4)' > sim.log 2>&1");
  gSystem->ChangeDirectory("../");
  gSystem->Exec("mkdir BackgroundSDigits");
  gSystem->Exec("cp Background/*SDigits.root BackgroundSDigits");
  gSystem->Exec("cp BackgroundFull/galice.root BackgroundSDigits/");
  gSystem->Exec("cp BackgroundFull/AliESDs.root BackgroundSDigits/");

  // Merged simulation
  gSystem->Setenv("CONFIG_EMBEDDING","kMerged");

  cout<<">>>>> MERGED SIMULATION <<<<<<"<< endl;
  gSystem->Exec("mkdir Merged");
  gSystem->Exec("cp Config.C Merged/");
  gSystem->Exec("cp sim.C Merged/");
  gSystem->Exec("cp rec.C Merged/");
  gSystem->ChangeDirectory("Merged/");
  gSystem->Exec("aliroot -b -q 'sim.C(1)' > sim.log 2>&1");
  cout<<">>>>> MERGED RECONSTRUCTION <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q 'rec.C(1)' > rec.log 2>&1");
  gSystem->ChangeDirectory("../");

  // Pure signal re-reconstruction
  gSystem->Setenv("CONFIG_EMBEDDING","kSignal");
  
  cout<<">>>>> SIGNAL SIMULATION <<<<<<"<< endl;
  gSystem->Exec("mkdir Signal");
  gSystem->Exec("cp Config.C Signal/");
  gSystem->Exec("cp sim.C Signal/");
  gSystem->Exec("cp rec.C Signal/");
  gSystem->Exec("cp Merged/*SDigits*.root Signal/");
  gSystem->Exec("cp Merged/galice.root Signal/");
  gSystem->Exec("cp Merged/Kinematics.root Signal/");
  gSystem->ChangeDirectory("Signal/");
  gSystem->Exec("aliroot -b -q 'sim.C(2)' > sim.log 2>&1");
  cout<<">>>>> SIGNAL RECONSTRUCTION <<<<<"<<endl;
  gSystem->Exec("aliroot -b -q 'rec.C(2)' > rec.log 2>&1");
  gSystem->ChangeDirectory("../");
  
  //  cout<<">>>>> TAG <<<<<"<<endl;
  //  gSystem->Exec("aliroot -b -q tag.C > tag.log 2>&1");
  //  cout<<">>>>> CHECK ESD <<<<<"<<endl;
  //  gSystem->Exec("aliroot -b -q CheckESD.C > check.log 2>&1");
  //  cout<<">>>>> AOD <<<<<"<<endl;
  //  gSystem->Exec("aliroot -b -q CreateAODfromESD.C > aod.log 2>&1");
}
