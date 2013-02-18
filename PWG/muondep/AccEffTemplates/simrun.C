// #define VERBOSEARGS
// simrun.C
{
// set job and simulation variables as :
// root.exe -b -q simrun.C  --run <x> --chunk <y> --event <n> 

  int nrun = 0;
  int nchunk = 0;
  int nevent = 0;
  int seed = 0;
  Bool_t snapshot(kFALSE);
  
  char sseed[1024];
  char srun[1024];
  char schunk[1024];

  sprintf(srun,"");
  sprintf(schunk,"");

  for (int i=0; i< gApplication->Argc();i++){
#ifdef VERBOSEARGS
    printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
#endif
    if (!(strcmp(gApplication->Argv(i),"--run")))
      nrun = atoi(gApplication->Argv(i+1));
    sprintf(srun,"%d",nrun);

    if (!(strcmp(gApplication->Argv(i),"--chunk")))
      nchunk = atoi(gApplication->Argv(i+1));
    sprintf(schunk,"%d",nchunk);

    if (!(strcmp(gApplication->Argv(i),"--event")))
      nevent = atoi(gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--snapshot")))
      snapshot = kTRUE;
  }

  seed = nrun * 10000 + nchunk;
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
  gSystem->Setenv("DC_RUN",srun); // used in AliSimulation.cxx
  gSystem->Setenv("DC_EVENT",schunk); // Not used in Config.C (used for setting seed)
  
// Needed to produce simulated RAW data
  gSystem->Setenv("ALIMDC_RAWDB1","./mdc1");
  gSystem->Setenv("ALIMDC_RAWDB2","./mdc2");
  gSystem->Setenv("ALIMDC_TAGDB","./mdc1/tag");
  gSystem->Setenv("ALIMDC_RUNDB","./mdc1/meta");
  cout<< "SIMRUN:: Run " << gSystem->Getenv("DC_RUN")
          << " Chunk " << gSystem->Getenv("DC_EVENT")
	  << endl;

  cout<<">>>>> SIMULATION <<<<<"<<endl;
  if ( snapshot )
  {
    gSystem->Setenv("OCDB_SNAPSHOT_CREATE","kTRUE");
    gSystem->Setenv("OCDB_SNAPSHOT_FILENAME","OCDB_sim.root");
  }
  
  gSystem->Exec(Form("aliroot -b -q sim.C\\(%d\\) > sim.log 2>&1",nevent));
  gSystem->Exec("mv syswatch.log simwatch.log");

  if ( snapshot )
  {
    gSystem->Setenv("OCDB_SNAPSHOT_FILENAME","OCDB_rec.root");
  }
  
  cout<<">>>>> RECONSTRUCTION <<<<<"<<endl;
  
  if ( snapshot )
  {
    // for some reason must include ITS objects in the snapshot
    // (to be able to instantiante the vertexer later on ?)
    gSystem->Exec("aliroot -b -q rec.C\\(1\\) > rec.log 2>&1");
  }
  else
  {
    gSystem->Exec("aliroot -b -q rec.C > rec.log 2>&1");
  }
  
  gSystem->Exec("mv syswatch.log recwatch.log");

  if ( snapshot )
  {
    gSystem->Exec(Form("mkdir -p OCDB/%s",srun));
    gSystem->Exec(Form("mv OCDB_*.root OCDB/%s/",srun));
  }
  else
  {
    //  cout<<">>>>> TAG <<<<<"<<endl;
    //  gSystem->Exec("aliroot -b -q tag.C > tag.log 2>&1");
  
    cout<<">>>>> CHECK ESD <<<<<"<<endl;
    gSystem->Exec("aliroot -b -q CheckESD.C > checkesd.log 2>&1");
    cout<<">>>>> MAKE AOD <<<<<"<<endl;
    gSystem->Exec("aliroot -b -q AODtrain.C > aod.log 2>&1");
    cout<<">>>>> CHECK AOD <<<<<"<<endl;
    gSystem->Exec("aliroot -b -q CheckAOD.C > checkaod.log 2>&1");
  }
  
  return 0;
}
