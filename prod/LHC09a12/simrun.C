// #define VERBOSEARGS
// simrun.C
{
  // extract the run and event variables given with 
  // --run <x> --event <y> --type <t>
  // where "type" can be "ppMBias" or "pptrg2mu"
  int nrun = 0;
  int nevent = 0;
  int seed = 0;
  char sseed[1024];
  char srun[1024];
  char sevent[1024];
  char type[1024];
  sprintf(srun,"");
  sprintf(sevent,"");
  for (int i=0; i< gApplication->Argc();i++){
#ifdef VERBOSEARGS
    printf("Arg %d:  %s\n",i,gApplication->Argv(i));
#endif
    if (!(strcmp(gApplication->Argv(i),"--run")))
      nrun = atoi(gApplication->Argv(i+1));
    sprintf(srun,"%d",nrun);
    if (!(strcmp(gApplication->Argv(i),"--event")))
      nevent = atoi(gApplication->Argv(i+1));
    sprintf(sevent,"%d",nevent);
    if (!(strcmp(gApplication->Argv(i),"--type")))
      strcpy(type,gApplication->Argv(i+1));
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
  
  gSystem->Setenv("CONFIG_SEED",sseed);
  gSystem->Setenv("CONFIG_RUN_TYPE",type);
  gSystem->Setenv("DC_RUN",srun);

  gSystem->Setenv("ALIMDC_RAWDB1","./mdc1");
  gSystem->Setenv("ALIMDC_RAWDB2","./mdc2");
  gSystem->Setenv("ALIMDC_TAGDB","./mdc1/tag");
  gSystem->Setenv("ALIMDC_RUNDB","./mdc1/meta");
  
  gSystem->Exec("cp $ALICE_ROOT/.rootrc .rootrc");
  gSystem->Exec("aliroot -b -q sim.C 2>&1 | tee sim.log");
  gSystem->Exec("mkdir generated");
  gSystem->Exec("mv *.root generated");
  gSystem->Exec("ln -s generated/geometry.root");
  gSystem->Exec("ln -s generated/raw.root");
  gSystem->Exec("aliroot -b -q rec.C 2>&1 | tee rec.log");
  gSystem->Exec("mv galice.root galice_rec.root");
  gSystem->Exec("mv generated/*.root .");
  gSystem->Exec("aliroot -b -q tag.C 2>&1 | tee tag.log");
  gSystem->Exec("aliroot -b -q CheckESD.C 2>&1 | tee check.log");

}
