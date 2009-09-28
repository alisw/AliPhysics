//#define VERBOSEARGS
{
  //Process environmental variables from command line:
  char anainputdata[1024];
  char smode[1024];
  char sconfig1[1024];
  char sconfig2[1024];
  char sconfig3[1024];
  char sevent[1024];

  sprintf(anainputdata,"");
  sprintf(smode,"");
  sprintf(sconfig1,"");
  sprintf(sconfig2,"");
  sprintf(sconfig3,"");
  sprintf(sevent,"");


  for (int i=0; i< gApplication->Argc();i++){
#ifdef VERBOSEARGS
    printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
#endif

    if (!(strcmp(gApplication->Argv(i),"--anaInputData")))
       sprintf(anainputdata,"%s",gApplication->Argv(i+1));	
    if (!(strcmp(gApplication->Argv(i),"--mode")))
       sprintf(smode,"%s",gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i),"--config1")))
       sprintf(sconfig1,"%s",gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i),"--config2")))
       sprintf(sconfig2,"%s",gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i),"--config3")))
       sprintf(sconfig3,"%s",gApplication->Argv(i+1));
    if (!(strcmp(gApplication->Argv(i),"--sevent")))
       sprintf(sevent,"%s",gApplication->Argv(i+1));
  }
	
    gSystem->Setenv("anaInputData", anainputdata);
    gSystem->Setenv("MODE"  ,  smode);
    gSystem->Setenv("CONFIG1", sconfig1);
    gSystem->Setenv("CONFIG2", sconfig2);
    gSystem->Setenv("CONFIG3", sconfig3);
    gSystem->Setenv("SEVENT" , sevent);

	
    printf("Run.C: Variables: anaInputData %s, mode %s, config1 %s, config2 %s, config3 %s, first event %s\n", anainputdata,smode,sconfig1,sconfig2,sconfig3,sevent);

    gSystem->Exec("root -b -q -l anaJete.C > anaJete.log 2>&1");

}
