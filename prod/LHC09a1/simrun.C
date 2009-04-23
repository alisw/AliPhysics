//#define VERBOSEARGS
{
  // set job and simulation variables as : alienoutner needed for seed
  // root run.C  --run <x> --event <y> --aliencounter <counter> --process <proc> --minhard <min> --maxhard <max> --minpt <minpt>


  char* program = "aliroot";
  //  char* program = "alienaliroot";
  
  int nrun = 0;
  int nevent = 0;
  int seed = 0;
  int naliencounter = 0;
  
  char sseed[1024];
  char srun[1024];
  char sevent[1024];
  char saliencounter[1024];
  char sprocess[1024];
  char sminpthard[1024];
  char smaxpthard[1024];
  char sminptgammapi0[1024];
  char squench[1024];
  char sqhat[1024];

  Int_t npthardbin = 0;
  const Int_t nPtHardBins = 3;
  const Int_t ptHardLo[nPtHardBins] = {15,50,100};
  const Int_t ptHardHi[nPtHardBins] = {50,100,-1};




  sprintf(srun,"");
  sprintf(sevent,"");
  sprintf(saliencounter,"");
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
    if (!(strcmp(gApplication->Argv(i),"--aliencounter")))
      naliencounter = atoi(gApplication->Argv(i+1));
    sprintf(saliencounter,"%d",naliencounter);

    if (!(strcmp(gApplication->Argv(i),"--process")))
      sprintf(sprocess, gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--pthardbin")))
      npthardbin = atoi(gApplication->Argv(i+1));

    
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

  if(!(strcmp(sminpthard,""))){
    sprintf(sminpthard,"%d",ptHardLo[npthardbin%nPtHardBins]);
    Printf(">>> MinPtHard %s %d",sminpthard,npthardbin%nPtHardBins);
  }
  if(!(strcmp(smaxpthard,""))){
    sprintf(smaxpthard,"%d",ptHardHi[npthardbin%nPtHardBins]);
    Printf(">>> MaxPtHard %s %d",smaxpthard,npthardbin%nPtHardBins);
  }
 

  
  /*
  char cmd[200] ; 
  sprintf(cmd, ".! tar zxvf DBpp.tgz") ; 
  gROOT->ProcessLine(cmd) ; 
  sprintf(cmd, ".! tar zxvf TPCCalib_v4-16-Rev-01.tgz") ; 
  gROOT->ProcessLine(cmd) ; 
  */
  gSystem->Exec("echo $PWD ");
  gSystem->Exec("ls -l ");

  seed = nrun * 100000 + naliencounter * 1000 + nevent;
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
//  gSystem->Setenv("DC_RUN_TYPE","kPyJetJet");
//  gSystem->Setenv("PTHARDMIN","40");
//  gSystem->Setenv("PTHARDMAX","-1");

  gSystem->Setenv("PTGAMMAPI0MIN",sminptgammapi0);//"1");
  gSystem->Setenv("QUENCHING",squench);
  gSystem->Setenv("QHAT",sqhat);
//  gSystem->Setenv("QUENCHING","0");
//  gSystem->Setenv("QHAT","0");

  gSystem->Setenv("ECMS","10000");
//  gSystem->Setenv("ECMS","14000");
  gSystem->Setenv("ALIMDC_RAWDB1","./mdc1");
  gSystem->Setenv("ALIMDC_RAWDB2","./mdc2");
  gSystem->Setenv("ALIMDC_TAGDB","./mdc1/tag");
  gSystem->Setenv("ALIMDC_RUNDB","./mdc1/meta");
  cout<< "SIMRUN:: Run " << gSystem->Getenv("DC_RUN") << " Event " << gSystem->Getenv("DC_EVENT") 
	  << " Process "   << gSystem->Getenv("DC_RUN_TYPE") 
	  << " minpthard " << gSystem->Getenv("PTHARDMIN") 
	  << " maxpthard " << gSystem->Getenv("PTHARDMAX") 
	  << " minpt "     << gSystem->Getenv("PTGAMMAPI0MIN") 
	  << endl;
  //gSystem->Exec("cp $ROOTSYS/etc/system.rootrc .rootrc");
  cout<<">>>>> SIMULATION <<<<<"<<endl;

  gSystem->Exec("echo ALICE_ROOT $ALICE_ROOT");
  gSystem->Exec("echo ROOTSYS $ROOTSYS");
  gSystem->Exec("echo LD_LIB $LD_LIBRARY_PATH");
  gSystem->Exec("echo PATH $PATH");
  gSystem->Exec("ls -l $ALICE_ROOT/bin/");
  gSystem->Exec("ls -l $ALICE_ROOT/bin/*");

  gSystem->Exec(Form("%s -b -q \"sim.C(%d)\" > sim.log 2>&1",program,nevent));
  

  // residual


  gSystem->MakeDirectory("residual");
  gSystem->Exec("rm residual/*.root");
  gSystem->Exec("ln -s ${PWD}/*.root residual/");

  // residual_trd
  gSystem->MakeDirectory("residual_trd");
  gSystem->Exec("rm residual_trd/*.root");
  gSystem->Exec("ln -s ${PWD}/*.root residual_trd/");

  gSystem->ChangeDirectory("residual");
  cout<<">>>>> RECONSTRUCTION <<<<<"<<endl;
  Int_t iStatus  = 0;
  iStatus = gSystem->Exec(Form("%s -b -q ../rec.C > rec.log 2>&1",program));
  Printf("%d",iStatus);
  cout<<">>>>> TAG <<<<<"<<endl;
  //  iStatus = gSystem->Exec("aliroot -b -q ../tag.C 2>&1 tag.log");
  iStatus = gSystem->Exec(Form("%s -b -q ../tag.C > tag.log 2>&1",program));
  Printf("%d",iStatus);
  cout<<">>>>> CHECK ESD <<<<<"<<endl;
  iStatus = gSystem->Exec(Form("%s -b -q ../CheckESD.C > check.log 2>&1",program));
  cout<<">>>>> CREATE AOD <<<<<"<<endl;
  iStatus = gSystem->Exec(Form("%s -b -q ../runAODFilterJets.C > aod.log 2>&1",program));
  Printf("%d",iStatus);

  gSystem->Exec(Form("mv Run*tag.root ../"));
  gSystem->Exec(Form("mv *RecPoints.root ../")); // for debugging
  TString addString("");
  // logs
  gSystem->Exec(Form("mv rec.log ../rec%s.log",addString.Data()));
  gSystem->Exec(Form("mv tag.log ../tag%s.log",addString.Data()));
  gSystem->Exec(Form("mv check.log ../check%s.log",addString.Data()));
  gSystem->Exec(Form("mv aod.log ../aod%s.log",addString.Data()));
  // roots
  gSystem->Exec(Form("mv AliESDs.root ../AliESDs%s.root",addString.Data()));
  gSystem->Exec(Form("mv AliESDfriends.root ../AliESDfriends%s.root",addString.Data()));
  gSystem->Exec(Form("mv AliAOD.root ../AliAOD%s.root",addString.Data()));
  gSystem->ChangeDirectory("../");


  

  // residual_trd

  gSystem->ChangeDirectory("residual_trd");
  /*
  sprintf(cmd, ".! tar zxvf ../DBpp.tgz") ; 
  gROOT->ProcessLine(cmd) ; 
  sprintf(cmd, ".! tar zxvf ../TPCCalib_v4-16-Rev-01.tgz") ; 
  gROOT->ProcessLine(cmd) ; 
  */
  gSystem->Exec("ls -l");
  cout<<">>>>> REC <<<<<"<<endl;
  iStatus = gSystem->Exec(Form("%s -b -q ../rec2.C > rec.log 2>&1",program));
  
  Printf("%d",iStatus);
  cout<<">>>>> CHECK ESD <<<<<"<<endl;
  iStatus = gSystem->Exec(Form("%s -b -q ../CheckESD.C > check.log 2>&1",program));
  Printf("%d",iStatus);
  cout<<">>>>> CREATE AOD <<<<<"<<endl;
  iStatus = gSystem->Exec(Form("%s -b -q ../runAODFilterJets.C > aod.log 2>&1",program));
  Printf("%d",iStatus);
  TString addString("TRD");
  // logs
  gSystem->Exec(Form("mv rec.log ../rec%s.log",addString.Data()));
  gSystem->Exec(Form("mv check.log ../check%s.log",addString.Data()));
  gSystem->Exec(Form("mv aod.log ../aod%s.log",addString.Data()));
  // roots
  gSystem->Exec(Form("mv AliESDs.root ../AliESDs%s.root",addString.Data()));
  gSystem->Exec(Form("mv AliESDfriends.root ../AliESDfriends%s.root",addString.Data()));
  gSystem->Exec(Form("mv AliAOD.root ../AliAOD%s.root",addString.Data()));
  gSystem->ChangeDirectory("../");
  gSystem->Exec("ls -l");
}
