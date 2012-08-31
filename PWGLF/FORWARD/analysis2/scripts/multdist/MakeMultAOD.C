/*
void Setup900GeV(MakeAODTrain p);
void Setup900GeVMC(MakeAODTrain p);
void Add900GeVruns(MakeAODTrain p);
*/
class MakeAODTrain;

void MakeMultAOD()
{
  bool        usePar  =  true;
  bool        mc      =  true;
  Int_t       nEvents = -1;
  UShort_t    proof   = 0;

  const char* name    = "test_flatMult_withbg";
  UShort_t    type    = 1;  // pp==1, PbPb==2
  UShort_t    cms     = 900; // 2750 for 2760GeV
  Short_t     field   = 5; // -5 for 2760GeV!!!
  Bool_t      useCent = false;
  const char* builder = "$(ALICE_ROOT)/PWG2/FORWARD/analysis2/trains/BuildTrain.C";
  gROOT->LoadMacro(builder);
  BuildTrain("MakeAODTrain");

  std::cout << "making AOD train " << std::endl;
  MakeAODTrain t(name, type, cms, field, useCent);
  t.SetDataSet("");
  t.SetNReplica(2);
  t.SetAllowOverwrite(true);
  t.SetProofServer(Form("workers=%d",proof));
 
  
  t.SetROOTVersion("v5-28-00f");
  t.SetAliROOTVersion("v4-21-33-AN");

  //if pp coll
  if(type==1){
    if (cms == 900) { 
      if (!mc) Setup900GeV(t);
      else     Setup900GeVMC(t);
    }
    if (cms == 2750) {
      if (!mc) Setup2760GeV(t);
      else Setup2760GeVMC(t);
    }
    
    if (cms == 7000) { 
      if (!mc) Setup7000GeV(t);
      else     Setup7000GeVMC(t);
    }
  }
  else
    SetupPbPb(t);
   
    
  t.Run("GRID", "FULL", nEvents, mc, usePar);
}
//________________________________________________________

void SetupPbPb(MakeAODTrain& p){
  p.SetDataDir("/alice/data/2010/LHC10h") ;
  p.SetESDPass(2); 


  p.AddRun(138190);
  p.AddRun(138534);
  p.AddRun(138364);
  p.AddRun(138442);
  //p.AddRun(139465);
  p.AddRun(138396);
  //p.AddRun(137722);
  //p.AddRun(139107);
  //p.AddRun(139437);
  p.AddRun(138653);
  //p.AddRun(139038);
  //p.AddRun(137608);
}

void Setup900GeV(MakeAODTrain& p){
  p.SetDataDir("/alice/data/2010/LHC10c");
  p.SetESDPass(3); 

  Add900GeVruns(p);

}

void Setup900GeVMC(MakeAODTrain& p){
  //Flat mult
  p.SetDataDir("/alice/sim/LHC10f1");

  //nomal MC (Pythia)
  //p.SetDataDir("/alice/sim/LHC11b1a");

  //900 Gev normal Phojet MC
  //p.SetDataDir("/alice/sim/LHC11c1");
  
  Add900GeVruns(p);

}

void Add900GeVruns(MakeAODTrain& p){
  p.AddRun(118506);
  p.AddRun(118507);
  p.AddRun(118512);
  p.AddRun(118556);
  p.AddRun(118558);
  p.AddRun(118560);
  p.AddRun(118561);
  p.AddRun(121039);
  p.AddRun(121040);
  
}
void Setup2760GeV(MakeAODTrain& p){
  p.SetDataDir("/alice/data/2011/LHC11a");
  p.SetESDPass(2);
  p.SetPassPostfix("_with_SDD");
  Add2760GeVruns(p);
}
void Setup2760GeVMC(MakeAODTrain& p){
  p.SetDataDir("/alice/sim/LHC11b10c");
  Add2760GeVruns(p);
}

void Add2760GeVruns(MakeAODTrain& p){
  /*
  p.AddRun(146860); 
  p.AddRun(146859);	
  p.AddRun(146858);	
  p.AddRun(146857);	
  p.AddRun(146856);
  */
  p.AddRun(146824); 
  p.AddRun(146817);
  p.AddRun(146807);
  p.AddRun(146806);
  p.AddRun(146805);
  p.AddRun(146804);
  p.AddRun(146803);
  p.AddRun(146802);
  p.AddRun(146801);
  p.AddRun(146748);
  p.AddRun(146747);
  p.AddRun(146746);
  p.AddRun(146689);
  p.AddRun(146688);
  p.AddRun(146686);
  
}

void Setup7000GeV(MakeAODTrain& p){
  p.SetDataDir("/alice/data/2010/LHC10d");
  p.AddRun(126424);
  p.AddRun(126422);  
  p.AddRun(126409);
  p.AddRun(126408);
  p.AddRun(126407);
  p.AddRun(126406);
  p.AddRun(126404);
  p.AddRun(126359);
  p.AddRun(126352);
  p.AddRun(126351);
  p.AddRun(126097);
  p.AddRun(125855);
  p.AddRun(125851);
  p.AddRun(125850);
  p.AddRun(125849);
  p.AddRun(125848);
  p.AddRun(125847);
  p.AddRun(125101);
  p.AddRun(125097);
  p.AddRun(125085);
  
  p.SetESDPass(2); 
}
void Setup7000GeVMC(MakeAODTrain& p){
  p.SetDataDir("/alice/sim/LHC10h16");
  p.AddRun(125186);
  p.AddRun(125633);
  p.AddRun(125842);
  p.AddRun(125851);
  p.AddRun(126008);
  p.AddRun(126081);
  p.AddRun(126082);
  p.AddRun(126090);
  p.AddRun(126097);
  p.AddRun(126167);
  p.AddRun(126283);
  p.AddRun(126359);
  p.AddRun(126408);
  p.AddRun(126409);
  p.AddRun(126422);
  p.AddRun(126425);
  p.AddRun(126437);
}


//EOF
