void runCorrectL1Phase()
{
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I$HOME/RootUtils/MyClass");//for my class

  gROOT->LoadMacro("GetCorrectL1Phase.C+g");
  printf("GetCorrectL1Phase.C will run.\n");

  TString trigger = "kINT7";
  cout << "selected trigger is " << trigger << endl;

  Int_t runarray[] = {244340, 244343, 244351, 244355, 244359, 244364, 244377, 244411, 244416, 244418,
                 244421, 244453, 244480, 244481, 244482, 244483, 244484, 244531, 244540, 244542,
                 244617, 244618, 244627, 244628};

  const Int_t Nrun = sizeof(runarray)/sizeof(runarray[0]);
  cout << "Nrun = " << Nrun << endl;

  Int_t run=0;
  const Int_t Nddl=21; //20 SRUs + 1 STU
  Int_t CorrectL1phase[21]={};
  for(Int_t i=0;i<21;i++) CorrectL1phase[i] = -999;

  TFile *outfile = new TFile("CorrectL1Phase_LHC15n.root","RECREATE");
  TTree *tree = new TTree("CorrectL1Phasetree","Tree of L1 Phase based on BC%4");
  tree->Branch("run",&run,"run/I");
  tree->Branch("CorrectL1phase",CorrectL1phase,"CorrectL1phase[21]/I"); 

  for(Int_t i=0;i<Nrun;i++){
    run = runarray[i];
    GetCorrectL1Phase(run[i],trigger,CorrectL1phase);
    for(Int_t iddl=0;iddl<21;iddl++) cout << "run = " << run << " , ddl = " << iddl << " , L1phase = " << CorrectL1phase[iddl] << endl;
    tree->Fill();
  }

  outfile->WriteTObject(tree);
  outfile->Close();

}

