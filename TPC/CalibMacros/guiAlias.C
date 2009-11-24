/*
  //Make default Aliases for guiTime:
  // 1. Run the guiTime
  //  guiTime
  // 2. Define aliases
  .L $ALICE_ROOT/TPC/CalibMacros/guiAlias.C
  Init();
  // 3. Use aliases inside the guiTime
  //    You can use them as custom cuts
  //
  // browse special streamers
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chainDCS = tool.MakeChainRandom("time.txt","dcs",0,10000);
  TChain * chainCTP = tool.MakeChainRandom("time.txt","ctp",0,10000);
  TChain * chainAlign = tool.MakeChainRandom("time.txt","align",0,10000);
  
*/

TTree * guiTree =  guiTime->GetChain();

void guiAlias(){
  guiTree =  guiTime->GetChain();
  MakeAliasCE(4);
  MakeAliasLT(4);
  MakeAliasCosmic(4);   
  guiTree->SetAlias("goCut","abs(goofie.fElements[3]-2.677)<0.1");  
}

void SetStyle(){ 
  Float_t mx0=0.2, mx1=0.1, my0=0.15, my1=0.1;
  guiTime->GetCanvas()->SetTicks(1,1);
  guiTime->GetCanvas()->SetMargin(mx0,mx1,my0,my1);
  gStyle->SetTitleYSize(0.03);
  gStyle->SetTitleXSize(0.03);
  gStyle->SetTitleXOffset(2);
  gStyle->SetTitleYOffset(6);
}

void MakeAliasCE(Double_t deltaT){
  //
  // Aliases cuts for CE
  //
  guiTree->SetAlias("ceCut0", "tdriftCE.fElements[72]>100 && tdriftCE.fElements[73]>100");
  guiTree->SetAlias("dceCutTime", Form("sqrt(dcea^2+dcec^2)<%f",deltaT*3600));
  guiTree->SetAlias("ceCut","dceCutTime&&ceCut0");
};

void MakeAliasLT(Double_t deltaT){
  guiTree->SetAlias("ltCut", Form("sqrt(dla^2+dlc^2)<%f", deltaT*3600)); 
}

void MakeAliasCosmic(Double_t deltaT){
  guiTree->SetAlias("cosmicCut", Form("abs(dcosmic)<%f", deltaT*3600));
  guiTree->SetAlias("itsCut", Form("((dits!=0)&&abs(dits)<%f)", deltaT*3600));
}


void SetGoofieAlias(){
  //
  // goofie aliases
  //
  guiTree->SetAlias("ptrelG","(goofie.fElements[17]/0.3426-1)");
  guiTree->SetAlias("vdriftGN","goofie.fElements[3]/(1+ptrelG)");
  //
}


