void PrintCorrectL1Phase(const TString period="LHC15n")
{
  const TString infile = Form("CorrectL1Phase_%s.root",period.Data());

  TFile *rootfile = TFile::Open(infile,"READ");

  TTree *tree = (TTree*)rootfile->Get("CorrectL1Phasetree");
  Int_t run=0;
  const Int_t Nddl=21;
  Int_t CorrectL1phase[Nddl]={};
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("CorrectL1phase",CorrectL1phase);

  const Int_t Nrun = tree->GetEntriesFast();
  printf("There are %d runs in %s where PHOS was in.\n",Nrun,period.Data());

  //there are 21 DDLs for PHOS. 6-19 are for SRUs, 20 is for STU in Run2.
  //0-5 are reserved for M0 and a top half of M1 which does not exist.

  for(Int_t irun=0;irun<Nrun;irun++){


    tree->GetEntry(irun);

  //  if(run != 244628) continue;


    printf("irun:%d == official run number:%d\n",irun,run);

    for(Int_t iddl=0;iddl<Nddl;iddl++){
      printf("  DDL:%d | correct L1phase:%d\n",iddl,CorrectL1phase[iddl]);
    }
 
  }
 
  rootfile->Clone();
}

