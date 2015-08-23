void runPedestal()
{
  //Daiki Sekihata (Hiroshima University)
  //run with "aliroot -l -b -q runPedestal.C"

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I$HOME/SRU/APD/");

  //gROOT->LoadMacro("Pedestals.C+g");
  //printf("Pedestals.C will run.\n");

  const Int_t runnumber = 233043;

  printf("Analyzing Run : %d\n",runnumber);

  //TString textfile = Form("%d.txt",runnumber);
  //Pedestals(textfile);

  gROOT->LoadMacro("AliPHOSFEEMapRun2.cxx+g");
  gROOT->LoadMacro("CreatePedestalTable.C+g");
  printf("CreatePedestalTable.C will run.\n");

  CreatePedestalTable(runnumber);

  gSystem->Exec(Form("tar cf PedestalTable_%d.tar Pedestal*.txt",runnumber));
  gSystem->Exec("rm Pedestal*.txt");

  printf("Done!\n");

}

