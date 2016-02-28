void runPedestal(const Int_t runnumber = 233043)
{
  //Daiki Sekihata (Hiroshima University)
  //run with "aliroot -l -b -q runPedestal.C"

  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I$HOME/SRU/APD/");

  //gROOT->LoadMacro("Pedestals.C+g");
  //printf("Pedestals.C will run.\n");

  printf("Analyzing Run : %d\n",runnumber);

  //TString textfile = Form("%d.txt",runnumber);
  //Pedestals(textfile);

  gROOT->LoadMacro("AliPHOSFEEMapRun2.cxx+g");
  gROOT->LoadMacro("CreatePedestalTable.C+g");
  printf("CreatePedestalTable.C will run\n");

  if (CreatePedestalTable(runnumber)) {
    gSystem->Exec(Form("tar czf PedestalTable_%d.tgz %d/Pedestal*.txt",
		       runnumber,runnumber));
    gSystem->Exec(Form("rm -fr %d",runnumber));
    printf("Done! Pedestals for loading to ALTRO stored in PedestalTable_%d.tgz\n",
	   runnumber);
  }
  else
    return;
}
