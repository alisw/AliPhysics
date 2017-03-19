AliGenerator* AddMCGenEPOS_3(double energy = 7000.0)
{

  // prepare environment
  gROOT->ProcessLine(".! source /cvmfs/alice.cern.ch/etc/login.sh");
  gROOT->ProcessLine(".! eval $(alienv printenv GCC-Toolchain/v4.9.3-7)");
  gROOT->ProcessLine(".! eval $(alienv printenv CRMC::v1.5.4-2)");

  // force path to requested output
  gROOT->ProcessLine(".! sed -e s#__CRMC_BASEDIR__#"$CRMC_ROOT"");



  
  //TString comment = comment.Append(" pp: EPOS-LHC");
  // //    EPOS LHC
   printf("--- LAUNCHING CRMC ---\n");




   
   TString cmd = Form("$CRMC_BASEDIR/bin/crmc -t -c crmc.local.param -f crmceventfifo -o hepmc -p%d -P-%d -n201 -m0 crmceventfifo", (Int_t)energy/2, (Int_t)energy/2);
   printf("%s\n", cmd.Data());
   printf("----- CRMC PARAM -----\n");
   gROOT->ProcessLine(".! cp crmc.param crmc.local.param");
   gROOT->ProcessLine(".! sed -ibak 's,BASEDIR,'\"$CRMC_BASEDIR\"',' crmc.local.param");
   gROOT->ProcessLine(".! cat crmc.local.param");
   printf("----------------------\n");
   gROOT->ProcessLine(Form(".! %s &", cmd.Data()));

  //gROOT->ProcessLine(".! echo \"pass1\"");
  //gROOT->ProcessLine(".! mkfifo crmceventfifo");
  //gROOT->ProcessLine(".! ls ./");
  //gROOT->ProcessLine(".! echo \"pass2\"");
  //gROOT->ProcessLine(".! sh gen_eposlhc.sh crmceventfifo &");
  //gROOT->ProcessLine(".! echo \"pass3\"");

  AliGenReaderHepMC *reader = new AliGenReaderHepMC();
  reader->SetFileName("crmceventfifo");
  AliGenExtFile *gener = new AliGenExtFile(-1);
  gener->SetReader(reader);

  return gener;
}
