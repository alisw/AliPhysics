
AliGenerator* AddMCGenppEPOSLHC(double energy = 13000.)
{

  gROOT->ProcessLine(".! source /cvmfs/alice.cern.ch/etc/login.sh");
  gROOT->ProcessLine(".! eval $(alienv enter VO_ALICE@GCC-Toolchain::v4.9.3-alice3-5)");
  gROOT->ProcessLine(".! eval $(alienv enter VO_ALICE@CRMC::v1.5.4-8)");

  printf("----- CRMC PARAM -----\n");
  gROOT->ProcessLine(".! cp ${ALICE_PHYSICS}/PWGLF/CommonUtils/MConthefly/macros/crmc_template.param crmc.local.param");
  //  gROOT->ProcessLine(".! sed -ibak 's,BASEDIR,'\"$CRMC_BASEDIR\"',' crmc.local.param");
  gROOT->ProcessLine(".! cat crmc.local.param");
  printf("----------------------\n");

  printf("--- LAUNCHING CRMC ---\n");
  TString cmd = Form("crmc -t -c crmc.local.param -f crmceventfifo -o hepmc -p%d -P-%d -n201 -m0", (Int_t)energy / 2, (Int_t)energy / 2);
  printf("%s\n", cmd.Data());

  gROOT->ProcessLine(Form(".! %s &", cmd.Data()));

  AliGenReaderHepMC *reader = new AliGenReaderHepMC();
  reader->SetFileName("crmceventfifo");
  AliGenExtFile *gener = new AliGenExtFile(-1);
  gener->SetReader(reader);

  return gener;
}

