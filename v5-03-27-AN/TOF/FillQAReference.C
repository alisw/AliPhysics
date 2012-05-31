/*
 * this macro will be used to setup QA reference repository.
 * reference histos are retrieved from input file and stored 
 * in the QA reference repository.
 *
 * autor: Roberto Preghenella (R+)
 * contact: preghenella@bo.infn.it
 */

FillQAReference(const Char_t *fileName = "$ALICE_ROOT/TOF/data/TOFQA.LHC09a4.80050.root")
{

  TFile *file = TFile::Open(fileName);
  TH2F *hTOFDig_ClusMap = (TH2F *)file->Get("TOFDig_ClusMap");
  TH1F *hTOFDig_ClusTime = (TH1F *)file->Get("TOFDig_ClusTime");
  TH1F *hTOFDig_ClusToT = (TH1F *)file->Get("TOFDig_ClusToT");
  TH1F *hTOFDig_NClus = (TH1F *)file->Get("TOFDig_NClus");

  if (!hTOFDig_ClusMap || !hTOFDig_ClusTime || !hTOFDig_ClusToT || !hTOFDig_NClus) {
    printf("couldn't retrieve all reference histos from file. abort.\n");
    return;
  }

  printf("reference histo successfully retrieved.\n");

  printf("filling QA reference repository...\n");
  printf("nothing done\n");
  
}
