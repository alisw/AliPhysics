void batch()
{
  ifstream ifs("goodruns.txt");
  int run = 0;
  const char* dir = "../rootfiles/res_LHC10e_09122010/mergedruns";
  char* inFileName = "";
  char* outFileName = "";

  while (ifs >> run){
    inFileName = Form("%s/merged_run%i.root", dir, run);

    outFileName = Form("output/histos_%i.root", run);

    char* cmd = Form("root -b -q 'runAutoCorr.C+g(\"%s\", \"%s\")'", 
		     inFileName, outFileName); 

    cout << cmd << endl;

    gSystem->Exec(cmd);

  }
  return;
}
