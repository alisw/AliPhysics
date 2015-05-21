void runRun()
{
    printf("RUN.C");
    
    printf("current run:%s\n",gSystem->Getenv("CURRENT_RUN_NUMBER"));
    
    gROOT->ProcessLine(TString::Format(".x visscan_init.C(\"local:///daqSFS/ONLRECO/cdb\", \"/daqSFS/ONLRECO/reco/run%s\", 0, 1, 1)",gSystem->Getenv("CURRENT_RUN_NUMBER")));
}