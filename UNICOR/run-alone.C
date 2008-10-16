{
gROOT->LoadMacro("makechain.C"); 
tr = makechain("esdTree","/u/sma/data/mc/v4-13-Rev-01/fulltrd_pdc_0.txt",10);

gSystem->Load("libEG.so");
gSystem->Load("libTree.so");
gSystem->Load("libVMC.so"); 
gSystem->Load("libSTEERBase.so");
gSystem->Load("libESD.so");
gSystem->Load("libANALYSIS");
gSystem->Load("libUNICOR.so");

d0=new AliDEventAliceESD();
lo = new AliDLoop(tr,d0,"unicor-result.root");
lo->Run();
}
