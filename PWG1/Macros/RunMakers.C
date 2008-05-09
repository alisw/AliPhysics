// Run this macro to correlate MC and reconstruction information (PWG1 library).
// Macro must be run in directory containing MC and ESD trees
//
RunMakers()
{
// load AliRoot libraries
gSystem->Load("libANALYSIS.so");
gSystem->Load("libANALYSISalice.so");
gSystem->Load("libPWG1.so");

// collect MC information
AliGenInfoMaker *infoMC = new AliGenInfoMaker("galice.root","genTracks.root",0,0);
infoMC->Exec();

// correlate MC and reconstruction information
AliRecInfoMaker *infoMCR =
new AliRecInfoMaker("genTracks.root","cmpESDTracks.root","galice.root",0,0);
infoMCR->Exec();
}
