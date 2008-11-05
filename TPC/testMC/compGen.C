void compGen()
{

  //
  // Run gen info maker  
  //
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG1.so");
  AliGenInfoMaker *t = new AliGenInfoMaker("galice.root","genTracks.root",0,0);
  t->Exec();
}
