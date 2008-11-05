void compRec()
{

  //
  // Run gen info maker  
  //
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPWG0base.so");
  gSystem->Load("libPWG0dep.so");
  gSystem->Load("libPWG1.so");
  AliRecInfoMaker *t2 = new AliRecInfoMaker("genTracks.root","cmpESDTracks.root","galice.root",0,0);
  t2->Exec();


}
