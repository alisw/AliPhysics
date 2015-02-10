/// \file compRec.C

void compRec()
{

  /// Run gen info maker

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWGPP");
  AliRecInfoMaker *t2 = new AliRecInfoMaker("genTracks.root","cmpESDTracks.root","galice.root",0,0);
  t2->Exec();


}
