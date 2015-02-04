/// \file compGen.C

void compGen()
{

  /// Run gen info maker

  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");
  gSystem->Load("libPWGPP");
  AliGenInfoMaker *t = new AliGenInfoMaker("galice.root","genTracks.root",0,0);
  t->Exec();
}
