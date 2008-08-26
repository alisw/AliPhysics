/*

Simple macro to test EMCAL Reconstruction

J.L. Klay
LLNL

*/


void TestEMCALReconstruction(Int_t nev =-1) {

  //  AliLog::SetModuleDebugLevel("EMCAL",100);

  AliReconstruction rec;

  rec.SetRunTracking("");
  rec.SetRunVertexFinder(kFALSE);

  //calls local reconstruction of EMCAL and filling of ESD
  rec.SetRunLocalReconstruction("EMCAL");  //only do emcal
  rec.SetFillESD("EMCAL");
  rec.SetEventRange(0,nev);
  //rec.SetInput("raw.root");
  //rec.SetRunQA(":");

  // **** The field map settings must be the same as in Config.C !
  AliMagFMaps *field=new AliMagFMaps("Maps","Maps",2,1.,10.,AliMagFMaps::k5kG);
  Bool_t uniform=kFALSE;
  AliTracker::SetFieldMap(field,uniform);

  TStopwatch timer;
  timer.Start();

  rec.Run();
  timer.Stop();
  timer.Print();
  gObjectTable->Print();

}
