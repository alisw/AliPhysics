/*

Simple macro to test EMCAL Reconstruction

J.L. Klay
LLNL

*/


void TestEMCALReconstruction(Int_t nev =3) {

  AliReconstruction rec;

  TStopwatch timer;
  timer.Start();
  rec.SetRunTracking("");
  rec.SetRunVertexFinder(kFALSE);
  //calls local reconstruction of EMCAL and filling of ESD
  rec.SetRunLocalReconstruction("EMCAL");  //only do emcal
  rec.SetFillESD("EMCAL");
  //Start at event 1, not event 0, and do all the rest...
  rec.Run(0,-1);
  timer.Stop();
  timer.Print();

}
