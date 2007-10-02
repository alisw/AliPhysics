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
  rec.SetEventRange(0,nev);
  rec.Run();
  timer.Stop();
  timer.Print();

}
