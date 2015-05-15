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
  rec.SetRunQA(":");

  // Decomment this line in case of real data,
  // add the proper name of the file
  //rec.SetInput("raw.root");
  
  //OCDB settings
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  // Decommetn this line in case of anchored MC runs or data,
  // with the appropriate year
  //rec.SetDefaultStorage("alien://Folder=/alice/data/2011/OCDB");

  rec.SetSpecificStorage("GRP/GRP/Data",
                         Form("local://%s",gSystem->pwd()));

  TStopwatch timer;
  timer.Start();

  rec.Run();
  timer.Stop();
  timer.Print();
  gObjectTable->Print();

}
