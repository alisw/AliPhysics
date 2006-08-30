void RawReconstruction(TString file="2006run2211.root")
{
  // Script to reconstruct raw data from the PHOS beam test July-August 2006
  // Raw data are available in /castor/cern.ch/alice/testbeam/phos/2006/2006run*.raw
  // and has to be converted to root format
  // Author: Boris Polichtchouk, August 2006

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("PHOS","local://CalibDB");

  AliReconstruction rec ;
  rec.SetOption("PHOS","OldRCUFormat");
  rec.SetRunTracking("PHOS") ;
  rec.SetRunVertexFinder(kFALSE) ; 
  rec.SetRunLocalReconstruction("PHOS") ;
  rec.SetFillESD("PHOS") ; 

  rec.SetInput(file.Data());  // read RAW data
  rec.Run();

}
