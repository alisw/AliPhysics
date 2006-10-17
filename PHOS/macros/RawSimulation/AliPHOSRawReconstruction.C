void AliPHOSRawReconstruction(TString file="raw.root")
{
  // Reconstruction of RAW data from the input file raw.root
  // Boris Polichtchouk, 13 October 2006

  AliReconstruction rec ;
//   rec.SetOption("PHOS","OldRCUFormat");
  rec.SetRunTracking("PHOS") ;
  rec.SetRunVertexFinder(kFALSE) ; 
  rec.SetRunLocalReconstruction("PHOS") ;
  rec.SetFillESD("") ;

  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  AliTracker::SetFieldMap(field,kFALSE); 

  rec.SetInput(file.Data());  // read RAW data
  rec.Run();

}
