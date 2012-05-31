void AliPHOSRawReconstruction(TString file="raw.root")
{
  // Reconstruction of RAW data from the input file raw.root
  // Boris Polichtchouk, 31 Aug 2007


  //AliLog::SetGlobalDebugLevel(1);

  AliReconstruction rec ;
//   rec.SetOption("PHOS","OldRCUFormat");
  rec.SetRunTracking("") ;
  rec.SetRunVertexFinder(kFALSE) ; 
  rec.SetRunLocalReconstruction("PHOS") ;
  rec.SetFillESD("PHOS") ;

  //Set rec. parameters different from the default ones.
  AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  recEmc->SetSubtractPedestals(kFALSE); // do not sibtract pedestals!

  AliPHOSReconstructor::SetRecoParamEmc(recEmc);

  rec.SetInput(file.Data());  // read RAW data
  rec.Run();


}
