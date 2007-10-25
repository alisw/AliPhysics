void Rec(const Int_t debLevel=0)
{

  AliLog::SetGlobalDebugLevel(debLevel);
  AliReconstruction rec ;
  rec.SetRunReconstruction("PHOS") ;
  rec.SetRunTracking("PHOS") ;
  rec.SetRunVertexFinder(kFALSE) ; 
  rec.SetFillESD("PHOS") ; 

  //Uncomment the following lines to use rec. parameters 
  //other than default ones.
  // AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  //  recEmc->SetSubtractPedestals(kTRUE);
  //  recEmc->SetMinE(0.01);
  //  recEmc->SetClusteringThreshold(0.01);
  //  AliPHOSReconstructor::SetRecoParamEmc(recEmc);

  rec.Run();

}
