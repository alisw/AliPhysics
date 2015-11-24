
//AliAnalysisHadEtMonteCarloLocal * ConfigHadEtMonteCarloLocal(){
AliAnalysisHadEtMonteCarlo * ConfigHadEtMonteCarlo(){
  //cout<<"Hello I am configuring you"<<endl;
  cout<<"You are analyzing HIJING Pb+Pb simulations"<<endl;
  //AliAnalysisHadEtMonteCarloLocal *hadEtMC = new AliAnalysisHadEtMonteCarloLocal();
  AliAnalysisHadEtMonteCarlo *hadEtMC = new AliAnalysisHadEtMonteCarlo();
  //Whether or not to investigate the effects of efficiency, momentum resolution, PID, etc.
  hadEtMC->InvestigateSmearing(kTRUE);

  //Turns off O(100) histograms that we do not normally use UNLESS we need to calculate new corrections
  hadEtMC->RunLightweight(kFALSE);

  //Whether or not to look at Et(sim)-Et(reco) for full acceptance
  hadEtMC->InvestigateFull(kTRUE);

  //Whether or not to look at Et(sim)-Et(reco) for EMCAL acceptance
  hadEtMC->InvestigateEMCAL(kFALSE);

  //Whether or not to look at Et(sim)-Et(reco) for PHOS acceptance
  hadEtMC->InvestigatePHOS(kFALSE);

  //Whether or not to look at Et(sim)-Et(reco) for Pi/K/p in full acceptance (full acceptance must be turned on)
  hadEtMC->InvestigatePiKP(kFALSE);

  //Look at ITS+TPC tracks
  hadEtMC->RequireITSHits(kTRUE);

  //Look at the 2010 p+p data...
  hadEtMC->SetDataSet(2015);
  hadEtMC->SetV0ScaleDataSet(20100);

  //Turn baryon enhancement on and off
  hadEtMC->EnhanceBaryons(kTRUE);

  //For these data we only have MC w/added signals available so cut those out
  hadEtMC->CheckForHIJINGLabel();

  //hadEtMC->SetNumberOfCentralityBins(41);

  hadEtMC->Init();
  return hadEtMC;
}
