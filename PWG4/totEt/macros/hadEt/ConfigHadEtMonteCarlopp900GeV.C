
AliAnalysisHadEtMonteCarlo * ConfigHadEtMonteCarlo(){
  //cout<<"Hello I am configuring you"<<endl;
  cout<<"You are analyzing 900 GeV p+p simulations using the 2010 track cuts and 2009 V0 scaling"<<endl;
  AliAnalysisHadEtMonteCarlo *hadEtMC = new AliAnalysisHadEtMonteCarlo();
  //Whether or not to investigate the effects of efficiency, momentum resolution, PID, etc.
  hadEtMC->InvestigateSmearing(kFALSE);

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
  hadEtMC->SetDataSet(2009);
  hadEtMC->SetV0ScaleDataSet(2009);

  //Turn baryon enhancement on and off
  hadEtMC->EnhanceBaryons(kTRUE);

  hadEtMC->Init();
  return hadEtMC;
}
