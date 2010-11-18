
AliAnalysisHadEtMonteCarlo * ConfigHadEtMonteCarlo(){
  //cout<<"Hello I am configuring you"<<endl;
  AliAnalysisHadEtMonteCarlo *hadEtMC = new AliAnalysisHadEtMonteCarlo();
  //Whether or not to investigate the effects of efficiency, momentum resolution, PID, etc.
  hadEtMC->InvestigateSmearing(kFALSE);

  //Whether or not to look at Et(sim)-Et(reco) for full acceptance
  hadEtMC->InvestigateFull(kFALSE);

  //Whether or not to look at Et(sim)-Et(reco) for EMCAL acceptance
  hadEtMC->InvestigateEMCAL(kTRUE);

  //Whether or not to look at Et(sim)-Et(reco) for PHOS acceptance
  hadEtMC->InvestigatePHOS(kTRUE);

  //Whether or not to look at Et(sim)-Et(reco) for Pi/K/p in full acceptance (full acceptance must be turned on)
  hadEtMC->InvestigatePiKP(kFALSE);

  //Look at ITS+TPC tracks
  hadEtMC->RequireITSHits(kTRUE);

  hadEtMC->Init();
  return hadEtMC;
}
