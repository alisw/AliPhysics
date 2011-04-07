
AliAnalysisEtMonteCarlo * ConfigEtMonteCarlo(Bool_t isEmcal = true){
  Bool_t EMCAL = isEmcal;
  if(EMCAL){
    AliAnalysisEtMonteCarloEmcal *totEtMC = new AliAnalysisEtMonteCarloEmcal();
    //Look at the 2010 p+p data...
    totEtMC->SetDataSet(2010);
    
    totEtMC->Init();
    return totEtMC;
  }
  else{
    AliAnalysisEtMonteCarloPhos *totEtMC = new AliAnalysisEtMonteCarloPhos();
    //Look at the 2010 p+p data...
    totEtMC->SetDataSet(2010);
    
    totEtMC->Init();
    return totEtMC;
  }
}
