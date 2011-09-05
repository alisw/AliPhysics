
AliAnalysisEtMonteCarlo * ConfigEtMonteCarlo(Bool_t EMCAL = true, Bool_t DETAIL = false){
  //Bool_t EMCAL = true;
  if(EMCAL){
    if(DETAIL){
      AliAnalysisEmEtMonteCarlo *totEtMC = new AliAnalysisEmEtMonteCarlo();
      //Look at the 2010 p+p data...
      totEtMC->SetDataSet(2010);
      totEtMC->Init();
      return totEtMC;
    }
    else{
      AliAnalysisEtMonteCarloEmcal *totEtMC = new AliAnalysisEtMonteCarloEmcal();
      //Look at the 2010 p+p data...
      totEtMC->SetDataSet(2010);
      totEtMC->Init();
      return totEtMC;
    }
  }
  else{
    AliAnalysisEtMonteCarloPhos *totEtMC = new AliAnalysisEtMonteCarloPhos();
    //Look at the 2010 p+p data...
    totEtMC->SetDataSet(2010);
    totEtMC->Init();
    return totEtMC;
  }
}
