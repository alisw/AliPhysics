//These functions are needed for the plugin to work.  Values have to be set by hand.
Bool_t GetIsEMCAL(){
  return kTRUE;
}

Bool_t GetIsMC(){
  return kTRUE;
}


AliAnalysisEtMonteCarlo * ConfigEtMonteCarlo(Bool_t EMCAL = true, Bool_t DETAIL = false){
  //Bool_t EMCAL = true;
  if(EMCAL){
    if(DETAIL){
      AliAnalysisEmEtMonteCarlo *totEtMC = new AliAnalysisEmEtMonteCarlo();
      //Look at the 2010 Pb+Pb data...
      totEtMC->SetDataSet(20100);
      totEtMC->Init();
      return totEtMC;
    }
    else{
      AliAnalysisEtMonteCarloEmcal *totEtMC = new AliAnalysisEtMonteCarloEmcal();
      //Look at the 2010 Pb+Pb data...
      totEtMC->SetDataSet(20100);
      totEtMC->CalcTrackMatchVsMult();
      totEtMC->CalcForKaonCorrection();
      totEtMC->Init();
      return totEtMC;
    }
  }
  else{
    AliAnalysisEtMonteCarloPhos *totEtMC = new AliAnalysisEtMonteCarloPhos();
    //Look at the 2010 Pb+Pb data...
    totEtMC->SetDataSet(20100);
      totEtMC->CalcTrackMatchVsMult();
      totEtMC->CalcForKaonCorrection();
    totEtMC->Init();
    return totEtMC;
  }
}
