
AliAnalysisEtReconstructed * ConfigEtReconstructed(Bool_t EMCAL = true, Bool_t DETAIL = false){
  gInterpreter->GenerateDictionary("std::map<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;map")  ;
  gInterpreter->GenerateDictionary("std::pair<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;utility");
  //Bool_t EMCAL = true;
  TFile *infile = new TFile("corrections.root");
  corrections = (AliAnalysisHadEtCorrections *)infile->Get("hadCorrectionEMCAL");
  cout<<"Setting the AliAnalysisHadEtCorrections to "<<corrections->GetName()<<endl;
  cout<<"eta cut is "<<corrections->GetEtaCut()<<endl;

  if(EMCAL){
    if(DETAIL){
      AliAnalysisEmEtReconstructed *totEtReco = new AliAnalysisEmEtReconstructed();
      totEtReco->SetDataSet(2010);
      //Set corrections...
      totEtReco->SetCorrections(corrections);    
      totEtReco->Init();
      return totEtReco;
    }
    else{
      AliAnalysisEtReconstructedEmcal *totEtReco = new AliAnalysisEtReconstructedEmcal();
      totEtReco->SetDataSet(2010);
      //Set corrections...
      totEtReco->SetCorrections(corrections);    
      totEtReco->Init();
      return totEtReco;
    }
  }
  else{
    AliAnalysisEtReconstructedPhos *totEtReco = new AliAnalysisEtReconstructedPhos();
    totEtReco->SetDataSet(2010);
    //Set corrections...
    totEtReco->SetCorrections(corrections);    
    totEtReco->Init();
    return totEtReco;
  }
}

