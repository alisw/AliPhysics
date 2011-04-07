
AliAnalysisEtReconstructed * ConfigEtReconstructed(Bool_t isEmcal = true){
  Bool_t EMCAL = isEmcal;
  if(EMCAL){
    AliAnalysisEtReconstructedEmcal *totEtReco = new AliAnalysisEtReconstructedEmcal();
    totEtReco->SetDataSet(2010);
    //Set corrections...
    
    TFile *infile = new TFile("corrections.root");
    corrections = (AliAnalysisHadEtCorrections *)infile->Get("hadCorrectionEMCAL");
    cout<<"Setting the AliAnalysisHadEtCorrections to "<<corrections->GetName()<<endl;
    cout<<"eta cut is "<<corrections->GetEtaCut()<<endl;
    totEtReco->SetCorrections(corrections);
    
    totEtReco->Init();
    return totEtReco;
  }
  else{
    AliAnalysisEtReconstructedPhos *totEtReco = new AliAnalysisEtReconstructedPhos();
    totEtReco->SetDataSet(2010);
    //Set corrections...
    
    TFile *infile = new TFile("corrections.root");
    corrections = (AliAnalysisHadEtCorrections *)infile->Get("hadCorrectionPHOS");
    cout<<"Setting the AliAnalysisHadEtCorrections to "<<corrections->GetName()<<endl;
    cout<<"eta cut is "<<corrections->GetEtaCut()<<endl;
    totEtReco->SetCorrections(corrections);
    
    totEtReco->Init();
    return totEtReco;
  }
  gInterpreter->GenerateDictionary("std::map<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;map")  ;
  gInterpreter->GenerateDictionary("std::pair<int, AliPhysicsSelection*>", "AliPhysicsSelection.h;utility");
}

