
AliAnalysisHadEtReconstructed * ConfigHadEtReconstructed(){
  //cout<<"Hello I am configuring you"<<endl;
  AliAnalysisHadEtReconstructed *hadEtReco = new AliAnalysisHadEtReconstructed();
  hadEtReco->SetDataSet(20111);
  //Set corrections...

  TFile *infile = new TFile("corrections.root");
  corrections = (AliAnalysisHadEtCorrections *)infile->Get("hadCorrectionEMCAL");
  cout<<"Setting the AliAnalysisHadEtCorrections to "<<corrections->GetName()<<endl;
  cout<<"eta cut is "<<corrections->GetEtaCut()<<endl;
  hadEtReco->SetCorrections(corrections);

  hadEtReco->Init();
  return hadEtReco;
}

