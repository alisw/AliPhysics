
AliAnalysisHadEtReconstructed * ConfigHadEtReconstructed(){
  //cout<<"Hello I am configuring you"<<endl;
  AliAnalysisHadEtReconstructed *hadEtReco = new AliAnalysisHadEtReconstructed();
  hadEtReco->Init();
  return hadEtReco;
}

