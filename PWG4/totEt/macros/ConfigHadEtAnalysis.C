
AliAnalysisHadEtCorrections * ConfigHadEtAnalysis(){
  //cout<<"Hello I am configuring you"<<endl;

  TFile *infile = new TFile("corrections.root");
  corrections = (AliAnalysisHadEtCorrections *)infile->Get("hadCorrectionEMCAL");
  cout<<"Setting the AliAnalysisHadEtCorrections to "<<corrections->GetName()<<endl;
  cout<<"eta cut is "<<corrections->GetEtaCut()<<endl;
  return corrections;
}
