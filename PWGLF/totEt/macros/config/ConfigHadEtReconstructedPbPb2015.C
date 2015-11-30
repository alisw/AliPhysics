
//AliAnalysisHadEtReconstructedLocal * ConfigHadEtReconstructedLocal(){
AliAnalysisHadEtReconstructed * ConfigHadEtReconstructed(){
  //cout<<"Hello I am configuring you"<<endl;
  //AliAnalysisHadEtReconstructedLocal *hadEtReco = new AliAnalysisHadEtReconstructedLocal();
  AliAnalysisHadEtReconstructed *hadEtReco = new AliAnalysisHadEtReconstructed();
  hadEtReco->SetDataSet(2015);
  //Set corrections...

  TGrid::Connect("alien://") ;
  TString infilename = "alien:///alice/cern.ch/user/c/cnattras/rootFiles/corrections/corrections.LHC13e1abc.PbPb.168464.ForData.2015.root";
  TString outfilename = "corrections.root";
  TFile::Cp(infilename.Data(),outfilename.Data());


  TFile *infile = new TFile("corrections.root");
  //corrections = (AliAnalysisHadEtCorrectionsLocal *)infile->Get("hadCorrectionEMCAL");
  //cout<<"Setting the AliAnalysisHadEtCorrectionsLocal to "<<corrections->GetName()<<endl;
  corrections = (AliAnalysisHadEtCorrections *)infile->Get("hadCorrectionEMCAL");
  cout<<"Setting the AliAnalysisHadEtCorrections to "<<corrections->GetName()<<endl;
  cout<<"eta cut is "<<corrections->GetEtaCut()<<endl;
  cout<<"Data set is 2015"<<endl;
  //hadEtReco->SetCorrectionsLocal(corrections);
  hadEtReco->SetCorrections(corrections);

  //hadEtReco->SetNumberOfCentralityBins(41);
  hadEtReco->Init();
  return hadEtReco;
}

