AliAnalysisTaskCaloFilter * ConfigCaloFilter(){
  
  AliAnalysisTaskCaloFilter * filter = new AliAnalysisTaskCaloFilter("EMCALFilter");
  filter->SelectCollisionCandidates(); 
  filter->SetCaloFilter(AliAnalysisTaskCaloFilter::kBoth); //kPHOS or kBoth
  filter->SwitchOnClusterCorrection();
  //filter->SetDebugLevel(10);
  AliEMCALRecoUtils * reco = filter->GetEMCALRecoUtils();
  reco->SetParticleType(AliEMCALRecoUtils::kPhoton);
  reco->SetW0(4.5);
  
  reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);
  
  TGeoHMatrix *matrix[4];
  
  double rotationMatrix[4][9] = {-0.014587, -0.999892, -0.002031, 0.999892, -0.014591,  0.001979, -0.002009, -0.002002,  0.999996,
    -0.014587,  0.999892,  0.002031, 0.999892,  0.014591, -0.001979, -0.002009,  0.002002, -0.999996,
    -0.345864, -0.938278, -0.003412, 0.938276, -0.345874,  0.003010, -0.004004, -0.002161,  0.999990,
    -0.345861,  0.938280,  0.003412, 0.938276,  0.345874, -0.003010, -0.004004,  0.002161, -0.999990};
  
  double translationMatrix[4][3] = {0.351659,    447.576446,  176.269742,
    1.062577,    446.893974, -173.728870,
    -154.213287, 419.306156,  176.753692,
    -153.018950, 418.623681, -173.243605};
  for(int j=0; j<4; j++)
  {
    matrix[j] = new TGeoHMatrix();
    matrix[j]->SetRotation(rotationMatrix[j]);
    matrix[j]->SetTranslation(translationMatrix[j]);
    matrix[j]->Print();
    filter->SetEMCALGeometryMatrixInSM(matrix[j],j);
  }
  
  filter->SwitchOnLoadOwnEMCALGeometryMatrices();
  
  reco->SetNonLinearityFunction(AliEMCALRecoUtils::kNoCorrection);
  
  //Time dependent corrections    
  //Recover file from alien  /alice/cern.ch/user/g/gconesab/TimeDepCorrectionDB
//  reco->SwitchOnTimeDepCorrection();
//  char cmd[200] ;
//  sprintf(cmd, ".!tar xvfz CorrectionFiles.tgz") ;
//  gROOT->ProcessLine(cmd) ;
  
  //Recalibration factors
  //Recover the file from alien  /alice/cern.ch/user/g/gconesab/RecalDB
//  reco->SwitchOnRecalibration();
//  TFile * f = new TFile("RecalibrationFactors.root","read");
//  TH2F * h0 = (TH2F*)f->Get("EMCALRecalFactors_SM0");
//  TH2F * h1 = (TH2F*)f->Get("EMCALRecalFactors_SM1");
//  TH2F * h2 = (TH2F*)f->Get("EMCALRecalFactors_SM2");
//  TH2F * h3 = (TH2F*)f->Get("EMCALRecalFactors_SM3");
//  
//  reco->SetEMCALChannelRecalibrationFactors(0,h0);
//  reco->SetEMCALChannelRecalibrationFactors(1,h1);
//  reco->SetEMCALChannelRecalibrationFactors(2,h2);
//  reco->SetEMCALChannelRecalibrationFactors(3,h3);
//  
//  //Bad channels
//  //Recover the file from alien  /alice/cern.ch/user/g/gconesab/BadChannelsDB
//  reco->SwitchOnBadChannelsRemoval();
//  reco->SwitchOnDistToBadChannelRecalculation();
//  TFile * fbad = new TFile("BadChannels.root","read");
//  TH2I * hbad0 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod0");
//  TH2I * hbad1 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod1");
//  TH2I * hbad2 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod2");
//  TH2I * hbad3 = (TH2I*)fbad->Get("EMCALBadChannelMap_Mod3");
//  reco->SetEMCALChannelStatusMap(0,hbad0);
//  reco->SetEMCALChannelStatusMap(1,hbad1);
//  reco->SetEMCALChannelStatusMap(2,hbad2);
//  reco->SetEMCALChannelStatusMap(3,hbad3);
  
  //reco->Print("");
  filter->PrintInfo();   

  return filter;

}