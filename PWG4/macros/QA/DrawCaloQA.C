DrawCaloQA(TString calo = "EMCAL", Bool_t kine = kFALSE)
{
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPWG4PartCorrBase.so");
  gSystem->Load("libPWG4PartCorrDep.so");
  
  TFile *f = new TFile("AnalysisResults.root","read");
  TDirectoryFile *dir = (TDirectoryFile *)f->Get("CaloQA");
  printf("Get list\n");
  TList* outputList = (TList*)dir->Get("CaloQA");
  
  AliAnaCalorimeterQA* qa = new AliAnaCalorimeterQA();
  qa->AddToHistogramsName(calo+"_");
  qa->SwitchOnPlotsMaking();
  qa->SetStyleMacro("style.C");
  if(kine) qa->SwitchOnDataMC() ;
  else     qa->SwitchOffDataMC() ;
  if(calo == "EMCAL"){
    qa->SetNumberOfModules(4);
    qa->SwitchOnCalorimetersCorrelation();

  }
  else if (calo=="PHOS"){
    qa->SetNumberOfModules(3);
    qa->SwitchOffCalorimetersCorrelation();

  }

  //Set Histrograms bins and ranges
  qa->SetHistoPtRangeAndNBins(0, 5, 100) ;
  if(calo=="EMCAL"){
    qa->SetHistoPhiRangeAndNBins(1.37, 2.23, 25) ;
    qa->SetHistoEtaRangeAndNBins(-0.8, 0.8, 20) ;
  }
  else if (calo=="PHOS"){
    qa->SetHistoPhiRangeAndNBins(255*TMath::DegToRad(), 325*TMath::DegToRad(), 200) ;
    qa->SetHistoEtaRangeAndNBins(-0.13, 0.13, 160) ;
  }
  qa->SetHistoMassRangeAndNBins(0., 0.6, 300) ;
  qa->SetHistoAsymmetryRangeAndNBins(0., 1. , 25) ;
  qa->SetHistoPOverERangeAndNBins(0,10.,100);
  qa->SetHistodEdxRangeAndNBins(0.,400.,200);
  qa->SetHistodRRangeAndNBins(0.,TMath::Pi(),300);
  qa->SetHistoTimeRangeAndNBins(0.,1000,1000);
  qa->SetHistoRatioRangeAndNBins(0.,2.,100);
  qa->SetHistoVertexDistRangeAndNBins(0.,500.,100);
  qa->SetHistoNClusterCellRangeAndNBins(0,300,300);
  qa->SetHistoXRangeAndNBins(-250,100,25);
  qa->SetHistoYRangeAndNBins(370,450,100);
  qa->SetHistoZRangeAndNBins(-350,350,50);
  qa->SetHistoRRangeAndNBins(420,460,300);

  //Make the histograms
  qa->Terminate(outputList);
  
  
}
