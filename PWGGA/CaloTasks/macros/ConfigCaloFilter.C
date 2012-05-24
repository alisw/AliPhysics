AliAnalysisTaskCaloFilter * ConfigCaloFilter()
{
  // Configure calorimeter filter analysis
  AliAnalysisTaskCaloFilter * filter = new AliAnalysisTaskCaloFilter("EMCALFilter");
  filter->SelectCollisionCandidates(); 
  filter->SetCaloFilter(AliAnalysisTaskCaloFilter::kBoth); //kPHOS or kBoth
  filter->SwitchOffClusterCorrection();
  //filter->SetDebugLevel(10);
  
  filter->SetVzCut(10.);
  filter->SetEnergyCut(10.);
  filter->SetNcellsCut(2);
  
  filter->SwitchOnFillTracks();
  filter->SwitchOnFillAODFile();
  
  filter->PrintInfo();   

  return filter;

}