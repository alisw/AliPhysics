AliAnalysisTaskCaloFilter * ConfigCaloFilter()
{
  // Configure calorimeter filter analysis
  AliAnalysisTaskCaloFilter * filter = new AliAnalysisTaskCaloFilter("EMCALFilter");
  filter->SelectCollisionCandidates(); 
  filter->SetCaloFilter(AliAnalysisTaskCaloFilter::kEMCAL); //kPHOS, kEMCAL or kBoth
  filter->SwitchOffClusterCorrection();
  //filter->SetDebugLevel(10);
  
  filter->SetVzCut(10.); // kEMCAL to have effect
  filter->SetEnergyCut(10.); // kEMCAL to have effect
  filter->SetNcellsCut(2);
  
  filter->SwitchOnFillTracks();
  filter->SwitchOnFillAODFile();
  
  filter->PrintInfo();   

  return filter;

}