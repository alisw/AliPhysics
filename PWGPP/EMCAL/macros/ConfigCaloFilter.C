///
/// \file ConfigCaloFilter.C
/// \ingroup EMCALPerformanceMacros
/// \brief Configuration analysis task filtering events and calorimeter data into AOD format.
///
/// Basic example of configuration, using alternative configuration way.
/// Configuration is applied at the initialization time of the task.
/// The suggested macro is though AddTaskCaloFilter.C
///
/// \author : Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///
AliAnalysisTaskCaloFilter * ConfigCaloFilter()
{
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