/*************************************************************************
   Configure Macro to set all cuts for Fluctuation Study 
   AliEbyEEventSelector | Implimented
   Origin: Satyajit Jena
   Date:Mon Nov 22 22:15:40 CET 2010
 *************************************************************************/

AliEbyEEventSelector *GetEbyEAnalysisBaseObject(const char* analysisLevel = "ESD",
						const char* esdAnalysisType = "kGlobal",
						Bool_t kUseOnlineTrigger = kFALSE,
						Bool_t kUseOfflineTrigger = kFALSE,
						const char* CentralityType = "HardType") 

{
  
  AliEbyEEventSelector *base = new AliEbyEEventSelector();
  base->SetAnalysisLevel(analysisLevel);
  
  base->SetCentralityBin(50);
  if(CentralityType == "HardType") {
    base->SetCentralityType(AliEbyEEventSelector::kHardFlat);
    base->SetCentralityEstimator("V0M");

    printf(" I am Inside centrality Flat \n");
  }

  if (CentralityType == "Central") {
      base->SetCentralityType(AliEbyEEventSelector::kFlat);
     base->SetCentralityEstimator("V0M");
     base->SetCentralityInputFiles("$ALICE_ROOT/ANALYSIS/macros/AliCentralityBy1D_LHC10g2a_100.root",
				   "$ALICE_ROOT/ANALYSIS/macros/AliCentralityByFunction_LHC10g2a_100.root");
  }

  if(analysisLevel == "ESD") {
    /*  if(kAnalyzeMC)
        base->SetTriggerMode(AliEbyEEventSelector::kMB2);
        if(kUseOnlineTrigger) base->UseOnlineTrigger();
        if(kUseOfflineTrigger) base->OfflineTriggerInit();
      
    */

    switch(esdAnalysisType) {
    case "TPC":
      base->SetAnalysisMode(AliEbyEEventSelector::kTPC);
      base->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
      break;
    
    case "ITS":
      base->SetAnalysisMode(AliEbyEEventSelector::kITS);
      base->SetPhaseSpace(9, -0.9, 0.9, 6, 0.45, 1.05);
      break;
      
    case "TPCnITS":
      base->SetAnalysisMode(AliEbyEEventSelector::kTPCnITS);
      base->SetPhaseSpace(9, -0.9, 0.9, 6, 0.45, 1.05);
      break;
      
    case "Global":
      base->SetAnalysisMode(AliEbyEEventSelector::kGlobal);
      base->SetPhaseSpace(20, -1.0, 1.0, 48, 0.3, 1.5);
      break;

    case "Forward":
      base->SetAnalysisMode(AliEbyEEventSelector::kForward);
      base->SetPhaseSpace(20, -1.0, 1.0, 48, 0.3, 1.5);
      break;
      
    default:
      break;
    }
   base->SetAcceptedVertexDiamond(1.,1.,10.);
     
  }
  if(analysisLevel == "MC") 
    base->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
  if(analysisLevel == "AOD")
    base->SetPhaseSpace(10, -0.5, 0.5, 16, 0.5, 0.9);
  
  return base;
}
