///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskMcKnoUe Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskMcKnoUeSyst* AddTaskMcKnoUeSyst(const Char_t* taskname="McKnoUe", Bool_t  useMC  = kTRUE, Bool_t performMCclosuretest = kFALSE, Bool_t isPythia=kTRUE, Double_t minpT=0.5,Bool_t TPCclustersLow = kFALSE, Bool_t TPCclustersHigh = kFALSE, Bool_t NcrLow = kFALSE, Bool_t NcrHigh = kFALSE, Bool_t ChisqTPCLow = kFALSE, Bool_t ChisqTPCHigh = kFALSE, Bool_t ChisqITSLow = kFALSE, Bool_t ChisqITSHigh = kFALSE, Bool_t ChisqITSmTPCLow = kFALSE, Bool_t ChisqITSmTPCHigh = kFALSE, Bool_t DcazLow = kFALSE, Bool_t DcazHigh = kFALSE, Bool_t GeoTPCLow1 = kFALSE, Bool_t GeoTPCLow2 = kFALSE, Bool_t GeoTPCHigh1 = kFALSE, Bool_t GeoTPCHigh2 = kFALSE, Bool_t SPDreqVar1 = kFALSE,Bool_t VertexZCutLow = kFALSE,Bool_t VertexZCutHigh = kFALSE)
{
    // get the manager via the static access member. since it's static, you don't need
    // an instance of the class to call the function


    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler this handler is part of the managing system and feeds events to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // now you create an instance of your task
    AliAnalysisTaskMcKnoUeSyst* taskUE = new AliAnalysisTaskMcKnoUeSyst("taskKno");
    if(!taskUE) return 0x0;
    taskUE->SetUseMC(useMC);
    taskUE->SetMCclosureTest(performMCclosuretest);
    taskUE->SetParametrizationEfficiency(isPythia);
    // add your task to the manager
    taskUE->SetPtMin(minpT);
    
    // Systematic -------------------------------
   // taskKno->SetNchTScut(IsTPConly);
    
    taskUE->SetTPCclustersLow(TPCclustersLow);
    taskUE->SetTPCclustersHigh(TPCclustersHigh);
    taskUE->SetNcrLow(NcrLow);
    taskUE->SetNcrHigh(NcrHigh);
    taskUE->SetChisqTPCLow(ChisqTPCLow);
    taskUE->SetChisqTPCHigh(ChisqTPCHigh);
    taskUE->SetChisqITSLow(ChisqITSLow);
    taskUE->SetChisqITSHigh(ChisqITSHigh);
    taskUE->SetChisqITSmTPCLow(ChisqITSmTPCLow);
    taskUE->SetChisqITSmTPCHigh(ChisqITSmTPCHigh);
    taskUE->SetDcazLow(DcazLow);
    taskUE->SetDcazHigh(DcazHigh);
    taskUE->SetGeoTPCLow1(GeoTPCLow1);
    taskUE->SetGeoTPCLow2(GeoTPCLow2);
    taskUE->SetGeoTPCHigh1(GeoTPCHigh1);
    taskUE->SetGeoTPCHigh2(GeoTPCHigh2);
    taskUE->SetSPDreqVar1(SPDreqVar1);
    taskUE->SetVertexZCutLow(VertexZCutLow);
    taskUE->SetVertexZCutHigh(VertexZCutHigh);
    // Systematic -------------------------------
    
    char const* TrackCutName;
    
    if(TPCclustersLow) {TrackCutName ="TPCclustersLow02";
        Printf("&&&&&&&&&&&&&&&&&&& FracSharedTPCClus02  is 0.2 ################");
    }//1
    else if(TPCclustersHigh) {TrackCutName ="TPCclustersHigh1";}//2
    else if(NcrLow) {TrackCutName ="NcrLow07";}//3
    else if(NcrHigh) {TrackCutName ="NcrHigh09";}//4
    else if(ChisqTPCLow) {TrackCutName ="ChisqTPCLow3";}//5
    else if(ChisqTPCHigh) {TrackCutName ="ChisqTPCHigh5";}//6
    
    else if(ChisqITSLow) {TrackCutName ="ChisqITSLow25";}//7
    else if(ChisqITSHigh) {TrackCutName ="ChisqITSHigh49";}//8
    
    else if(ChisqITSmTPCLow) {TrackCutName ="ChisqITSmTPCLow25";}//9
    else if(ChisqITSmTPCHigh) {TrackCutName ="ChisqITSmTPCHigh49";}//10
    
    else if(DcazLow) {TrackCutName ="DcazLow1";}//11
    else if(DcazHigh) {TrackCutName ="DcazHigh5";}//12
    
    else if(GeoTPCLow1) {TrackCutName ="GeoTPCLow12";}//13
    else if(GeoTPCHigh1) {TrackCutName ="GeoTPCHigh14";}//14
    
    else if(GeoTPCLow2) {TrackCutName ="GeoTPCLow2120";}//15
    else if(GeoTPCHigh2) {TrackCutName ="GeoTPCHigh2140";}//16
    
    else if (SPDreqVar1){TrackCutName ="fSPDreqVar1";}//17
    
    else if(VertexZCutLow) {TrackCutName ="VertexZCutLow5";}//18
    else if(VertexZCutHigh) {TrackCutName ="VertexZCutHigh15";}//19
    else{TrackCutName ="Default";}
    
    
    mgr->AddTask(taskUE);
   
    mgr->ConnectInput(taskUE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskUE,1,mgr->CreateContainer(Form("cList%s_%s_%1.2f",taskname,TrackCutName,minpT), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

    return taskUE;
}
