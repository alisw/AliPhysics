///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskMcKnoUeSystem Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
using namespace std;

AliAnalysisTaskMcKnoUeSystem* AddTaskMcKnoUeSystem(
    const Char_t *taskname="McKnoUe", Bool_t fUseMC=kFALSE, Bool_t IsMCclosureTest=kFALSE, Bool_t IsPythia=kFALSE, Bool_t IsppData=kFALSE,
    Bool_t IspPbData=kFALSE, Double_t minpT=0.5, Double_t PtLmin=5.0, Double_t PtLmax=40.0,
    Float_t fCutMaxFractionSharedTPCClusters=0.4, Float_t fCutMinRatioCrossedRowsOverFindableClustersTPC=0.8, Float_t fCutGeoNcrNclZone=3.0, 
    Float_t fCutGeoNcrNclLength=130.0, Bool_t fIsRequirementSPD=kTRUE, Float_t fCutMaxChi2PerClusterITS=36, Float_t fCutMaxChi2PerClusterTPC=4, 
    Float_t fCutMaxChi2TPCConstrainedVsGlobal=36, Float_t fCutMaxDCAToVertexZ=2,
    Int_t fCutMinNClusterTPC=50, Float_t fNchCutMaxChi2PerClusterTPC=4, Float_t fNchCutMaxDCAToVertexZ=3.2, Float_t fCutMaxDCAToVertexXY=2.4, 
    Bool_t fIsITSRefit=kTRUE, Bool_t fIsTPCRefit=kTRUE, Float_t fCutVertexZposition = 10)
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
    AliAnalysisTaskMcKnoUeSystem* taskUE = new AliAnalysisTaskMcKnoUeSystem("taskKno");
    if(!taskUE) return 0x0;
    taskUE->SetUseMC(fUseMC);
    taskUE->SetMCclosureTest(IsMCclosureTest);
    taskUE->SetParametrizationEfficiency(IsPythia);
    taskUE->SetParametrizationEfficiencyppdata(IsppData);
    taskUE->SetParametrizationEfficiencypPbdata(IspPbData);
    // add your task to the manager
    taskUE->SetPtMin(minpT);         //0.5  GeV/c
    taskUE->SetLeadingPtMin(PtLmin); //5.0  GeV/c
    taskUE->SetLeadingPtMax(PtLmax); //40.0 GeV/c
    // systematic
    // event selection for vertex position
    taskUE->SetUpVertexZposition(fCutVertexZposition); //10(default), 5, 15
    // track selection for leading particle
    taskUE->SetUpMaxFractionSharedTPCClusters(fCutMaxFractionSharedTPCClusters);  //0.4(default), 0.2, 1.0
    taskUE->SetUpMinRatioCrossedRowsOverFindableClustersTPC(fCutMinRatioCrossedRowsOverFindableClustersTPC);  //0.8(default), 0.7, 0.9
    taskUE->SetUpCutGeoNcrNcl(fCutGeoNcrNclZone, fCutGeoNcrNclLength);  //3.0cm(default), 2.0cm, 4.0cm;     130.0(default), 120.0, 140.0
    taskUE->SetUpClusterRequirementITS(fIsRequirementSPD);    //kTRUE--kAny(default), kFALSE--kNone
    taskUE->SetUpMaxChi2PerClusterITS(fCutMaxChi2PerClusterITS); //36(default), 25, 49
    taskUE->SetUpMaxChi2PerClusterTPC(fCutMaxChi2PerClusterTPC); //4(default), 3, 5
    taskUE->SetUpMaxChi2TPCConstrainedGlobal(fCutMaxChi2TPCConstrainedVsGlobal); //36(default), 25, 49
    taskUE->SetUpMaxDCAToVertexZ(fCutMaxDCAToVertexZ);  //2(default) 1cm, 5cm
    // track selection for Nch
    taskUE->SetUpMinNClustersTPC(fCutMinNClusterTPC);   //50(default), 40, 60
    taskUE->SetUpMaxChi2PerClusterTPC_nch(fNchCutMaxChi2PerClusterTPC);  //4(default),  3, 5
    taskUE->SetUpMaxDCAToVertexZ_nch(fNchCutMaxDCAToVertexZ);      //3.2(default), 2.2, 4.2
    taskUE->SetUpMaxDCAToVertexXY(fCutMaxDCAToVertexXY);          //2.4(default), 1.4, 3.4
    taskUE->SetUpRequireITSRefit(fIsITSRefit);  //kTRUE(default), kFALSE
    taskUE->SetUpRequireTPCRefit(fIsTPCRefit);  //kTRUE(default), kFALSE 

    cout<<"\n====================================================================================="<<endl;
    cout<<"pT min  = ("<<minpT<<")"<<endl;
    cout<<"leading pT min  = ("<<PtLmin<<")"<<endl;
    cout<<"leading pT max  = ("<<PtLmax<<")"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"*** event selection for vertex position ***"<<endl;
    cout<<"|Zvertex| <= ("<<fCutVertexZposition<<")"<<endl;
    cout<<"*** tarck selection for leading particle ***"<<endl;
    cout<<"cut-1-1  --- MaxFractionSharedTPCClusters = ("<<fCutMaxFractionSharedTPCClusters<<")"<<endl;
    cout<<"cut-1-2  --- MinRatioCrossedRowsOverFindableClustersTPC = ("<<fCutMinRatioCrossedRowsOverFindableClustersTPC<<")"<<endl;
    cout<<"cut-1-3  --- CutGeoNcrNcl = ("<<fCutGeoNcrNclZone<<", "<<fCutGeoNcrNclLength<<")"<<endl;
    cout<<"cut-1-4  --- ClusterRequirementITS = ("<<fIsRequirementSPD<<")"<<endl;
    cout<<"cut-1-5  --- MaxChi2PerClusterITS = ("<<fCutMaxChi2PerClusterITS<<")"<<endl;
    cout<<"cut-1-6  --- MaxChi2PerClusterTPC = ("<<fCutMaxChi2PerClusterTPC<<")"<<endl;
    cout<<"cut-1-7  --- MaxChi2TPCConstrainedGlobal = ("<<fCutMaxChi2TPCConstrainedVsGlobal<<")"<<endl;
    cout<<"cut-1-8  --- MaxDCAToVertexZ = ("<<fCutMaxDCAToVertexZ<<")"<<endl;
    cout<<"\n*** tarck selection for Nch ***"<<endl;
    cout<<"cut-2-1  --- MinNClustersTPC = ("<<fCutMinNClusterTPC<<")"<<endl;
    cout<<"cut-2-2  --- MaxChi2PerClusterTPC = ("<<fNchCutMaxChi2PerClusterTPC<<")"<<endl;
    cout<<"cut-2-3  --- MaxDCAToVertexZNch = ("<<fNchCutMaxDCAToVertexZ<<")"<<endl;
    cout<<"cut-2-4  --- MaxDCAToVertexXY = ("<<fCutMaxDCAToVertexXY<<")"<<endl;
    cout<<"cut-2-5  --- RequireITSRefit = ("<<fIsITSRefit<<")"<<endl;
    cout<<"cut-2-6  --- RequireTPCRefit = ("<<fIsTPCRefit<<")"<<endl;
    cout<<"=====================================================================================\n"<<endl;
    
    mgr->AddTask(taskUE);


    mgr->ConnectInput(taskUE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskUE,1,mgr->CreateContainer(Form("cList%s_%1.2f",taskname,minpT), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

    return taskUE;
}


