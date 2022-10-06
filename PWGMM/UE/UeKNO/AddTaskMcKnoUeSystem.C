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
    Bool_t IspPbData=kFALSE, Double_t minpT=0.5, Double_t PtLmin=5, Double_t PtLmax=40,
    Float_t fCutMaxFractionSharedTPCClusters=0.4, Float_t fCutMinRatioCrossedRowsOverFindableClustersTPC=0.8, Float_t fCutGeoNcrNclZone=3, 
    Float_t fCutGeoNcrNclLength=130, Bool_t fIsRequirementSPD=kTRUE, Float_t fCutMaxChi2PerClusterITS=36, Float_t fCutMaxChi2PerClusterTPC=4, 
    Float_t fCutMaxChi2TPCConstrainedVsGlobal=36, Float_t fCutMaxDCAToVertexZ=2,
    Int_t fCutMinNClusterTPC=50, Float_t fCutMaxChi2PerClusterTPC_nch=4, Float_t fCutMaxDCAToVertexZ_nch=3.2, Float_t fCutMaxDCAToVertexXY=2.4, 
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
    taskUE->SetUpCutGeoNcrNcl(fCutGeoNcrNclZone, fCutGeoNcrNclLength);  //3cm(default), 2cm, 4cm;     130(default), 120, 140
    taskUE->SetUpClusterRequirementITS(fIsRequirementSPD);    //kTRUE--kAny(default), kFALSE--don't use that cut
    taskUE->SetUpMaxChi2PerClusterITS(fCutMaxChi2PerClusterITS); //36(default), 25, 49
    taskUE->SetUpMaxChi2PerClusterTPC(fCutMaxChi2PerClusterTPC); //4(default), 3, 5
    taskUE->SetUpMaxChi2TPCConstrainedGlobal(fCutMaxChi2TPCConstrainedVsGlobal); //36(default), 25, 49
    taskUE->SetUpMaxDCAToVertexZ(fCutMaxDCAToVertexZ);  //2(default) 1cm, 5cm
    // track selection for Nch
    taskUE->SetUpMinNClustersTPC(fCutMinNClusterTPC);   //50(default), 30, 70
    taskUE->SetUpMaxChi2PerClusterTPC_nch(fCutMaxChi2PerClusterTPC_nch);  //4(default),  3, 5
    taskUE->SetUpMaxDCAToVertexZ_nch(fCutMaxDCAToVertexZ_nch);      //3.2(default), 2.0, 4.0
    taskUE->SetUpMaxDCAToVertexXY(fCutMaxDCAToVertexXY);          //2.4(default), 1.0, 4.0
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
    cout<<"cut-trk-1  --- MaxFractionSharedTPCClusters = ("<<fCutMaxFractionSharedTPCClusters<<")"<<endl;
    cout<<"cut-trk-2  --- MinRatioCrossedRowsOverFindableClustersTPC = ("<<fCutMinRatioCrossedRowsOverFindableClustersTPC<<")"<<endl;
    cout<<"cut-trk-3  --- CutGeoNcrNcl = ("<<fCutGeoNcrNclZone<<", "<<fCutGeoNcrNclLength<<")"<<endl;
    cout<<"cut-trk-4  --- ClusterRequirementITS = ("<<fIsRequirementSPD<<")"<<endl;
    cout<<"cut-trk-5  --- MaxChi2PerClusterITS = ("<<fCutMaxChi2PerClusterITS<<")"<<endl;
    cout<<"cut-trk-6  --- MaxChi2PerClusterTPC = ("<<fCutMaxChi2PerClusterTPC<<")"<<endl;
    cout<<"cut-trk-7  --- MaxChi2TPCConstrainedGlobal = ("<<fCutMaxChi2TPCConstrainedVsGlobal<<")"<<endl;
    cout<<"cut-trk-8  --- MaxDCAToVertexZ = ("<<fCutMaxDCAToVertexZ<<")"<<endl;
    cout<<"\n*** tarck selection for Nch ***"<<endl;
    cout<<"cut-nch-1  --- MinNClustersTPC = ("<<fCutMinNClusterTPC<<")"<<endl;
    cout<<"cut-nch-2  --- MaxChi2PerClusterTPC = ("<<fCutMaxChi2PerClusterTPC_nch<<")"<<endl;
    cout<<"cut-nch-3  --- MaxDCAToVertexZNch = ("<<fCutMaxDCAToVertexZ_nch<<")"<<endl;
    cout<<"cut-nch-4  --- MaxDCAToVertexXY = ("<<fCutMaxDCAToVertexXY<<")"<<endl;
    cout<<"cut-nch-5  --- RequireITSRefit = ("<<fIsITSRefit<<")"<<endl;
    cout<<"cut-nch-6  --- RequireTPCRefit = ("<<fIsTPCRefit<<")"<<endl;
    cout<<"=====================================================================================\n"<<endl;

 
    Float_t  vtx  = fCutVertexZposition; //10(default), 5, 15
    Float_t  trk1 = fCutMaxFractionSharedTPCClusters;                //0.4(default), 0.2, 1.0
    Float_t  trk2 = fCutMinRatioCrossedRowsOverFindableClustersTPC;  //0.8(default), 0.7, 0.9
    Float_t  trk3_1 = fCutGeoNcrNclZone;       //3cm(default), 2cm, 4cm;
    Float_t  trk3_2 = fCutGeoNcrNclLength;     //130(default), 120, 140
    Bool_t   trk4 = fIsRequirementSPD;         //kTRUE--kAny(default), kFALSE--don't use that cut
    Float_t  trk5 = fCutMaxChi2PerClusterITS;  //36(default), 25, 49
    Float_t  trk6 = fCutMaxChi2PerClusterTPC;  //4(default), 3, 5
    Float_t  trk7 = fCutMaxChi2TPCConstrainedVsGlobal;  //36(default), 25, 49
    Float_t  trk8 = fCutMaxDCAToVertexZ;                //2(default) 1cm, 5cm
    Int_t    nch1 = fCutMinNClusterTPC;          //50(default), 40, 60
    Float_t  nch2 = fCutMaxChi2PerClusterTPC_nch; //4(default),  3, 5
    Float_t  nch3 = fCutMaxDCAToVertexZ_nch;      //3.2(default), 2, 4
    Float_t  nch4 = fCutMaxDCAToVertexXY;        //2.4(default), 1, 4
    Bool_t   nch5 = fIsITSRefit; //kTRUE(default), kFALSE    
    Bool_t   nch6 = fIsTPCRefit; //kTRUE(default), kFALSE
 
    const Char_t *cutName = "default";
    //vtx cut
    if (vtx<10) cutName = "vtx_low";
    else if (vtx>10) cutName = "vtx_high";
    
    //track cut for leading particle
    if (trk1<0.39) cutName = "maxFracTPC_low";
    else if (trk1>0.41) cutName = "maxFracTPC_high";
    
    if (trk2<0.79) cutName = "minRatioTPC_low";
    else if (trk2>0.81) cutName = "minRatioTPC_high";
    
    if (trk3_1<3) cutName = "geoZone_low";
    else if (trk3_1>3) cutName = "geoZone_high";
    
    if (trk3_2<130) cutName = "geoLength_low";
    else if (trk3_2>130) cutName = "geoLength_high";
    
    if (!trk4) cutName = "reqSPD_off";
    
    if (trk5<36) cutName = "maxChiqITS_low";
    else if (trk5>36) cutName = "maxChiqITS_high";
    
    if (trk6<4) cutName = "maxChiqTPC_low";
    else if (trk6>4) cutName = "maxChiqTPC_high";
    
    if (trk7<36) cutName = "maxChiqTPCglobal_low";
    else if (trk7>36) cutName = "maxChiqTPCglobal_high";
    
    if (trk8<2) cutName = "maxDCAz_low";
    else if (trk8>2) cutName = "maxDCAz_high";
    
    //track cut for nch
    if (nch1<50) cutName = "minClusterTPC_low";
    else if (nch1>50) cutName = "minClusterTPC_high";
    
    if (nch2<4) cutName = "maxChiqTPC_nch_low";
    else if (nch2>4) cutName = "maxChiqTPC_nch_high";
       
    if (nch3<3.1) cutName = "maxDCAz_nch_low";
    else if (nch3>3.3) cutName = "maxDCAz_nch_high";
    
    if (nch4<2.3) cutName = "maxDCAxy_low";
    else if (nch4>2.5) cutName = "maxDACxy_high";
       
    if (!nch5) cutName = "ITSrefit_off";
    if (!nch6) cutName = "TPCrefit_off";
    
    cout<<"\ncutName = "<<cutName<<"\n"<<endl;
    

    mgr->AddTask(taskUE);

    mgr->ConnectInput(taskUE,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskUE,1,mgr->CreateContainer(Form("cList%s_%s_%1.2f", taskname, cutName, minpT), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

    return taskUE;
}


