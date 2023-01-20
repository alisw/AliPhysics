///////////////////////////////////////////////////////////////////
//                                                               //
//            AddTaskChargedVsRT Macro to run on grids               //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskChargedVsRT* AddTaskChargedVsRT(const Char_t* taskname="ChargedVsRT",
                                   Bool_t  useMC  = kTRUE,
                                   Bool_t performMCclosuretest = kFALSE,
                                   Bool_t IsHyAna = kFALSE,
                                   Bool_t MultV0 = kFALSE,
                                   Double_t minpT=0.5,
                                   Double_t PtLmin = 1.0,
                                   Double_t PtLmax = 15.0,
                                   Int_t NchNbins = 200,
                                   Double_t maxbinNch = 199.5,
                                   Bool_t TPCclustersVar1 = kFALSE,
                                   Bool_t TPCclustersVar2 = kFALSE,
                                   Bool_t NcrVar1 = kFALSE,
                                   Bool_t NcrVar2 = kFALSE,
                                   Bool_t ChisqTPCVar1 = kFALSE,
                                   Bool_t ChisqTPCVar2 = kFALSE,
                                   Bool_t ChisqITSVar1 = kFALSE,
                                   Bool_t ChisqITSVar2 = kFALSE,
                                   Bool_t ChisqITSmTPCVar1 = kFALSE,
                                   Bool_t ChisqITSmTPCVar2 = kFALSE,
                                   Bool_t DcazVar1 = kFALSE,
                                   Bool_t DcazVar2 = kFALSE,
                                   Bool_t GeoTPCVar1 = kFALSE,
                                   Bool_t GeoTPCVar2 = kFALSE,
                                   Bool_t GeoTPCVar3 = kFALSE,
                                   Bool_t GeoTPCVar4 = kFALSE,
                                   Bool_t SPDreqVar1 = kFALSE)
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
    AliAnalysisTaskChargedVsRT* taskKno = new AliAnalysisTaskChargedVsRT("taskKno");
    if(!taskKno) return 0x0;
    taskKno->SetUseMC(useMC);
    taskKno->SetMCclosureTest(performMCclosuretest);
    // add your task to the manager
    taskKno->SetPtMin(minpT);
    taskKno->SetIsHybridAnalysis(IsHyAna);
    taskKno->SetMultPercenV0(MultV0);
    taskKno->SetLeadingPtMin(PtLmin);
    taskKno->SetLeadingPtMax(PtLmax);
    taskKno->SetNchNbin(NchNbins);
    taskKno->SetNchBinMax(maxbinNch);
    // Systematic -------------------------------
    taskKno->SetTPCclustersVar1(TPCclustersVar1);
    taskKno->SetTPCclustersVar2(TPCclustersVar2);
    taskKno->SetNcrVar1(NcrVar1);
    taskKno->SetNcrVar2(NcrVar2);
    taskKno->SetChisqTPCVar1(ChisqTPCVar1);
    taskKno->SetChisqTPCVar2(ChisqTPCVar2);
    taskKno->SetChisqITSVar1(ChisqITSVar1);
    taskKno->SetChisqITSVar2(ChisqITSVar2);
    //taskKno->SetChisqITSmTPCVar1(ChisqITSmTPCVar1);
    //taskKno->SetChisqITSmTPCVar2(ChisqITSmTPCVar2);
    taskKno->SetDcazVar1(DcazVar1);
    taskKno->SetDcazVar2(DcazVar2);
    //taskKno->SetGeoTPCVar1(GeoTPCVar1);
    //taskKno->SetGeoTPCVar2(GeoTPCVar2);
    //taskKno->SetGeoTPCVar3(GeoTPCVar3);
    //taskKno->SetGeoTPCVar4(GeoTPCVar4);
    //taskKno->SetSPDreqVar1(SPDreqVar1);
    // Systematic -------------------------------
    mgr->AddTask(taskKno);

    mgr->ConnectInput(taskKno,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(taskKno,1,mgr->CreateContainer(Form("cList%s_%1.2f_%1.1f_%1.1f",taskname,minpT,PtLmin,PtLmax), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:%s", AliAnalysisManager::GetCommonFileName(),taskname)));

    return taskKno;
}
