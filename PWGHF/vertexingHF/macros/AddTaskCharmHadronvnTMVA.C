#include "AliHFQnVectorHandler.h"

AliAnalysisTaskSECharmHadronvnTMVA *AddTaskCharmHadronvnTMVA(TString BDTfilename="",
                                                    TString tenderTaskName = "HFTenderQnVectors",
                                                    int harm=2, 
                                                    TString filename="alien:///alice/cern.ch/user/f/fgrosa/DplustoKpipiCuts_3050_central_topod0cut_kINT7.root",
                                                    int decCh=AliAnalysisTaskSECharmHadronvnTMVA::kDplustoKpipi,
                                                    TString cutsobjname="AnalysisCuts", 
                                                    TString suffix="", 
                                                    int flagep=AliAnalysisTaskSECharmHadronvnTMVA::kVZERO/*kTPC,kTPCVZERO,kVZEROA,kVZEROC*/,
                                                    float minC=30.,
                                                    float maxC=50.,
                                                    int calibType=AliHFQnVectorHandler::kQnFrameworkCalib/*kQnCalib*/, 
                                                    TString OADBfilename="", 
                                                    int normMethod=AliHFQnVectorHandler::kQoverM/*kQoverQlength,kQoverSqrtM,kNone*/,
                                                    AliAnalysisTaskSECharmHadronvnTMVA::FlowMethod meth=AliAnalysisTaskSECharmHadronvnTMVA::kEP/*kSP,kEvShapeEP,kEPVsMass,kEvShapeEPVsMass*/,
                                                    AliAnalysisTaskSECharmHadronvnTMVA::q2Method q2meth=AliAnalysisTaskSECharmHadronvnTMVA::kq2TPC/*kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC*/,
                                                    int useAODProtection=1)
{
    //
    // Test macro for the AliAnalysisTaskSE for D-mesons vn analyses

    // Get the pointer to the existing analysis manager via the static access method.
    //============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AliAnalysisTaskSECharmHadronvnTMVA", "No analysis manager to connect to.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AliAnalysisTaskSECharmHadronvnTMVA", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(type.Contains("ESD")){
        ::Error("AliAnalysisTaskSECharmHadronvnTMVA", "This task requires to run on AOD");
        return NULL;
    }

    bool stdcuts=kFALSE;
    TFile* filecuts;
    if( filename.EqualTo("") ) {
        stdcuts=kTRUE;
    } 
    else {
        filecuts=TFile::Open(filename.Data());
        if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
            Printf("FATAL: Input file not found : check your cut object");
            return NULL;
        }
    }
    
    AliRDHFCuts *analysiscuts=0x0;
    int pdgmes=-1;
    //Analysis cuts
    
    if(decCh==AliAnalysisTaskSECharmHadronvnTMVA::kDplustoKpipi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsDplustoKpipi();
        analysiscuts->SetStandardCutsPbPb2011();
        } else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Dplus");
        pdgmes=411;
    } else if(decCh==AliAnalysisTaskSECharmHadronvnTMVA::kD0toKpi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsD0toKpi();
        analysiscuts->SetStandardCutsPbPb2011();
        } else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Dzero");
        pdgmes=421;
    } else if(decCh==AliAnalysisTaskSECharmHadronvnTMVA::kDstartoKpipi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsDStartoKpipi();
        analysiscuts->SetStandardCutsPbPb2011();
        } else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Dstar");
        pdgmes=413;
    }
    else if(decCh==AliAnalysisTaskSECharmHadronvnTMVA::kDstoKKpi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsDstoKKpi();
        //analysiscuts->SetStandardCutsPbPb2011(); //to be implemented in AliRDHFCutsDstoKKpi
        } else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Ds");
        pdgmes=431;
    }
    if(pdgmes==-1){
        Printf("FATAL: Wrong meson setting");
        return NULL;
    }
    if(!analysiscuts){
        Printf("FATAL: Specific AliRDHFCuts not found");
        return NULL;
    }
        
    // Analysis task
    AliAnalysisTaskSECharmHadronvnTMVA *vnTask = new AliAnalysisTaskSECharmHadronvnTMVA("HFvnAnalysis",analysiscuts,decCh);
    vnTask->SetHarmonic(harm);
    if(decCh == AliAnalysisTaskSECharmHadronvnTMVA::kDstartoKpipi) {
        vnTask->SetNMassBins(200);
    } else if(decCh == AliAnalysisTaskSECharmHadronvnTMVA::kDplustoKpipi || decCh == AliAnalysisTaskSECharmHadronvnTMVA::kD0toKpi) {
        vnTask->SetNMassBins(104);
        vnTask->SetMassLimits(0.2,pdgmes);
    } else if(decCh == AliAnalysisTaskSECharmHadronvnTMVA::kDstoKKpi) {
        vnTask->SetNMassBins(200);  // to be decided
    }
    vnTask->SetMinCentrality(minC);
    vnTask->SetMaxCentrality(maxC);
    vnTask->SetQnVectorDetConf(flagep);
    vnTask->SetNormMethod(normMethod);
    vnTask->SetCalibrationType(calibType);
    vnTask->SetAODMismatchProtection(useAODProtection);
    vnTask->SetOADBFileName(OADBfilename);
    vnTask->SetTenderTaskName(tenderTaskName);
    if(meth==AliAnalysisTaskSECharmHadronvnTMVA::kEvShapeEP || meth==AliAnalysisTaskSECharmHadronvnTMVA::kEvShapeSP)
        vnTask->SetqnMethod(q2meth);

    if(flagep==AliAnalysisTaskSECharmHadronvnTMVA::kTPC) {
        suffix+="TPC";
    } else if(flagep==AliAnalysisTaskSECharmHadronvnTMVA::kTPCVZERO) {
        suffix+="TPCVZERO";
    } else if(flagep==AliAnalysisTaskSECharmHadronvnTMVA::kVZERO) {
        suffix+="VZERO";
    } else if(flagep==AliAnalysisTaskSECharmHadronvnTMVA::kVZEROA) {
        suffix+="VZEROA";
    } else if(flagep==AliAnalysisTaskSECharmHadronvnTMVA::kVZEROC) {
        suffix+="VZEROC";
    } else if(flagep==AliAnalysisTaskSECharmHadronvnTMVA::kPosTPCVZERO) {
        suffix+="POSTPCVZERO";
    } else if(flagep==AliAnalysisTaskSECharmHadronvnTMVA::kNegTPCVZERO) {
        suffix+="NEGTPCVZERO";
    }
    vnTask->SetFlowMethod(meth);
    if(meth==AliAnalysisTaskSECharmHadronvnTMVA::kEP) {
        suffix+="_EP";
    } else if(meth==AliAnalysisTaskSECharmHadronvnTMVA::kSP) {
        suffix+="_SP";
    } else if(meth==AliAnalysisTaskSECharmHadronvnTMVA::kEvShapeEP) {
        suffix+="_EvShapeEP";
    } else if(meth==AliAnalysisTaskSECharmHadronvnTMVA::kEvShapeSP) {
        suffix+="_EvShapeSP";
    } else if(meth==AliAnalysisTaskSECharmHadronvnTMVA::kEPVsMass) {
        suffix+="_EPVsMass";
    } else if(meth==AliAnalysisTaskSECharmHadronvnTMVA::kEvShapeEPVsMass) {
        suffix+="_EvShapeEPVsMass";
    }
    
    Int_t Nptbins = analysiscuts->GetNPtBins();
    Float_t *ptbin = analysiscuts->GetPtBinLimits();
    TFile *fileBDT = TFile::Open(BDTfilename);
    if(!fileBDT ||(fileBDT&& !fileBDT->IsOpen())) ::Fatal("AddTaskCharmHadronvnTMVA", "BDT file not found : check your BDT object");
    
    TList *bdtlist = new TList();
    for(Int_t i=0;i<Nptbins;i++){
        TString BDTobjname = "BDT";
        BDTobjname += Form("1_%.0f_%.0f",ptbin[i],ptbin[i+1]);
        AliRDHFBDT *thisbdt = (AliRDHFBDT*)(fileBDT->Get(BDTobjname)->Clone(Form("_%s",BDTobjname.Data())));
        if(!thisbdt) ::Fatal("AddTaskCharmHadronvnTMVAnew", "Failed to find BDT named BDT1");
        //   std::cout<<thisbdt->GetDesc()<<endl;
        bdtlist->Add(thisbdt);
        
        TString BDT2objname1 = "BDT";
        TString BDT2objname2 = "BDT";
        TString BDT2objname3 = "BDT";
        BDT2objname1 += Form("2_%.0f_%.0f_0",ptbin[i],ptbin[i+1]);
        BDT2objname2 += Form("2_%.0f_%.0f_1",ptbin[i],ptbin[i+1]);
        BDT2objname3 += Form("2_%.0f_%.0f_2",ptbin[i],ptbin[i+1]);
        
        AliRDHFBDT *thisbdt2_0 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname1)->Clone(Form("_%s",BDT2objname1.Data())));
        AliRDHFBDT *thisbdt2_1 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname2)->Clone(Form("_%s",BDT2objname2.Data())));
        AliRDHFBDT *thisbdt2_2 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname3)->Clone(Form("_%s",BDT2objname3.Data())));
        if(!thisbdt2_0) ::Fatal("AddTaskCharmHadronvnTMVAnew", "Failed to find BDT named BDT2");
        if(!thisbdt2_1) ::Fatal("AddTaskCharmHadronvnTMVAnew", "Failed to find BDT named BDT2");
        if(!thisbdt2_2) ::Fatal("AddTaskCharmHadronvnTMVAnew", "Failed to find BDT named BDT2");
        //  std::cout<<thisbdt2_0->GetDesc()<<endl;
        bdtlist->Add(thisbdt2_0);
        bdtlist->Add(thisbdt2_1);
        bdtlist->Add(thisbdt2_2);
        
        
    }
    fileBDT->Close();
    vnTask->SetBDTList(bdtlist);
    
    
    mgr->AddTask(vnTask);
        
    // Create containers for input/output
        
    TString contname=Form("cinputvn%s",suffix.Data());
    AliAnalysisDataContainer *cinputvn = mgr->CreateContainer(contname.Data(),TChain::Class(),AliAnalysisManager::kInputContainer);
        
    TString outputfile = AliAnalysisManager::GetCommonFileName();
    TString outputhistos = Form("%s:PWGHF_D2H_HFvn_%s",outputfile.Data(),suffix.Data());
            
    contname=Form("coutputvn%s",suffix.Data());
    AliAnalysisDataContainer *coutputvn = mgr->CreateContainer(contname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputhistos.Data());
        
    contname=Form("cutobj%s",suffix.Data());
    AliAnalysisDataContainer *cutobj = mgr->CreateContainer(contname.Data(),AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer,outputhistos.Data());
    

        
    mgr->ConnectInput(vnTask,0,mgr->GetCommonInputContainer());
            
    mgr->ConnectOutput(vnTask,1,coutputvn);
        
    mgr->ConnectOutput(vnTask,2,cutobj);
    

        
    return vnTask;
}
