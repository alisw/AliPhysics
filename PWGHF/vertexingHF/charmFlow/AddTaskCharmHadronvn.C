#include "AliHFQnVectorHandler.h"

AliAnalysisTaskSECharmHadronvn *AddTaskCharmHadronvn(TString tenderTaskName = "HFTenderQnVectors",
                                                    int harm=2, 
                                                    TString filename="alien:///alice/cern.ch/user/f/fgrosa/DplustoKpipiCuts_3050_central_topod0cut_kINT7.root",
                                                    int decCh=AliAnalysisTaskSECharmHadronvn::kDplustoKpipi,
                                                    TString cutsobjname="AnalysisCuts", 
                                                    TString suffix="", 
                                                    int flagep=AliAnalysisTaskSECharmHadronvn::kVZERO/*kTPC,kTPCVZERO,kVZEROA,kVZEROC*/,
                                                    float minC=30.,
                                                    float maxC=50.,
                                                    int calibType=AliHFQnVectorHandler::kQnFrameworkCalib/*kQnCalib*/, 
                                                    TString OADBfilename="", 
                                                    int normMethod=AliHFQnVectorHandler::kQoverM/*kQoverQlength,kQoverSqrtM,kNone*/,
                                                    AliAnalysisTaskSECharmHadronvn::FlowMethod meth=AliAnalysisTaskSECharmHadronvn::kEP/*kSP,kEvShapeEP,kEPVsMass,kEvShapeEPVsMass*/, 
                                                    AliAnalysisTaskSECharmHadronvn::q2Method q2meth=AliAnalysisTaskSECharmHadronvn::kq2TPC/*kq2PosTPC,kq2NegTPC,kq2VZERO,kq2VZEROA,kq2VZEROC*/, 
                                                    int useAODProtection=1)
{
    //
    // Test macro for the AliAnalysisTaskSE for D-mesons vn analyses

    // Get the pointer to the existing analysis manager via the static access method.
    //============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        ::Error("AliAnalysisTaskSECharmHadronvn", "No analysis manager to connect to.");
        return NULL;
    }
    if (!mgr->GetInputEventHandler()) {
        ::Error("AliAnalysisTaskSECharmHadronvn", "This task requires an input event handler");
        return NULL;
    }
    TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if(type.Contains("ESD")){
        ::Error("AliAnalysisTaskSECharmHadronvn", "This task requires to run on AOD");
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
    
    if(decCh==AliAnalysisTaskSECharmHadronvn::kDplustoKpipi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsDplustoKpipi();
        analysiscuts->SetStandardCutsPbPb2011();
        } else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Dplus");
        pdgmes=411;
    } else if(decCh==AliAnalysisTaskSECharmHadronvn::kD0toKpi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsD0toKpi();
        analysiscuts->SetStandardCutsPbPb2011();
        } else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Dzero");
        pdgmes=421;
    } else if(decCh==AliAnalysisTaskSECharmHadronvn::kDstartoKpipi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsDStartoKpipi();
        analysiscuts->SetStandardCutsPbPb2011();
        } else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Dstar");
        pdgmes=413;
    }
    else if(decCh==AliAnalysisTaskSECharmHadronvn::kDstoKKpi) {
        if(stdcuts) {
        analysiscuts = new AliRDHFCutsDstoKKpi();
        //analysiscuts->SetStandardCutsPbPb2011(); //to be implemented in AliRDHFCutsDstoKKpi
        } else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutsobjname);
        suffix.Prepend("Ds");
        pdgmes=431;
    }
    else if(decCh==AliAnalysisTaskSECharmHadronvn::kLctopK0S) {
        if(stdcuts) {
         analysiscuts = new AliRDHFCutsLctoV0();
         analysiscuts->SetStandardCutsPbPb2011();
        }else analysiscuts = (AliRDHFCutsLctoV0*)filecuts->Get(cutsobjname);
        suffix.Prepend("Lc2V0");
        pdgmes=4122;
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
    AliAnalysisTaskSECharmHadronvn *vnTask = new AliAnalysisTaskSECharmHadronvn("HFvnAnalysis",analysiscuts,decCh);
    vnTask->SetHarmonic(harm);
    if(decCh == AliAnalysisTaskSECharmHadronvn::kDstartoKpipi) {
        vnTask->SetNMassBins(200);
    } else if(decCh == AliAnalysisTaskSECharmHadronvn::kDplustoKpipi || decCh == AliAnalysisTaskSECharmHadronvn::kD0toKpi) {
        vnTask->SetNMassBins(104);
        vnTask->SetMassLimits(0.2,pdgmes);
    } else if(decCh == AliAnalysisTaskSECharmHadronvn::kDstoKKpi) {
        vnTask->SetNMassBins(200);  // to be decided
    } else if(decCh == AliAnalysisTaskSECharmHadronvn::kLctopK0S) {
        vnTask->SetNMassBins(250);
        vnTask->SetMassLimits(0.25,pdgmes);
    }
    vnTask->SetMinCentrality(minC);
    vnTask->SetMaxCentrality(maxC);
    vnTask->SetQnVectorDetConf(flagep);
    vnTask->SetNormMethod(normMethod);
    vnTask->SetCalibrationType(calibType);
    vnTask->SetAODMismatchProtection(useAODProtection);
    vnTask->SetOADBFileName(OADBfilename);
    vnTask->SetTenderTaskName(tenderTaskName);
    if(meth==AliAnalysisTaskSECharmHadronvn::kEvShapeEP || meth==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) 
        vnTask->SetqnMethod(q2meth);

    if(flagep==AliAnalysisTaskSECharmHadronvn::kTPC) {
        suffix+="TPC";
    } else if(flagep==AliAnalysisTaskSECharmHadronvn::kTPCVZERO) {
        suffix+="TPCVZERO";
    } else if(flagep==AliAnalysisTaskSECharmHadronvn::kVZERO) {
        suffix+="VZERO";
    } else if(flagep==AliAnalysisTaskSECharmHadronvn::kVZEROA) {
        suffix+="VZEROA";
    } else if(flagep==AliAnalysisTaskSECharmHadronvn::kVZEROC) {
        suffix+="VZEROC";
    } else if(flagep==AliAnalysisTaskSECharmHadronvn::kPosTPCVZERO) {
        suffix+="POSTPCVZERO";
    } else if(flagep==AliAnalysisTaskSECharmHadronvn::kNegTPCVZERO) {
        suffix+="NEGTPCVZERO";
    }
    vnTask->SetFlowMethod(meth);
    if(meth==AliAnalysisTaskSECharmHadronvn::kEP) {
        suffix+="_EP";
    } else if(meth==AliAnalysisTaskSECharmHadronvn::kSP) {
        suffix+="_SP";
    } else if(meth==AliAnalysisTaskSECharmHadronvn::kEvShapeEP) {
        suffix+="_EvShapeEP";
    } else if(meth==AliAnalysisTaskSECharmHadronvn::kEvShapeSP) {
        suffix+="_EvShapeSP";
    } else if(meth==AliAnalysisTaskSECharmHadronvn::kEPVsMass) {
        suffix+="_EPVsMass";
    } else if(meth==AliAnalysisTaskSECharmHadronvn::kEvShapeEPVsMass) {
        suffix+="_EvShapeEPVsMass";
    }
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
