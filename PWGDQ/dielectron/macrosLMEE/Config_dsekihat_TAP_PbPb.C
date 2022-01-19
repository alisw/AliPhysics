#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisTaskTagAndProbe.h>

TString cutnames[] = {
  "anaFilter_pileup0_track0_Nsc01_pid0"
  ,"anaFilter_pileup0_track0_Nsc01_pid9"
  ,"anaFilter_pileup0_track0_Nsc01_TPCTOForTPCITS1"
  ,"anaFilter_pileup0_track0_Nsc01_TPCTOForTPCITS5"
  ,"anaFilter_pileup0_track0_Nsc01_TPCTOForTPCITS6"
};

const Int_t n = sizeof(cutnames)/sizeof(cutnames[0]);
Int_t GetN(){return n;}

void Config_dsekihat_TAP_PbPb(
    AliAnalysisTaskTagAndProbe *task,
    const TString type,
    const Bool_t isMC
    )
{
    //tag criteria should be tight.
    AliESDtrackCuts *esdTrackCuts     = LMEECutLib::SetupESDtrackCuts();
    AliAnalysisCuts *filter_track_tag = LMEECutLib::SetupTrackCutsForTag();
    AliAnalysisCuts *filter_pid_tag   = LMEECutLib::SetupPIDCutsForTag();

    AliAnalysisCuts *prefilter_pid_probe   = LMEECutLib::SetupPrePIDCutsForProbe();//loose cut, eff = 1, only to reduce meaningless track loop
    AliAnalysisCuts *prefilter_track_probe = LMEECutLib::SetupPreTrackCutsForProbe();//loose cut, eff = 1, only to reduce meaningless track loop

    task->GetTagFilter()->AddCuts(esdTrackCuts);
    task->GetTagFilter()->AddCuts(prefilter_track_probe);
    task->GetTagFilter()->AddCuts(filter_pid_tag);//tight electron

    task->GetProbeFilter()->AddCuts(esdTrackCuts);
    task->GetProbeFilter()->AddCuts(prefilter_track_probe);
    task->GetProbeFilter()->AddCuts(prefilter_pid_probe);

    if(type.Contains("PID",TString::kIgnoreCase)){
        printf("PID efficiency will be evaluated\n");
        for(Int_t i=0;i<n;i++){
            LMEECutLib *lib = new LMEECutLib(cutnames[i]);
            AliAnalysisFilter *ppfilter = new AliAnalysisFilter(cutnames[i],cutnames[i]);//passing probe filter
            AliAnalysisCuts *filter_pid_passingprobe = lib->SetupPIDCutsForPassingProbe();//same as real anslysis for PID efficiency
            ppfilter->AddCuts(esdTrackCuts);
            ppfilter->AddCuts(prefilter_track_probe);//same as probe
            ppfilter->AddCuts(filter_pid_passingprobe);//difference between probe and passingprobe should be only PID
            task->AddPassingProbeFilter(ppfilter);
            //delete ppfileter;
        }
    }
//    else if(type.Contains("Track",TString::kIgnoreCase)){
//        printf("Track efficiency will be evaluated\n");
//        AliAnalysisCuts *filter_track_passingprobe = lib->SetupTrackCutsForPassingProbe();//same as real anslysis for tracking efficiency
//        task->GetProbeFilter()->AddCuts(prefilter_track_probe);
//        //task->GetProbeFilter()->AddCuts(filter_pid_passingprobe);
//
//        task->GetPassingProbeFilter()->AddCuts(filter_track_passingprobe);//difference between probe and passingprobe should be only track cuts
//        //task->GetPassingProbeFilter()->AddCuts(filter_pid_passingprobe);
//    }
    if(!isMC){
        if(cutnames[0].Contains("pileup0")){
            printf("apply pileup cut!\n");
            TF1 *f1min = new TF1("f1min","pol2(0)",0,1e+8);
            f1min->SetNpx(1000);
            f1min->FixParameter(0,-3000);
            f1min->FixParameter(1,0.0099);
            f1min->FixParameter(2,9.42e-10);

            AliDielectronEventCuts*  pileupcuts = new AliDielectronEventCuts("pileupcuts","pileupcuts");
            pileupcuts->SetMinCorrCutFunction(f1min, AliDielectronVarManager::kNTPCclsEvent, AliDielectronVarManager::kNSDDSSDclsEvent);
            pileupcuts->Print();
            task->GetEventFilter()->AddCuts(pileupcuts);
        }
    }


}
//______________________________________________________________________________________
//______________________________________________________________________________________

