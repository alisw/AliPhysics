#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisTaskTagAndProbe.h>

void Config_dsekihat_TAP_PbPb(
		AliAnalysisTaskTagAndProbe *task,
		const TString type,
		const TString cutname
		)
{

  LMEECutLib *lib = new LMEECutLib(cutname);

	//tag criteria should be tight.
  AliESDtrackCuts *esdTrackCuts     = lib->SetupESDtrackCuts();
  AliAnalysisCuts *filter_pid_tag   = lib->SetupPIDCutsForTag();
  AliAnalysisCuts *filter_track_tag = lib->SetupTrackCutsForTag();

	//tagged criteria should be tight. //i.e. high purity electron
	task->GetTagFilter()->AddCuts(esdTrackCuts);
	task->GetTagFilter()->AddCuts(filter_pid_tag);
	task->GetTagFilter()->AddCuts(filter_track_tag);

  AliAnalysisCuts *prefilter_pid_probe   = lib->SetupPrePIDCutsForProbe();//loose cut, eff = 1, only to reduce meaningless track loop
  AliAnalysisCuts *prefilter_track_probe = lib->SetupPreTrackCutsForProbe();//loose cut, eff = 1, only to reduce meaningless track loop

	AliAnalysisCuts *filter_pid_passingprobe   = lib->SetupPIDCutsForPassingProbe();//same as real anslysis
	AliAnalysisCuts *filter_track_passingprobe = lib->SetupTrackCutsForPassingProbe();//same as real anslysis

	if(type.Contains("PID",TString::kIgnoreCase)){
	    printf("PID efficiency will be evaluated");
		task->GetProbeFilter()->AddCuts(esdTrackCuts);
		task->GetProbeFilter()->AddCuts(filter_track_passingprobe);
		task->GetProbeFilter()->AddCuts(prefilter_pid_probe);

		task->GetPassingProbeFilter()->AddCuts(esdTrackCuts);
		task->GetPassingProbeFilter()->AddCuts(filter_track_passingprobe);//same as probe
		task->GetPassingProbeFilter()->AddCuts(filter_pid_passingprobe);//difference between probe and passingprobe should be only PID
	}
	else if(type.Contains("Track",TString::kIgnoreCase)){
	    printf("Track efficiency will be evaluated");
		task->GetProbeFilter()->AddCuts(prefilter_track_probe);
		task->GetProbeFilter()->AddCuts(filter_pid_passingprobe);

		task->GetPassingProbeFilter()->AddCuts(filter_track_passingprobe);//difference between probe and passingprobe should be only track cuts
		task->GetPassingProbeFilter()->AddCuts(filter_pid_passingprobe);
	}


	if(cutname.Contains("woPU")){
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
//______________________________________________________________________________________
//______________________________________________________________________________________

