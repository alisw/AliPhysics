// $Id$

AliAnalysisTaskGammaHadron* AddTaskGammaHadron(
		Int_t       InputGammaOrPi0 = 0,                 //..gamma analysis=0, pi0 analyis=1, pi0 SB1=2, pi0 SB2=3,
		Bool_t      InputSeMe       = 0,                 //..same event=0 mixed event =1
		Bool_t      InputMCorData   = 0,                 // 0->MC, 1->Data
		UInt_t      evtTriggerType  = AliVEvent::kEMCEGA,//..use this type of events to combine gammas(trigger) with hadrons
		UInt_t      evtMixingType   = AliVEvent::kAnyINT,//..use only this type of events to fill your mixed event pool with tracks
		Bool_t      isRun2          = 1,                 //..changes some settigs and cuts depending on 2013 or 2015/2016 data
		Double_t    trackptcut      = 0.15,              //..
		Double_t    clusEcut        = 0.30,              //..
		Bool_t      SavePool        = 0,                 //..saves a mixed event pool to the output event
		const char *trackName       = "usedefault",
		const char *clusName        = "usedefault",
		const char *taskname        = "AliAnalysisTask",
		const char *suffix          = ""
)
{  
	  return AliAnalysisTaskGammaHadron::AddTaskGammaHadron(InputGammaOrPi0,InputSeMe,InputMCorData,
			  	  	  	  	  	  	  	  	  	  	  	   evtTriggerType,evtMixingType,isRun2,
			  	  	  	  	  	  	  	  	  	  	  	   trackptcut,clusEcut,SavePool,
			  	  	  	  	  	  	  	  	  	  	  	   trackName,clusName,taskname,suffix);
}
