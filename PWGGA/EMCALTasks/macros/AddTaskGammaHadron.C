// $Id$

AliAnalysisTaskGammaHadron* AddTaskGammaHadron(
		Int_t       InputGammaOrPi0 = 0,                 //..gamma analysis=0, pi0 analyis=1, pi0 SB1=2, pi0 SB2=3,
		Int_t       InputSeMe       = 0,                 //..same event=0 mixed event =1 mixed trigger = 2
		Bool_t      InputMCorData   = 0,                 // 0->Data, 1->MC
		UInt_t      evtTriggerType  = AliVEvent::kEMCEGA,//..use this type of events to combine gammas(trigger) with hadrons
		UInt_t      evtMixingType   = AliVEvent::kAnyINT,//..use only this type of events to fill your mixed event pool with tracks
		Bool_t      isRun2          = 1,                 //..changes some settigs and cuts depending on 2013 or 2015/2016 data
		Double_t    trackptcut      = 0.15,              //..
		Double_t    clusEcut        = 0.30,              //..
		Bool_t      SavePool        = 0,                 //..saves a mixed event pool to the output event
		const char *trackName       = "usedefault",
		const char *clusName        = "usedefault",
		const char *taskname        = "AliAnalysisTask",
		const char *poolFilePath    = "",                 //..Path to file with pool manager. e.g. /alice/cern.ch/user/ ... /TriggerPool.root
		const char *poolFileName    = "TriggerPool.root",
		const char *suffix          = ""
)
{
		AliAnalysisTaskGammaHadron * task = AliAnalysisTaskGammaHadron::AddTaskGammaHadron(InputGammaOrPi0,InputSeMe,InputMCorData,
			  	  	  	  	  	  	  	  	  	  	  	   evtTriggerType,evtMixingType,isRun2,
			  	  	  	  	  	  	  	  	  	  	  	   trackptcut,clusEcut,SavePool,
			  	  	  	  	  	  	  	  	  	  	  	   trackName,clusName,taskname,suffix);

		if (InputSeMe == 2) {
			printf("AddTask: Adding mixed event pool\n");
			// For Mixed Trigger Mode, copy alien file
			printf("Copying file from alien:/%s/%s\n",poolFilePath,poolFileName);
			gSystem->Exec(Form("alien_cp alien:/%s/%s .",poolFilePath,poolFileName));
			TFile * fPoolFile = TFile::Open(poolFileName);
			if (fPoolFile) {
				AliEventPoolManager * fPool = 0x0; 
				fPool = (AliEventPoolManager *) fPoolFile->Get("AliEventPoolManager");
				if (fPool) {
					task->SetExternalEventPoolManager(fPool);
				} else printf("AddTask: Could not find pool manager\n");
				fPoolFile->Close();
			}
		}

		return task;
//	  return AliAnalysisTaskGammaHadron::AddTaskGammaHadron(InputGammaOrPi0,InputSeMe,InputMCorData,
}
