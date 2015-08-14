void AddTaskSpectraBothManyTasks(Bool_t mc=kFALSE,
					     Double_t CentCutMin=0,
					     Double_t CentCutMax=100,
					     Double_t QvecCutMin=0,
					     Double_t QvecCutMax=100,
					     Double_t EtaMin=-0.8,
					     Double_t EtaMax=0.8,
					     Double_t Nsigmapid=3.,
					     Double_t pt=5.,
					     Double_t p=5.,
					     Double_t ymin=-0.5,
					     Double_t ymax=.5, 	
					     Double_t ptTofMatch=.6,
					     UInt_t trkbit=1,
					     UInt_t trkbitQVector=1,
					     Bool_t UseCentPatchAOD049=kFALSE,
					     Double_t DCA=100000,
					     UInt_t minNclsTPC=70,
					     Int_t nrebin=0,
					     TString centestimator="V0M",
					     Int_t pidmethod=3,
					     TString taskname="",
					     Int_t minmul=-1,
	                                     Int_t maxmul=-1,
					     Double_t etamulcut=0.8,	
					     Int_t minrun=-1,
					     Int_t maxrun=-1,
					     Bool_t usedAdditionalCuts=kFALSE,	
					     Bool_t makeQAhisto=kFALSE,
					     Bool_t dotheMCLoopAfterEventCuts=kFALSE,
					     Bool_t dotheeventcutsinmultselection=kFALSE,
					     Double_t ptTofMatchpi=0.6,
					     Double_t ptTofMatchka=0.5,
					     Double_t ptTofMatchpr=0.8)
{
	gROOT->LoadMacro("$ALICE_PHYSICS/../src/PWGLF/SPECTRA/PiKaPr/TestAOD/AddTaskSpectraBoth.C");
	

	for (int i=0;i<7;i++)
	{
		TString opt=Form("Cent%.3fto%.3f_QVec%.1fto%.1f_Eta%.1fto%.1f_%.1fSigmaPID_TrBit%dEst_%s_Pid_%d_Y%.1fto%.1f_tasknumber%d",CentCutMin,CentCutMax,QvecCutMin,QvecCutMax,EtaMin,EtaMax,Nsigmapid,trkbit,centestimator.Data(),pidmethod,ymin,ymax,i);
		AliAnalysisTaskSpectraBoth* task=AddTaskSpectraBoth(mc,CentCutMin,CentCutMax,QvecCutMin,QvecCutMax,EtaMin,EtaMax,Nsigmapid,pt,p,ymin,ymax,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,nrebin,centestimator,pidmethod,opt);
		if(minmul>-1&&maxmul>-1)
		{
			task->GetEventCuts()->SetMultiplicityCut(minmul,maxmul);
			task->GetEventCuts()->SetEtaRangeforMultiplictyCut(etamulcut);
		}
		AliESDtrackCuts* cut=0x0;
		if(i!=1)
			cut=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(0,1);
		else
			cut=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(0,0);

		if(i==2)
			cut->SetMaxDCAToVertexZ(1.0);
		if(i==3)
			cut->SetMaxDCAToVertexZ(3.0);
		if(i==4)
			cut->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);	
		
		task->SetAliESDtrackCuts(cut);
		task->GetTrackCuts()->SetUsedAdditionalCuts(usedAdditionalCuts);	
		
		if(i==0)
			task->SetMakePIDQAHisto(makeQAhisto);
		else
			task->SetMakePIDQAHisto(kFALSE);
		
		task->SetdotheMCLoopAfterEventCuts(dotheMCLoopAfterEventCuts);
		task->GetEventCuts()->SetDotheeventcutsinmultselection(dotheeventcutsinmultselection);
	
		if(minrun>0&&maxrun>0)
			task->GetEventCuts()->SetRunNumberRange(minrun,maxrun);
		if(ptTofMatchpi>0.0&&ptTofMatchka>0.0&&ptTofMatchpr>0.0)
			task->GetTrackCuts()->SetPtTOFMatchingPartDepended(ptTofMatchpi,ptTofMatchka,ptTofMatchpr);
	



	}
	
}
