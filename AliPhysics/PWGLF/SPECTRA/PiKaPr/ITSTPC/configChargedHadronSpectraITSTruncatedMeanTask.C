AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* GetAliAnalysisChargedHadronSpectraITSTruncatedMeanTask(int usemc=0)
{
	
	AliAnalysisChargedHadronSpectraITSTruncatedMeanTask* task1 = new AliAnalysisChargedHadronSpectraITSTruncatedMeanTask("AliAnalysisChargedHadronSpectraITSTruncatedMeanTask");
	AliESDtrackCuts*  cuts1 = new AliESDtrackCuts("cuts1","cuts1");
 	cuts1->SetRequireTPCRefit(kTRUE); // but only for pass4 or later, pass2 without requiring refit !!!!
	cuts1->SetRequireITSRefit(kTRUE);
 	cuts1->SetAcceptKinkDaughters(kFALSE);
 	cuts1->SetMinNClustersTPC(70);
 	cuts1->SetMaxChi2PerClusterTPC(4);
	cuts1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);	
	cuts1->SetDCAToVertex2D(kFALSE);
  	cuts1->SetRequireSigmaToVertex(kFALSE);
	cuts1->SetMaxDCAToVertexZ(2);
	AliESDpidCuts* tpcpidcut=new AliESDpidCuts();
	//tpcpidcut->SetTPCnSigmaCut(2,3);//pion
//	tpcpidcut->SetTPCnSigmaCut(3,3);//kaon
	//tpcpidcut->SetTPCnSigmaCut(4,3);//proton
	//task1->SetTPCPIDCUT(tpcpidcut);
	
	task1->SetAliESDtrackCuts(cuts1);
	task1->Setsigmacut(1.0);
	task1->SetYcut(0.5);
	task1->SetChargeCut(50.0);
 	task1->SetDCA2010();
	
	
	//task1->SetNsigmaDCAcut(5.0,5.0);
	if(usemc)
  	{	
		task1->SetMCOn();	
		/*TGraph* k0=new TGraph("/home/marek/Analysis/Spectra/feeddown/K0ratio.txt");
		TGraph* lambda=new TGraph("/home/marek/Analysis/Spectra/feeddown/Lambdaratio.txt");
		TGraph* antilambda=new TGraph("/home/marek/Analysis/Spectra/feeddown/AntiLambdaratio.txt");
		task1->SetWeights(k0,lambda,antilambda);*/
		Double_t par[5]={1.2,32.55,0.965,0.953,1.8};
		task1->SetFunctionParam(par);
					
 	}
  	else
  	{
		//Float_t par[5]={1.12,17.72,1.00,0.906,3.47};
		Double_t par[5]={0.91,20.44,1.11,0.948,3.68};
		task1->SetFunctionParam(par);	
  	}
	return task1;
}
