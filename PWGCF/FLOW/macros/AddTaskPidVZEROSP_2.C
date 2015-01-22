AddTaskPidVZEROSP_2(Int_t centralityselection=AliVEvent::kAny,Float_t etamin=-0.8,Float_t etamax=0.8,Int_t side=0,Int_t filterbit=1,Bool_t TOFbeta=kFALSE){
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGCF/FLOW/macros/AddTaskFlowCentralityPIDSP.C");

  const Int_t ncentr = 4;
  Int_t cmin[ncentr]={30,50,60,70};
  Int_t cmax[ncentr]={40,60,70,80};

  for(Int_t i=0;i < ncentr;i++){
    if(!TOFbeta){
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kFALSE,AliPID::kPion,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax,"",side,filterbit); // no pid
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kPion,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax,"",side,filterbit);
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kKaon,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax,"",side,filterbit);
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTOFbayesian,-1,2,0,etamin,etamax,"",side,filterbit);
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax,"",side,filterbit);
    }
    else{
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kFALSE,AliPID::kPion,AliFlowTrackCuts::kTOFbeta,0,2,0,etamin,etamax,"",side,filterbit); // no pid
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kPion,AliFlowTrackCuts::kTOFbeta,0,2,0,etamin,etamax,"",side,filterbit);
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kKaon,AliFlowTrackCuts::kTOFbeta,0,2,0,etamin,etamax,"",side,filterbit);
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTOFbeta,-1,2,0,etamin,etamax,"",side,filterbit);
      AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTOFbeta,0,2,0,etamin,etamax,"",side,filterbit);
    }
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kPion,AliFlowTrackCuts::kTPCbayesian,0,2,0,etamin,etamax,"",side,filterbit);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kKaon,AliFlowTrackCuts::kTPCbayesian,0,2,0,etamin,etamax,"",side,filterbit);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTPCbayesian,-1,2,0,etamin,etamax,"",side,filterbit);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTPCbayesian,0,2,0,etamin,etamax,"",side,filterbit);
  }
}

createSPres(){
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libNetx");
  gSystem->Load("libPWGflowBase");
  
  char name[200];
  char *spe[4]={"pion","kaon","antipr","proton"};
  char *tech[2]={"TOF","TPC"};
  const Int_t ncentr = 4;
  Int_t cmin[ncentr]={30,50,60,70};
  Int_t cmax[ncentr]={40,60,70,80};
  
  TFile *f = new TFile("AnalysisResults.root");
  TFile *fo = new TFile("results.root","RECREATE");
  TDirectory* directory = dynamic_cast<TDirectory*>(f->Get("outputSPanalysisTPCstandalone"));
  TList* listTemp = directory->GetListOfKeys();
  for(Int_t i=0;i < ncentr;i++){
    TList* list2 = dynamic_cast<TList*>(directory->Get(listTemp->At(i*ncentr)->GetName()));
    AliFlowAnalysisWithScalarProduct* sp2 = new AliFlowAnalysisWithScalarProduct();
    sp2->GetOutputHistograms(list2);
    sp2->Finish();
    AliFlowCommonHistResults* res2=sp2->GetCommonHistsRes();
    TH1D *h2=res2->GetHistDiffFlowPtPOI();
    sprintf(name,"v2SP_%s_%i_%i","AllCharged",cmin[i],cmax[i]);
    h2->SetName(name);
    fo->cd();
    h2->Write();
    for(Int_t j=0;j < 2;j++){
      for(Int_t k=0;k < 4;k++){
	TList* list = dynamic_cast<TList*>(directory->Get(listTemp->At(i*ncentr+j*4+k+1)->GetName()));
	AliFlowAnalysisWithScalarProduct* sp = new AliFlowAnalysisWithScalarProduct();
	sp->GetOutputHistograms(list);
	sp->Finish();
	AliFlowCommonHistResults* res=sp->GetCommonHistsRes();
	TH1D *h=res->GetHistDiffFlowPtPOI();
	sprintf(name,"v2SP_%s_%i_%i%s",spe[k],cmin[i],cmax[i],tech[j]);
	h->SetName(name);
	fo->cd();
	h->Write();
      }
    }
  }
  fo->Close();
}




