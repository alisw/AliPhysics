AddTaskPidVZEROSP(Int_t centralityselection=AliVEvent::kAny,Float_t etamin=-0.8,Float_t etamax=0.8){
  gROOT->LoadMacro("$ALICE_ROOT/PWGCF/FLOW/macros/AddTaskFlowCentralityPIDSP.C");

  const Int_t ncentr = 9;
  Int_t cmin[ncentr]={0,5,10,20,30,40,50,60,70};
  Int_t cmax[ncentr]={5,10,20,30,40,50,60,70,80};

  for(Int_t i=0;i < ncentr;i++){
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kFALSE,AliPID::kPion,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax); // no pid
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kPion,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kKaon,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTOFbayesian,-1,2,0,etamin,etamax);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTOFbayesian,0,2,0,etamin,etamax);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kPion,AliFlowTrackCuts::kTPCbayesian,0,2,0,etamin,etamax);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kKaon,AliFlowTrackCuts::kTPCbayesian,0,2,0,etamin,etamax);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTPCbayesian,-1,2,0,etamin,etamax);
    AddTaskFlowCentralityPIDSP(centralityselection,cmin[i],cmax[i],"AnalysisResults",kTRUE,AliPID::kProton,AliFlowTrackCuts::kTPCbayesian,0,2,0,etamin,etamax);
  }
}

createSPres(){
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libNetx.so");
  gSystem->Load("libPWGflowBase.so");
  
  char name[200];
  char *spe[4]={"pion","kaon","antipr","proton"};
  char *tech[2]={"TOF","TPC"};
  Int_t cmin[9]={0,5,10,20,30,40,50,60,70};
  Int_t cmax[9]={5,10,20,30,40,50,60,70,80};
  
  TFile *f = new TFile("AnalysisResults.root");
  TFile *fo = new TFile("results.root","RECREATE");
  TDirectory* directory = dynamic_cast<TDirectory*>(f->Get("outputSPanalysisTPCstandalone"));
  TList* listTemp = directory->GetListOfKeys();
  for(Int_t i=0;i < 9;i++){
    TList* list2 = dynamic_cast<TList*>(directory->Get(listTemp->At(i*9)->GetName()));
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
	TList* list = dynamic_cast<TList*>(directory->Get(listTemp->At(i*9+j*4+k+1)->GetName()));
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




