AliAnalysisTaskHFE* ConfigHFECalSys_PbPb(Bool_t useMC, int TPCclust, int Nits, int ITSstat, double nSigMim, double Mimeop, double Maxeop, int PIDorder){
  //
  // HFE standard task configuration
  //

  cout << "==== Summary of cuts  ===== " << endl;
  cout << "TPC ; " << TPCclust << endl;
  cout << "ITS ; " << Nits << endl;
  cout << "PID ; " << nSigMim << " ; " << Mimeop << " ; " << Maxeop << endl;

  AliHFEcuts *hfecuts = new AliHFEcuts("hfeCuts","HFE Standard Cuts for EMCal");
  hfecuts->CreateStandardCuts();
  //hfecuts->SetMinNClustersTPC(120);
  hfecuts->SetMinNClustersTPC(TPCclust);
  hfecuts->SetMinRatioTPCclusters(0.6);
  hfecuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);

  //hfecuts->SetMinNClustersITS(3);
  hfecuts->SetMinNClustersITS(Nits);
  hfecuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  if(ITSstat==1)hfecuts->SetCutITSpixel(AliHFEextraCuts::kFirst);
  hfecuts->SetCheckITSLayerStatus(kFALSE);

  //hfecuts->UnsetVertexRequirement();
  hfecuts->SetMaxImpactParam(3.,3.);
  hfecuts->SetPtRange(2.0,60.0);

  hfecuts->SetVertexRange(10.);
  //hfecuts->SetMaxChi2perClusterITS(36);
  //hfecuts->SetSigmaToVertex(10);
  //hfecuts->SetTOFPIDStep(kTRUE);
  //hfecuts->SetTOFMISMATCHStep(kTRUE);
  //hfecuts->SetTPCPIDCleanUpStep(kTRUE);
  hfecuts->SetQAOn();
  //hfecuts->SetMinNTrackletsTRD(5);

  AliAnalysisTaskHFE *task = new AliAnalysisTaskHFE("HFEanalysisStandardEMCal");
  task->SetHFECuts(hfecuts);
  task->SetPbPbAnalysis(kTRUE);
  //task->SetRemovePileUp(kTRUE);
  task->GetPIDQAManager()->SetHighResolutionHistos();

  // Define Variables
  AliHFEvarManager *vm = task->GetVarManager();
  //vm->AddVariable("pt");
  //vm->AddVariable("eta");
  
  const Double_t ptbinning[50] = {1., 1.1, 1.2, 1.3, 1.4, 
                                  1.5, 1.75, 2., 2.25, 2.5, 
                                  2.75, 3., 3.5, 4., 4.5, 
                                  5., 5.5, 6., 7., 8., 
                                  9., 10., 11., 12., 13., 
                                  14., 15., 16., 17., 18.,
                                  20., 22., 24., 26., 28.,
                                  30., 32., 34., 36., 38.,
				  40., 45., 50., 55., 60.,
			          65., 70., 80., 90., 100}; 

  const Double_t etabinning[17] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8}; 

  vm->AddVariable("pt", 49, ptbinning);
  vm->AddVariable("eta", 16, etabinning);
  vm->AddVariable("phi");
  vm->AddVariable("charge");
  vm->AddVariable("source");
  vm->AddVariable("centrality");
  /*
  if(!useMC){

      for(Int_t a=0;a<12;a++)
      {
	  TF1 *hBackground = new TF1("hadronicBackgroundFunction","TMath::Exp([0]/x + [1])", 0., 20.);
	  hBackground->SetParameter(0, -43.87);
	  hBackground->SetParameter(1, 2.85);
	  task->SetBackGroundFactorsFunction(hBackground,a);
      }


  }
  */
  // Define PID
  AliHFEpid *pid = task->GetPID();
  if(useMC) pid->SetHasMCData(kTRUE);

  if(PIDorder == 0)
     {
      pid->AddDetector("EMCAL", 1);
      pid->AddDetector("TPC", 0);
     }
  else if
    {
     pid->AddDetector("EMCAL", 0);
     pid->AddDetector("TPC", 1);
   }

//  pid->ConfigureTPCrejection();

  if(!useMC){

      Double_t params_centr_0_5[1];
      Double_t params_centr_5_10[1];
      Double_t params_centr_10_20[1];
      Double_t params_centr_20_30[1];
      Double_t params_centr_per[1];
      //params_centr_0_5[0]=0.16;  // cut tuned for 0-10%
      //params_centr_5_10[0]=0.16; // cut tuned for 0-10%
      //params_centr_10_20[0]=0.29;
      //params_centr_20_30[0]=0.38;
      //params_centr_per[0]=0.44;
      /*
      params_centr_0_5[0]=-1.5;  // cut tuned for 0-10%
      params_centr_5_10[0]=-1.5; // cut tuned for 0-10%
      params_centr_10_20[0]=-1.5;
      params_centr_20_30[0]=-1.5;
      params_centr_per[0]=-1.5;
      */
      params_centr_0_5[0] = nSigMim;  // cut tuned for 0-10%
      params_centr_5_10[0] = nSigMim; // cut tuned for 0-10%
      params_centr_10_20[0] = nSigMim;
      params_centr_20_30[0] = nSigMim;
      params_centr_per[0] = nSigMim;


      char *cutmodel;
      cutmodel="pol0";


      for(Int_t a=0;a<11;a++)
      {
	  if(a>3)  pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_per,3.0);
	  if(a==0) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_0_5,3.0);    //  0-5%
	  if(a==1) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_5_10,3.0);    //  5-10%
	  if(a==2) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_10_20,3.0);    //  10-20%
	  if(a==3) pid->ConfigureTPCcentralityCut(a,cutmodel,params_centr_20_30,3.0);    //  20-30%
      }


  }


  // V0 tagged tracks
  AliHFEcuts *v0trackCuts = new AliHFEcuts("V0trackCuts", "Track Cuts for tagged track Analysis");
  v0trackCuts->CreateStandardCuts();
  //v0trackCuts->SetMinNClustersTPC(120);
  v0trackCuts->SetMinNClustersTPC(TPCclust);
  v0trackCuts->SetMinRatioTPCclusters(0.6);
  v0trackCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
  v0trackCuts->SetMinNClustersITS(1);
  v0trackCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
  v0trackCuts->SetCheckITSLayerStatus(kFALSE);
  v0trackCuts->UnsetVertexRequirement();
  //v0trackCuts->SetMaxChi2perClusterITS(36);
  //hfecuts->SetSigmaToVertex(10);
  v0trackCuts->SetTOFPIDStep(kTRUE);
//  v0trackCuts->SetTOFMISMATCHStep(kTRUE);
  //v0trackCuts->SetTPCPIDCleanUpStep(kTRUE);
  v0trackCuts->SetQAOn();
  
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kTaggedTrackAnalysis);
  task->SetTaggedTrackCuts(v0trackCuts);
  task->SetCleanTaggedTrack(kFALSE);
  
  // change E/p cuts
  AliHFEpidEMCAL *emcpid = pid->AliHFEpid::GetDetPID(AliHFEpid::kEMCALpid);
  emcpid->SetEoPMax(Maxeop);
  emcpid->SetEoPMim(Mimeop);

  // QA
  task->SetQAOn(AliAnalysisTaskHFE::kPIDqa);
//  task->SetFillSignalOnly(kFALSE);    // for DE pluging for MC
  task->SetQAOn(AliAnalysisTaskHFE::kMCqa);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kIsElecBackGround);
  //task->SwitchOnPlugin(AliAnalysisTaskHFE::kSecVtx);
  task->SwitchOnPlugin(AliAnalysisTaskHFE::kDEstep);



  printf("*************************************\n");
  printf("Configuring standard Task:\n");
  task->Print();
  pid->PrintStatus();
  printf("*************************************\n"); 
  return task;
}
