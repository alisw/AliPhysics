void get_dedx(void){

  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/build/include -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG2/FLOW/AliFlowCommon -I$ALICE_ROOT/PWG2/FLOW/AliFlowTasks -I$ALICE_ROOT/PWG3/dielectron/ -g");

  gSystem->Load("libCore");// no
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libXMLIO");// no
  gSystem->Load("libPhysics");
  gSystem->Load("libXMLParser");
  gSystem->Load("libProof");
  gSystem->Load("libMinuit");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFrec");
  gSystem->Load("libT0base");
  gSystem->Load("libT0rec");
  gSystem->Load("libPWG2flowCommon");
  gSystem->Load("libPWG2flowTasks");

  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");



  gSystem->Load("libCORRFW.so");
  gSystem->Load("libPWG3base.so");
  gSystem->Load("libPWG3dielectron.so");
  gSystem->Load("libPWG3hfe.so");


  TChain *chain = new TChain("esdTree");
  for(int i=1; i!=4; ++i)
    chain->Add( Form("/home/gunji/softwares/dielectron/data/esd137549035/%d0/AliESDs.root",i) );

  AliAnalysisManager *mgr = new AliAnalysisManager("DielectronAnalysisManager");
  AliESDInputHandler *esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
  AddTaskPhysicsSelection(kFALSE);
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *taskCentrality =AddTaskCentrality();
  //  taskCentrality->SetPass(2);

  /*  gROOT->LoadMacro("AliDielectronDebugTreeTaku.cxx++");
  gROOT->LoadMacro("AliDielectronHistosTaku.cxx++");
  gROOT->LoadMacro("AliDielectronTaku.cxx++");
  gROOT->LoadMacro("AliAnalysisTaskMultiDielectronNewTaku.cxx++");
  */

  gSystem->Load("./AliDielectronHistosTaku_cxx.so");
  gSystem->Load("./AliDielectronDebugTreeTaku_cxx.so");
  gSystem->Load("./AliDielectronTaku_cxx.so");


  float mom[5][200];
  float dedx[5][200];



  /*
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();

  if ( inputHandler->GetPIDResponse() ){
    AliDielectronVarManager::SetPIDResponse( inputHandler->GetPIDResponse() );
  } else {
    if (isESD){
      if (!AliDielectronVarManager::GetESDpid()){
        if (AliDielectronMC::Instance()->HasMC()) {
	  AliDielectronVarManager::InitESDpid();
        } else {
	  cout<<" set pid as 1"<<endl;
	  AliDielectronVarManager::InitESDpid(1);
        }
      }
    }
  }

  AliPIDResponse * f= (AliPIDResponse*) AliDielectronVarManager::GetESDpid() ;
  AliTPCPIDResponse *ff = (AliTPCPIDResponse*)f->GetTPCResponse();

  for (Int_t j=0; j<AliPID::kSPECIES; j++) {
    //AliPID::EParticleType type=AliPID::EParticleType(j);
    for(int ii=0;ii<100;ii++){
      float p = 0.05+0.1*ii;
      //Double_t bethe=fTPCResponse.GetExpectedSignal(mom,type); 
      Double_t bethe=ff->GetExpectedSignal(p,j); 
      mom[j][ii] = p;
      dedx[j][ii] = bethe;

      //cout<<j<<" "<<p<<" "<<bethe<<endl;
    }
  }
  */

  Double_t fAlephParam[5]={2.11543/57,
                           20.3394,
                           1.0411e-12,
                           2.25543,
                           3.18663
  };

  /*
  Double_t fAlephParam[5]={2.11543/122,
                           42.3394,
                           2.0411e-22,
			   2.25543,
			   6.89
  };
  */

  AliESDpid *fESDpid = new AliESDpid();
  fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParam[0],
                                                    fAlephParam[1],
                                                    fAlephParam[2],
                                                    fAlephParam[3],
                                                    fAlephParam[4]);


  AliTPCPIDResponse *ff = (AliTPCPIDResponse*)fESDpid->GetTPCResponse();

  for (Int_t j=0; j<AliPID::kSPECIES; j++) {                                                                                           
    for(int ii=0;ii<200;ii++){                                                                                                                                                     
      float p = 0.025+0.05*ii;                                                                                                                                                   
      Double_t bethe=ff->GetExpectedSignal(p,j);                                                                                                                                   
      mom[j][ii] = p;                                                                                                                                                               
      dedx[j][ii] = bethe;                                                                                                                                                          
      cout<<j<<" "<<p<<" "<<bethe<<endl;
    }                                                                                                                                                                              
  }                   


  TGraph *g[5];
  char name[100];
  for(int i=0;i<5;i++){
    g[i] = new TGraph(200, mom[i], dedx[i]);
    sprintf(name,"g%d",i);
    g[i]->SetName(name);    
    g[i]->SetLineColor(i+2);
  }
  /*
  TFile *fin=new TFile("ana/tmp.root");
  fin->cd();
  hdedx_pt->SetAxisRange(0,4);
  hdedx_pt->Draw();
  g[0]->Draw("pc");
  g[2]->Draw("pc");
  g[3]->Draw("pc");
  g[4]->Draw("pc");
  */
  TFile *fout = new TFile("dedx.root","recreate");
  fout->cd();
  for(int i=0;i<5;i++){
    g[i]->Write();
  }
  


}




