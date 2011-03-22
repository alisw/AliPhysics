//gSystem->AddIncludePath("-I/$ALICE_ROOT/PWG3/vertexingHF");
  
#include <Riostream.h>
#include <TFile.h>
//#include <AliRDHFCutsD0toKpipipi.h>
#include <TClonesArray.h>
#include <TParameter.h>

//Author: Chiara Bianchin, cbianchi@pd.infn.it, Fabio Colamaria, fabio.colamaria@ba.infn.it

//macro to make a .root file which contains an AliRDHFCutsD0toKpipipi for AliAnalysisTaskSESelectHF4Prong task

void MakeCuts4Charm4Prong(){

gSystem->Load("libTree.so");
gSystem->Load("libGeom.so");
gSystem->Load("libPhysics.so");
gSystem->Load("libVMC.so");
gSystem->Load("libMinuit.so");
gSystem->Load("libSTEERBase.so");
gSystem->Load("libESD.so");
gSystem->Load("libAOD.so"); 
gSystem->Load("libANALYSIS.so");
gSystem->Load("libANALYSISalice.so");
gSystem->Load("libCORRFW.so");
gSystem->Load("libPWG3base.so");
gSystem->Load("libPWG3vertexingHF.so");
gSystem->Load("libPWG3muon.so");

  AliRDHFCutsD0toKpipipi* RDHFCharm4Prong=new AliRDHFCutsD0toKpipipi();
  RDHFCharm4Prong->SetName("Charm4ProngCuts");
  RDHFCharm4Prong->SetTitle("Cuts for D0 analysis");
 
  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4); // default is 5
  //esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);

  RDHFCharm4Prong->AddTrackCuts(esdTrackCuts);


  
  const Int_t nvars=9;

  const Int_t nptbins=5;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.; 
  ptbins[1]=4.5;
  ptbins[2]=6.;
  ptbins[3]=8.;
  ptbins[4]=12.;
  ptbins[5]=16.;
  
  RDHFCharm4Prong->SetPtBins(nptbins+1,ptbins);
  

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }

  //setting my cut values
    //cuts order
    //     printf("    |M-MD0| [GeV]");
    //     printf("    DCA [cm]");
    //     printf("    Dist 2-trk Vtx to PrimVtx [cm]");
    //     printf("    Dist 3-trk Vtx to PrimVtx [cm]");
    //     printf("    Dist 4-trk Vtx to PrimVtx [cm]");
    //     printf("    cosThetaPoint");
    //     printf("    pt [GeV/c]");
    //     printf("    rho mass");
    //     printf("    PID cut");

  //setting my cut values
  rdcutsvalmine[0][0]=0.2;
  rdcutsvalmine[1][0]=0.025;
  rdcutsvalmine[2][0]=0.0900;
  rdcutsvalmine[3][0]=0.0900;
  rdcutsvalmine[4][0]=0.0600;
  rdcutsvalmine[5][0]=0.995;
  rdcutsvalmine[6][0]=0.;
  rdcutsvalmine[7][0]=0.1;
  rdcutsvalmine[8][0]=0.;

  rdcutsvalmine[0][1]=0.2;
  rdcutsvalmine[1][1]=0.015;
  rdcutsvalmine[2][1]=0.1200;
  rdcutsvalmine[3][1]=0.1000;
  rdcutsvalmine[4][1]=0.0850;
  rdcutsvalmine[5][1]=0.99;
  rdcutsvalmine[6][1]=0.;
  rdcutsvalmine[7][1]=0.1;
  rdcutsvalmine[8][1]=0.;

  rdcutsvalmine[0][2]=0.2;
  rdcutsvalmine[1][2]=0.025;
  rdcutsvalmine[2][2]=0.0975;
  rdcutsvalmine[3][2]=0.0900;
  rdcutsvalmine[4][2]=0.0800;
  rdcutsvalmine[5][2]=0.99;
  rdcutsvalmine[6][2]=0.;
  rdcutsvalmine[7][2]=0.1;
  rdcutsvalmine[8][2]=0.;

  rdcutsvalmine[0][3]=0.2;
  rdcutsvalmine[1][3]=0.015;
  rdcutsvalmine[2][3]=0.0975;
  rdcutsvalmine[3][3]=0.0700;
  rdcutsvalmine[4][3]=0.0750;
  rdcutsvalmine[5][3]=0.9825;
  rdcutsvalmine[6][3]=0.;
  rdcutsvalmine[7][3]=0.1;
  rdcutsvalmine[8][3]=0.;

  rdcutsvalmine[0][4]=rdcutsvalmine[0][3];
  rdcutsvalmine[1][4]=rdcutsvalmine[1][3];
  rdcutsvalmine[2][4]=rdcutsvalmine[2][3];
  rdcutsvalmine[3][4]=rdcutsvalmine[3][3];
  rdcutsvalmine[4][4]=rdcutsvalmine[4][3];
  rdcutsvalmine[5][4]=rdcutsvalmine[5][3];
  rdcutsvalmine[6][4]=rdcutsvalmine[6][3];
  rdcutsvalmine[7][4]=rdcutsvalmine[7][3];
  rdcutsvalmine[8][4]=rdcutsvalmine[8][3];

  RDHFCharm4Prong->SetCuts(nvars,nptbins,rdcutsvalmine);
  cout<<"This is the odject I'm going to save:"<<endl;
  RDHFCharm4Prong->PrintAll();
  TFile* fout=new TFile("Charm4ProngCutsDef_16GeV.root","recreate");   //set this!! 
  fout->cd();
  RDHFCharm4Prong->Write();
  fout->Close();

}

void MakeCuts4Charm4ProngForMaxim(){

gSystem->Clear();
gSystem->Load("libTree.so");
gSystem->Load("libGeom.so");
gSystem->Load("libPhysics.so");
gSystem->Load("libVMC.so");
gSystem->Load("libMinuit.so");
gSystem->Load("libSTEERBase.so");
gSystem->Load("libESD.so");
gSystem->Load("libAOD.so"); 
gSystem->Load("libANALYSIS.so");
gSystem->Load("libANALYSISalice.so");
gSystem->Load("libCORRFW.so");
gSystem->Load("libPWG3base.so");
gSystem->Load("libPWG3vertexingHF.so");
gSystem->Load("libPWG3muon.so");

  AliRDHFCutsD0toKpipipi* RDHFCharm4Prong=new AliRDHFCutsD0toKpipipi();
  RDHFCharm4Prong->SetName("loosercuts");
  RDHFCharm4Prong->SetTitle("Cuts for significance maximization");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4);
  //esdTrackCuts->SetMinNClustersTPC(120);

  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10); 
  
  RDHFCharm4Prong->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=9;

  const Int_t nptbins=5; //change this when adding pt bins!
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.; 
  ptbins[1]=4.5;
  ptbins[2]=6.;
  ptbins[3]=8.;
  ptbins[4]=12.;
  ptbins[5]=16.; //o 25!

  RDHFCharm4Prong->SetPtBins(nptbins+1,ptbins);

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }

  //setting my cut values
    //cuts order
    //     printf("    |M-MD0| [GeV]");
    //     printf("    DCA [cm]");
    //     printf("    Dist 2-trk Vtx to PrimVtx [cm]");
    //     printf("    Dist 3-trk Vtx to PrimVtx [cm]");
    //     printf("    Dist 4-trk Vtx to PrimVtx [cm]");
    //     printf("    cosThetaPoint");
    //     printf("    pt [GeV/c]");
    //     printf("    rho mass");
    //     printf("    PID cut");

  Float_t cutsMatrixD0toKpipipiStand[nptbins][nvars]={{0.2,0.04,0.04,0.04,0.04,0.90,0.,0.1,0.},

{0.2,0.04,0.04,0.04,0.04,0.90,0.,0.1,0.},

{0.2,0.04,0.04,0.04,0.04,0.90,0.,0.1,0.}, 

{0.2,0.04,0.04,0.04,0.04,0.90,0.,0.1,0.},

{0.2,0.04,0.04,0.04,0.04,0.90,0.,0.1,0.}};

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpipipiStand[ibin][ivar];
    }
  }
  RDHFCharm4Prong->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  Int_t nvarsforopt=RDHFCharm4Prong->GetNVarsForOpt();
  Int_t dim=5; //set this!!
  Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(dim>nvarsforopt){
    cout<<"Number of variables for optimization has probably changed, check and edit accordingly"<<endl;
    return;
  } else {
    if(dim==nvarsforopt){
      boolforopt=RDHFCharm4Prong->GetVarsForOpt();
    }else{
      TString *names;
      names=new TString[nvars];
      TString answer="";
      Int_t checktrue=0;
      names=RDHFCharm4Prong->GetVarNames();
      for(Int_t i=0;i<nvars;i++){
	cout<<names[i]<<" for opt? (y/n)"<<endl;
	cin>>answer;
	if(answer=="y") {
	  boolforopt[i]=kTRUE;
	  checktrue++;
	}
	else boolforopt[i]=kFALSE;
      }
      if (checktrue!=dim) {
	cout<<"Error! You set "<<checktrue<<" kTRUE instead of "<<dim<<endl;
	return;
      }
      RDHFCharm4Prong->SetVarsForOpt(dim,boolforopt);
    }
  }

// Tight values for variables to be maximized
// (DCA, Dist12, Dist3, Dist4, CosThetaPoint)
// Indexes: variable, ptbin

  Float_t tighterval[5][5];
 
  //number of steps for each variable is set in the AddTask arguments (default=8)

  tighterval[0][0]=0.025;
  tighterval[1][0]=0.09;
  tighterval[2][0]=0.09;
  tighterval[3][0]=0.06;
  tighterval[4][0]=0.955;

  tighterval[0][1]=0.015;
  tighterval[1][1]=0.12;
  tighterval[2][1]=0.10;
  tighterval[3][1]=0.085;
  tighterval[4][1]=0.99;

  tighterval[0][2]=0.025;
  tighterval[1][2]=0.0975;
  tighterval[2][2]=0.09;
  tighterval[3][2]=0.08;
  tighterval[4][2]=0.99;

  tighterval[0][3]=0.015;
  tighterval[1][3]=0.0975;
  tighterval[2][3]=0.07;
  tighterval[3][3]=0.075;
  tighterval[4][3]=0.9825;

  tighterval[0][4]=0.015;
  tighterval[1][4]=0.0975;
  tighterval[2][4]=0.07;
  tighterval[3][4]=0.075;
  tighterval[4][4]=0.9825;

  TString name=""; 
  Int_t arrdim=dim*nptbins;
  cout<<"Will save "<<arrdim<<" TParameter<float>"<<endl;
  TClonesArray max("TParameter<float>",arrdim);
  for(Int_t ival=0;ival<dim;ival++){
    for(Int_t jpt=0;jpt<nptbins;jpt++){
      name=Form("par%dptbin%d",ival,jpt);
      cout<<"Setting "<<name.Data()<<" to "<<tighterval[ival][jpt]<<endl;
      new(max[jpt*dim+ival])TParameter<float>(name.Data(),tighterval[ival][jpt]);
    }
  }

  Bool_t flagPID=kFALSE; 
  RDHFCharm4Prong->SetUsePID(flagPID);

  RDHFCharm4Prong->PrintAll();
  printf("Use PID? %s\n",flagPID ? "yes" : "no");

/*
  //pid settings
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4D0");
  Int_t mode=1;
  const Int_t nlims=2;
  Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
  Bool_t compat=kTRUE; //effective only for this mode
  Bool_t asym=kTRUE;
  Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
  pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
  pidObj->SetMatch(mode);
  pidObj->SetPLimit(plims,nlims);
  pidObj->SetSigma(sigmas);
  pidObj->SetCompat(compat);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  RDHFCharm4Prong->SetPidHF(pidObj);

  RDHFCharm4Prong->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF
*/
  
  //activate pileup rejection (for pp)
  RDHFCharm4Prong->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  //Do not recalculate the vertex
  RDHFCharm4Prong->SetRemoveDaughtersFromPrim(kTRUE); //activate for pp

/* Per il pp levo la centralità!
//Metto kCentOff! Per il Pb c'era kCentV0M
  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=20,maxc=80;
  RDHFCharm4Prong->SetMinCentrality(minc);
  RDHFCharm4Prong->SetMaxCentrality(maxc);
  cent=Form("%.0f%.0f",minc,maxc);
  RDHFCharm4Prong->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
*/

  //temporary
//  RDHFCharm4Prong->SetFixRefs();

  TFile* fout=new TFile("Charm4ProngCutsForMaxim.root","recreate");
  fout->cd();
  RDHFCharm4Prong->Write();
  max.Write();
  fout->Close();
}
