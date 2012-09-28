#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsLctoV0.h>
#include <AliAODPidHF.h>
#include <TClonesArray.h>
#include <TParameter.h>
#include <TF1.h>


//Use:
//Set hard coded commentet with //set this!!
// root[] .L makeInput...C++
// root[] makeInputAliAnalysisTaskSE...()
//similar macros for D mesons as well as for Lc->3prongs

//Author: Annalisa De Caro - decaro@sa.infn.it


//macro to make a .root file which contains an AliRDHFCutsLctoV0 for AliAnalysisTaskSELc2pK0S task

void makeInputAliAnalysisTaskSELctoV0bachelor(){

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(0);//(4); // default is 5
  esdTrackCuts->SetMinNClustersTPC(70);
  //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  //	 				   AliESDtrackCuts::kAny); 
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  esdTrackCuts->SetAcceptKinkDaughters(kFALSE);



  AliESDtrackCuts* esdTrackCutsV0daughters=new AliESDtrackCuts();
  esdTrackCutsV0daughters->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCutsV0daughters->SetRequireTPCRefit(kTRUE);
  esdTrackCutsV0daughters->SetRequireITSRefit(kFALSE);//kTRUE);
  esdTrackCutsV0daughters->SetMinNClustersITS(0);//(4); // default is 5
  esdTrackCutsV0daughters->SetMinNClustersTPC(70);
  //esdTrackCutsV0daughters->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  //	      				              AliESDtrackCuts::kAny); 
  // default is kBoth, otherwise kAny
  esdTrackCutsV0daughters->SetMinDCAToVertexXY(0.);
  esdTrackCutsV0daughters->SetPtRange(0.,1.e10);
  esdTrackCutsV0daughters->SetAcceptKinkDaughters(kFALSE);



  AliRDHFCutsLctoV0* RDHFLctoV0An=new AliRDHFCutsLctoV0();
  RDHFLctoV0An->SetName("LctoV0AnalysisCuts");
  RDHFLctoV0An->SetTitle("Analysis cuts for Lc analysis");

  AliRDHFCutsLctoV0* RDHFLctoV0Prod=new AliRDHFCutsLctoV0();
  RDHFLctoV0Prod->SetName("LctoV0ProductionCuts");
  RDHFLctoV0Prod->SetTitle("Production cuts for Lc analysis");

  RDHFLctoV0Prod->AddTrackCuts(esdTrackCuts);
  RDHFLctoV0An->AddTrackCuts(esdTrackCuts);

  RDHFLctoV0Prod->AddTrackCutsV0daughters(esdTrackCutsV0daughters);
  RDHFLctoV0An->AddTrackCutsV0daughters(esdTrackCutsV0daughters);

  RDHFLctoV0Prod->SetPidSelectionFlag(1); // 0 -> TOF AND TPC
                                          // 1 -> if (TOF) TOF else TPC w veto
  RDHFLctoV0An->SetPidSelectionFlag(1); // 0 -> TOF AND TPC
                                        // 1 -> if (TOF) TOF else TPC w veto

  const Int_t nptbins=1;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=99999999.;
  RDHFLctoV0Prod->SetPtBins(nptbins+1,ptbins);
  RDHFLctoV0An->SetPtBins(nptbins+1,ptbins);

  const Int_t nvars=9 ;

  Float_t** prodcutsval;
  prodcutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){prodcutsval[ic]=new Float_t[nptbins];}
  for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
   prodcutsval[0][ipt2]=1.;    // inv. mass if K0S [GeV/c2]
   prodcutsval[1][ipt2]=1.;    // inv. mass if Lambda [GeV/c2]
   prodcutsval[2][ipt2]=0.05;  // inv. mass V0 if K0S [GeV/c2]
   prodcutsval[3][ipt2]=0.05;  // inv. mass V0 if Lambda [GeV/c2]
   prodcutsval[4][ipt2]=0.3;   // pT min bachelor track [GeV/c] // AOD by construction
   prodcutsval[5][ipt2]=0.;    // pT min V0-positive track [GeV/c]
   prodcutsval[6][ipt2]=0.;    // pT min V0-negative track [GeV/c]
   prodcutsval[7][ipt2]=1000.; // dca cascade cut [cm]
   prodcutsval[8][ipt2]=1000.; // dca V0 cut [nSigma] // it's 1.5 x offline V0s
  }
  RDHFLctoV0Prod->SetCuts(nvars,nptbins,prodcutsval);

  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
   anacutsval[0][ipt2]=0.25;   // inv. mass if K0S [GeV/c2]
   anacutsval[1][ipt2]=0.25;   // inv. mass if Lambda [GeV/c2]
   anacutsval[2][ipt2]=0.0075; // inv. mass V0 if K0S [GeV/c2]
   anacutsval[3][ipt2]=0.0030; // inv. mass V0 if Lambda [GeV/c2]
   anacutsval[4][ipt2]=0.3;    // pT min bachelor track [GeV/c] // AOD by construction
   anacutsval[5][ipt2]=0.;     // pT min V0-positive track [GeV/c]
   anacutsval[6][ipt2]=0.;     // pT min V0-negative track [GeV/c]
   anacutsval[7][ipt2]=1000.; // dca cascade cut [cm]
   anacutsval[8][ipt2]=1.5;   // dca V0 cut [nSigma] // it's 1.5 x offline V0s
  }
  RDHFLctoV0An->SetCuts(nvars,nptbins,anacutsval);


  //RDHFLc->SetRecoKF(); //set this if you want to recompute the secondary vertex with the KF package

  //pid settings
  //1. bachelor: default one
  AliAODPidHF* pidObjBachelor = new AliAODPidHF();
  Double_t sigmasBac[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjBachelor->SetSigma(sigmasBac);
  pidObjBachelor->SetAsym(kFALSE);
  pidObjBachelor->SetMatch(1);
  pidObjBachelor->SetTPC(kTRUE);
  pidObjBachelor->SetTOF(kTRUE);
  pidObjBachelor->SetTOFdecide(kFALSE);

  RDHFLctoV0An->SetPidHF(pidObjBachelor);
  RDHFLctoV0Prod->SetPidHF(pidObjBachelor);

  //2. V0pos
  AliAODPidHF* pidObjV0pos = new AliAODPidHF();
  Double_t sigmasV0pos[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjV0pos->SetSigma(sigmasV0pos);
  pidObjV0pos->SetAsym(kFALSE);
  pidObjV0pos->SetMatch(1);
  pidObjV0pos->SetTPC(kTRUE);
  pidObjV0pos->SetTOF(kTRUE);
  pidObjV0pos->SetTOFdecide(kFALSE);

  RDHFLctoV0An->SetPidV0pos(pidObjV0pos);
  RDHFLctoV0Prod->SetPidV0pos(pidObjV0pos);

  //2. V0neg
  AliAODPidHF* pidObjV0neg = new AliAODPidHF();
  Double_t sigmasV0neg[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjV0neg->SetSigma(sigmasV0neg);
  pidObjV0neg->SetAsym(kFALSE);
  pidObjV0neg->SetMatch(1);
  pidObjV0neg->SetTPC(kTRUE);
  pidObjV0neg->SetTOF(kTRUE);
  pidObjV0neg->SetTOFdecide(kFALSE);

  RDHFLctoV0An->SetPidV0neg(pidObjV0neg);
  RDHFLctoV0Prod->SetPidV0neg(pidObjV0neg);


  // uncomment these lines for Baysian PID:
  // Double_t threshold=0.3;
  // SetupCombinedPID(RDHFLctoV0An  ,threshold);
  // RDHFLctoV0An  ->SetPIDStrategy(AliRDHFCutsLctoV0::kCombined);
  //


  //uncomment these lines to apply cuts with the KF package
  //RDHFLctoV0An->SetCutsStrategy(AliRDHFCutsLctoV0::kKF);
  //for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
  //   anacutsval[0][ipt2]=1.;  //if <0., no topological constraint
  //   anacutsval[1][ipt2]=2.;  //cut on the Chi2/Ndf
  // }

  Bool_t pidflag=kTRUE;
  RDHFLctoV0An->SetUsePID(pidflag);
  RDHFLctoV0Prod->SetUsePID(kFALSE);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

  cout<<"This is the object I'm going to save:"<<endl;
  RDHFLctoV0Prod->PrintAll();
  cout<<"This is the object I'm going to save:"<<endl;
  RDHFLctoV0An->PrintAll();
  TFile* fout=new TFile("Lc2pK0SCuts.root","RECREATE"); 
  fout->cd();
  RDHFLctoV0Prod->Write();
  RDHFLctoV0An->Write();
  fout->Close();
  delete fout;

  delete pidObjBachelor;
  delete pidObjV0neg;
  delete pidObjV0pos;
  delete RDHFLctoV0Prod;
  delete RDHFLctoV0An;

}


//macro to make a .root file (for significance maximization) which contains an AliRDHFCutsLctoV0 with loose set of cuts  and TParameter with the tighest value of these cuts

void makeInputAliAnalysisTaskSESignificanceMaximization(){

  AliRDHFCutsLctoV0* RDHFLctoV0=new AliRDHFCutsLctoV0();
  RDHFLctoV0->SetName("loosercuts");
  RDHFLctoV0->SetTitle("Cuts for significance maximization");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4);

  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  
  AliESDtrackCuts* esdTrackCutsV0daughters=new AliESDtrackCuts();
  esdTrackCutsV0daughters->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCutsV0daughters->SetRequireTPCRefit(kTRUE);
  esdTrackCutsV0daughters->SetRequireITSRefit(kTRUE);
  esdTrackCutsV0daughters->SetMinNClustersITS(4); // default is 5
  esdTrackCutsV0daughters->SetMinNClustersTPC(70);
  esdTrackCutsV0daughters->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
  // default is kBoth, otherwise kAny
  esdTrackCutsV0daughters->SetMinDCAToVertexXY(0.);
  esdTrackCutsV0daughters->SetPtRange(0.,1.e10);

  RDHFLctoV0->AddTrackCuts(esdTrackCuts);

  RDHFLctoV0->AddTrackCutsV0daughters(esdTrackCutsV0daughters);

  RDHFLctoV0->SetPidSelectionFlag(1); // 0 -> TOF AND TPC
                                      // 1 -> if (TOF) TOF else TPC w veto

  const Int_t nvars=9;

  const Int_t nptbins=1; //change this when adding pt bins!
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=999999999.;

  RDHFLctoV0->SetPtBins(nptbins+1,ptbins);

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }

  //setting my cut values
  // inv. mass if K0s [GeV]
  // inv. mass if Lambda [GeV]
  // inv. mass V0 if K0S [GeV]
  // inv. mass V0 if Lambda [GeV]
  // pT min bachelor track [GeV/c]
  // pT min V0-positive track [GeV/c]
  // pT min V0-negative track [GeV/c]
  // dca cascade cut (cm)
  // dca V0 cut (cm)
  Float_t cutsMatrixLctoV0Stand[nptbins][nvars]=
    { {0.25,0.25,0.0075,0.0075,0.,0.,0.,1.5,1.5} };

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar=0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixLctoV0Stand[ibin][ivar];
    }
  }
  RDHFLctoV0->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);


  Int_t nvarsforopt=RDHFLctoV0->GetNVarsForOpt();
  Int_t dim=7; //set this!!
  Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(dim>nvarsforopt){
    cout<<"Number of variables for optimization has probably changed, check and edit accordingly"<<endl;
    return;
  } else {
    if(dim==nvarsforopt){
      boolforopt=RDHFLctoV0->GetVarsForOpt();
    }else{
      TString *names;
      names=new TString[nvars];
      TString answer="";
      Int_t checktrue=0;
      names=RDHFLctoV0->GetVarNames();
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
      RDHFLctoV0->SetVarsForOpt(dim,boolforopt);
    }
  }


  Float_t tighterval[dim][nptbins];
  // 0(2): inv. mass V0 if K0S [GeV]
  // 1(3): inv. mass V0 if Lambda [GeV]
  // 2(4): pT min bachelor track [GeV/c]
  // 3(5): pT min V0-positive track [GeV/c]
  // 4(6): pT min V0-negative track [GeV/c]
  // 5(7): dca cascade cut (cm)
  // 6(8): dca V0 cut (cm)

  // number of steps for each variable is set in the AddTask arguments (default=8)
  // set this!!
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    tighterval[0][ipt]=0.075; // inv. mass V0 if K0S [GeV]
    tighterval[1][ipt]=0.040; // inv. mass V0 if Lambda [GeV]
    tighterval[2][ipt]=0.1; // pT min bachelor track [GeV/c]
    tighterval[3][ipt]=0.1; // pT min V0-positive track [GeV/c]
    tighterval[4][ipt]=0.1; // pT min V0-negative track [GeV/c]
    tighterval[5][ipt]=10.; // dca cascade cut (cm)
    tighterval[6][ipt]=10.; // dca v0 cut (cm)
  }

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
  RDHFLctoV0->SetUsePID(flagPID);

  RDHFLctoV0->PrintAll();
  printf("Use PID? %s\n",flagPID ? "yes" : "no");

  //pid settings
  //1. bachelor: default one
  AliAODPidHF* pidObjBachelor = new AliAODPidHF();
  Double_t sigmasBac[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjBachelor->SetSigma(sigmasBac);
  pidObjBachelor->SetAsym(kFALSE);
  pidObjBachelor->SetMatch(1);
  pidObjBachelor->SetTPC(kTRUE);
  pidObjBachelor->SetTOF(kTRUE);
  pidObjBachelor->SetTOFdecide(kFALSE);
  RDHFLctoV0->SetPidHF(pidObjBachelor);

  //2. V0pos
  AliAODPidHF* pidObjV0pos = new AliAODPidHF();
  Double_t sigmasV0pos[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjV0pos->SetSigma(sigmasV0pos);
  pidObjV0pos->SetAsym(kFALSE);
  pidObjV0pos->SetMatch(1);
  pidObjV0pos->SetTPC(kTRUE);
  pidObjV0pos->SetTOF(kTRUE);
  pidObjV0pos->SetTOFdecide(kFALSE);
  RDHFLctoV0->SetPidV0pos(pidObjV0pos);

  //2. V0neg
  AliAODPidHF* pidObjV0neg = new AliAODPidHF();
  Double_t sigmasV0neg[5]={3.,1.,1.,3.,3.}; // 0, 1(A), 2(A) -> TPC; 3 -> TOF; 4 -> ITS
  pidObjV0neg->SetSigma(sigmasV0neg);
  pidObjV0neg->SetAsym(kFALSE);
  pidObjV0neg->SetMatch(1);
  pidObjV0neg->SetTPC(kTRUE);
  pidObjV0neg->SetTOF(kTRUE);
  pidObjV0neg->SetTOFdecide(kFALSE);
  RDHFLctoV0->SetPidV0neg(pidObjV0neg);

  //activate pileup rejection (for pp)
  //RDHFLctoV0->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  //Do not recalculate the vertex
  RDHFLctoV0->SetRemoveDaughtersFromPrim(kFALSE); //activate for pp

  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=20,maxc=80;
  RDHFLctoV0->SetMinCentrality(minc);
  RDHFLctoV0->SetMaxCentrality(maxc);
  cent=Form("%.0f%.0f",minc,maxc);
  //RDHFLctoV0->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid
  RDHFLctoV0->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

  //temporary
  //RDHFLctoV0->SetFixRefs();

  TFile* fout=new TFile(Form("cuts4SignifMaxim%s%s%sRecVtx%sPileupRej.root",
			     RDHFLctoV0->GetUseCentrality()==0 ? "pp" : "PbPb",
			     cent.Data(),
			     RDHFLctoV0->GetIsPrimaryWithoutDaughters() ? "" : "No",
			     RDHFLctoV0->GetOptPileUp() ? "" : "No"),"recreate"); //set this!! 

  fout->cd();
  RDHFLctoV0->Write();
  max.Write();
  fout->Close();
 
}
