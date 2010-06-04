#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsDplustoKpipi.h>
#include <TClonesArray.h>
#include <TParameter.h>

//macro to make a .root file which contains an AliRDHFCutsDplustoKpipi with loose set of cuts (for significance maximization) and TParameter with the tighest value of these cuts
//Needed for AliAnalysisTaskSESignificance

//Use:
//Set hard coded commented with //set this!!
//.x makeTFile4CutsDplustoKpipi.C++

//similar macros for the other D mesons

//Author: Chiara Bianchin, cbianchi@pd.infn.it
//        Giacomo Ortona, ortona@to.infn.it

void makeTFile4CutsDplustoKpipi(){

  AliRDHFCutsDplustoKpipi* RDHFDplustoKpipi=new AliRDHFCutsDplustoKpipi();
  RDHFDplustoKpipi->SetName("loosercuts");
  RDHFDplustoKpipi->SetTitle("Cuts for significance maximization");

//   AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
//   esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
//   //default
//   esdTrackCuts->SetRequireTPCRefit(kTRUE);
//   esdTrackCuts->SetRequireITSRefit(kTRUE);
//   esdTrackCuts->SetMinNClustersITS(4); // default is 5
//   esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
// 					 AliESDtrackCuts::kAny); 
//  // default is kBoth, otherwise kAny
//   esdTrackCuts->SetMinDCAToVertexXY(0.);
//   esdTrackCuts->SetPtRange(0.3,1.e10);


//   RDHFDplustoKpipi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=12;

  const Int_t nptbins=4;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=2.;
  ptbins[2]=3.;
  ptbins[3]=5.;
  ptbins[4]=99999.;
  
  RDHFDplustoKpipi->SetPtBins(nptbins+1,ptbins);
    
  //setting my cut values
    
  Float_t** prodcutsval;
  prodcutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){prodcutsval[ic]=new Float_t[nptbins];}  
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    prodcutsval[0][ipt]=0.2;
    prodcutsval[1][ipt]=0.4;
    prodcutsval[2][ipt]=0.4;
    prodcutsval[3][ipt]=0.;
    prodcutsval[4][ipt]=0.;
    prodcutsval[5][ipt]=0.01;
    prodcutsval[6][ipt]=0.06;
    prodcutsval[7][ipt]=0.02;
    prodcutsval[8][ipt]=0.;
    prodcutsval[9][ipt]=0.85;
    prodcutsval[10][ipt]=0.;
    prodcutsval[11][ipt]=10000000.0;
    
  }

  RDHFDplustoKpipi->SetCuts(nvars,nptbins,prodcutsval);

  Int_t nvarsforopt=RDHFDplustoKpipi->GetNVarsForOpt();
  const Int_t dim=5; //set this!!
  Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(dim>nvarsforopt){
    cout<<"Number of variables for optimization has probably changed, check and edit accordingly"<<endl;
    return;
  } else {
    if(dim==nvarsforopt){
      boolforopt=RDHFDplustoKpipi->GetVarsForOpt();
    }else{
      TString *names;
      names=new TString[nvars];
      TString answer="";
      Int_t checktrue=0;
      names=RDHFDplustoKpipi->GetVarNames();
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
      RDHFDplustoKpipi->SetVarsForOpt(dim,boolforopt);
    }
  }


  Float_t tighterval[dim][nptbins];
  //sigmavert
  //declength
  //Pmax
  //costhetapoint
  //sumd02

  //number of steps for each variable is 4 now
  //set this!!

  // tighterval[var][ptbin]
  tighterval[0][0]=0.022100;
  tighterval[0][1]=0.034;
  tighterval[0][2]=0.020667;
  tighterval[0][3]=0.023333;

  tighterval[1][0]=0.08;
  tighterval[1][1]=0.09;
  tighterval[1][2]=0.095;
  tighterval[1][3]=0.115;

  tighterval[2][0]=1.;
  tighterval[2][1]=1.;
  tighterval[2][2]=1.;
  tighterval[2][3]=1.;
  
  tighterval[3][0]=0.979;
  tighterval[3][1]=0.9975;
  tighterval[3][2]=0.995;
  tighterval[3][3]=0.9975;
  
  tighterval[4][0]=0.0055;
  tighterval[4][1]=0.0028;
  tighterval[4][2]=0.000883;
  tighterval[4][3]=0.000883;

  TString name=""; 
  Int_t arrdim=dim*nptbins;
  TClonesArray max("TParameter<float>",arrdim);
  for (Int_t ipt=0;ipt<nptbins;ipt++){
    for(Int_t ival=0;ival<dim;ival++){
      name=Form("par%dptbin%d",ival,ipt);
      new(max[ipt*dim+ival])TParameter<float>(name.Data(),tighterval[ival][ipt]);
    }
  }

  TFile* fout=new TFile("cuts4SignifMaximDplus.root","recreate");   //set this!! 
  fout->cd();
  RDHFDplustoKpipi->Write();
  max.Write();
  fout->Close();

}
