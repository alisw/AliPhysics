#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsD0toKpi.h>
#include <TClonesArray.h>
#include <TParameter.h>

//macro to make a .root file which contains an AliRDHFCutsD0toKpi with loose set of cuts (for significance maximization) and TParameter with the tighest value of these cuts
//Needed for AliAnalysisTaskSESignificance

//Use:
//Set hard coded commentet with //set this!!
//.x makeTFile4CutsD0toKpi.C++

//similar macros for the other D mesons

//Author: Chiara Bianchin, cbianchi@pd.infn.it

void makeTFile4CutsD0toKpi(){

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetName("loosercuts");
  RDHFD0toKpi->SetTitle("Cuts for significance maximization");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4);
  
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.9,0.9);
  esdTrackCuts->SetPtRange(0.1,1.e10);
  
  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=9;

  const Int_t nptbins=5;
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=5.;
  ptbins[5]=10.;
  
  RDHFD0toKpi->SetPtBins(nptbins+1,ptbins);
  

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }

  //setting my cut values
    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);

  /*
  //setting PPR cut values
  rdcutsvalPPR[0][0]=0.7;
  rdcutsvalPPR[1][0]=0.04;
  rdcutsvalPPR[2][0]=0.8;
  rdcutsvalPPR[3][0]=0.5;
  rdcutsvalPPR[4][0]=0.5;
  rdcutsvalPPR[5][0]=0.05;
  rdcutsvalPPR[6][0]=0.05;
  rdcutsvalPPR[7][0]=-0.0002;
  rdcutsvalPPR[8][0]=0.5;

  rdcutsvalPPR[0][1]=rdcutsvalPPR[0][2]=0.7;
  rdcutsvalPPR[1][1]=rdcutsvalPPR[1][2]=0.02;
  rdcutsvalPPR[2][1]=rdcutsvalPPR[2][2]=0.8;
  rdcutsvalPPR[3][1]=rdcutsvalPPR[3][2]=0.7;
  rdcutsvalPPR[4][1]=rdcutsvalPPR[4][2]=0.7;
  rdcutsvalPPR[5][1]=rdcutsvalPPR[5][2]=0.05;
  rdcutsvalPPR[6][1]=rdcutsvalPPR[6][2]=0.05;
  rdcutsvalPPR[7][1]=rdcutsvalPPR[7][2]=-0.0002;
  rdcutsvalPPR[8][1]=rdcutsvalPPR[8][2]=0.6;

  rdcutsvalPPR[0][3]=0.7;
  rdcutsvalPPR[1][3]=0.02;
  rdcutsvalPPR[2][3]=0.8;
  rdcutsvalPPR[3][3]=0.7;
  rdcutsvalPPR[4][3]=0.7;
  rdcutsvalPPR[5][3]=0.05;
  rdcutsvalPPR[6][3]=0.05;
  rdcutsvalPPR[7][3]=-0.0001;
  rdcutsvalPPR[8][3]=0.8;

  rdcutsvalPPR[0][4]=0.7;
  rdcutsvalPPR[1][4]=0.02;
  rdcutsvalPPR[2][4]=0.8;
  rdcutsvalPPR[3][4]=0.7;
  rdcutsvalPPR[4][4]=0.7;
  rdcutsvalPPR[5][4]=0.05;
  rdcutsvalPPR[6][4]=0.05;
  rdcutsvalPPR[7][4]=-0.00005;
  rdcutsvalPPR[8][4]=0.8;
  */
  //setting my cut values

  rdcutsvalmine[0][0]=0.065;
  rdcutsvalmine[1][0]=0.04;
  rdcutsvalmine[2][0]=0.8;
  rdcutsvalmine[3][0]=0.7;
  rdcutsvalmine[4][0]=0.7;
  rdcutsvalmine[5][0]=1000.*1E-4;
  rdcutsvalmine[6][0]=1000.*1E-4;
  rdcutsvalmine[7][0]=-0.0004;
  rdcutsvalmine[8][0]=0.71;

  rdcutsvalmine[0][1]=0.065;
  rdcutsvalmine[1][1]=0.04;
  rdcutsvalmine[2][1]=0.8;
  rdcutsvalmine[3][1]=0.7;
  rdcutsvalmine[4][1]=0.7;
  rdcutsvalmine[5][1]=1.;
  rdcutsvalmine[6][1]=1.;
  rdcutsvalmine[7][1]=-0.0003;
  rdcutsvalmine[8][1]=0.79;

  rdcutsvalmine[0][2]=0.65;
  rdcutsvalmine[1][2]=200.*1E-4;
  rdcutsvalmine[2][2]=0.8;
  rdcutsvalmine[3][2]=0.7;
  rdcutsvalmine[4][2]=0.7;
  rdcutsvalmine[5][2]=1000.*1E-4;
  rdcutsvalmine[6][2]=1000.*1E-4;
  rdcutsvalmine[7][2]=-0.0001;
  rdcutsvalmine[8][2]=0.83;

  rdcutsvalmine[0][3]=0.65;
  rdcutsvalmine[1][3]=200.*1E-4;
  rdcutsvalmine[2][3]=0.8;
  rdcutsvalmine[3][3]=0.7;
  rdcutsvalmine[4][3]=0.7;
  rdcutsvalmine[5][3]=1000.*1E-4;
  rdcutsvalmine[6][3]=1000.*1E-4;
  rdcutsvalmine[7][3]=-0.00005;
  rdcutsvalmine[8][3]=0.78;

  rdcutsvalmine[0][4]=0.65;
  rdcutsvalmine[1][4]=200.*1E-4;
  rdcutsvalmine[2][4]=0.8;
  rdcutsvalmine[3][4]=0.7;
  rdcutsvalmine[4][4]=0.7;
  rdcutsvalmine[5][4]=1000.*1E-4;
  rdcutsvalmine[6][4]=1000.*1E-4;
  rdcutsvalmine[7][4]=-0.00001;
  rdcutsvalmine[8][4]=0.79;


  RDHFD0toKpi->SetCuts(nvars,nptbins,rdcutsvalmine);

  Int_t nvarsforopt=RDHFD0toKpi->GetNVarsForOpt();
  Int_t dim=2; //set this!!
  Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(dim>nvarsforopt){
    cout<<"Number of variables for optimization has probably changed, check and edit accordingly"<<endl;
    return;
  } else {
    if(dim==nvarsforopt){
      boolforopt=RDHFD0toKpi->GetVarsForOpt();
    }else{
      TString *names;
      names=new TString[nvars];
      TString answer="";
      Int_t checktrue=0;
      names=RDHFD0toKpi->GetVarNames();
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
      RDHFD0toKpi->SetVarsForOpt(dim,boolforopt);
    }
  }


  Float_t tighterval[dim][nptbins];
  //dca
  //costhetastar
  //d0d0 <-this 
  //costhetapoint <-this 

  
  //number of steps for each variable is 4 now
  //set this!!
  // tighterval[0][0]=0.01;
  // tighterval[1][0]=0.8;
  tighterval[0][0]=-0.0007;
  tighterval[1][0]=0.99;

  // tighterval[0][1]=0.01;
  // tighterval[1][1]=0.8;
  tighterval[0][1]=-0.0006;
  tighterval[1][1]=0.99;
 
 // tighterval[0][2]=0.01;
  // tighterval[1][2]=0.8;
  tighterval[0][2]=-0.0004;
  tighterval[1][2]=0.99;
 
  // tighterval[0][3]=0.01;
  // tighterval[1][3]=0.8;
  tighterval[0][3]=-0.00035;
  tighterval[1][3]=0.98;

  // tighterval[0][4]=0.01;
  // tighterval[1][4]=0.8;
  tighterval[0][4]=-0.0003;
  tighterval[1][4]=0.99;


  TString name=""; 
  Int_t arrdim=dim*nptbins;
  cout<<"Will save "<<arrdim<<" TParameter<float>"<<endl;
  TClonesArray max("TParameter<float>",arrdim);
  for(Int_t i=0;i<dim;i++){
    for(Int_t j=0;j<nptbins;j++){
      name=Form("par%dptbin%d",i,j);
      cout<<"Setting "<<name.Data()<<" to "<<tighterval[i][j]<<endl;
      new(max[i*nptbins+j])TParameter<float>(name.Data(),tighterval[i][j]);
    }
  }
 
  TFile* fout=new TFile("cuts4SignifMaxim.root","recreate");   //set this!! 
  fout->cd();
  RDHFD0toKpi->Write();
  max.Write();
  fout->Close();
 
}
