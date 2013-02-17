#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsD0toKpi.h>
#include <TClonesArray.h>
#include <TParameter.h>


//Use:
//Set hard coded commentet with //set this!!
// root[] .L makeInputD0tasks.C++
// root[] makeInputAliAnalysisTaskSED0Mass()
// root[] makeInputAliAnalysisTaskSESignificanceMaximization()
//similar macros for the other D mesons

//Author: Chiara Bianchin, cbianchi@pd.infn.it


//macro to make a .root file which contains an AliRDHFCutsD0toKpi for AliAnalysisTaskSED0Mass task
void ModifyFromStandardCuts(Int_t system=1 /*0=pp, 1=PbPb, 2=pp 2.76TeV*/){
  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  TString info="",cent="";

  if(system==0) {
    RDHFD0toKpi->SetStandardCutsPP2010();
    info+="ppStd";
  }
  if(system==1) {
    RDHFD0toKpi->SetStandardCutsPbPb2010();
    //Change centrality if needed
    /*
    Float_t minc=20,maxc=80;
    RDHFD0toKpi->SetMinCentrality(minc);
    RDHFD0toKpi->SetMaxCentrality(maxc);
    */
    RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

    cent=Form("%.0f%.0f",RDHFD0toKpi->GetMinCentrality(),RDHFD0toKpi->GetMaxCentrality());
    info+="PbPbStd";

  }
  if(system==2) {
    RDHFD0toKpi->SetStandardCutsPP2011_276TeV();
    info+="pp276Std";
  }

  RDHFD0toKpi->SetName("D0toKpiCuts");
  RDHFD0toKpi->SetTitle("Cuts for D0 analysis");
  

  //here add what you need! Find examples below

  //Trigger mask
  /*
  RDHFD0toKpi->SetTriggerMask(0);
  RDHFD0toKpi->SetTriggerMask(AliVEvent::kEMC1 | AliVEvent::kEMC7);
  RDHFD0toKpi->SetTriggerClass("CEMC");

  info+="EMCTr";
  */

  //Event selection
  RDHFD0toKpi->SetUsePhysicsSelection(kFALSE);
  info+="noPhysSel";

  cout<<"This is the odject I'm going to save:"<<endl;
  cout<<info<<endl;
  RDHFD0toKpi->PrintAll();
  TFile* fout=new TFile(Form("D0toKpiCuts%s%sRecVtx%sPileupRej%s.root", RDHFD0toKpi->GetUseCentrality()==0 ? "" : cent.Data(),RDHFD0toKpi->GetIsPrimaryWithoutDaughters() ? "" : "No",RDHFD0toKpi->GetOptPileUp() ? "" : "No",info.Data()),"recreate");   //set this!! 

  fout->cd();
  RDHFD0toKpi->Write();
  fout->Close();

}

void makeInputAliAnalysisTaskSED0Mass(){

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetName("D0toKpiCuts");
  RDHFD0toKpi->SetTitle("Cuts for D0 analysis");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  //esdTrackCuts->SetMinNClustersTPC(120);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.8,1.e10);


  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=11;

  const Int_t nptbins=13;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=8.;
  ptbins[9]=12.;
  ptbins[10]=16.;
  ptbins[11]=20.;
  ptbins[12]=24.;
  ptbins[13]=9999.;

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
    //     printf("    |cosThetaPointXY| < %f\n",fD0toKpiCuts[9]);
    //     printf("    NormDecayLenghtXY    > %f\n",fD0toKpiCuts[10]);


  Double_t arrcuts[11]={0.3,0.03,0.8,0.8,0.8,0.1,0.1,-0.0004,0.9,0.998,5.}; //put the last 2 values at 0. for pp

  //setting my cut values
  //0-0.5 GeV
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][0]=arrcuts[ic];

  //0.5-1 GeV/c
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][1]=arrcuts[ic];

  //1-2 GeV 
  arrcuts[1]=0.025;  arrcuts[7]=-0.0003;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][2]=arrcuts[ic];

  //2-3 GeV
  arrcuts[7]=-0.00026;
  arrcuts[9]=0.998;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][3]=arrcuts[ic];

  //3-4 GeV
  arrcuts[7]=-0.00015; arrcuts[8]=0.85;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][4]=arrcuts[ic];

  //4-5 GeV 
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][5]=arrcuts[ic];

  //5-6 GeV
  arrcuts[7]=-0.0001; arrcuts[8]=0.85;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][6]=arrcuts[ic];

  //6-8 GeV
  arrcuts[2]=1.;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][7]=arrcuts[ic];

  //8-12 GeV
  arrcuts[8]=0.8;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][8]=arrcuts[ic];

  //12-16 GeV
  arrcuts[1]=0.03;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][9]=arrcuts[ic];

  //16-20 GeV
  arrcuts[1]=0.035;
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][10]=arrcuts[ic];

  //20-24 GeV
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][11]=arrcuts[ic];

  //24-9999 GeV
  for(Int_t ic=0;ic<nvars;ic++) rdcutsvalmine[ic][12]=arrcuts[ic];

  RDHFD0toKpi->SetCuts(nvars,nptbins,rdcutsvalmine);

  Bool_t pidflag=kTRUE;
  RDHFD0toKpi->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

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
  RDHFD0toKpi->SetPidHF(pidObj);

  RDHFD0toKpi->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF

  //activate pileup rejection (for pp)
  //RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  //Do not recalculate the vertex
  RDHFD0toKpi->SetRemoveDaughtersFromPrim(kFALSE); //activate for pp

  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=20,maxc=80;
  RDHFD0toKpi->SetMinCentrality(minc);
  RDHFD0toKpi->SetMaxCentrality(maxc);
  cent=Form("%.0f%.0f",minc,maxc);
  RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

  //temporary
  //RDHFD0toKpi->SetFixRefs();

  cout<<"This is the odject I'm going to save:"<<endl;
  RDHFD0toKpi->PrintAll();
  TFile* fout=new TFile(Form("D0toKpiCuts%s%s%sRecVtx%sPileupRej.root", RDHFD0toKpi->GetUseCentrality()==0 ? "pp" : "PbPb",cent.Data(),RDHFD0toKpi->GetIsPrimaryWithoutDaughters() ? "" : "No",RDHFD0toKpi->GetOptPileUp() ? "" : "No"),"recreate");   //set this!! 

  fout->cd();
  RDHFD0toKpi->Write();
  fout->Close();

}
 
//macro to make a .root file (for significance maximization) which contains an AliRDHFCutsD0toKpi with loose set of cuts  and TParameter with the tighest value of these cuts

void makeInputAliAnalysisTaskSESignificanceMaximization(){

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetName("loosercuts");
  RDHFD0toKpi->SetTitle("Cuts for significance maximization");

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
  esdTrackCuts->SetPtRange(0.8,1.e10);
  
  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=11;

  const Int_t nptbins=14; //change this when adding pt bins!
  Float_t ptbins[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;
  ptbins[2]=1.;
  ptbins[3]=2.;
  ptbins[4]=3.;
  ptbins[5]=4.;
  ptbins[6]=5.;
  ptbins[7]=6.;
  ptbins[8]=8.;
  ptbins[9]=10.;
  ptbins[10]=12.;
  ptbins[11]=16.;
  ptbins[12]=20.;
  ptbins[13]=24.;
  ptbins[14]=9999.;

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
    //     printf("    |cosThetaPointXY| < %f\n",fD0toKpiCuts[9]);
    //     printf("    NormDecayLenghtXY    > %f\n",fD0toKpiCuts[10]);


    Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,400.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,-40000.*1E-8,0.75,0.,2.},/* pt<0.5*/
						  {0.400,400.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,-40000.*1E-8,0.75,0.,2.},/* 0.5<pt<1*/
						  {0.400,400.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,-33000.*1E-8,0.75,0.,2.},/* 1<pt<2 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-30000.*1E-8,0.85,0.994,2.},/* 2<pt<3 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-26000.*1E-8,0.85,0.994,2.},/* 3<pt<4 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-15000.*1E-8,0.85,0.994,2.},/* 4<pt<5 */
						  {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-11000.*1E-8,0.82,0.994,2.},/* 5<pt<6 */
						  {0.400,270.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-10000.*1E-8,0.78,0.994,2.},/* 6<pt<8 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-1000.*1E-8,0.7,0.994,2.},/* 8<pt<10 */
						  {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,-1000.*1E-8,0.7,0.994,2.},/* 10<pt<12 */
						  {0.400,350.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-500.*1E-8,0.7,0.994,2.},/* 12<pt<16 */
						  {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-500.*1E-8,0.7,0.994,2.},/* 16<pt<20 */
						  {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-500.*1E-8,0.7,0.994,2.},/* 20<pt<24 */
						  {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,-500.*1E-8,0.7,0.994,2.}};/* pt>24 */

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }
  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);


  Bool_t stdvaropt=kFALSE;
  Int_t dim=4; //set this!!
  Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(stdvaropt){
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


  Float_t tighterval[dim][nptbins];
  //dca  
  //costhetastar
  //d0d0 <-this 
  //costhetapoint <-this 
  //cosThetaPointXY <-this 
  //NormDecayLength <-this 
  
  //number of steps for each variable is set in the AddTask arguments (default=8)
  //set this!!
  //0-0.5
  tighterval[0][0]=-0.00065;
  tighterval[1][0]=1.;
  tighterval[2][0]=0.4;
  tighterval[3][0]=6.;

  //0.5-1.
  tighterval[0][1]=-0.00065;
  tighterval[1][1]=1.;
  tighterval[2][1]=0.4;
  tighterval[3][1]=6.;

  //1-2
  tighterval[0][2]=-0.00065;
  tighterval[1][2]=1.;
  tighterval[2][2]=0.4;
  tighterval[3][2]=6.;
 
  //2-3
  tighterval[0][3]=-0.0006;
  tighterval[1][3]=1.;
  tighterval[2][3]=1.;
  tighterval[3][3]=6.;

  //3-4
  tighterval[0][4]=-0.00046;
  tighterval[1][4]=1.;
  tighterval[2][4]=1.;
  tighterval[3][4]=6.;
 
  //4-5
  tighterval[0][5]=-0.00045;
  tighterval[1][5]=1.;
  tighterval[2][5]=1.;
  tighterval[3][5]=6.;
 
  //5-6
  tighterval[0][6]=-0.00031;
  tighterval[1][6]=1.;
  tighterval[2][6]=1.;
  tighterval[3][6]=6.;

  //6-8
  tighterval[0][7]=-0.00021;
  tighterval[1][7]=0.98;
  tighterval[2][7]=1.;
  tighterval[3][7]=6.;

  //8-10
  tighterval[0][8]=-0.0001;
  tighterval[1][8]=0.98;
  tighterval[2][8]=1.;
  tighterval[3][8]=6.;

  //10-12
  tighterval[0][9]=-0.0001;
  tighterval[1][9]=0.9;
  tighterval[2][9]=1.;
  tighterval[3][9]=6.;
 
  //12-16
  tighterval[0][10]=-0.00005;
  tighterval[1][10]=0.9;
  tighterval[2][10]=1.;
  tighterval[3][10]=6.;

  //16-20
  tighterval[0][11]=-0.00005;
  tighterval[1][11]=0.9;
  tighterval[2][11]=1.;
  tighterval[3][11]=6.;

  //20-24
  tighterval[0][12]=-0.00005;
  tighterval[1][12]=0.9;
  tighterval[2][12]=1.;
  tighterval[3][12]=6.;

  //>24
  tighterval[0][13]=-0.00005;
  tighterval[1][13]=0.9;
  tighterval[2][13]=1.;
  tighterval[3][13]=6.;


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

  Bool_t flagPID=kTRUE;
  RDHFD0toKpi->SetUsePID(flagPID);

  RDHFD0toKpi->PrintAll();
  printf("Use PID? %s\n",flagPID ? "yes" : "no");

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
  RDHFD0toKpi->SetPidHF(pidObj);

  RDHFD0toKpi->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF
  
  //activate pileup rejection (for pp)
  //RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  //Do not recalculate the vertex
  RDHFD0toKpi->SetRemoveDaughtersFromPrim(kFALSE); //activate for pp

  TString cent="";
  //centrality selection (Pb-Pb)
  Float_t minc=20,maxc=80;
  RDHFD0toKpi->SetMinCentrality(minc);
  RDHFD0toKpi->SetMaxCentrality(maxc);
  cent=Form("%.0f%.0f",minc,maxc);
  RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

  //temporary
  RDHFD0toKpi->SetFixRefs();

  TFile* fout=new TFile(Form("cuts4SignifMaxim%s%s%sRecVtx%sPileupRej.root", RDHFD0toKpi->GetUseCentrality()==0 ? "pp" : "PbPb",cent.Data(),RDHFD0toKpi->GetIsPrimaryWithoutDaughters() ? "" : "No",RDHFD0toKpi->GetOptPileUp() ? "" : "No"),"recreate");   //set this!! 

  fout->cd();
  RDHFD0toKpi->Write();
  max.Write();
  fout->Close();
 
}

