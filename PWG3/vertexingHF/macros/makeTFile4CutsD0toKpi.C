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

void makeInputAliAnalysisTaskSED0Mass(){

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetName("D0toKpiCuts");
  RDHFD0toKpi->SetTitle("Cuts for D0 analysis");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  //esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);


  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=9;

  const Int_t nptbins=11;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=1.;
  ptbins[2]=2.;
  ptbins[3]=3.;
  ptbins[4]=4.;
  ptbins[5]=5.;
  ptbins[6]=6.;
  ptbins[7]=8.;
  ptbins[8]=12.;
  ptbins[9]=16.;
  ptbins[10]=20.;
  ptbins[11]=24.;
  //ptbins[12]=99999.;
  
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
  Double_t arrcuts[9]={0.2,0.03,0.8,0.3,0.3,0.1,0.1,-0.0004,0.7};
  //cout<<"here"<<endl;
  //setting my cut values

  //0-1 GeV
  rdcutsvalmine[0][0]=arrcuts[0];
  rdcutsvalmine[1][0]=arrcuts[1];
  rdcutsvalmine[2][0]=arrcuts[2];
  rdcutsvalmine[3][0]=arrcuts[3];
  rdcutsvalmine[4][0]=arrcuts[4];
  rdcutsvalmine[5][0]=arrcuts[5];
  rdcutsvalmine[6][0]=arrcuts[6];
  rdcutsvalmine[7][0]=arrcuts[7];
  rdcutsvalmine[8][0]=arrcuts[8];

  //1-2 GeV
  arrcuts[1]=0.02; arrcuts[3]=arrcuts[4]=0.4; arrcuts[7]=-0.00032; arrcuts[8]=0.8;
  rdcutsvalmine[0][1]=arrcuts[0];
  rdcutsvalmine[1][1]=arrcuts[1];
  rdcutsvalmine[2][1]=arrcuts[2];
  rdcutsvalmine[3][1]=arrcuts[3];
  rdcutsvalmine[4][1]=arrcuts[4];
  rdcutsvalmine[5][1]=arrcuts[5];
  rdcutsvalmine[6][1]=arrcuts[6];
  rdcutsvalmine[7][1]=arrcuts[7];
  rdcutsvalmine[8][1]=arrcuts[8];

  //2-3 GeV
  arrcuts[3]=arrcuts[4]=0.7; arrcuts[7]=-0.00026;
  rdcutsvalmine[0][2]=arrcuts[0];
  rdcutsvalmine[1][2]=arrcuts[1];
  rdcutsvalmine[2][2]=arrcuts[2];
  rdcutsvalmine[3][2]=arrcuts[3];
  rdcutsvalmine[4][2]=arrcuts[4];
  rdcutsvalmine[5][2]=arrcuts[5];
  rdcutsvalmine[6][2]=arrcuts[6];
  rdcutsvalmine[7][2]=arrcuts[7];
  rdcutsvalmine[8][2]=arrcuts[8];


  //3-4 GeV - 4-5 GeV 
  arrcuts[7]=-0.00015;
  rdcutsvalmine[0][3]=rdcutsvalmine[0][4]=arrcuts[0];
  rdcutsvalmine[1][3]=rdcutsvalmine[1][4]=arrcuts[1];
  rdcutsvalmine[2][3]=rdcutsvalmine[2][4]=arrcuts[2];
  rdcutsvalmine[3][3]=rdcutsvalmine[3][4]=arrcuts[3];
  rdcutsvalmine[4][3]=rdcutsvalmine[4][4]=arrcuts[4];
  rdcutsvalmine[5][3]=rdcutsvalmine[5][4]=arrcuts[5];
  rdcutsvalmine[6][3]=rdcutsvalmine[6][4]=arrcuts[6];
  rdcutsvalmine[7][3]=rdcutsvalmine[7][4]=arrcuts[7];
  rdcutsvalmine[8][3]=rdcutsvalmine[8][4]=arrcuts[8];


  //5-6 GeV - 6-8 GeV - 8-12 GeV - 12-16 GeV
  arrcuts[1]=0.015;  arrcuts[7]=-0.0001;

  rdcutsvalmine[0][5]=rdcutsvalmine[0][6]=rdcutsvalmine[0][7]=rdcutsvalmine[0][8]=arrcuts[0];
  rdcutsvalmine[1][5]=rdcutsvalmine[1][6]=rdcutsvalmine[1][7]=rdcutsvalmine[1][8]=arrcuts[1];
  rdcutsvalmine[2][5]=rdcutsvalmine[2][6]=rdcutsvalmine[2][7]=rdcutsvalmine[2][8]=arrcuts[2];
  rdcutsvalmine[3][5]=rdcutsvalmine[3][6]=rdcutsvalmine[3][7]=rdcutsvalmine[3][8]=arrcuts[3];
  rdcutsvalmine[4][5]=rdcutsvalmine[4][6]=rdcutsvalmine[4][7]=rdcutsvalmine[4][8]=arrcuts[4];
  rdcutsvalmine[5][5]=rdcutsvalmine[5][6]=rdcutsvalmine[5][7]=rdcutsvalmine[5][8]=arrcuts[5];
  rdcutsvalmine[6][5]=rdcutsvalmine[6][6]=rdcutsvalmine[6][7]=rdcutsvalmine[6][8]=arrcuts[6];
  rdcutsvalmine[7][5]=rdcutsvalmine[7][6]=rdcutsvalmine[7][7]=rdcutsvalmine[7][8]=arrcuts[7];
  rdcutsvalmine[8][5]=rdcutsvalmine[8][6]=rdcutsvalmine[8][7]=rdcutsvalmine[8][8]=arrcuts[8];

  //16-20 GeV - 20-24 GeV --to be optimized
  arrcuts[7]=0.0001; arrcuts[8]=0.7;
  rdcutsvalmine[0][9]=rdcutsvalmine[0][10]=arrcuts[0];
  rdcutsvalmine[1][9]=rdcutsvalmine[1][10]=arrcuts[1];
  rdcutsvalmine[2][9]=rdcutsvalmine[2][10]=arrcuts[2];
  rdcutsvalmine[3][9]=rdcutsvalmine[3][10]=arrcuts[3];
  rdcutsvalmine[4][9]=rdcutsvalmine[4][10]=arrcuts[4];
  rdcutsvalmine[5][9]=rdcutsvalmine[5][10]=arrcuts[5];
  rdcutsvalmine[6][9]=rdcutsvalmine[6][10]=arrcuts[6];
  rdcutsvalmine[7][9]=rdcutsvalmine[7][10]=arrcuts[7];
  rdcutsvalmine[8][9]=rdcutsvalmine[8][10]=arrcuts[8];

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

  //activate pileup rejection
  RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  cout<<"This is the odject I'm going to save:"<<endl;
  RDHFD0toKpi->PrintAll();
  TFile* fout=new TFile("D0toKpiCutsStdpileup.root","recreate");   //set this!! 
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
  
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.9,0.9);
  esdTrackCuts->SetPtRange(0.1,1.e10);
  
  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  const Int_t nvars=9;

  const Int_t nptbins=10; //change this when adding pt bins!
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
  ptbins[9]=12.;
  //ptbins[10]=16.;
  ptbins[10]=9999.;

  
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
  /*
  //setting my cut values

  Float_t passcut[9]={0.7,0.04,0.8,0.7,0.7,0.1,0.1,-0.00015,0.5};

  //setting my cut values
  //0-1
  rdcutsvalmine[0][0]=passcut[0];//0.7;
  rdcutsvalmine[1][0]=passcut[1];//0.04;
  rdcutsvalmine[2][0]=passcut[2];//0.8;
  rdcutsvalmine[3][0]=passcut[3];//0.5;
  rdcutsvalmine[4][0]=passcut[4];//0.5;
  rdcutsvalmine[5][0]=passcut[5];//0.05;
  rdcutsvalmine[6][0]=passcut[6];//0.05;
  rdcutsvalmine[7][0]=passcut[7];//-0.00025;
  rdcutsvalmine[8][0]=passcut[8];//0.7;

  //1-2;2-3
  passcut[1]=0.02;
  passcut[8]=0.7;
  rdcutsvalmine[0][1]=rdcutsvalmine[0][2]=passcut[0];//0.7;
  rdcutsvalmine[1][1]=rdcutsvalmine[1][2]=passcut[1];//0.02;
  rdcutsvalmine[2][1]=rdcutsvalmine[2][2]=passcut[2];//0.8;
  rdcutsvalmine[3][1]=rdcutsvalmine[3][2]=passcut[3];//0.7;
  rdcutsvalmine[4][1]=rdcutsvalmine[4][2]=passcut[4];//0.7;
  rdcutsvalmine[5][1]=rdcutsvalmine[5][2]=passcut[5];//1.;
  rdcutsvalmine[6][1]=rdcutsvalmine[6][2]=passcut[6];//1.;
  rdcutsvalmine[7][1]=rdcutsvalmine[7][2]=passcut[7];//-0.00025;
  rdcutsvalmine[8][1]=rdcutsvalmine[8][2]=passcut[8];//0.8;

  //3-5
  passcut[7]=-0.00005;
  passcut[8]=0.6;
  rdcutsvalmine[0][3]=passcut[0];//0.7;
  rdcutsvalmine[1][3]=passcut[1];//0.02;
  rdcutsvalmine[2][3]=passcut[2];//0.8;
  rdcutsvalmine[3][3]=passcut[3];//0.7;
  rdcutsvalmine[4][3]=passcut[4];//0.7;
  rdcutsvalmine[5][3]=passcut[5];//0.05;
  rdcutsvalmine[6][3]=passcut[6];//0.05;
  rdcutsvalmine[7][3]=passcut[7];//-0.00015;
  rdcutsvalmine[8][3]=passcut[8];//0.8;

  //5-8
  passcut[7]=-0.00001;
  passcut[8]=0.6;
  rdcutsvalmine[0][4]=passcut[0];//0.7;
  rdcutsvalmine[1][4]=passcut[1];//0.02;
  rdcutsvalmine[2][4]=passcut[2];//0.8;
  rdcutsvalmine[3][4]=passcut[3];//0.7;
  rdcutsvalmine[4][4]=passcut[4];//0.7;
  rdcutsvalmine[5][4]=passcut[5];//0.05;
  rdcutsvalmine[6][4]=passcut[6];//0.05;
  rdcutsvalmine[7][4]=passcut[7];//-0.00015;
  rdcutsvalmine[8][4]=passcut[8];//0.9;

  //8-12
  passcut[7]=-0.00000;
  rdcutsvalmine[0][5]=passcut[0];
  rdcutsvalmine[1][5]=passcut[1];
  rdcutsvalmine[2][5]=passcut[2];
  rdcutsvalmine[3][5]=passcut[3];
  rdcutsvalmine[4][5]=passcut[4];
  rdcutsvalmine[5][5]=passcut[5];
  rdcutsvalmine[6][5]=passcut[6];
  rdcutsvalmine[7][5]=passcut[7];
  rdcutsvalmine[8][5]=passcut[8];

  //>12
  passcut[7]=-0.00000;
  rdcutsvalmine[0][6]=passcut[0];
  rdcutsvalmine[1][6]=passcut[1];
  rdcutsvalmine[2][6]=passcut[2];
  rdcutsvalmine[3][6]=passcut[3];
  rdcutsvalmine[4][6]=passcut[4];
  rdcutsvalmine[5][6]=passcut[5];
  rdcutsvalmine[6][6]=passcut[6];
  rdcutsvalmine[7][6]=passcut[7];
  rdcutsvalmine[8][6]=passcut[8];
  */

  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.2,400.*1E-4,0.7,0.3,0.3,1000.*1E-4,1000.*1E-4,-0.00005,0.85},/* pt<0.5*/
						  {0.2,400.*1E-4,0.7,0.3,0.3,1000.*1E-4,1000.*1E-4,-0.00005,0.85},/* 0.5<pt<1*/ 
						  {0.2,400.*1E-4,0.7,0.4,0.4,1000.*1E-4,1000.*1E-4,-0.00005,0.85},/* 1<pt<2 */
						  {0.2,400.*1E-4,0.7,0.7,0.7,1000.*1E-4,1000.*1E-4,-0.00005,0.7},/* 2<pt<3 */
						  {0.2,400.*1E-4,0.7,0.7,0.7,1000.*1E-4,1000.*1E-4,-0.00005,0.7},/* 3<pt<4 */
						  {0.2,400.*1E-4,0.7,0.7,0.7,1000.*1E-4,1000.*1E-4,-0.00005,0.7},/* 4<pt<5 */
						  {0.2,400.*1E-4,0.7,0.7,0.7,1000.*1E-4,1000.*1E-4,-0.00005,0.7},/* 5<pt<6 */
						  {0.2,400.*1E-4,0.7,0.7,0.7,1000.*1E-4,1000.*1E-4,-0.00005,0.7},/* 6<pt<8 */
						  {0.2,400.*1E-4,0.7,0.7,0.7,1000.*1E-4,1000.*1E-4,0.0002,0.7},/* 8<pt<12 */ 
						  {0.2,400.*1E-4,0.7,0.7,0.7,1000.*1E-4,1000.*1E-4,0.0002,0.7}};/*pt>12 */ 

  
  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];
  for (Int_t ibin=0;ibin<nptbins;ibin++){
    for (Int_t ivar = 0; ivar<nvars; ivar++){
      cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }
  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);


  Int_t nvarsforopt=RDHFD0toKpi->GetNVarsForOpt();
  Int_t dim=3; //set this!!
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
  //dca  <- this
  //costhetastar
  //d0d0 <-this 
  //costhetapoint <-this 

  
  //number of steps for each variable is 4 now
  //set this!!
  // tighterval[0][0]=0.01;
  tighterval[0][0]=100e-4;
  tighterval[1][0]=-0.0006;
  tighterval[2][0]=1.;

  // tighterval[0][1]=0.01;
  tighterval[0][1]=100e-4;
  tighterval[1][1]=-0.0007; //try with tighter in this bin
  tighterval[2][1]=1.;
 
 // tighterval[0][2]=0.01;
  tighterval[0][2]=100e-4;
  tighterval[1][2]=-0.0006;
  tighterval[2][2]=1.;
 
  // tighterval[0][3]=0.01;
  tighterval[0][3]=100e-4;
  tighterval[1][3]=-0.0006;
  tighterval[2][3]=0.95;

  // tighterval[0][4]=0.01;
  tighterval[0][4]=100e-4;
  tighterval[1][4]=-0.0006;
  tighterval[2][4]=0.95;

  // tighterval[0][5]=0.01;
  tighterval[0][5]=100e-4;
  tighterval[1][5]=-0.0006;
  tighterval[2][5]=0.95;

  // tighterval[0][6]=0.01;
  tighterval[0][6]=100e-4;
  tighterval[1][6]=-0.0006;
  tighterval[2][6]=0.95;

  // tighterval[0][6]=0.01;
  tighterval[0][7]=100e-4;
  tighterval[1][7]=-0.0006;
  tighterval[2][7]=0.95;

  // tighterval[0][6]=0.01;
  tighterval[0][8]=100e-4;
  tighterval[1][8]=-0.0006;
  tighterval[2][8]=0.95;

  // tighterval[0][6]=0.01;
  tighterval[0][9]=100e-4;
  tighterval[1][9]=-0.0006;
  tighterval[2][9]=0.95;


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
  
  //activate pileup rejection
  RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);


  TFile* fout=new TFile("cuts4SignifMaxim10binspileup.root","recreate");   //set this!! 
  fout->cd();
  RDHFD0toKpi->Write();
  max.Write();
  fout->Close();
 
}

