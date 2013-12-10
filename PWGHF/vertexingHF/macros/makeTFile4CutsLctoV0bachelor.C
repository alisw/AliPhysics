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
  esdTrackCuts->SetPtRange(0.,1.e10);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
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
  esdTrackCutsV0daughters->SetEtaRange(-0.8,0.8);
  esdTrackCutsV0daughters->SetAcceptKinkDaughters(kFALSE);

  AliRDHFCutsLctoV0* RDHFLctoV0An=new AliRDHFCutsLctoV0();
  RDHFLctoV0An->SetName("LctoV0AnalysisCuts");
  RDHFLctoV0An->SetTitle("Analysis cuts for Lc analysis");
  RDHFLctoV0An->SetKinkRejection(!esdTrackCuts->GetAcceptKinkDaughters());
  RDHFLctoV0An->AddTrackCuts(esdTrackCuts);
  RDHFLctoV0An->AddTrackCutsV0daughters(esdTrackCutsV0daughters);
  RDHFLctoV0An->SetUseTrackSelectionWithFilterBits(kTRUE);
  RDHFLctoV0An->SetPidSelectionFlag(2); // 0 -> TOF AND TPC
                                        // 1 -> if (TOF) TOF else TPC w veto
                                        // 2 -> if (p>1) TPC@3s else if (1<=p<2.5) {if (TOF) TOF@3s AND TPC@3s} else (p>=2.5) {if (TOF) -2s<TOF<3s AND TPC@3s}
                                        // 3 -> if (p>1) TPC@3s else if (1<=p<2.5) {if (TOF) TOF@3s AND TPC@3s} else if (2.5<=p<3) {if (TOF) -2s<TOF<3s AND TPC@3s} else (p>=3) {if (TOF) -2s<TOF<3s AND -3s<TPC<2s}

  const Int_t nptbins=14;
  RDHFLctoV0An->SetNPtBins(nptbins);

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 0.;
  ptbins[1]= 1.;
  ptbins[2]= 2.;
  ptbins[3]= 3.;
  ptbins[4]= 4.;
  ptbins[5]= 5.;
  ptbins[6]= 6.;
  ptbins[7]= 7.;
  ptbins[8]= 8.;
  ptbins[9]= 9.;
  ptbins[10]=10.;
  ptbins[11]=12.;
  ptbins[12]=17.;
  ptbins[13]=25.;
  ptbins[14]=35.;
  RDHFLctoV0An->SetPtBins(nptbins+1,ptbins);

  const Int_t nvars=17;

  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
   anacutsval[0][ipt2]=0.25;   // inv. mass if K0S [GeV/c2]
   anacutsval[1][ipt2]=0.25;   // inv. mass if Lambda [GeV/c2]
   anacutsval[3][ipt2]=0.05;   // inv. mass V0 if Lambda [GeV/c2] ---> WE ARE SEARCHING Lc -> p+K0S, so cut on m(Lambda) has to be leave as it was at filtering level!!!
   anacutsval[5][ipt2]=0.2;    // pT min V0-positive track [GeV/c]
   anacutsval[6][ipt2]=0.2;    // pT min V0-negative track [GeV/c]
   anacutsval[7][ipt2]=1000.;  // dca cascade cut [cm]
   anacutsval[8][ipt2]=1.5;    // dca V0 cut [nSigma] // it's 1.5 x offline V0s
   anacutsval[11][ipt2]=0.09;  // d0 max V0 wrt PV [cm]
   anacutsval[12][ipt2]=0.;    // mass K0S veto [GeV/c2]
   anacutsval[13][ipt2]=0.005; // mass Lambda/LambdaBar veto [GeV/c2]
   anacutsval[14][ipt2]=0.300; // mass Gamma veto [GeV/c2]
   anacutsval[16][ipt2]=0.;    // V0 type cut
  }

  anacutsval[2][0]=0.00346*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][1]=0.00360*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][2]=0.00385*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][3]=0.00406*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][4]=0.00429*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][5]=0.00450*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][6]=0.00450*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][7]=0.00450*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][8]=0.00520*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][9]=0.00520*2.; // inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][10]=0.00520*2.;// inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][11]=0.00520*2.;// inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][12]=0.00520*2.;// inv. mass V0 if K0S [GeV/c2]
  anacutsval[2][13]=0.00520*2.;// inv. mass V0 if K0S [GeV/c2]

  anacutsval[4][0]=0.6;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][1]=0.5;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][2]=0.6;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][3]=0.7;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][4]=0.8;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][5]=0.9;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][6]=1.0;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][7]=1.2;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][8]=1.4;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][9]=1.7;    // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][10]=2.0;   // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][11]=2.0;   // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][12]=2.0;   // pT min bachelor track [GeV/c] // AOD by construction
  anacutsval[4][13]=2.0;   // pT min bachelor track [GeV/c] // AOD by construction

  anacutsval[9][0]=0.994;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][1]=0.994;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][2]=0.998;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][3]=0.999;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][4]=0.999;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][5]=0.999;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][6]=0.999;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][7]=0.999;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][8]=0.999;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][9]=0.999;  // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][10]=0.999; // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][11]=0.999; // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][12]=0.999; // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
  anacutsval[9][13]=0.999; // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level

  anacutsval[10][0]=0.05;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][1]=0.04;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][2]=0.03;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][3]=0.05;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][4]=0.05;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][5]=0.02;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][6]=0.02;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][7]=0.02;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][8]=0.02;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][9]=0.02;  // d0 max bachelor wrt PV [cm]
  anacutsval[10][10]=0.02; // d0 max bachelor wrt PV [cm]
  anacutsval[10][11]=0.02; // d0 max bachelor wrt PV [cm]
  anacutsval[10][12]=0.02; // d0 max bachelor wrt PV [cm]
  anacutsval[10][13]=0.02; // d0 max bachelor wrt PV [cm]


  anacutsval[15][0]=0.5; // pT min V0 track [GeV/c]
  anacutsval[15][1]=0.6; // pT min V0 track [GeV/c]
  anacutsval[15][2]=0.7; // pT min V0 track [GeV/c]
  anacutsval[15][3]=0.7; // pT min V0 track [GeV/c]
  anacutsval[15][4]=0.8; // pT min V0 track [GeV/c]
  anacutsval[15][5]=0.9; // pT min V0 track [GeV/c]
  anacutsval[15][6]=1.0; // pT min V0 track [GeV/c]
  anacutsval[15][7]=1.0; // pT min V0 track [GeV/c]
  anacutsval[15][8]=1.2; // pT min V0 track [GeV/c]
  anacutsval[15][9]=1.2; // pT min V0 track [GeV/c]
  anacutsval[15][10]=1.2; // pT min V0 track [GeV/c]
  anacutsval[15][11]=1.2; // pT min V0 track [GeV/c]
  anacutsval[15][12]=1.2; // pT min V0 track [GeV/c]
  anacutsval[15][13]=1.2; // pT min V0 track [GeV/c]



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


  // uncomment these lines for Baysian PID:
  // Double_t threshold=0.3;
  // SetupCombinedPID(RDHFLctoV0An  ,threshold);
  // RDHFLctoV0An  ->SetPIDStrategy(AliRDHFCutsLctoV0::kCombined);
  //


  //uncomment these lines to apply cuts with the KF package
  //RDHFLctoV0An->SetCutsStrategy(AliRDHFCutsLctoV0::kKF);
  //for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
  //   anacutsval[0]=1.;  //if <0., no topological constraint
  //   anacutsval[1]=2.;  //cut on the Chi2/Ndf
  // }

  Bool_t pidflag=kTRUE;
  RDHFLctoV0An->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

  cout<<"This is the (anal) object I'm going to save:"<<endl;
  RDHFLctoV0An->PrintAll();
  TFile* fout=new TFile("Lc2pK0SCuts.root","RECREATE"); 
  fout->cd();
  RDHFLctoV0An->Write();
  fout->Close();
  delete fout;

  delete pidObjBachelor;
  delete pidObjV0neg;
  delete pidObjV0pos;
  //delete RDHFLctoV0An;

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
  esdTrackCuts->SetRequireITSRefit(kFALSE);
  esdTrackCuts->SetMinNClustersITS(0);
  //esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.,1.e10);
  
  AliESDtrackCuts* esdTrackCutsV0daughters=new AliESDtrackCuts();
  esdTrackCutsV0daughters->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCutsV0daughters->SetRequireTPCRefit(kTRUE);
  esdTrackCutsV0daughters->SetMinNClustersTPC(70);
  esdTrackCutsV0daughters->SetRequireITSRefit(kFALSE);
  esdTrackCutsV0daughters->SetMinNClustersITS(0);
  //esdTrackCutsV0daughters->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny); // default is kBoth, otherwise kAny
  esdTrackCutsV0daughters->SetMinDCAToVertexXY(0.);
  esdTrackCutsV0daughters->SetPtRange(0.,1.e10);

  RDHFLctoV0->AddTrackCuts(esdTrackCuts);

  RDHFLctoV0->AddTrackCutsV0daughters(esdTrackCutsV0daughters);

  RDHFLctoV0->SetPidSelectionFlag(2); // 0 -> TOF AND TPC
                                      // 1 -> if (TOF) TOF else TPC w veto
                                      // 2 -> if (p>1) TPC@3s else if (1<=p<2.5) {if (TOF) TOF@3s AND TPC@3s} else (p>=2.5) {if (TOF) -2s<TOF<3s AND TPC@3s}
                                      // 3 -> if (p>1) TPC@3s else if (1<=p<2.5) {if (TOF) TOF@3s AND TPC@3s} else if (2.5<=p<3) {if (TOF) -2s<TOF<3s AND TPC@3s} else (p>=3) {if (TOF) -2s<TOF<3s AND -3s<TPC<2s}

  const Int_t nvars=14;

  const Int_t nptbins=9;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 0.;
  ptbins[1]= 1.;
  ptbins[2]= 2.;
  ptbins[3]= 3.;
  ptbins[4]= 4.;
  ptbins[5]= 5.;
  ptbins[6]= 6.;
  ptbins[7]= 8.;
  ptbins[8]=12.;
  ptbins[9]=99999999.;
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
  // dca V0 cut (nSigmas)
  // cosPA V0 cut
  // d0 max bachelor wrt PV [cm]
  // d0 max V0 wrt PV [cm]
  // mass K0S veto [GeV/c2]
  // mass Lambda/LambdaBar veto [GeV/c2]
  // mass Gamma veto [GeV/c2]
  // pT min V0 track [GeV/c]
  Float_t **cutsMatrixLctoV0Stand = new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++)
    cutsMatrixLctoV0Stand[ic]=new Float_t[nptbins];
  for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
    cutsMatrixLctoV0Stand[0][ipt2] =0.0075;
    cutsMatrixLctoV0Stand[1][ipt2] =0.0030;
    cutsMatrixLctoV0Stand[2][ipt2] =0.4;
    cutsMatrixLctoV0Stand[3][ipt2] =0.2;
    cutsMatrixLctoV0Stand[4][ipt2] =0.2;
    cutsMatrixLctoV0Stand[5][ipt2] =1000.;
    cutsMatrixLctoV0Stand[6][ipt2] =1.5;
    cutsMatrixLctoV0Stand[7][ipt2] =0.99;
    cutsMatrixLctoV0Stand[8][ipt2] =0.05;
    cutsMatrixLctoV0Stand[9][ipt2] =0.1;
    cutsMatrixLctoV0Stand[10][ipt2] =0.0;
    cutsMatrixLctoV0Stand[11][ipt2] =0.005;
    cutsMatrixLctoV0Stand[12][ipt2] =0.100;
  }
  cutsMatrixLctoV0Stand[0][13]=0.0; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[1][13]=0.6; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[2][13]=0.8; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[3][13]=0.8; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[4][13]=0.8; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[5][13]=1.0; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[6][13]=1.0; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[7][13]=1.0; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[8][13]=0.0; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[9][13]=0.0; // pT min V0 track [GeV/c]


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
  Int_t dim=14; //set this!!
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
  // 6(8): dca V0 cut (nSigmas)
  // 7(9): cosPA V0 cut
  // 8(10): d0 max bachelor wrt PV [cm]
  // 9(11): d0 max V0 wrt PV [cm]
  // 10(12): mass K0S veto [GeV/c2]
  // 11(13): mass Lambda/LambdaBar veto [GeV/c2]
  // 12(14): mass Gamma veto [GeV/c2]
  // 13(15): pT min V0 track [GeV/c]

  // number of steps for each variable is set in the AddTask arguments (default=8)
  // set this!!
  for(Int_t ipt=0;ipt<nptbins;ipt++){
    tighterval[0][ipt] =0.075;  // inv. mass V0 if K0S [GeV]
    tighterval[1][ipt] =0.040;  // inv. mass V0 if Lambda [GeV]
    tighterval[2][ipt] =0.4;    // pT min bachelor track [GeV/c]
    tighterval[3][ipt] =0.2;    // pT min V0-positive track [GeV/c]
    tighterval[4][ipt] =0.2;    // pT min V0-negative track [GeV/c]
    tighterval[5][ipt] =100.;   // dca cascade cut (cm)
    tighterval[6][ipt] =1.5;    // dca v0 cut
    tighterval[7][ipt] =0.99;   // cosPA v0 cut
    tighterval[8][ipt] =0.05;   // d0 max bachelor wrt PV [cm]
    tighterval[9][ipt] =0.1;    // d0 max V0 wrt PV [cm]
    tighterval[10][ipt] =0.0;   // mass K0S veto [GeV/c2]
    tighterval[11][ipt] =0.005; // mass Lambda/LambdaBar veto [GeV/c2]
    tighterval[12][ipt] =0.100; // mass Gamma veto [GeV/c2]
  }
  tighterval[0][13]=0.0; // pT min V0 track [GeV/c]
  tighterval[1][13]=0.6; // pT min V0 track [GeV/c]
  tighterval[2][13]=0.8; // pT min V0 track [GeV/c]
  tighterval[3][13]=0.8; // pT min V0 track [GeV/c]
  tighterval[4][13]=0.8; // pT min V0 track [GeV/c]
  tighterval[5][13]=1.0; // pT min V0 track [GeV/c]
  tighterval[6][13]=1.0; // pT min V0 track [GeV/c]
  tighterval[7][13]=1.0; // pT min V0 track [GeV/c]
  tighterval[8][13]=0.0; // pT min V0 track [GeV/c]
  tighterval[9][13]=0.0; // pT min V0 track [GeV/c]

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
