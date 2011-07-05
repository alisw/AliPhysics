// gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/PWG3 -I$ALICE_ROOT/PWG3/vertexingHF -g");

#include <Riostream.h>
#include <TFile.h>
#include "AliRDHFCutsDStartoKpipi.h"
#include <TClonesArray.h>
#include <TParameter.h>


//Use:
//Set hard coded commentet with //set this!!
// root[] .L makeTFile4CutsDStar.....C++
// root[] makeInputAliAnalysisTaskSED0Mass()
// root[] makeInputAliAnalysisTaskSESignificanceMaximization()
//similar macros for the other D mesons

//Author: Alessandro Grelli, a.grelli@uu.nl


//macro to make a .root file which contains an AliRDHFCutsDStartoKpipi for AliAnalysisTaskSEDStarSpectra task and CF task
void makeInputAliAnalysisTaskSEDStarSpectra(const char *set_cuts="utrecht"){

  AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi();
  RDHFDStartoKpipi->SetName("DStartoKpipiCuts");
  RDHFDStartoKpipi->SetTitle("Cuts for D* analysis");

  //Centrality selection
  RDHFDStartoKpipi->SetUseCentrality(kFALSE);
  RDHFDStartoKpipi->SetMinCentrality(40);
  RDHFDStartoKpipi->SetMaxCentrality(80);

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  //esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny); 
 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);

 // soft pion pre-selections 
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdSoftPicuts->SetRequireTPCRefit(kFALSE);
  esdSoftPicuts->SetRequireITSRefit(kFALSE);
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny); //test d0 asimmetry
  esdSoftPicuts->SetPtRange(0.0,1.e10);

  // set pre selections
  RDHFDStartoKpipi->AddTrackCuts(esdTrackCuts);
  RDHFDStartoKpipi->AddTrackCutsSoftPi(esdSoftPicuts);

  const Int_t nvars=16;
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
  ptbins[8]=7.;
  ptbins[9]=8.;
  ptbins[10]=12.;
  ptbins[11]=16.;
  ptbins[12]=24.;
  ptbins[13]=999.;

  RDHFDStartoKpipi->SetPtBins(nptbins+1,ptbins);
  
  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }

  if(set_cuts=="utrecht"){
    //0-0.5
    rdcutsvalmine[0][0]=0.7;
    rdcutsvalmine[1][0]=0.03;
    rdcutsvalmine[2][0]=0.8;
    rdcutsvalmine[3][0]=0.3;
    rdcutsvalmine[4][0]=0.3;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=0.00035;
    rdcutsvalmine[8][0]=0.73;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;
    //0.5-1
    rdcutsvalmine[0][1]=0.7;
    rdcutsvalmine[1][1]=0.03;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.3;
    rdcutsvalmine[4][1]=0.3;
    rdcutsvalmine[5][1]=0.1;
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00020;
    rdcutsvalmine[8][1]=0.73;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;
    //1-2
    rdcutsvalmine[0][2]=0.7;
    rdcutsvalmine[1][2]=0.02;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.4;
    rdcutsvalmine[4][2]=0.4;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.00036;
    rdcutsvalmine[8][2]=0.82;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.7;
    rdcutsvalmine[1][3]=0.02;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=-0.00016;
    rdcutsvalmine[8][3]=0.90;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;
    //3-4
    rdcutsvalmine[0][4]=0.7;
    rdcutsvalmine[1][4]=0.05;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.;
    rdcutsvalmine[4][4]=1.;
    rdcutsvalmine[5][4]=0.042;
    rdcutsvalmine[6][4]=0.056;
    rdcutsvalmine[7][4]=-0.000065;
    rdcutsvalmine[8][4]=0.9;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.1;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=100.;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=-1.;
    rdcutsvalmine[15][4]=0.;
    //4-5
    rdcutsvalmine[0][5]=0.7;
    rdcutsvalmine[1][5]=0.08;
    rdcutsvalmine[2][5]=0.9;
    rdcutsvalmine[3][5]=1.2;
    rdcutsvalmine[4][5]=1.2;
    rdcutsvalmine[5][5]=0.07;
    rdcutsvalmine[6][5]=0.07;
    rdcutsvalmine[7][5]=0.0001;
    rdcutsvalmine[8][5]=0.9;
    rdcutsvalmine[9][5]=0.3;
    rdcutsvalmine[10][5]=0.1;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=100.;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=-1.;
    rdcutsvalmine[15][5]=0.;
    //5-6
    rdcutsvalmine[0][6]=0.7;
    rdcutsvalmine[1][6]=0.1;
    rdcutsvalmine[2][6]=1.0;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.08;
    rdcutsvalmine[6][6]=0.08;
    rdcutsvalmine[7][6]=0.0005;
    rdcutsvalmine[8][6]=0.8;
    rdcutsvalmine[9][6]=0.3;
    rdcutsvalmine[10][6]=0.1;
    rdcutsvalmine[11][6]=0.05;
    rdcutsvalmine[12][6]=100000.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=-1.;
    rdcutsvalmine[15][6]=0.;
    //6-7
    rdcutsvalmine[0][7]=0.7;
    rdcutsvalmine[1][7]=0.1;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=0.001;
    rdcutsvalmine[8][7]=0.7;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.05;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=-1.;
    rdcutsvalmine[15][7]=0.;
    //7-8
    rdcutsvalmine[0][8]=0.7;
    rdcutsvalmine[1][8]=0.1;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=0.001;
    rdcutsvalmine[8][8]=0.7;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.1;
    rdcutsvalmine[11][8]=0.05;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=-1.;
    rdcutsvalmine[15][8]=0.;
    //8-12
    rdcutsvalmine[0][9]=0.7;
    rdcutsvalmine[1][9]=0.1;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.;
    rdcutsvalmine[4][9]=1.;
    rdcutsvalmine[5][9]=0.1;
    rdcutsvalmine[6][9]=0.1;
    rdcutsvalmine[7][9]=0.006;
    rdcutsvalmine[8][9]=0.7;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.1;
    rdcutsvalmine[11][9]=0.05;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=-1.;
    rdcutsvalmine[15][9]=0.;
    //12-16
    rdcutsvalmine[0][10]=0.7;
    rdcutsvalmine[1][10]=0.1;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.3;
    rdcutsvalmine[4][10]=0.3;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=0.01;
    rdcutsvalmine[8][10]=0.7;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.1;
    rdcutsvalmine[11][10]=0.05;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=-1.;
    rdcutsvalmine[15][10]=0.;
    //16-24
    rdcutsvalmine[0][11]=0.7;
    rdcutsvalmine[1][11]=0.2;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=.3;
    rdcutsvalmine[4][11]=.3;
    rdcutsvalmine[5][11]=0.15;
    rdcutsvalmine[6][11]=0.15;
    rdcutsvalmine[7][11]=0.01;
    rdcutsvalmine[8][11]=0.7;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.05;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=0.5;
    rdcutsvalmine[14][11]=-1.;
    rdcutsvalmine[15][11]=0.;
    //>24
    rdcutsvalmine[0][12]=0.7;
    rdcutsvalmine[1][12]=0.6;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=.4;
    rdcutsvalmine[4][12]=.4;
    rdcutsvalmine[5][12]=0.5;
    rdcutsvalmine[6][12]=0.5;
    rdcutsvalmine[7][12]=0.1;
    rdcutsvalmine[8][12]=0.7;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=0.5;
    rdcutsvalmine[14][12]=-1.;
    rdcutsvalmine[15][12]=0.;
  }
  if(set_cuts=="heidelberg"){

    //0-0.5
    rdcutsvalmine[0][0]=0.7;
    rdcutsvalmine[1][0]=0.03;
    rdcutsvalmine[2][0]=0.7;
    rdcutsvalmine[3][0]=0.8;
    rdcutsvalmine[4][0]=0.8;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=-0.00002;
    rdcutsvalmine[8][0]=0.9;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;
    //0.5-1
    rdcutsvalmine[0][1]=0.7;
    rdcutsvalmine[1][1]=0.03;
    rdcutsvalmine[2][1]=0.7;
    rdcutsvalmine[3][1]=0.8;
    rdcutsvalmine[4][1]=0.8;
    rdcutsvalmine[5][1]=0.1;
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00002;
    rdcutsvalmine[8][1]=0.9;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;

    //1-2
    rdcutsvalmine[0][2]=0.7;
    rdcutsvalmine[1][2]=0.03;
    rdcutsvalmine[2][2]=0.7;
    rdcutsvalmine[3][2]=0.8;
    rdcutsvalmine[4][2]=0.8;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.00002;
    rdcutsvalmine[8][2]=0.9;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.7;
    rdcutsvalmine[1][3]=0.03;
    rdcutsvalmine[2][3]=0.7;
    rdcutsvalmine[3][3]=0.8;
    rdcutsvalmine[4][3]=0.8;
    rdcutsvalmine[5][3]=0.1;
    rdcutsvalmine[6][3]=0.1;
    rdcutsvalmine[7][3]=-0.00002;
    rdcutsvalmine[8][3]=0.9;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;
    //3-4
    rdcutsvalmine[0][4]=0.7;
    rdcutsvalmine[1][4]=0.03;
    rdcutsvalmine[2][4]=0.7;
    rdcutsvalmine[3][4]=0.9;
    rdcutsvalmine[4][4]=0.9;
    rdcutsvalmine[5][4]=0.1;
    rdcutsvalmine[6][4]=0.1;
    rdcutsvalmine[7][4]=0.000002;
    rdcutsvalmine[8][4]=0.8;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.1;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=100.;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=-1.;
    rdcutsvalmine[15][4]=0.;
    //4-5
    rdcutsvalmine[0][5]=0.7;
    rdcutsvalmine[1][5]=0.03;
    rdcutsvalmine[2][5]=0.7;
    rdcutsvalmine[3][5]=0.9;
    rdcutsvalmine[4][5]=0.9;
    rdcutsvalmine[5][5]=0.1;
    rdcutsvalmine[6][5]=0.1;
    rdcutsvalmine[7][5]=0.000002;
    rdcutsvalmine[8][5]=0.8;
    rdcutsvalmine[9][5]=0.3;
    rdcutsvalmine[10][5]=0.1;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=100.;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=-1.;
    rdcutsvalmine[15][5]=0.;
    //5-6
    rdcutsvalmine[0][6]=0.7;
    rdcutsvalmine[1][6]=0.03;
    rdcutsvalmine[2][6]=0.7;
    rdcutsvalmine[3][6]=1.0;
    rdcutsvalmine[4][6]=1.0;
    rdcutsvalmine[5][6]=0.1;
    rdcutsvalmine[6][6]=0.1;
    rdcutsvalmine[7][6]=0.000002;
    rdcutsvalmine[8][6]=0.8;
    rdcutsvalmine[9][6]=0.3;
    rdcutsvalmine[10][6]=0.1;
    rdcutsvalmine[11][6]=0.05;
    rdcutsvalmine[12][6]=100.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=-1.;
    rdcutsvalmine[15][6]=0.;
    //6-8
    rdcutsvalmine[0][7]=0.7;
    rdcutsvalmine[1][7]=0.03;
    rdcutsvalmine[2][7]=0.7;
    rdcutsvalmine[3][7]=1.0;
    rdcutsvalmine[4][7]=1.0;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=0.000002;
    rdcutsvalmine[8][7]=0.8;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.05;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=-1.;
    rdcutsvalmine[15][7]=0.;
    //8-12
    rdcutsvalmine[0][8]=0.7;
    rdcutsvalmine[1][8]=0.03;
    rdcutsvalmine[2][8]=0.7;
    rdcutsvalmine[3][8]=1.0;
    rdcutsvalmine[4][8]=1.0;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=0.000002;
    rdcutsvalmine[8][8]=0.8;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.1;
    rdcutsvalmine[11][8]=0.05;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=-1.;
    rdcutsvalmine[15][8]=0.;
    //12-16
    rdcutsvalmine[0][9]=0.7;
    rdcutsvalmine[1][9]=0.03;
    rdcutsvalmine[2][9]=0.7;
    rdcutsvalmine[3][9]=1.0;
    rdcutsvalmine[4][9]=1.0;
    rdcutsvalmine[5][9]=0.1;
    rdcutsvalmine[6][9]=0.1;
    rdcutsvalmine[7][9]=0.000002;
    rdcutsvalmine[8][9]=0.8;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.1;
    rdcutsvalmine[11][9]=0.05;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=-1.;
    rdcutsvalmine[15][9]=0.;
    //16-20
    rdcutsvalmine[0][10]=0.7;
    rdcutsvalmine[1][10]=0.03;
    rdcutsvalmine[2][10]=0.7;
    rdcutsvalmine[3][10]=1.0;
    rdcutsvalmine[4][10]=1.0;
    rdcutsvalmine[5][10]=0.1;
    rdcutsvalmine[6][10]=0.1;
    rdcutsvalmine[7][10]=0.000002;
    rdcutsvalmine[8][10]=0.8;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.1;
    rdcutsvalmine[11][10]=0.05;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=-1.;
    rdcutsvalmine[15][10]=0.;
    //20-24
    rdcutsvalmine[0][11]=0.7;
    rdcutsvalmine[1][11]=0.03;
    rdcutsvalmine[2][11]=0.7;
    rdcutsvalmine[3][11]=1.0;
    rdcutsvalmine[4][11]=1.0;
    rdcutsvalmine[5][11]=0.1;
    rdcutsvalmine[6][11]=0.1;
    rdcutsvalmine[7][11]=0.000002;
    rdcutsvalmine[8][11]=0.8;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.05;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=0.5;
    rdcutsvalmine[14][11]=-1.;
    rdcutsvalmine[15][11]=0.;
    //>24
    rdcutsvalmine[0][12]=0.7;
    rdcutsvalmine[1][12]=0.03;
    rdcutsvalmine[2][12]=0.7;
    rdcutsvalmine[3][12]=1.0;
    rdcutsvalmine[4][12]=1.0;
    rdcutsvalmine[5][12]=0.1;
    rdcutsvalmine[6][12]=0.1;
    rdcutsvalmine[7][12]=0.000002;
    rdcutsvalmine[8][12]=0.8;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=0.5;
    rdcutsvalmine[14][12]=-1.;
    rdcutsvalmine[15][12]=0.;
  }

  RDHFDStartoKpipi->SetCuts(nvars,nptbins,rdcutsvalmine);

  Bool_t pidflag=kFALSE;
  RDHFDStartoKpipi->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

  // PID SETTINGS
  AliAODPidHF* pidObj=new AliAODPidHF();
  // pidObj->SetName("pid4DSatr");
  Int_t mode=1;
  Double_t priors[5]={0.01,0.001,0.3,0.3,0.3};
  pidObj->SetPriors(priors);
  pidObj->SetMatch(mode);
  pidObj->SetSigma(0,2); // TPC
  pidObj->SetSigma(3,3); // TOF
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);
  RDHFDStartoKpipi->SetPidHF(pidObj);

  //activate pileup rejection
  RDHFDStartoKpipi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  RDHFDStartoKpipi->PrintAll();

  TFile* fout=new TFile("DStartoKpipiCuts.root","recreate");   //set this!! 
  fout->cd();
  RDHFDStartoKpipi->Write();
  fout->Close();

}

//macro to make a .root file (for significance maximization) which contains an AliRDHFCutsDStartoKpipi with loose set of cuts  and TParameter with the tighest value of these cuts
// copy form D0 ... NOT TESTED YIET ... to be tested!!

void makeInputAliAnalysisTaskSEDstarSignificanceMaximization(){

  AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi();
  RDHFDStartoKpipi->SetName("loosercuts");
  RDHFDStartoKpipi->SetTitle("Cuts for significance maximization");

  //Centrality selection
  RDHFDStartoKpipi->SetUseCentrality(kFALSE);
  RDHFDStartoKpipi->SetMinCentrality(40);
  RDHFDStartoKpipi->SetMaxCentrality(80);

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetMinNClustersITS(4);
  
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetPtRange(0.1,1.e10);
  
 // soft pion pre-selections 
  // 
  AliESDtrackCuts* esdSoftPicuts=new AliESDtrackCuts();
  esdSoftPicuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdSoftPicuts->SetRequireTPCRefit(kFALSE);
  esdSoftPicuts->SetRequireITSRefit(kFALSE);
  esdSoftPicuts->SetMinNClustersITS(4); // default is 4
  esdSoftPicuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					  AliESDtrackCuts::kAny); //test d0 asimmetry
  esdSoftPicuts->SetPtRange(0.0,1.e10);

  // set pre selections
  RDHFDStartoKpipi->AddTrackCuts(esdTrackCuts);
  RDHFDStartoKpipi->AddTrackCutsSoftPi(esdSoftPicuts);

  const Int_t nvars=16;
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
  ptbins[8]=7.;
  ptbins[9]=8.;
  ptbins[10]=12.;
  ptbins[11]=16.;
  ptbins[12]=24.;
  ptbins[13]=999.;
  RDHFDStartoKpipi->SetPtBins(nptbins+1,ptbins);
  

  Float_t** rdcutsvalmine;
  rdcutsvalmine=new Float_t*[nvars];
  for(Int_t iv=0;iv<nvars;iv++){
    rdcutsvalmine[iv]=new Float_t[nptbins];
  }
  
    //0-0.5
    rdcutsvalmine[0][0]=0.7;
    rdcutsvalmine[1][0]=0.03;
    rdcutsvalmine[2][0]=0.8;
    rdcutsvalmine[3][0]=0.3;
    rdcutsvalmine[4][0]=0.3;
    rdcutsvalmine[5][0]=0.1;
    rdcutsvalmine[6][0]=0.1;
    rdcutsvalmine[7][0]=0.00035;
    rdcutsvalmine[8][0]=0.73;
    rdcutsvalmine[9][0]=0.3;
    rdcutsvalmine[10][0]=0.1;
    rdcutsvalmine[11][0]=0.05;
    rdcutsvalmine[12][0]=100.;
    rdcutsvalmine[13][0]=0.5;
    rdcutsvalmine[14][0]=-1.;
    rdcutsvalmine[15][0]=0.;
    //0.5-1
    rdcutsvalmine[0][1]=0.7;
    rdcutsvalmine[1][1]=0.03;
    rdcutsvalmine[2][1]=0.8;
    rdcutsvalmine[3][1]=0.3;
    rdcutsvalmine[4][1]=0.3;
    rdcutsvalmine[5][1]=0.1;
    rdcutsvalmine[6][1]=0.1;
    rdcutsvalmine[7][1]=-0.00020;
    rdcutsvalmine[8][1]=0.73;
    rdcutsvalmine[9][1]=0.3;
    rdcutsvalmine[10][1]=0.1;
    rdcutsvalmine[11][1]=0.05;
    rdcutsvalmine[12][1]=100.;
    rdcutsvalmine[13][1]=0.5;
    rdcutsvalmine[14][1]=-1.;
    rdcutsvalmine[15][1]=0.;
    //1-2
    rdcutsvalmine[0][2]=0.7;
    rdcutsvalmine[1][2]=0.02;
    rdcutsvalmine[2][2]=0.8;
    rdcutsvalmine[3][2]=0.4;
    rdcutsvalmine[4][2]=0.4;
    rdcutsvalmine[5][2]=0.1;
    rdcutsvalmine[6][2]=0.1;
    rdcutsvalmine[7][2]=-0.00036;
    rdcutsvalmine[8][2]=0.82;
    rdcutsvalmine[9][2]=0.3;
    rdcutsvalmine[10][2]=0.1;
    rdcutsvalmine[11][2]=0.05;
    rdcutsvalmine[12][2]=100.;
    rdcutsvalmine[13][2]=0.5;
    rdcutsvalmine[14][2]=-1.;
    rdcutsvalmine[15][2]=0.;
    //2-3
    rdcutsvalmine[0][3]=0.7;
    rdcutsvalmine[1][3]=0.02;
    rdcutsvalmine[2][3]=0.8;
    rdcutsvalmine[3][3]=0.7;
    rdcutsvalmine[4][3]=0.7;
    rdcutsvalmine[5][3]=0.08;
    rdcutsvalmine[6][3]=0.08;
    rdcutsvalmine[7][3]=-0.00016;
    rdcutsvalmine[8][3]=0.90;
    rdcutsvalmine[9][3]=0.3;
    rdcutsvalmine[10][3]=0.1;
    rdcutsvalmine[11][3]=0.05;
    rdcutsvalmine[12][3]=100.;
    rdcutsvalmine[13][3]=0.5;
    rdcutsvalmine[14][3]=-1.;
    rdcutsvalmine[15][3]=0.;
    //3-4
    rdcutsvalmine[0][4]=0.7;
    rdcutsvalmine[1][4]=0.05;
    rdcutsvalmine[2][4]=0.8;
    rdcutsvalmine[3][4]=1.;
    rdcutsvalmine[4][4]=1.;
    rdcutsvalmine[5][4]=0.042;
    rdcutsvalmine[6][4]=0.056;
    rdcutsvalmine[7][4]=-0.000065;
    rdcutsvalmine[8][4]=0.9;
    rdcutsvalmine[9][4]=0.3;
    rdcutsvalmine[10][4]=0.1;
    rdcutsvalmine[11][4]=0.05;
    rdcutsvalmine[12][4]=100.;
    rdcutsvalmine[13][4]=0.5;
    rdcutsvalmine[14][4]=0.;
    rdcutsvalmine[15][4]=0.;
    //4-5
    rdcutsvalmine[0][5]=0.7;
    rdcutsvalmine[1][5]=0.08;
    rdcutsvalmine[2][5]=0.9;
    rdcutsvalmine[3][5]=1.2;
    rdcutsvalmine[4][5]=1.2;
    rdcutsvalmine[5][5]=0.07;
    rdcutsvalmine[6][5]=0.07;
    rdcutsvalmine[7][5]=0.0001;
    rdcutsvalmine[8][5]=0.9;
    rdcutsvalmine[9][5]=0.3;
    rdcutsvalmine[10][5]=0.1;
    rdcutsvalmine[11][5]=0.05;
    rdcutsvalmine[12][5]=100.;
    rdcutsvalmine[13][5]=0.5;
    rdcutsvalmine[14][5]=-1.;
    rdcutsvalmine[15][5]=0.;
    //5-6
    rdcutsvalmine[0][6]=0.7;
    rdcutsvalmine[1][6]=0.1;
    rdcutsvalmine[2][6]=1.0;
    rdcutsvalmine[3][6]=1.;
    rdcutsvalmine[4][6]=1.;
    rdcutsvalmine[5][6]=0.08;
    rdcutsvalmine[6][6]=0.08;
    rdcutsvalmine[7][6]=0.0005;
    rdcutsvalmine[8][6]=0.8;
    rdcutsvalmine[9][6]=0.3;
    rdcutsvalmine[10][6]=0.1;
    rdcutsvalmine[11][6]=0.05;
    rdcutsvalmine[12][6]=100000.;
    rdcutsvalmine[13][6]=0.5;
    rdcutsvalmine[14][6]=-1.;
    rdcutsvalmine[15][6]=0.;
    //6-7
    rdcutsvalmine[0][7]=0.7;
    rdcutsvalmine[1][7]=0.1;
    rdcutsvalmine[2][7]=1.0;
    rdcutsvalmine[3][7]=1.;
    rdcutsvalmine[4][7]=1.;
    rdcutsvalmine[5][7]=0.1;
    rdcutsvalmine[6][7]=0.1;
    rdcutsvalmine[7][7]=0.001;
    rdcutsvalmine[8][7]=0.7;
    rdcutsvalmine[9][7]=0.3;
    rdcutsvalmine[10][7]=0.1;
    rdcutsvalmine[11][7]=0.05;
    rdcutsvalmine[12][7]=100.;
    rdcutsvalmine[13][7]=0.5;
    rdcutsvalmine[14][7]=-1.;
    rdcutsvalmine[15][7]=0.;
    //7-8
    rdcutsvalmine[0][8]=0.7;
    rdcutsvalmine[1][8]=0.1;
    rdcutsvalmine[2][8]=1.0;
    rdcutsvalmine[3][8]=1.;
    rdcutsvalmine[4][8]=1.;
    rdcutsvalmine[5][8]=0.1;
    rdcutsvalmine[6][8]=0.1;
    rdcutsvalmine[7][8]=0.001;
    rdcutsvalmine[8][8]=0.7;
    rdcutsvalmine[9][8]=0.3;
    rdcutsvalmine[10][8]=0.1;
    rdcutsvalmine[11][8]=0.05;
    rdcutsvalmine[12][8]=100.;
    rdcutsvalmine[13][8]=0.5;
    rdcutsvalmine[14][8]=-1.;
    rdcutsvalmine[15][8]=0.;
    //8-12
    rdcutsvalmine[0][9]=0.7;
    rdcutsvalmine[1][9]=0.1;
    rdcutsvalmine[2][9]=1.0;
    rdcutsvalmine[3][9]=1.;
    rdcutsvalmine[4][9]=1.;
    rdcutsvalmine[5][9]=0.1;
    rdcutsvalmine[6][9]=0.1;
    rdcutsvalmine[7][9]=0.006;
    rdcutsvalmine[8][9]=0.7;
    rdcutsvalmine[9][9]=0.3;
    rdcutsvalmine[10][9]=0.1;
    rdcutsvalmine[11][9]=0.05;
    rdcutsvalmine[12][9]=100.;
    rdcutsvalmine[13][9]=0.5;
    rdcutsvalmine[14][9]=-1.;
    rdcutsvalmine[15][9]=0.;
    //12-16
    rdcutsvalmine[0][10]=0.7;
    rdcutsvalmine[1][10]=0.1;
    rdcutsvalmine[2][10]=1.0;
    rdcutsvalmine[3][10]=0.3;
    rdcutsvalmine[4][10]=0.3;
    rdcutsvalmine[5][10]=0.15;
    rdcutsvalmine[6][10]=0.15;
    rdcutsvalmine[7][10]=0.01;
    rdcutsvalmine[8][10]=0.7;
    rdcutsvalmine[9][10]=0.3;
    rdcutsvalmine[10][10]=0.1;
    rdcutsvalmine[11][10]=0.05;
    rdcutsvalmine[12][10]=100.;
    rdcutsvalmine[13][10]=0.5;
    rdcutsvalmine[14][10]=-1.;
    rdcutsvalmine[15][10]=0.;
    //16-24
    rdcutsvalmine[0][11]=0.7;
    rdcutsvalmine[1][11]=0.2;
    rdcutsvalmine[2][11]=1.0;
    rdcutsvalmine[3][11]=.3;
    rdcutsvalmine[4][11]=.3;
    rdcutsvalmine[5][11]=0.15;
    rdcutsvalmine[6][11]=0.15;
    rdcutsvalmine[7][11]=0.01;
    rdcutsvalmine[8][11]=0.7;
    rdcutsvalmine[9][11]=0.3;
    rdcutsvalmine[10][11]=0.1;
    rdcutsvalmine[11][11]=0.05;
    rdcutsvalmine[12][11]=100.;
    rdcutsvalmine[13][11]=0.5;
    rdcutsvalmine[14][11]=-1.;
    rdcutsvalmine[15][11]=0.;
    //>24
    rdcutsvalmine[0][12]=0.7;
    rdcutsvalmine[1][12]=0.6;
    rdcutsvalmine[2][12]=1.0;
    rdcutsvalmine[3][12]=.4;
    rdcutsvalmine[4][12]=.4;
    rdcutsvalmine[5][12]=0.5;
    rdcutsvalmine[6][12]=0.5;
    rdcutsvalmine[7][12]=0.1;
    rdcutsvalmine[8][12]=0.7;
    rdcutsvalmine[9][12]=0.3;
    rdcutsvalmine[10][12]=0.1;
    rdcutsvalmine[11][12]=0.05;
    rdcutsvalmine[12][12]=100.;
    rdcutsvalmine[13][12]=0.5;
    rdcutsvalmine[14][12]=-1.;
    rdcutsvalmine[15][12]=0.;

  RDHFDStartoKpipi->SetCuts(nvars,nptbins,rdcutsvalmine);

  Int_t nvarsforopt=RDHFDStartoKpipi->GetNVarsForOpt();
  Int_t dim=2; //set this!!
  Bool_t *boolforopt;
  boolforopt=new Bool_t[nvars];
  if(dim>nvarsforopt){
    cout<<"Number of variables for optimization has probably changed, check and edit accordingly"<<endl;
    return;
  } else {
    if(dim==nvarsforopt){
      boolforopt=RDHFDStartoKpipi->GetVarsForOpt();
    }else{
      TString *names;
      names=new TString[nvars];
      TString answer="";
      Int_t checktrue=0;
      names=RDHFDStartoKpipi->GetVarNames();
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
      RDHFDStartoKpipi->SetVarsForOpt(dim,boolforopt);
    }
  }


  Float_t tighterval[dim][nptbins];
  //dca
  //costhetastar
  //d0d0 <-this 
  //costhetapoint <-this 

  
  //number of steps for each variable is 4 now
  //set this!!
  tighterval[0][0]=-0.0007;
  tighterval[1][0]=0.99;

  tighterval[0][1]=-0.0006;
  tighterval[1][1]=0.99;
 
  tighterval[0][2]=-0.0004;
  tighterval[1][2]=0.99;
 
  tighterval[0][3]=-0.00035;
  tighterval[1][3]=0.98;

  tighterval[0][4]=-0.0003;
  tighterval[1][4]=0.99;


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
 
  TFile* fout=new TFile("cuts4SignifMaxim.root","recreate");   //set this!! 
  fout->cd();
  RDHFDStartoKpipi->Write();
  max.Write();
  fout->Close();
 
}

