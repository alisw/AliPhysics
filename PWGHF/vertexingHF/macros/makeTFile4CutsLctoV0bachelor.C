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

enum ECollisionSystem{kpp7TeV,kpPb,kPbPb};
enum EDataKind{kData,kLHC10f7a,kLHC11b2,kLHC15a2a,kLHC13d3};
enum EReconstructionPass{kPass2,kPass4};

// Configuration
const Int_t nptbins=6;
Int_t fCollisionSystem=kpp7TeV;
Int_t fDataType=kData;
Int_t fPass=kPass2;
Bool_t fChangeK0ScutOnMC=kTRUE;
Bool_t fChangeAlsoLambdaVetoOnMC=kTRUE;

const Float_t kmInvPipPim_pp7TeV_pass2[nptbins] = {0.00432330, 0.00454560, 0.00481579, 0.00505379, 0.00542286, 0.00602017};
const Float_t kmInvPipPim_pp7TeV_pass4[nptbins] = {0.00398116, 0.00412519, 0.00432494, 0.00454682, 0.00486837, 0.00537445};
const Float_t kmInvPipPim_pPb[nptbins]          = {0.00383156, 0.00394778, 0.00413934, 0.00438431, 0.00480869, 0.00544472};
const Float_t kmInvPipPim_LHC10f7a[nptbins]     = {0.00398520, 0.00419456, 0.00445022, 0.00462474, 0.00491519, 0.00544176};
const Float_t kmInvPipPim_LHC11b2[nptbins]      = {0.00385273, 0.00394676, 0.00412335, 0.00420992, 0.00445422, 0.00493645};
const Float_t kmInvPipPim_LHC15a2a[nptbins]     = {0.00376002, 0.00386896, 0.00399314, 0.00412475, 0.00435170, 0.00466817};
const Float_t kmInvPipPim_LHC13d3[nptbins]      = {0.00358776, 0.00370684, 0.00383341, 0.00402011, 0.00439230, 0.00489222};

const Float_t kmInvPPi_pp7TeV_pass2[nptbins]    = {0.00203633, 0.00196967, 0.00196595, 0.00199441, 0.00202691, 0.00211474};
const Float_t kmInvPPi_pp7TeV_pass4[nptbins]    = {0.00194356, 0.00185442, 0.00184579, 0.00187065, 0.00193578, 0.00210853};
const Float_t kmInvPPi_pPb[nptbins]             = {0.00214349, 0.00199650, 0.00192680, 0.00192083, 0.00199756, 0.00213533};
const Float_t kmInvPPi_LHC10f7a[nptbins]        = {0.00190966, 0.00182307, 0.00177527, 0.00183747, 0.00181431, 0.00204811};
const Float_t kmInvPPi_LHC11b2[nptbins]         = {0.00185630, 0.00180454, 0.00179920, 0.00182468, 0.00189292, 0.00204497};
const Float_t kmInvPPi_LHC15a2a[nptbins]        = {0.00187642, 0.00175389, 0.00171987, 0.00171276, 0.00180925, 0.00194709};
const Float_t kmInvPPi_LHC13d3[nptbins]         = {0.00205320, 0.00191146, 0.00179558, 0.00182894, 0.00185072, 0.00208818};

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
  RDHFLctoV0An->SetUseTrackSelectionWithFilterBits(kFALSE);//(kTRUE);
  if (fCollisionSystem==kpPb) {
    // //
    // // Trigger selection
    // //
    RDHFLctoV0An->SetTriggerClass("");     
    if (fDataType==kData) {
      RDHFLctoV0An->SetTriggerMask(AliVEvent::kINT7);
    } else if (fDataType==kLHC13d3) {
      RDHFLctoV0An->SetTriggerMask(AliVEvent::kINT7 | AliVEvent::kMB);
    }
    // multi vertexer pileup rejection
    RDHFLctoV0An-> SetOptPileup(AliRDHFCuts::kRejectMVPileupEvent); 
  }
  RDHFLctoV0An->SetLowPtCut(1.0); // default value 1.0 GeV/c
  if (fCollisionSystem==kpp7TeV) {
    RDHFLctoV0An->SetHighPtCut(3.0); // default value 2.5 GeV/c
    RDHFLctoV0An->SetPidSelectionFlag(4); // 0 -> TOF AND TPC
  } else if (fCollisionSystem==kpPb) {
    RDHFLctoV0An->SetHighPtCut(999.); // default value 2.5 GeV/c
    RDHFLctoV0An->SetPidSelectionFlag(2); // 0 -> TOF AND TPC
                                          // 1 -> if (TOF) TOF else TPC w veto
                                          // 2 -> if (p<1) TPC@2s else if (1<=p<2.5) {if (TOF) TOF@3s AND TPC@3s} else (p>=2.5) {if (TOF) -2s<TOF<3s AND TPC@3s}
                                          // 3 -> if (p<1) TPC@2s else if (1<=p<2.5) {if (TOF) TOF@3s AND TPC@3s} else (p>=2.5) {if (TOF) -2s<TOF<3s AND -3s<TPC<2s}
                                          // 4 -> if (p<1) TPC@2s else if (1<=p<2.5) {if (TOF) TOF@3s} else if (p>=2.5) {if (TOF) -2s<TOF<3s}
                                          // 5 -> if (p<1) TPC@2s else if (1<=p<2.5) {if (TOF) TOF@3s} else if (p>=2.5) {if (TOF) -2s<TOF<3s else TPC@3s}
                                          // 6 -> if (p<1) TPC@2s else if (1<=p<2.5) {if (TOF) TOF@3s} else if (p>=2.5) {if (TOF) -2s<TOF<3s else -3s<TPC<2s}
                                          // 7 -> if (p<1) TPC@2s else if (1<=p<2.5) {if (TOF) TOF@3s else TPC@3s} else (p>=2.5) {if (TOF) -2s<TOF<3s else TPC@3s}
                                          // 8 -> if (p<1) TPC@2s else if (1<=p<2.5) {if (TOF) TOF@3s else TPC@3s} else (p>=2.5) {if (TOF) -2s<TOF<3s else -2<TPC<3s}
                                          // 9 -> combined: TOF, TPC
  } else if (fCollisionSystem==kPbPb) {
    RDHFLctoV0An->SetHighPtCut(999.); // default value 2.5 GeV/c
    RDHFLctoV0An->SetPidSelectionFlag(10); // 10: -2<tof<3 for pt>3 GeV/c, -3<tpc<2 for pt>4 GeV/c 
  }
  RDHFLctoV0An->SetNPtBins(nptbins);

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]= 2.;
  ptbins[1]= 3.;
  ptbins[2]= 4.;
  ptbins[3]= 5.;
  ptbins[4]= 6.;
  ptbins[5]= 8.;
  ptbins[6]= 12.;
  RDHFLctoV0An->SetPtBins(nptbins+1,ptbins);

  const Int_t nvars=21;

  Float_t** anacutsval;
  anacutsval=new Float_t*[nvars];
  for(Int_t ic=0;ic<nvars;ic++){anacutsval[ic]=new Float_t[nptbins];}
  for(Int_t ipt2=0;ipt2<nptbins;ipt2++){
   anacutsval[0][ipt2]=0.25;   // inv. mass if K0S [GeV/c2]
   anacutsval[1][ipt2]=1.00;   // inv. mass if Lambda [GeV/c2]
   anacutsval[3][ipt2]=0.05;   // inv. mass V0 if Lambda [GeV/c2] ---> WE ARE SEARCHING Lc -> p+K0S, so cut on m(Lambda) has to be leave as it was at filtering level!!!
   anacutsval[7][ipt2]=1000.;  // dca cascade cut [cm]
   anacutsval[8][ipt2]=1.5;    // dca V0 cut [nSigma] // it's 1.5 x offline V0s
   anacutsval[9][ipt2]=0.99;   // cosPA V0 cut // it's 0.90 x offline V0s at reconstruction level, 0.99 at filtering level
   anacutsval[12][ipt2]=0.;    // mass K0S veto [GeV/c2]
   anacutsval[13][ipt2]=0.005; // mass Lambda/LambdaBar veto [GeV/c2]
   anacutsval[16][ipt2]=9999; // max cos proton emission angle in Lc CMS max
   anacutsval[17][ipt2]=-9999; // min cos proton emission angle in Lc CMS max
   anacutsval[18][ipt2]=-9999; // Re-signed d0 min [cm]
   anacutsval[19][ipt2]=-9999; // V0 qT/|alpha| min
   anacutsval[20][ipt2]=0.;    // V0 type cut
  }

  if (fCollisionSystem==kpp7TeV) {
    if (fPass==kPass2) {
      anacutsval[2][0]=0.00500; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][1]=0.00812; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][2]=0.00858; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][3]=0.00900; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][4]=0.01000; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][5]=0.01000; // inv. mass V0 if K0S [GeV/c2]
      if (fDataType==kLHC10f7a) {
	for (Int_t ii=0; ii<nptbins; ii++) {
	  if (fChangeK0ScutOnMC) {
	    anacutsval[2][ii]=anacutsval[2][ii]/kmInvPipPim_pp7TeV_pass2[ii]*kmInvPipPim_LHC10f7a[ii]; // inv. mass V0 if K0S [GeV/c2]
	    if (fChangeAlsoLambdaVetoOnMC) {
	      anacutsval[13][ii]=anacutsval[13][ii]/kmInvPPi_pp7TeV_pass2[ii]*kmInvPPi_LHC10f7a[ii]; // inv. mass Lambda/LambdaBar [GeV/c2]
	    }
	  }
	}
      } else if (fDataType==kLHC11b2) {
	for (Int_t ii=0; ii<nptbins; ii++) {
	  if (fChangeK0ScutOnMC) {
	    anacutsval[2][ii]=anacutsval[2][ii]/kmInvPipPim_pp7TeV_pass2[ii]*kmInvPipPim_LHC11b2[ii]; // inv. mass V0 if K0S [GeV/c2]
	    if (fChangeAlsoLambdaVetoOnMC) {
	      anacutsval[13][ii]=anacutsval[13][ii]/kmInvPPi_pp7TeV_pass2[ii]*kmInvPPi_LHC11b2[ii]; // inv. mass Lambda/LambdaBar [GeV/c2]
	    }
	  }
	}
      }

    } else if (fPass==kPass4) { // as for pass2 -> not optimized for pass4 data!
      anacutsval[2][0]=0.00500; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][1]=0.00812; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][2]=0.00858; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][3]=0.00900; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][4]=0.01000; // inv. mass V0 if K0S [GeV/c2]
      anacutsval[2][5]=0.01000; // inv. mass V0 if K0S [GeV/c2]
      // just to start to see pass4
      for (Int_t ii=0; ii<nptbins; ii++) {
	anacutsval[2][ii]=anacutsval[2][ii]/kmInvPipPim_pp7TeV_pass2[ii]*kmInvPipPim_pp7TeV_pass4[ii]; // inv. mass V0 if K0S [GeV/c2]
      }
      anacutsval[13][0]=0.005; // inv. mass V0 if Lambda/LambdaBar [GeV/c2]
      anacutsval[13][1]=0.005; // inv. mass V0 if Lambda/LambdaBar [GeV/c2]
      anacutsval[13][2]=0.005; // inv. mass V0 if Lambda/LambdaBar [GeV/c2]
      anacutsval[13][3]=0.005; // inv. mass V0 if Lambda/LambdaBar [GeV/c2]
      anacutsval[13][4]=0.005; // inv. mass V0 if Lambda/LambdaBar [GeV/c2]
      anacutsval[13][5]=0.005; // inv. mass V0 if Lambda/LambdaBar [GeV/c2]
      // just to start to see pass4
      for (Int_t ii=0; ii<nptbins; ii++) {
	anacutsval[13][ii]=anacutsval[13][ii]/kmInvPPi_pp7TeV_pass2[ii]*kmInvPPi_pp7TeV_pass4[ii]; // inv. mass V0 if Lambda/LambdaBar [GeV/c2]
      }
      if (fDataType==kLHC15a2a) {
	for (Int_t ii=0; ii<nptbins; ii++) {
	  if (fChangeK0ScutOnMC) {
	    anacutsval[2][ii]=anacutsval[2][ii]/kmInvPipPim_pp7TeV_pass4[ii]*kmInvPipPim_LHC15a2a[ii]; // inv. mass V0 if K0S [GeV/c2]
	    if (fChangeAlsoLambdaVetoOnMC) {
	      anacutsval[13][ii]=anacutsval[13][ii]/kmInvPPi_pp7TeV_pass4[ii]*kmInvPPi_LHC15a2a[ii]; // inv. mass Lambda/LambdaBar [GeV/c2]
	    }
	  }
	}
      }
    }

  } else if (fCollisionSystem==kpPb) {
    anacutsval[2][0]=0.006; // inv. mass V0 if K0S [GeV/c2]
    anacutsval[2][1]=0.007; // inv. mass V0 if K0S [GeV/c2]
    anacutsval[2][2]=0.008; // inv. mass V0 if K0S [GeV/c2]
    anacutsval[2][3]=0.009; // inv. mass V0 if K0S [GeV/c2]
    anacutsval[2][4]=0.010; // inv. mass V0 if K0S [GeV/c2]
    anacutsval[2][5]=0.011; // inv. mass V0 if K0S [GeV/c2]
    if (fDataType==kLHC13d3) {
      for (Int_t ii=0; ii<nptbins; ii++) {
	if (fChangeK0ScutOnMC) {
	  anacutsval[2][ii]=anacutsval[2][ii]/kmInvPipPim_pPb[ii]*kmInvPipPim_LHC13d3[ii]; // inv. mass V0 if K0S [GeV/c2]
	  if (fChangeAlsoLambdaVetoOnMC) {
	    anacutsval[13][ii]=anacutsval[13][ii]/kmInvPPi_pPb[ii]*kmInvPPi_LHC13d3[ii]; // inv. mass Lambda/LambdaBar [GeV/c2]
	  }
	}
      }
    }
  }

  if (fCollisionSystem==kpp7TeV) {

    anacutsval[4][0]=0.6;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][1]=0.6;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][2]=0.7;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][3]=0.7;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][4]=0.9;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][5]=1.1;    // pT min bachelor track [GeV/c] // AOD by construction

    anacutsval[5][0]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][1]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][2]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][3]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][4]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][5]=0.3;    // pT min V0-positive track [GeV/c]

    anacutsval[6][0]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][1]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][2]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][3]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][4]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][5]=0.3;    // pT min V0-negative track [GeV/c]

    anacutsval[10][0]=0.05;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][1]=0.05;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][2]=0.10;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][3]=0.10;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][4]=0.10;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][5]=0.10;  // d0 max bachelor wrt PV [cm]

    anacutsval[11][0]=0.05;  // d0 max V0 wrt PV [cm]
    anacutsval[11][1]=0.05;  // d0 max V0 wrt PV [cm]
    anacutsval[11][2]=0.09;  // d0 max V0 wrt PV [cm]
    anacutsval[11][3]=0.09;  // d0 max V0 wrt PV [cm]
    anacutsval[11][4]=999.;  // d0 max V0 wrt PV [cm]
    anacutsval[11][5]=999.;  // d0 max V0 wrt PV [cm]

    anacutsval[14][0]=0.100; // mass Gamma veto [GeV/c2]
    anacutsval[14][1]=0.100; // mass Gamma veto [GeV/c2]
    anacutsval[14][2]=0.300; // mass Gamma veto [GeV/c2]
    anacutsval[14][3]=0.300; // mass Gamma veto [GeV/c2]
    anacutsval[14][4]=0.300; // mass Gamma veto [GeV/c2]
    anacutsval[14][5]=0.300; // mass Gamma veto [GeV/c2]

    anacutsval[15][0]=0.5; // pT min V0 track [GeV/c]
    anacutsval[15][1]=0.6; // pT min V0 track [GeV/c]
    anacutsval[15][2]=0.7; // pT min V0 track [GeV/c]
    anacutsval[15][3]=1.0; // pT min V0 track [GeV/c]
    anacutsval[15][4]=1.1; // pT min V0 track [GeV/c]
    anacutsval[15][5]=1.2; // pT min V0 track [GeV/c]

  } else if (fCollisionSystem==kpPb) {

    anacutsval[4][0]=0.5;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][1]=0.6;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][2]=0.7;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][3]=0.9;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[4][4]=0.9;    // pT min bachelor track [GeV/c] // AOD by construction
    anacutsval[5][5]=0.9;    // pT min bachelor track [GeV/c] // AOD by construction

    anacutsval[5][0]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][1]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][2]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][3]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][4]=0.2;    // pT min V0-positive track [GeV/c]
    anacutsval[5][5]=0.3;    // pT min V0-positive track [GeV/c]

    anacutsval[6][0]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][1]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][2]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][3]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][4]=0.2;    // pT min V0-negative track [GeV/c]
    anacutsval[6][5]=0.3;    // pT min V0-negative track [GeV/c]

    anacutsval[10][0]=0.05;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][1]=0.05;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][2]=0.08;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][3]=0.10;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][4]=0.10;  // d0 max bachelor wrt PV [cm]
    anacutsval[10][5]=0.10;  // d0 max bachelor wrt PV [cm]

    for (Int_t ipt2=0; ipt2<nptbins; ipt2++) anacutsval[11][ipt2]=999.;  // d0 max V0 wrt PV [cm]

    anacutsval[14][0]=0.100; // mass Gamma veto [GeV/c2]
    anacutsval[14][1]=0.100; // mass Gamma veto [GeV/c2]
    anacutsval[14][2]=0.100; // mass Gamma veto [GeV/c2]
    anacutsval[14][3]=0.100; // mass Gamma veto [GeV/c2]
    anacutsval[14][4]=0.100; // mass Gamma veto [GeV/c2]
    anacutsval[14][5]=0.100; // mass Gamma veto [GeV/c2]

    anacutsval[15][0]=0.8; // pT min V0 track [GeV/c]
    anacutsval[15][1]=0.8; // pT min V0 track [GeV/c]
    anacutsval[15][2]=0.8; // pT min V0 track [GeV/c]
    anacutsval[15][3]=1.2; // pT min V0 track [GeV/c]
    anacutsval[15][4]=1.0; // pT min V0 track [GeV/c]
    anacutsval[15][5]=1.2; // pT min V0 track [GeV/c]

  }

  RDHFLctoV0An->SetCuts(nvars,nptbins,anacutsval);

  if (RDHFLctoV0An->GetPidSelectionFlag()==9) {
    Float_t minCombProb[nptbins];
    for (Int_t iBin=0; iBin<nptbins; iBin++) minCombProb[iBin]=0.;
    RDHFLctoV0An->SetMinCombinedProbability(nptbins,minCombProb);
  }

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
  if (RDHFLctoV0An->GetPidSelectionFlag()==9) {
    pidObjBachelor->SetUseCombined();
    pidObjBachelor->SetUseDefaultPriors(kTRUE);
    pidObjBachelor->SetCombDetectors(AliAODPidHF::kTPCTOF);
  }

  RDHFLctoV0An->SetPidHF(pidObjBachelor);

  // uncomment these lines for Baysian PID:
  // Double_t threshold=0.3;
  // SetupCombinedPID(RDHFLctoV0An  ,threshold);
  // RDHFLctoV0An  ->SetPIDStrategy(AliRDHFCutsLctoV0::kCombined);
  //

  //RDHFLc->SetRecoKF(); //set this if you want to recompute the secondary vertex with the KF package
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
  TFile* fout;
  if (fDataType==kData) {
    fout=new TFile("Lc2pK0SCuts.root","RECREATE");
  } else {
    fout=new TFile("Lc2pK0SCuts.root","RECREATE");
  }
  fout->cd();
  RDHFLctoV0An->Write();
  fout->Close();
  delete fout;

  delete pidObjBachelor;
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

  const Int_t nvars=14;//17-3

  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=2.;
  ptbins[1]=3.;
  ptbins[2]=4.;
  ptbins[3]=5.;
  ptbins[4]=6.;
  ptbins[5]=8.;
  ptbins[6]=12.;
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
    cutsMatrixLctoV0Stand[ 0][ipt2] =0.0075;
    cutsMatrixLctoV0Stand[ 1][ipt2] =0.0500;
    cutsMatrixLctoV0Stand[ 2][ipt2] =0.4;
    cutsMatrixLctoV0Stand[ 3][ipt2] =0.2;
    cutsMatrixLctoV0Stand[ 4][ipt2] =0.2;
    cutsMatrixLctoV0Stand[ 5][ipt2] =1000.;
    cutsMatrixLctoV0Stand[ 6][ipt2] =1.5;
    cutsMatrixLctoV0Stand[ 7][ipt2] =0.99;
    cutsMatrixLctoV0Stand[ 8][ipt2] =0.05;
    cutsMatrixLctoV0Stand[ 9][ipt2] =0.1;
    cutsMatrixLctoV0Stand[10][ipt2] =0.0;
    cutsMatrixLctoV0Stand[11][ipt2] =0.005;
    cutsMatrixLctoV0Stand[12][ipt2] =0.100;
  }
  cutsMatrixLctoV0Stand[13][0]=0.0; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[13][1]=0.6; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[13][2]=0.8; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[13][3]=0.8; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[13][4]=0.8; // pT min V0 track [GeV/c]
  cutsMatrixLctoV0Stand[13][5]=1.0; // pT min V0 track [GeV/c]


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
    tighterval[ 0][ipt] =0.075;  // inv. mass V0 if K0S [GeV]
    tighterval[ 1][ipt] =0.040;  // inv. mass V0 if Lambda [GeV]
    tighterval[ 2][ipt] =0.4;    // pT min bachelor track [GeV/c]
    tighterval[ 3][ipt] =0.2;    // pT min V0-positive track [GeV/c]
    tighterval[ 4][ipt] =0.2;    // pT min V0-negative track [GeV/c]
    tighterval[ 5][ipt] =100.;   // dca cascade cut (cm)
    tighterval[ 6][ipt] =1.5;    // dca v0 cut
    tighterval[ 7][ipt] =0.99;   // cosPA v0 cut
    tighterval[ 8][ipt] =0.05;   // d0 max bachelor wrt PV [cm]
    tighterval[ 9][ipt] =0.1;    // d0 max V0 wrt PV [cm]
    tighterval[10][ipt] =0.0;   // mass K0S veto [GeV/c2]
    tighterval[11][ipt] =0.005; // mass Lambda/LambdaBar veto [GeV/c2]
    tighterval[12][ipt] =0.100; // mass Gamma veto [GeV/c2]
  }
  tighterval[13][0]=0.0; // pT min V0 track [GeV/c]
  tighterval[13][1]=0.6; // pT min V0 track [GeV/c]
  tighterval[13][2]=0.8; // pT min V0 track [GeV/c]
  tighterval[13][3]=0.8; // pT min V0 track [GeV/c]
  tighterval[13][4]=0.8; // pT min V0 track [GeV/c]
  tighterval[13][5]=0.8; // pT min V0 track [GeV/c]

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
