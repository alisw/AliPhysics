/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

/// \ingroup macros
/// \file MUONRecoCheck.C
/// \brief Utility macro to check the muon reconstruction. 
///
/// Reconstructed tracks are compared to reference tracks. The reference tracks 
/// are built from AliTrackReference for the hit in chamber (0..9) and from 
/// kinematics (TreeK) for the vertex parameters.  
///
/// \author Jean-Pierre Cussonneau, Philippe Pillot, Subatech  

// ROOT includes
#include <Riostream.h>
#include "TMath.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGeoManager.h"

// STEER includes
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliLog.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONTrackParam.h"
#include "AliMUONRecoParam.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVCluster.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONTriggerTrack.h"

Double_t langaufun(Double_t *x, Double_t *par);
void     FitGausResVsMom(TH2* h, Int_t nBins, const Double_t mean0, const Double_t sigma0, const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma);
void     FitPDCAVsMom(TH2* h, Int_t nBins, const char* fitting, TGraphAsymmErrors* gSigma);
TCanvas* DrawVsAng(const char* name, const char* title, TH1* h1, TH2* h2);
TCanvas* DrawVsPos(const char* name, const char* title, TH2* h1, TH2* h2, TH2* h3);
TCanvas* DrawResMomVsMom(const char* name, const char* title, TH2* h, Int_t nBins, TF1* f2 = 0x0, const char* fitting = "");

//------------------------------------------------------------------------------------
void MUONRecoCheck (Int_t nEvent = -1, const char* pathSim="./generated/", const char* esdFileName="AliESDs.root",
		    const char* ocdbPath = "local://$ALICE_ROOT/OCDB", Int_t absorberRegion = -1)
{
  /// Associate the reconstructed tracks with the simulated ones and check the quality of the reconstruction
  /// (tracking/trigger efficiency; momentum, slope,... resolutions at first cluster and at vertex; cluster resolution).
  /// You can limit the calculation of track resolution at vertex to the tracks crossing the absorber in a given region
  /// with the flag "absorberRegion": -1=all, 1=[2,3]deg, 2=[3,10]deg.
  
  Double_t aAbsLimits[2];
  if (absorberRegion > -1) {
    if (absorberRegion == 1) {
      aAbsLimits[0] = 2.;
      aAbsLimits[1] = 3.;
    } else if (absorberRegion == 2) {
      aAbsLimits[0] = 3.;
      aAbsLimits[1] = 10.;
    } else {
      cout<<"Unknown absorber region. Valid choices are: -1=all, 1=[2,3]deg, 2=[3,10]deg"<<endl;
      return;
    }
  } else {
    aAbsLimits[0] = 0.;
    aAbsLimits[1] = 90.;
  }
  
  AliLog::SetClassDebugLevel("AliMCEvent",-1);
  
  // ###################################### define histograms ###################################### //
  // File for histograms and histogram booking
  TFile *histoFile = new TFile("MUONRecoCheck.root", "RECREATE");
  
  TH1F *hReconstructible = new TH1F("hReconstructible"," Nb of reconstructible tracks / evt",15,-0.5,14.5);
  TH1F *hReco = new TH1F("hReco"," Nb of reconstructed tracks / evt",15,-0.5,14.5);
  TH1F *hNClusterComp = new TH1F("hNClusterComp"," Nb of compatible clusters / track ",15,-0.5,14.5);
  TH1F *hTrackRefID = new TH1F("hTrackRefID"," track reference ID ",100,-0.5,99.5);
  TH1F *hTriggerable = new TH1F("hTriggerable"," Nb of triggerable tracks / evt",15,-0.5,14.5);
  TH1F *hTriggered = new TH1F("hTriggered"," Nb of triggered tracks / evt",15,-0.5,14.5);
  
  // momentum resolution at vertex
  histoFile->mkdir("momentumAtVertex","momentumAtVertex");
  histoFile->cd("momentumAtVertex");
  
  const Int_t pNBins = 30;
  const Double_t pEdges[2] = {0., 300.};
  const Int_t deltaPAtVtxNBins = 250;
  const Double_t deltaPAtVtxEdges[2] = {-35., 15.};
  
  TH1F *hResMomVertex = new TH1F("hResMomVertex"," delta P at vertex;#Delta_{p} (GeV/c)",deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TH2D *hResMomVertexVsMom = new TH2D("hResMomVertexVsMom","#Delta_{p} at vertex versus p;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMom_2_3_Deg = new TH2D("hResMomVertexVsMom_2_3_Deg","#Delta_{p} at vertex versus p for tracks between 2 and 3 degrees at absorber end;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMom_3_10_Deg = new TH2D("hResMomVertexVsMom_3_10_Deg","#Delta_{p} at vertex versus p for tracks between 3 and 10 degrees at absorber end;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMom_0_2_DegMC = new TH2D("hResMomVertexVsMom_0_2_DegMC","#Delta_{p} at vertex versus p for tracks with MC angle below 2 degrees;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtVtxNBins/10,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TH2D *hResMomVertexVsPosAbsEnd_0_2_DegMC = new TH2D("hResMomVertexVsPosAbsEnd_0_2_DegMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{p} (GeV/c)",1000,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsPosAbsEnd_2_3_DegMC = new TH2D("hResMomVertexVsPosAbsEnd_2_3_DegMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{p} (GeV/c)",1000,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsPosAbsEnd_3_10_DegMC = new TH2D("hResMomVertexVsPosAbsEnd_3_10_DegMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{p} (GeV/c)",1000,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TH2D *hResMomVertexVsAngle = new TH2D("hResMomVertexVsAngle","#Delta_{p} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{p} (GeV/c)",10,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH2D *hResMomVertexVsMCAngle = new TH2D("hResMomVertexVsMCAngle","#Delta_{p} at vertex versus MC angle;MC angle (Deg);#Delta_{p} (GeV/c)",10,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  TH3D *hResMomVertexVsAngleVsMom = new TH3D("hResMomVertexVsAngleVsMom","#Delta_{p} at vertex versus track position at absorber end converted to degrees versus momentum;p (GeV/c);angle (Deg);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],100,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  
  TGraphAsymmErrors* gMeanResMomVertexVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResMomVertexVsMom->SetName("gMeanResMomVertexVsMom");
  gMeanResMomVertexVsMom->SetTitle("<#Delta_{p}> at vertex versus p;p (GeV/c);<#Delta_{p}> (GeV/c)");
  TGraphAsymmErrors* gMostProbResMomVertexVsMom = new TGraphAsymmErrors(pNBins);
  gMostProbResMomVertexVsMom->SetName("gMostProbResMomVertexVsMom");
  gMostProbResMomVertexVsMom->SetTitle("Most probable #Delta_{p} at vertex versus p;p (GeV/c);Most prob. #Delta_{p} (GeV/c)");
  TGraphAsymmErrors* gSigmaResMomVertexVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResMomVertexVsMom->SetName("gSigmaResMomVertexVsMom");
  gSigmaResMomVertexVsMom->SetTitle("#sigma_{p}/p at vertex versus p;p (GeV/c);#sigma_{p}/p (%)");
  
  // momentum resolution at first cluster
  histoFile->mkdir("momentumAtFirstCluster","momentumAtFirstCluster");
  histoFile->cd("momentumAtFirstCluster");
  
  const Int_t deltaPAtFirstClNBins = 500;
  const Double_t deltaPAtFirstClEdges[2] = {-25., 25.};
  
  TH1F *hResMomFirstCluster = new TH1F("hResMomFirstCluster"," delta P at first cluster;#Delta_{p} (GeV/c)",deltaPAtFirstClNBins,deltaPAtFirstClEdges[0],deltaPAtFirstClEdges[1]);
  TH2D *hResMomFirstClusterVsMom = new TH2D("hResMomFirstClusterVsMom","#Delta_{p} at first cluster versus p;p (GeV/c);#Delta_{p} (GeV/c)",2*pNBins,pEdges[0],pEdges[1],deltaPAtFirstClNBins,deltaPAtFirstClEdges[0],deltaPAtFirstClEdges[1]);
  
  TGraphAsymmErrors* gMeanResMomFirstClusterVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResMomFirstClusterVsMom->SetName("gMeanResMomFirstClusterVsMom");
  gMeanResMomFirstClusterVsMom->SetTitle("<#Delta_{p}> at first cluster versus p;p (GeV/c);<#Delta_{p}> (GeV/c)");
  TGraphAsymmErrors* gSigmaResMomFirstClusterVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResMomFirstClusterVsMom->SetName("gSigmaResMomFirstClusterVsMom");
  gSigmaResMomFirstClusterVsMom->SetTitle("#sigma_{p}/p at first cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
  
  // angular resolution at vertex
  histoFile->mkdir("slopesAtVertex","slopesAtVertex");
  histoFile->cd("slopesAtVertex");
  
  const Int_t deltaSlopeAtVtxNBins = 500;
  const Double_t deltaSlopeAtVtxEdges[2] = {-0.05, 0.05};
  
  TH1F *hResSlopeXVertex = new TH1F("hResSlopeXVertex","#Delta_{slope_{X}} at vertex;#Delta_{slope_{X}}", deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH1F *hResSlopeYVertex = new TH1F("hResSlopeYVertex","#Delta_{slope_{Y}} at vertex;#Delta_{slope_{Y}}", deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeXVertexVsMom = new TH2D("hResSlopeXVertexVsMom","#Delta_{slope_{X}} at vertex versus p;p (GeV/c);#Delta_{slope_{X}}",2*pNBins,pEdges[0],pEdges[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeYVertexVsMom = new TH2D("hResSlopeYVertexVsMom","#Delta_{slope_{Y}} at vertex versus p;p (GeV/c);#Delta_{slope_{Y}}",2*pNBins,pEdges[0],pEdges[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  
  TH2D *hResSlopeXVertexVsPosAbsEnd_0_2_DegMC = new TH2D("hResSlopeXVertexVsPosAbsEnd_0_2_DegMC","#Delta_{slope_{X}} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{slope_{X}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeYVertexVsPosAbsEnd_0_2_DegMC = new TH2D("hResSlopeYVertexVsPosAbsEnd_0_2_DegMC","#Delta_{slope_{Y}} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{slope_{Y}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeXVertexVsPosAbsEnd_2_3_DegMC = new TH2D("hResSlopeXVertexVsPosAbsEnd_2_3_DegMC","#Delta_{slope_{X}} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{slope_{X}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeYVertexVsPosAbsEnd_2_3_DegMC = new TH2D("hResSlopeYVertexVsPosAbsEnd_2_3_DegMC","#Delta_{slope_{Y}} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{slope_{Y}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeXVertexVsPosAbsEnd_3_10_DegMC = new TH2D("hResSlopeXVertexVsPosAbsEnd_3_10_DegMC","#Delta_{slope_{X}} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{slope_{X}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeYVertexVsPosAbsEnd_3_10_DegMC = new TH2D("hResSlopeYVertexVsPosAbsEnd_3_10_DegMC","#Delta_{slope_{Y}} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{slope_{Y}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  
  TH2D *hResSlopeXVertexVsAngle = new TH2D("hResSlopeXVertexVsAngle","#Delta_{slope_{X}} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{slope_{X}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeYVertexVsAngle = new TH2D("hResSlopeYVertexVsAngle","#Delta_{slope_{Y}} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{slope_{Y}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeXVertexVsMCAngle = new TH2D("hResSlopeXVertexVsMCAngle","#Delta_{slope_{X}} at vertex versus MC angle;MC angle (Deg);#Delta_{slope_{X}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  TH2D *hResSlopeYVertexVsMCAngle = new TH2D("hResSlopeYVertexVsMCAngle","#Delta_{slope_{Y}} at vertex versus MC angle;MC angle (Deg);#Delta_{slope_{Y}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  
  TGraphAsymmErrors* gMeanResSlopeXVertexVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResSlopeXVertexVsMom->SetName("gMeanResSlopeXVertexVsMom");
  gMeanResSlopeXVertexVsMom->SetTitle("<#Delta_{slope_{X}}> at vertex versus p;p (GeV/c);<#Delta_{slope_{X}}>");
  TGraphAsymmErrors* gSigmaResSlopeXVertexVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResSlopeXVertexVsMom->SetName("gSigmaResSlopeXVertexVsMom");
  gSigmaResSlopeXVertexVsMom->SetTitle("#sigma_{slope_{X}} at vertex versus p;p (GeV/c);#sigma_{slope_{X}}");
  TGraphAsymmErrors* gMeanResSlopeYVertexVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResSlopeYVertexVsMom->SetName("gMeanResSlopeYVertexVsMom");
  gMeanResSlopeYVertexVsMom->SetTitle("<#Delta_{slope_{Y}}> at vertex versus p;p (GeV/c);<#Delta_{slope_{Y}}>");
  TGraphAsymmErrors* gSigmaResSlopeYVertexVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResSlopeYVertexVsMom->SetName("gSigmaResSlopeYVertexVsMom");
  gSigmaResSlopeYVertexVsMom->SetTitle("#sigma_{slope_{Y}} at vertex versus p;p (GeV/c);#sigma_{slope_{Y}}");
  
  // angular resolution at first cluster
  histoFile->mkdir("slopesAtFirstCluster","slopesAtFirstCluster");
  histoFile->cd("slopesAtFirstCluster");
  
  const Int_t deltaSlopeAtFirstClNBins = 500;
  const Double_t deltaSlopeAtFirstClEdges[2] = {-0.01, 0.01};
  
  TH1F *hResSlopeXFirstCluster = new TH1F("hResSlopeXFirstCluster","#Delta_{slope_{X}} at first cluster;#Delta_{slope_{X}}", deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  TH2D *hResSlopeXFirstClusterVsMom = new TH2D("hResSlopeXFirstClusterVsMom","#Delta_{slope_{X}} at first cluster versus p;p (GeV/c);#Delta_{slope_{X}}",2*pNBins,pEdges[0],pEdges[1], deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  TH1F *hResSlopeYFirstCluster = new TH1F("hResSlopeYFirstCluster","#Delta_{slope_{Y}} at first cluster;#Delta_{slope_{Y}}", deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  TH2D *hResSlopeYFirstClusterVsMom = new TH2D("hResSlopeYFirstClusterVsMom","#Delta_{slope_{Y}} at first cluster versus p;p (GeV/c);#Delta_{slope_{Y}}",2*pNBins,pEdges[0],pEdges[1], deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  
  TGraphAsymmErrors* gMeanResSlopeXFirstClusterVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResSlopeXFirstClusterVsMom->SetName("gMeanResSlopeXFirstClusterVsMom");
  gMeanResSlopeXFirstClusterVsMom->SetTitle("<#Delta_{slope_{X}}> at first cluster versus p;p (GeV/c);<#Delta_{slope_{X}}>");
  TGraphAsymmErrors* gSigmaResSlopeXFirstClusterVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResSlopeXFirstClusterVsMom->SetName("gSigmaResSlopeXFirstClusterVsMom");
  gSigmaResSlopeXFirstClusterVsMom->SetTitle("#sigma_{slope_{X}} at first cluster versus p;p (GeV/c);#sigma_{slope_{X}}");
  TGraphAsymmErrors* gMeanResSlopeYFirstClusterVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResSlopeYFirstClusterVsMom->SetName("gMeanResSlopeYFirstClusterVsMom");
  gMeanResSlopeYFirstClusterVsMom->SetTitle("<#Delta_{slope_{Y}}> at first cluster versus p;p (GeV/c);<#Delta_{slope_{Y}}>");
  TGraphAsymmErrors* gSigmaResSlopeYFirstClusterVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResSlopeYFirstClusterVsMom->SetName("gSigmaResSlopeYFirstClusterVsMom");
  gSigmaResSlopeYFirstClusterVsMom->SetTitle("#sigma_{slope_{Y}} at first cluster versus p;p (GeV/c);#sigma_{slope_{Y}}");
  
  // DCA resolution and MCS angular dispersion
  histoFile->mkdir("DCA","DCA");
  histoFile->cd("DCA");
  
  const Int_t deltaPDCANBins = 500;
  const Double_t deltaPDCAEdges[2] = {0., 1000.};
  const Double_t deltaPMCSAngEdges[2] = {-0.5, 0.5};
  
  TH1F *hPDCA = new TH1F("hPDCA","p #times DCA at vertex;p #times DCA (GeV #times cm)", deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  TH2D *hPDCAVsMom_2_3_Deg = new TH2D("hPDCAVsMom_2_3_Deg","p #times DCA versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DCA (GeV #times cm)",2*pNBins,pEdges[0],pEdges[1], deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  TH2D *hPDCAVsMom_3_10_Deg = new TH2D("hPDCAVsMom_3_10_Deg","p #times DCA versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DCA (GeV #times cm)",2*pNBins,pEdges[0],pEdges[1], deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  TH2D *hPMCSAngVsMom_2_3_Deg = new TH2D("hPMCSAngVsMom_2_3_Deg","p #times #Delta#theta_{MCS} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #Delta#theta_{MCS} (GeV)",2*pNBins,pEdges[0],pEdges[1], deltaPDCANBins, deltaPMCSAngEdges[0], deltaPMCSAngEdges[1]);
  TH2D *hPMCSAngVsMom_3_10_Deg = new TH2D("hPMCSAngVsMom_3_10_Deg","p #times #Delta#theta_{MCS} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #Delta#theta_{MCS} (GeV)",2*pNBins,pEdges[0],pEdges[1], deltaPDCANBins, deltaPMCSAngEdges[0], deltaPMCSAngEdges[1]);
  
  TH2D *hPDCAVsPosAbsEnd_0_2_DegMC = new TH2D("hPDCAVsPosAbsEnd_0_2_DegMC","p #times DCA versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);p #times DCA (GeV #times cm)",1000,0.,100.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  TH2D *hPDCAVsPosAbsEnd_2_3_DegMC = new TH2D("hPDCAVsPosAbsEnd_2_3_DegMC","p #times DCA}versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);p #times DCA (GeV #times cm)",1000,0.,100.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  TH2D *hPDCAVsPosAbsEnd_3_10_DegMC = new TH2D("hPDCAVsPosAbsEnd_3_10_DegMC","p #times DCA versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);p #times DCA (GeV #times cm)",1000,0.,100.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  
  TH2D *hPDCAVsAngle = new TH2D("hPDCAVsAngle","p #times DCA versus track position at absorber end converted to degrees;angle (Deg);p #times DCA (GeV #times cm)",10,0.,10.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  TH2D *hPDCAVsMCAngle = new TH2D("hPDCAVsMCAngle","p #times DCA versus MC angle;MC angle (Deg);p #times DCA (GeV #times cm)",10,0.,10.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  
  TGraphAsymmErrors* gSigmaPDCAVsMom_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAVsMom_2_3_Deg->SetName("gSigmaPDCAVsMom_2_3_Deg");
  gSigmaPDCAVsMom_2_3_Deg->SetTitle("#sigma_{p #times DCA} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCA} (GeV #times cm)");
  TGraphAsymmErrors* gSigmaPDCAVsMom_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPDCAVsMom_3_10_Deg->SetName("gSigmaPDCAVsMom_3_10_Deg");
  gSigmaPDCAVsMom_3_10_Deg->SetTitle("#sigma_{p #times DCA} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCA} (GeV #times cm)");
  TGraphAsymmErrors* gMeanPMCSAngVsMom_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPMCSAngVsMom_2_3_Deg->SetName("gMeanPMCSAngVsMom_2_3_Deg");
  gMeanPMCSAngVsMom_2_3_Deg->SetTitle("<p #times #Delta#theta_{MCS}> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times #Delta#theta_{MCS}> (GeV)");
  TGraphAsymmErrors* gSigmaPMCSAngVsMom_2_3_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPMCSAngVsMom_2_3_Deg->SetName("gSigmaPMCSAngVsMom_2_3_Deg");
  gSigmaPMCSAngVsMom_2_3_Deg->SetTitle("#sigma_{p #times #Delta#theta_{MCS}} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times #Delta#theta_{MCS}} (GeV)");
  TGraphAsymmErrors* gMeanPMCSAngVsMom_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gMeanPMCSAngVsMom_3_10_Deg->SetName("gMeanPMCSAngVsMom_3_10_Deg");
  gMeanPMCSAngVsMom_3_10_Deg->SetTitle("<p #times #Delta#theta_{MCS}> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times #Delta#theta_{MCS}> (GeV)");
  TGraphAsymmErrors* gSigmaPMCSAngVsMom_3_10_Deg = new TGraphAsymmErrors(pNBins);
  gSigmaPMCSAngVsMom_3_10_Deg->SetName("gSigmaPMCSAngVsMom_3_10_Deg");
  gSigmaPMCSAngVsMom_3_10_Deg->SetTitle("#sigma_{p #times #Delta#theta_{MCS}} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times #Delta#theta_{MCS}} (GeV)");
  
  // eta resolution at vertex
  histoFile->mkdir("etaAtVertex","etaAtVertex");
  histoFile->cd("etaAtVertex");
  
  const Int_t deltaEtaAtVtxNBins = 500;
  const Double_t deltaEtaAtVtxEdges[2] = {-0.5, 0.5};
  
  TH1F *hResEtaVertex = new TH1F("hResEtaVertex","#Delta_{eta} at vertex;#Delta_{eta}", deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  TH2D *hResEtaVertexVsMom = new TH2D("hResEtaVertexVsMom","#Delta_{eta} at vertex versus p;p (GeV/c);#Delta_{eta}",2*pNBins,pEdges[0],pEdges[1], deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  
  TH2D *hResEtaVertexVsPosAbsEnd_0_2_DegMC = new TH2D("hResEtaVertexVsPosAbsEnd_0_2_DegMC","#Delta_{eta} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{eta}",1000,0.,100.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  TH2D *hResEtaVertexVsPosAbsEnd_2_3_DegMC = new TH2D("hResEtaVertexVsPosAbsEnd_2_3_DegMC","#Delta_{eta} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{eta}",1000,0.,100.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  TH2D *hResEtaVertexVsPosAbsEnd_3_10_DegMC = new TH2D("hResEtaVertexVsPosAbsEnd_3_10_DegMC","#Delta_{eta} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{eta}",1000,0.,100.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  
  TH2D *hResEtaVertexVsAngle = new TH2D("hResEtaVertexVsAngle","#Delta_{eta} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{eta}",10,0.,10.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  TH2D *hResEtaVertexVsMCAngle = new TH2D("hResEtaVertexVsMCAngle","#Delta_{eta} at vertex versus MC angle;MC angle (Deg);#Delta_{eta}",10,0.,10.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  
  TGraphAsymmErrors* gMeanResEtaVertexVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResEtaVertexVsMom->SetName("gMeanResEtaVertexVsMom");
  gMeanResEtaVertexVsMom->SetTitle("<#Delta_{eta}> at vertex versus p;p (GeV/c);<#Delta_{eta}>");
  TGraphAsymmErrors* gSigmaResEtaVertexVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResEtaVertexVsMom->SetName("gSigmaResEtaVertexVsMom");
  gSigmaResEtaVertexVsMom->SetTitle("#sigma_{eta} at vertex versus p;p (GeV/c);#sigma_{eta}");
  
  // phi resolution at vertex
  histoFile->mkdir("phiAtVertex","phiAtVertex");
  histoFile->cd("phiAtVertex");
  
  const Int_t deltaPhiAtVtxNBins = 500;
  const Double_t deltaPhiAtVtxEdges[2] = {-0.5, 0.5};
  
  TH1F *hResPhiVertex = new TH1F("hResPhiVertex","#Delta_{phi} at vertex;#Delta_{phi}", deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  TH2D *hResPhiVertexVsMom = new TH2D("hResPhiVertexVsMom","#Delta_{phi} at vertex versus p;p (GeV/c);#Delta_{phi}",2*pNBins,pEdges[0],pEdges[1], deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  
  TH2D *hResPhiVertexVsPosAbsEnd_0_2_DegMC = new TH2D("hResPhiVertexVsPosAbsEnd_0_2_DegMC","#Delta_{phi} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{phi}",1000,0.,100.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  TH2D *hResPhiVertexVsPosAbsEnd_2_3_DegMC = new TH2D("hResPhiVertexVsPosAbsEnd_2_3_DegMC","#Delta_{phi} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{phi}",1000,0.,100.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  TH2D *hResPhiVertexVsPosAbsEnd_3_10_DegMC = new TH2D("hResPhiVertexVsPosAbsEnd_3_10_DegMC","#Delta_{phi} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{phi}",1000,0.,100.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  
  TH2D *hResPhiVertexVsAngle = new TH2D("hResPhiVertexVsAngle","#Delta_{phi} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{phi}",10,0.,10.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  TH2D *hResPhiVertexVsMCAngle = new TH2D("hResPhiVertexVsMCAngle","#Delta_{phi} at vertex versus MC angle;MC angle (Deg);#Delta_{phi}",10,0.,10.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  
  TGraphAsymmErrors* gMeanResPhiVertexVsMom = new TGraphAsymmErrors(pNBins);
  gMeanResPhiVertexVsMom->SetName("gMeanResPhiVertexVsMom");
  gMeanResPhiVertexVsMom->SetTitle("<#Delta_{phi}> at vertex versus p;p (GeV/c);<#Delta_{phi}>");
  TGraphAsymmErrors* gSigmaResPhiVertexVsMom = new TGraphAsymmErrors(pNBins);
  gSigmaResPhiVertexVsMom->SetName("gSigmaResPhiVertexVsMom");
  gSigmaResPhiVertexVsMom->SetTitle("#sigma_{phi} at vertex versus p;p (GeV/c);#sigma_{phi}");
  
  // cluster resolution
  histoFile->mkdir("clusters","clusters");
  histoFile->cd("clusters");
  
  TH1F* hResidualXInCh[AliMUONConstants::NTrackingCh()];
  TH1F* hResidualYInCh[AliMUONConstants::NTrackingCh()];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInCh[i] = new TH1F(Form("hResidualXInCh%d",i+1), Form("cluster-track residual-X distribution in chamber %d;#Delta_{X} (cm)",i+1), 1000, -1., 1.);
    hResidualYInCh[i] = new TH1F(Form("hResidualYInCh%d",i+1), Form("cluster-track residual-Y distribution in chamber %d;#Delta_{Y} (cm)",i+1), 1000, -0.5, 0.5);
  }
  
  TGraphErrors* gResidualXPerChMean = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualXPerChMean->SetName("gResidualXPerChMean");
  gResidualXPerChMean->SetTitle("cluster-trackRef residual-X per Ch: mean;chamber ID;<#Delta_{X}> (cm)");
  gResidualXPerChMean->SetMarkerStyle(kFullDotLarge);
  TGraphErrors* gResidualYPerChMean = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualYPerChMean->SetName("gResidualYPerChMean");
  gResidualYPerChMean->SetTitle("cluster-trackRef residual-Y per Ch: mean;chamber ID;<#Delta_{Y}> (cm)");
  gResidualYPerChMean->SetMarkerStyle(kFullDotLarge);
  TGraphErrors* gResidualXPerChSigma = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualXPerChSigma->SetName("gResidualXPerChSigma");
  gResidualXPerChSigma->SetTitle("cluster-trackRef residual-X per Ch: sigma;chamber ID;#sigma_{X} (cm)");
  gResidualXPerChSigma->SetMarkerStyle(kFullDotLarge);
  TGraphErrors* gResidualYPerChSigma = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gResidualYPerChSigma->SetName("gResidualYPerChSigma");
  gResidualYPerChSigma->SetTitle("cluster-trackRef residual-Y per Ch: sigma;chamber ID;#sigma_{Y} (cm)");
  gResidualYPerChSigma->SetMarkerStyle(kFullDotLarge);

  histoFile->mkdir("trigger");
  histoFile->cd("trigger");
  TH1F* hResidualTrigX11 = new TH1F("hResiudalTrigX11", "Residual X11", 100, -10., 10.);
  TH1F* hResidualTrigY11 = new TH1F("hResiudalTrigY11", "Residual Y11", 100, -10., 10.);
  TH1F* hResidualTrigSlopeY = new TH1F("hResiudalTrigSlopeY", "Residual Y slope", 100, -0.1, 0.1);
  TH1F* hTriggerableMatchFailed = new TH1F("hTriggerableMatchFailed", "Triggerable multiplicity for events with no match", 15, -0.5, 14.5);
  
  
  // ###################################### initialize ###################################### //
  AliMUONRecoCheck rc(esdFileName, pathSim);
  
  // load necessary data from OCDB
  AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
  AliCDBManager::Instance()->SetRun(rc.GetRunNumber());
  if (!AliMUONCDB::LoadField()) return;
  AliMUONTrackExtrap::SetField();
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry()) return;
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  AliMUONESDInterface::ResetTracker(recoParam);
  
  // get sigma cut from recoParam to associate clusters with TrackRefs in case the label are not used
  Double_t sigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking();
  // compute the mask of requested stations from recoParam
  UInt_t requestedStationMask = 0;
  for (Int_t i = 0; i < 5; i++) if (recoParam->RequestStation(i)) requestedStationMask |= ( 1 << i );
  // get from recoParam whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
  Bool_t request2ChInSameSt45 = !recoParam->MakeMoreTrackCandidates();
  
  Int_t nevents = rc.NumberOfEvents();
  
  if (nevents < nEvent || nEvent < 0) nEvent = nevents;
  
  Int_t ievent;
  Int_t nReconstructibleTracks = 0;
  Int_t nReconstructedTracks = 0;
  Int_t nReconstructibleTracksCheck = 0;
  AliMUONTrackParam *trackParam;
  Double_t x1,y1,z1,slopex1,slopey1,pX1,pY1,pZ1,p1,pT1,eta1,phi1;
  Double_t x2,y2,z2,slopex2,slopey2,pX2,pY2,pZ2,p2,pT2,eta2,phi2;
  Double_t xAbs,yAbs,dAbs,aAbs,aMCS,aMC;
  Double_t xDCA,yDCA,dca,pU;
  Double_t aMCSMoy = 0., aMCS2Moy = 0., dMCSMoy = 0., dMCS2Moy = 0., adMCSMoy = 0.;
  Int_t nMCS = 0;
  
  // ###################################### fill histograms ###################################### //
  for (ievent=0; ievent<nEvent; ievent++)
  {
    if ((ievent+1)%100 == 0) cout<<"\rEvent processing... "<<ievent+1<<flush;
    
    AliMUONVTrackStore* trackStore = rc.ReconstructedTracks(ievent, kFALSE);
    AliMUONVTrackStore* trackRefStore = rc.ReconstructibleTracks(ievent, requestedStationMask, request2ChInSameSt45);
    
    hReconstructible->Fill(trackRefStore->GetSize());
    hReco->Fill(trackStore->GetSize());
    
    nReconstructibleTracks += trackRefStore->GetSize();
    nReconstructedTracks += trackStore->GetSize();

    AliMUONVTriggerTrackStore* triggerTrackRefStore = rc.TriggerableTracks(ievent);
    AliMUONVTriggerTrackStore* triggerTrackStore = rc.TriggeredTracks(ievent);

    hTriggerable->Fill(triggerTrackRefStore->GetSize());
    hTriggered->Fill(triggerTrackStore->GetSize());

    // loop over trigger trackRef
    TIter nextTrig(triggerTrackRefStore->CreateIterator());
    AliMUONTriggerTrack* triggerTrackRef;
    Int_t nTriggerMatches = 0;
    while ( ( triggerTrackRef = static_cast<AliMUONTriggerTrack*>(nextTrig()) ) )
    {
      
      AliMUONTriggerTrack* triggerTrackMatched = 0x0;
      
      // loop over trackReco and look for compatible track
      TIter nextTrig2(triggerTrackStore->CreateIterator());
      AliMUONTriggerTrack* triggerTrackReco;
      while ( ( triggerTrackReco = static_cast<AliMUONTriggerTrack*>(nextTrig2()) ) )
      {
	
        // check if trackReco is compatible with trackRef
        if (triggerTrackReco->Match(*triggerTrackRef, sigmaCut)) {
          triggerTrackMatched = triggerTrackReco;
          nTriggerMatches++;
          break;
        }
      }
      
      if (triggerTrackMatched) { // tracking requirements verified, track is found
        hResidualTrigX11->Fill( triggerTrackMatched->GetX11() - triggerTrackRef->GetX11() );
        hResidualTrigY11->Fill( triggerTrackMatched->GetY11() - triggerTrackRef->GetY11() );
        hResidualTrigSlopeY->Fill( triggerTrackMatched->GetSlopeY() - triggerTrackRef->GetSlopeY() );
      }
    } // loop on trigger track ref
    
    if ( nTriggerMatches != triggerTrackStore->GetSize() )
      hTriggerableMatchFailed->Fill(triggerTrackRefStore->GetSize());
    
    // loop over trackRef
    TIter next(trackRefStore->CreateIterator());
    AliMUONTrack* trackRef;
    while ( ( trackRef = static_cast<AliMUONTrack*>(next()) ) )
    {
      
      hTrackRefID->Fill(trackRef->GetUniqueID());
      
      AliMUONTrack* trackMatched = 0x0;
      Int_t nMatchClusters = 0;
      
      // loop over trackReco and look for compatible track
      TIter next2(trackStore->CreateIterator());
      AliMUONTrack* trackReco;
      while ( ( trackReco = static_cast<AliMUONTrack*>(next2()) ) )
      {
	
	// check if trackReco is compatible with trackRef
	if (trackReco->Match(*trackRef, sigmaCut, nMatchClusters)) {
	  trackMatched = trackReco;
	  break;
	}
	
      }
      
      if (trackMatched) { // tracking requirements verified, track is found
        nReconstructibleTracksCheck++;
        hNClusterComp->Fill(nMatchClusters);
	
	// compute track position at the end of the absorber
        AliMUONTrackParam trackParamAtAbsEnd(*((AliMUONTrackParam*)trackMatched->GetTrackParamAtCluster()->First()));
	AliMUONTrackExtrap::ExtrapToZ(&trackParamAtAbsEnd, AliMUONConstants::AbsZEnd());
        xAbs = trackParamAtAbsEnd.GetNonBendingCoor();
        yAbs = trackParamAtAbsEnd.GetBendingCoor();
	dAbs = TMath::Sqrt(xAbs*xAbs + yAbs*yAbs);
	aAbs = TMath::ATan(-dAbs/AliMUONConstants::AbsZEnd()) * TMath::RadToDeg();
        pX2 = trackParamAtAbsEnd.Px();
        pY2 = trackParamAtAbsEnd.Py();
        pZ2 = trackParamAtAbsEnd.Pz();
        pT2 = TMath::Sqrt(pX2*pX2 + pY2*pY2);
	aMCS = TMath::ATan(-pT2/pZ2) * TMath::RadToDeg();
	
        trackParam = trackRef->GetTrackParamAtVertex();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        slopex1 = trackParam->GetNonBendingSlope();
        slopey1 = trackParam->GetBendingSlope();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        pT1 = TMath::Sqrt(pX1*pX1 + pY1*pY1);
	aMC = TMath::ATan(-pT1/pZ1) * TMath::RadToDeg();
	eta1 = TMath::Log(TMath::Tan(0.5*TMath::ATan(-pT1/pZ1)));
	phi1 = TMath::Pi()+TMath::ATan2(-pY1, -pX1);
	
	trackParam = trackMatched->GetTrackParamAtVertex();
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        slopex2 = trackParam->GetNonBendingSlope();
        slopey2 = trackParam->GetBendingSlope();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        pT2 = TMath::Sqrt(pX2*pX2 + pY2*pY2);
	eta2 = TMath::Log(TMath::Tan(0.5*TMath::ATan(-pT2/pZ2)));
	phi2 = TMath::Pi()+TMath::ATan2(-pY2, -pX2);
        
        AliMUONTrackParam trackParamAtDCA(*((AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First()));
	pU = trackParamAtDCA.P();
	AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamAtDCA, z2);
        xDCA = trackParamAtDCA.GetNonBendingCoor();
        yDCA = trackParamAtDCA.GetBendingCoor();
	dca = TMath::Sqrt(xDCA*xDCA + yDCA*yDCA);
	
        hResMomVertex->Fill(p2-p1);
	hResSlopeXVertex->Fill(slopex2-slopex1);
	hResSlopeYVertex->Fill(slopey2-slopey1);
	hPDCA->Fill(0.5*(p2+pU)*dca);
	hResEtaVertex->Fill(eta2-eta1);
	hResPhiVertex->Fill(phi2-phi1);
	if (aMC >= aAbsLimits[0] && aMC <= aAbsLimits[1]) {
	  hResMomVertexVsMom->Fill(p1,p2-p1);
	  hResSlopeXVertexVsMom->Fill(p1,slopex2-slopex1);
	  hResSlopeYVertexVsMom->Fill(p1,slopey2-slopey1);
	  hResEtaVertexVsMom->Fill(p1,eta2-eta1);
	  hResPhiVertexVsMom->Fill(p1,phi2-phi1);
	}
	hResMomVertexVsAngleVsMom->Fill(p1,aAbs,p2-p1);
	if (aAbs > 2. && aAbs < 3.) {
	  hResMomVertexVsMom_2_3_Deg->Fill(p1,p2-p1);
	  hPDCAVsMom_2_3_Deg->Fill(p1,0.5*(p2+pU)*dca);
	  hPMCSAngVsMom_2_3_Deg->Fill(p1,0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad());
	}
	else if (aAbs >= 3. && aAbs < 10.) {
	  hResMomVertexVsMom_3_10_Deg->Fill(p1,p2-p1);
	  hPDCAVsMom_3_10_Deg->Fill(p1,0.5*(p2+pU)*dca);
	  hPMCSAngVsMom_3_10_Deg->Fill(p1,0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad());
	  aMCSMoy += 0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad();
	  aMCS2Moy += (0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad()) * (0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad());
	  dMCSMoy += 0.5*(p2+pU)*(dAbs-pT1/pZ1*AliMUONConstants::AbsZEnd());
	  dMCS2Moy += (0.5*(p2+pU)*(dAbs-pT1/pZ1*AliMUONConstants::AbsZEnd())) * (0.5*(p2+pU)*(dAbs-pT1/pZ1*AliMUONConstants::AbsZEnd()));
	  adMCSMoy += (0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad()) * (0.5*(p2+pU)*(dAbs-pT1/pZ1*AliMUONConstants::AbsZEnd()));
	  nMCS++;
	}
	if (aMC < 2.) {
	  hResMomVertexVsMom_0_2_DegMC->Fill(p1,p2-p1);
	  hResMomVertexVsPosAbsEnd_0_2_DegMC->Fill(dAbs,p2-p1);
	  hResSlopeXVertexVsPosAbsEnd_0_2_DegMC->Fill(dAbs,slopex2-slopex1);
	  hResSlopeYVertexVsPosAbsEnd_0_2_DegMC->Fill(dAbs,slopey2-slopey1);
	  hPDCAVsPosAbsEnd_0_2_DegMC->Fill(dAbs,0.5*(p2+pU)*dca);
	  hResEtaVertexVsPosAbsEnd_0_2_DegMC->Fill(dAbs,eta2-eta1);
	  hResPhiVertexVsPosAbsEnd_0_2_DegMC->Fill(dAbs,phi2-phi1);
	}
	else if (aMC >= 2. && aMC < 3) {
	  hResMomVertexVsPosAbsEnd_2_3_DegMC->Fill(dAbs,p2-p1);
	  hResSlopeXVertexVsPosAbsEnd_2_3_DegMC->Fill(dAbs,slopex2-slopex1);
	  hResSlopeYVertexVsPosAbsEnd_2_3_DegMC->Fill(dAbs,slopey2-slopey1);
	  hPDCAVsPosAbsEnd_2_3_DegMC->Fill(dAbs,0.5*(p2+pU)*dca);
	  hResEtaVertexVsPosAbsEnd_2_3_DegMC->Fill(dAbs,eta2-eta1);
	  hResPhiVertexVsPosAbsEnd_2_3_DegMC->Fill(dAbs,phi2-phi1);
	}
	else if (aMC >= 3. && aMC < 10.) {
	  hResMomVertexVsPosAbsEnd_3_10_DegMC->Fill(dAbs,p2-p1);
	  hResSlopeXVertexVsPosAbsEnd_3_10_DegMC->Fill(dAbs,slopex2-slopex1);
	  hResSlopeYVertexVsPosAbsEnd_3_10_DegMC->Fill(dAbs,slopey2-slopey1);
	  hPDCAVsPosAbsEnd_3_10_DegMC->Fill(dAbs,0.5*(p2+pU)*dca);
	  hResEtaVertexVsPosAbsEnd_3_10_DegMC->Fill(dAbs,eta2-eta1);
	  hResPhiVertexVsPosAbsEnd_3_10_DegMC->Fill(dAbs,phi2-phi1);
	}
	hResMomVertexVsAngle->Fill(aAbs,p2-p1);
	hResSlopeXVertexVsAngle->Fill(aAbs,slopex2-slopex1);
	hResSlopeYVertexVsAngle->Fill(aAbs,slopey2-slopey1);
	hPDCAVsAngle->Fill(aAbs,0.5*(p2+pU)*dca);
	hResEtaVertexVsAngle->Fill(aAbs,eta2-eta1);
	hResPhiVertexVsAngle->Fill(aAbs,phi2-phi1);
	hResMomVertexVsMCAngle->Fill(aMC,p2-p1);
	hResSlopeXVertexVsMCAngle->Fill(aMC,slopex2-slopex1);
	hResSlopeYVertexVsMCAngle->Fill(aMC,slopey2-slopey1);
	hPDCAVsMCAngle->Fill(aMC,0.5*(p2+pU)*dca);
	hResEtaVertexVsMCAngle->Fill(aMC,eta2-eta1);
	hResPhiVertexVsMCAngle->Fill(aMC,phi2-phi1);
	
        trackParam = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->First();
        x1 = trackParam->GetNonBendingCoor();
        y1 = trackParam->GetBendingCoor();
        z1 = trackParam->GetZ();
        slopex1 = trackParam->GetNonBendingSlope();
        slopey1 = trackParam->GetBendingSlope();
        pX1 = trackParam->Px();
        pY1 = trackParam->Py();
        pZ1 = trackParam->Pz();
        p1  = trackParam->P();
        pT1 = TMath::Sqrt(pX1*pX1 + pY1*pY1);
	
        trackParam = (AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First();
        x2 = trackParam->GetNonBendingCoor();
        y2 = trackParam->GetBendingCoor();
        z2 = trackParam->GetZ();
        slopex2 = trackParam->GetNonBendingSlope();
        slopey2 = trackParam->GetBendingSlope();
        pX2 = trackParam->Px();
        pY2 = trackParam->Py();
        pZ2 = trackParam->Pz();
        p2  = trackParam->P();
        pT2 = TMath::Sqrt(pX2*pX2 + pY2*pY2);
        
        hResMomFirstCluster->Fill(p2-p1);
	hResMomFirstClusterVsMom->Fill(p1,p2-p1);
	
	hResSlopeXFirstCluster->Fill(slopex2-slopex1);
	hResSlopeYFirstCluster->Fill(slopey2-slopey1);
	hResSlopeXFirstClusterVsMom->Fill(p1,slopex2-slopex1);
	hResSlopeYFirstClusterVsMom->Fill(p1,slopey2-slopey1);
	
	// Fill residuals
	// Loop over clusters of first track
	AliMUONTrackParam* trackParamAtCluster1 = (AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->First();
	while (trackParamAtCluster1) {
	  AliMUONVCluster* cluster1 = trackParamAtCluster1->GetClusterPtr();
	  AliMUONTrackParam* trackParamAtCluster2 = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->First();
	  while (trackParamAtCluster2) {
	    AliMUONVCluster* cluster2 = trackParamAtCluster2->GetClusterPtr();
	    if (cluster1->GetDetElemId() == cluster2->GetDetElemId()) {
	      hResidualXInCh[cluster1->GetChamberId()]->Fill(cluster1->GetX() - cluster2->GetX());
	      hResidualYInCh[cluster1->GetChamberId()]->Fill(cluster1->GetY() - cluster2->GetY());
	      break;
	    }
	    trackParamAtCluster2 = (AliMUONTrackParam*) trackRef->GetTrackParamAtCluster()->After(trackParamAtCluster2);
	  }
	  trackParamAtCluster1 = (AliMUONTrackParam*) trackMatched->GetTrackParamAtCluster()->After(trackParamAtCluster1);
	}
	
      }
      
    } // end loop track ref.

  } // end loop on event  
  cout<<"\rEvent processing... "<<nevents<<" done"<<endl;
  
  // ###################################### compute stuff ###################################### //
  cout<<"\nWhen not specified, resolution at vertex is computed for ";
  if (absorberRegion == 1) cout<<"tracks in the absorber region [2,3] deg."<<endl;
  else if (absorberRegion == 2) cout<<"tracks in the absorber region [3,10] deg."<<endl;
  else cout<<"all tracks"<<endl;
  
  // compute momentum resolution at vertex versus p
  TF1 *f2 = new TF1("f2",langaufun,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1],4);
  Int_t rebinFactorX = TMath::Max(hResMomVertexVsMom->GetNbinsX()/pNBins, 1);
  for (Int_t i = rebinFactorX; i <= hResMomVertexVsMom->GetNbinsX(); i+=rebinFactorX) {
    cout<<"\rFitting momentum residuals at vertex... "<<i/rebinFactorX<<"/"<<pNBins<<flush;
    TH1D *tmp = hResMomVertexVsMom->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    f2->SetParameters(0.2,0.,(Double_t)tmp->GetEntries(),1.);
    tmp->Fit("f2","WWNQ");
    Double_t fwhm = f2->GetParameter(0);
    Double_t sigma = f2->GetParameter(3);
    Double_t sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Int_t rebin = TMath::Max(Int_t(0.5*sigmaP/tmp->GetBinWidth(1)),1);
    while (deltaPAtVtxNBins%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    tmp->Fit("f2","NQ");
    fwhm = f2->GetParameter(0);
    sigma = f2->GetParameter(3);
    sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Double_t fwhmErr = f2->GetParError(0);
    Double_t sigmaErr = f2->GetParError(3);
    Double_t sigmaPErr = TMath::Sqrt(sigma*sigma*sigmaErr*sigmaErr + fwhm*fwhm*fwhmErr*fwhmErr/(64.*log(2.)*log(2.))) / sigmaP;
    hResMomVertexVsMom->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = hResMomVertexVsMom->GetMean();
    hResMomVertexVsMom->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-hResMomVertexVsMom->GetBinLowEdge(i-rebinFactorX+1), hResMomVertexVsMom->GetBinLowEdge(i+1)-p};
    gMeanResMomVertexVsMom->SetPoint(i/rebinFactorX-1, p, tmp->GetMean());
    gMeanResMomVertexVsMom->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], tmp->GetMeanError(), tmp->GetMeanError());
    gMostProbResMomVertexVsMom->SetPoint(i/rebinFactorX-1, p, -f2->GetParameter(1));
    gMostProbResMomVertexVsMom->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], f2->GetParError(1), f2->GetParError(1));
    gSigmaResMomVertexVsMom->SetPoint(i/rebinFactorX-1, p, 100.*sigmaP/p);
    gSigmaResMomVertexVsMom->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], 100.*sigmaPErr/p, 100.*sigmaPErr/p);
    delete tmp;
  }
  cout<<"\rFitting momentum residuals at vertex... "<<pNBins<<"/"<<pNBins<<endl;
  
  // compute momentum relative resolution at first cluster versus p
  FitGausResVsMom(hResMomFirstClusterVsMom, pNBins, 0., 1., "momentum residuals at first cluster", gMeanResMomFirstClusterVsMom, gSigmaResMomFirstClusterVsMom);
  rebinFactorX = TMath::Max(hResMomFirstClusterVsMom->GetNbinsX()/pNBins, 1);
  for (Int_t i = rebinFactorX; i <= hResMomFirstClusterVsMom->GetNbinsX(); i+=rebinFactorX) {
    Double_t x,y;
    gSigmaResMomFirstClusterVsMom->GetPoint(i/rebinFactorX-1, x, y);
    gSigmaResMomFirstClusterVsMom->SetPoint(i/rebinFactorX-1, x, 100.*y/x);
    gSigmaResMomFirstClusterVsMom->SetPointEYlow(i/rebinFactorX-1, 100.*gSigmaResMomFirstClusterVsMom->GetErrorYlow(i/rebinFactorX-1)/x);
    gSigmaResMomFirstClusterVsMom->SetPointEYhigh(i/rebinFactorX-1, 100.*gSigmaResMomFirstClusterVsMom->GetErrorYhigh(i/rebinFactorX-1)/x);
  }
  
  // compute slopeX resolution at vertex versus p
  FitGausResVsMom(hResSlopeXVertexVsMom, pNBins, 0., 2.e-3, "slopeX residuals at vertex", gMeanResSlopeXVertexVsMom, gSigmaResSlopeXVertexVsMom);
  
  // compute slopeY resolution at vertex versus p
  FitGausResVsMom(hResSlopeYVertexVsMom, pNBins, 0., 2.e-3, "slopeY residuals at vertex", gMeanResSlopeYVertexVsMom, gSigmaResSlopeYVertexVsMom);
  
  // compute slopeX resolution at first cluster versus p
  FitGausResVsMom(hResSlopeXFirstClusterVsMom, pNBins, 0., 3.e-4, "slopeX residuals at first cluster", gMeanResSlopeXFirstClusterVsMom, gSigmaResSlopeXFirstClusterVsMom);
  
  // compute slopeY resolution at first cluster versus p
  FitGausResVsMom(hResSlopeYFirstClusterVsMom, pNBins, 0., 2.e-4, "slopeY residuals at first cluster", gMeanResSlopeYFirstClusterVsMom, gSigmaResSlopeYFirstClusterVsMom);
  
  // compute p*DCA resolution in the region [2,3] deg at absorber end
  FitPDCAVsMom(hPDCAVsMom_2_3_Deg, pNBins, "p*DCA (tracks in [2,3] deg.)", gSigmaPDCAVsMom_2_3_Deg);
  
  // compute p*DCA resolution in the region [3,10] deg at absorber end
  FitPDCAVsMom(hPDCAVsMom_3_10_Deg, pNBins, "p*DCA (tracks in [3,10] deg.)", gSigmaPDCAVsMom_3_10_Deg);
  
  // compute MCS angular dispersion in the region [2,3] deg at absorber end
  FitGausResVsMom(hPMCSAngVsMom_2_3_Deg, pNBins, 0., 2.e-3, "p*MCSAngle (tracks in [2,3] deg.)", gMeanPMCSAngVsMom_2_3_Deg, gSigmaPMCSAngVsMom_2_3_Deg);
  
  // compute MCS angular dispersion in the region [3,10] deg at absorber end
  FitGausResVsMom(hPMCSAngVsMom_3_10_Deg, pNBins, 0., 2.e-3, "p*MCSAngle (tracks in [3,10] deg.)", gMeanPMCSAngVsMom_3_10_Deg, gSigmaPMCSAngVsMom_3_10_Deg);
  
  // compute eta resolution at vertex versus p
  FitGausResVsMom(hResEtaVertexVsMom, pNBins, 0., 0.1, "eta residuals at vertex", gMeanResEtaVertexVsMom, gSigmaResEtaVertexVsMom);
  
  // compute phi resolution at vertex versus p
  FitGausResVsMom(hResPhiVertexVsMom, pNBins, 0., 0.01, "phi residuals at vertex", gMeanResPhiVertexVsMom, gSigmaResPhiVertexVsMom);
  
  // compute cluster-track residual mean and dispersion
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    hResidualXInCh[i]->GetXaxis()->SetRangeUser(-3.*hResidualXInCh[i]->GetRMS(), 3.*hResidualXInCh[i]->GetRMS());
    gResidualXPerChMean->SetPoint(i, i+1, hResidualXInCh[i]->GetMean());
    gResidualXPerChMean->SetPointError(i, 0., hResidualXInCh[i]->GetMeanError());
    gResidualXPerChSigma->SetPoint(i, i+1, hResidualXInCh[i]->GetRMS());
    gResidualXPerChSigma->SetPointError(i, 0., hResidualXInCh[i]->GetRMSError());
    hResidualXInCh[i]->GetXaxis()->SetRange(0,0);
    hResidualYInCh[i]->GetXaxis()->SetRangeUser(-3.*hResidualYInCh[i]->GetRMS(), 3.*hResidualYInCh[i]->GetRMS());
    gResidualYPerChMean->SetPoint(i, i+1, hResidualYInCh[i]->GetMean());
    gResidualYPerChMean->SetPointError(i, 0., hResidualYInCh[i]->GetMeanError());
    gResidualYPerChSigma->SetPoint(i, i+1, hResidualYInCh[i]->GetRMS());
    gResidualYPerChSigma->SetPointError(i, 0., hResidualYInCh[i]->GetRMSError());
    hResidualYInCh[i]->GetXaxis()->SetRange(0,0);
  }
  
  // ###################################### display histograms ###################################### //
  // diplay momentum residuals
  TCanvas* cResMom = DrawVsAng("cResMom", "momentum residual at vertex in 3 angular regions", hResMomVertex, hResMomVertexVsAngle);
  TCanvas* cResMomMC = DrawVsAng("cResMomMC", "momentum residual at vertex in 3 MC angular regions", hResMomVertex, hResMomVertexVsMCAngle);
  TCanvas* cResMomVsPos = DrawVsPos("cResMomVsPos", "momentum residual at vertex versus position at absorber end in 3 MC angular regions",
				    hResMomVertexVsPosAbsEnd_0_2_DegMC, hResMomVertexVsPosAbsEnd_2_3_DegMC, hResMomVertexVsPosAbsEnd_3_10_DegMC);
  TCanvas* cResMom_2_3_Deg = DrawResMomVsMom("cResMom_2_3_Deg", "momentum residual for tracks between 2 and 3 degrees",
					     hResMomVertexVsMom_2_3_Deg, 10, f2, "momentum residuals at vertex (tracks in [2,3] deg.)");
  TCanvas* cResMom_3_10_Deg = DrawResMomVsMom("cResMom_3_10_Deg", "momentum residual for tracks between 3 and 10 degrees",
					      hResMomVertexVsMom_3_10_Deg, 10, f2, "momentum residuals at vertex (tracks in [3,10] deg.)");
  TCanvas* cResMom_0_2_DegMC = DrawResMomVsMom("cResMom_0_2_DegMC", "momentum residuals for tracks with MC angle < 2 degrees", hResMomVertexVsMom_0_2_DegMC, 5);
  
  // diplay slopeX residuals
  TCanvas* cResSlopeX = DrawVsAng("cResSlopeX", "slope_{X} residual at vertex in 3 angular regions", hResSlopeXVertex, hResSlopeXVertexVsAngle);
  TCanvas* cResSlopeXMC = DrawVsAng("cResSlopeXMC", "slope_{X} residual at vertex in 3 MC angular regions", hResSlopeXVertex, hResSlopeXVertexVsMCAngle);
  TCanvas* cResSlopeXVsPos = DrawVsPos("cResSlopeXVsPos", "slope_{X} residual at vertex versus position at absorber end in 3 MC angular regions",
				       hResSlopeXVertexVsPosAbsEnd_0_2_DegMC, hResSlopeXVertexVsPosAbsEnd_2_3_DegMC, hResSlopeXVertexVsPosAbsEnd_3_10_DegMC);
  
  // diplay slopeY residuals
  TCanvas* cResSlopeY = DrawVsAng("cResSlopeY", "slope_{Y} residual at vertex in 3 angular regions", hResSlopeYVertex, hResSlopeYVertexVsAngle);
  TCanvas* cResSlopeYMC = DrawVsAng("cResSlopeYMC", "slope_{Y} residual at vertex in 3 MC angular regions", hResSlopeYVertex, hResSlopeYVertexVsMCAngle);
  TCanvas* cResSlopeYVsPos = DrawVsPos("cResSlopeYVsPos", "slope_{Y} residual at vertex versus position at absorber end in 3 MC angular regions",
				       hResSlopeYVertexVsPosAbsEnd_0_2_DegMC, hResSlopeYVertexVsPosAbsEnd_2_3_DegMC, hResSlopeYVertexVsPosAbsEnd_3_10_DegMC);
  
  // diplay P*DCA
  TCanvas* cPDCA = DrawVsAng("cPDCA", "p #times DCA in 3 angular regions", hPDCA, hPDCAVsAngle);
  TCanvas* cPDCAMC = DrawVsAng("cPDCAMC", "p #times DCA in 3 MC angular regions", hPDCA, hPDCAVsMCAngle);
  TCanvas* cPDCAVsPos = DrawVsPos("cPDCAVsPos", "p #times DCA versus position at absorber end in 3 MC angular regions",
				  hPDCAVsPosAbsEnd_0_2_DegMC, hPDCAVsPosAbsEnd_2_3_DegMC, hPDCAVsPosAbsEnd_3_10_DegMC);
  
  // diplay eta residuals
  TCanvas* cResEta = DrawVsAng("cResEta", "eta residual at vertex in 3 angular regions", hResEtaVertex, hResEtaVertexVsAngle);
  TCanvas* cResEtaMC = DrawVsAng("cResEtaMC", "eta residual at vertex in 3 MC angular regions", hResEtaVertex, hResEtaVertexVsMCAngle);
  TCanvas* cResEtaVsPos = DrawVsPos("cResEtaVsPos", "eta residual at vertex versus position at absorber end in 3 MC angular regions",
				    hResEtaVertexVsPosAbsEnd_0_2_DegMC, hResEtaVertexVsPosAbsEnd_2_3_DegMC, hResEtaVertexVsPosAbsEnd_3_10_DegMC);
  
  // diplay phi residuals
  TCanvas* cResPhi = DrawVsAng("cResPhi", "phi residual at vertex in 3 angular regions", hResPhiVertex, hResPhiVertexVsAngle);
  TCanvas* cResPhiMC = DrawVsAng("cResPhiMC", "phi residual at vertex in 3 MC angular regions", hResPhiVertex, hResPhiVertexVsMCAngle);
  TCanvas* cResPhiVsPos = DrawVsPos("cResPhiVsPos", "phi residual at vertex versus position at absorber end in 3 MC angular regions",
				    hResPhiVertexVsPosAbsEnd_0_2_DegMC, hResPhiVertexVsPosAbsEnd_2_3_DegMC, hResPhiVertexVsPosAbsEnd_3_10_DegMC);
  
  // ###################################### save histogram ###################################### //
  histoFile->Write();
  
  histoFile->cd("momentumAtVertex");
  gMeanResMomVertexVsMom->Write();
  gMostProbResMomVertexVsMom->Write();
  gSigmaResMomVertexVsMom->Write();
  cResMom->Write();
  cResMomMC->Write();
  cResMomVsPos->Write();
  cResMom_2_3_Deg->Write();
  cResMom_3_10_Deg->Write();
  cResMom_0_2_DegMC->Write();
  
  histoFile->cd("slopesAtVertex");
  gMeanResSlopeXVertexVsMom->Write();
  gMeanResSlopeYVertexVsMom->Write();
  gSigmaResSlopeXVertexVsMom->Write();
  gSigmaResSlopeYVertexVsMom->Write();
  cResSlopeX->Write();
  cResSlopeY->Write();
  cResSlopeXMC->Write();
  cResSlopeYMC->Write();
  cResSlopeXVsPos->Write();
  cResSlopeYVsPos->Write();
  
  histoFile->cd("DCA");
  gSigmaPDCAVsMom_2_3_Deg->Write();
  gSigmaPDCAVsMom_3_10_Deg->Write();
  gMeanPMCSAngVsMom_2_3_Deg->Write();
  gSigmaPMCSAngVsMom_2_3_Deg->Write();
  gMeanPMCSAngVsMom_3_10_Deg->Write();
  gSigmaPMCSAngVsMom_3_10_Deg->Write();
  cPDCA->Write();
  cPDCAMC->Write();
  cPDCAVsPos->Write();
  
  histoFile->cd("etaAtVertex");
  gMeanResEtaVertexVsMom->Write();
  gSigmaResEtaVertexVsMom->Write();
  cResEta->Write();
  cResEtaMC->Write();
  cResEtaVsPos->Write();
  
  histoFile->cd("phiAtVertex");
  gMeanResPhiVertexVsMom->Write();
  gSigmaResPhiVertexVsMom->Write();
  cResPhi->Write();
  cResPhiMC->Write();
  cResPhiVsPos->Write();
  
  histoFile->cd("momentumAtFirstCluster");
  gMeanResMomFirstClusterVsMom->Write();
  gSigmaResMomFirstClusterVsMom->Write();
  
  histoFile->cd("slopesAtFirstCluster");
  gMeanResSlopeXFirstClusterVsMom->Write();
  gMeanResSlopeYFirstClusterVsMom->Write();
  gSigmaResSlopeXFirstClusterVsMom->Write();
  gSigmaResSlopeYFirstClusterVsMom->Write();
  
  histoFile->cd("clusters");
  gResidualXPerChMean->Write();
  gResidualXPerChSigma->Write();
  gResidualYPerChMean->Write();
  gResidualYPerChSigma->Write();
  
  histoFile->Close();
  
  // ###################################### clean memory ###################################### //
  delete cResMom;
  delete cResMomMC;
  delete cResMomVsPos;
  delete cResMom_2_3_Deg;
  delete cResMom_3_10_Deg;
  delete cResMom_0_2_DegMC;
  delete cResSlopeX;
  delete cResSlopeY;
  delete cResSlopeXMC;
  delete cResSlopeYMC;
  delete cResSlopeXVsPos;
  delete cResSlopeYVsPos;
  delete cPDCA;
  delete cPDCAMC;
  delete cPDCAVsPos;
  delete cResEta;
  delete cResEtaMC;
  delete cResEtaVsPos;
  delete cResPhi;
  delete cResPhiMC;
  delete cResPhiVsPos;
  
  // ###################################### print statistics ###################################### //
  printf("\n");
  printf("nb of reconstructible tracks: %d \n", nReconstructibleTracks);
  printf("nb of reconstructed tracks: %d \n", nReconstructedTracks);
  printf("nb of reconstructible tracks which are reconstructed: %d \n", nReconstructibleTracksCheck);
  
  aMCSMoy /= (Double_t) nMCS;
  aMCS2Moy /= (Double_t) nMCS;
  dMCSMoy /= (Double_t) nMCS;
  dMCS2Moy /= (Double_t) nMCS;
  adMCSMoy /= (Double_t) nMCS;
  Double_t sigma2_ThetaMCS = aMCS2Moy - aMCSMoy*aMCSMoy;
  Double_t sigma2_PosMCS = dMCS2Moy - dMCSMoy*dMCSMoy;
  Double_t cov_ThetaPosMCS = - (adMCSMoy - aMCSMoy*dMCSMoy);
  printf("\nmultiple scattering of tracks between 3 and 10 deg. at absorber end:\n");
  printf(" sigma_ThetaMCS = %f\n", TMath::Sqrt(sigma2_ThetaMCS));
  printf(" sigma_PosMCS = %f\n", TMath::Sqrt(sigma2_PosMCS));
  printf(" cov_ThetaPosMCS = %f\n", cov_ThetaPosMCS);
  printf(" --> sigma_DCA = %f\n", TMath::Sqrt(AliMUONConstants::AbsZEnd()*AliMUONConstants::AbsZEnd()*sigma2_ThetaMCS
					      - 2.*AliMUONConstants::AbsZEnd()*cov_ThetaPosMCS + sigma2_PosMCS));
  printf("\n");
}

//------------------------------------------------------------------------------------
Double_t langaufun(Double_t *x, Double_t *par) {
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation), 
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0; // number of convolution steps
  Double_t sc = 5.0;   // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    //change x -> -x because the tail of the Landau is at the left here...
    fland = TMath::Landau(-xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    //change x -> -x because the tail of the Landau is at the left here...
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(-xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}

//------------------------------------------------------------------------------------
void FitGausResVsMom(TH2* h, Int_t nBins, const Double_t mean0, const Double_t sigma0,
		     const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma)
{
  /// generic function to fit residuals versus momentum with a gaussian
  static TF1* fGaus = 0x0;
  if (!fGaus) fGaus = new TF1("fGaus","gaus");
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/nBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,nBins)<<flush;
    TH1D *tmp = h->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    fGaus->SetParameters(tmp->GetEntries(), mean0, sigma0);
    tmp->Fit("fGaus","WWNQ");
    Int_t rebin = TMath::Max(Int_t(0.5*fGaus->GetParameter(2)/tmp->GetBinWidth(1)),1);
    while (tmp->GetNbinsX()%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    tmp->Fit("fGaus","NQ");
    h->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = h->GetMean();
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetBinLowEdge(i-rebinFactorX+1), h->GetBinLowEdge(i+1)-p};
    gMean->SetPoint(i/rebinFactorX-1, p, fGaus->GetParameter(1));
    gMean->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fGaus->GetParError(1), fGaus->GetParError(1));
    gSigma->SetPoint(i/rebinFactorX-1, p, fGaus->GetParameter(2));
    gSigma->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fGaus->GetParError(2), fGaus->GetParError(2));
    delete tmp;
  }
  cout<<Form("\rFitting %s... %d/%d",fitting,nBins,nBins)<<endl;
}

//------------------------------------------------------------------------------------
void FitPDCAVsMom(TH2* h, Int_t nBins, const char* fitting, TGraphAsymmErrors* gSigma)
{
  /// generic function to fit p*DCA distributions
  static TF1* fPGaus = 0x0;
  if (!fPGaus) fPGaus = new TF1("fPGaus","x*gaus");
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/nBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,nBins)<<flush;
    TH1D *tmp = h->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    fPGaus->SetParameters(1.,-100.,100.);
    Int_t rebin = 50.*(tmp->GetNbinsX()/(tmp->GetBinLowEdge(tmp->GetNbinsX()+1)-tmp->GetBinLowEdge(1)));
    while (tmp->GetNbinsX()%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    tmp->Fit("fPGaus","NQ");
    h->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = h->GetMean();
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetBinLowEdge(i-rebinFactorX+1), h->GetBinLowEdge(i+1)-p};
    gSigma->SetPoint(i/rebinFactorX-1, p, fPGaus->GetParameter(2));
    gSigma->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fPGaus->GetParError(2), fPGaus->GetParError(2));
    delete tmp;
  }
  cout<<Form("\rFitting %s... %d/%d",fitting,nBins,nBins)<<endl;
}

//------------------------------------------------------------------------------------
TCanvas* DrawVsAng(const char* name, const char* title, TH1* h1, TH2* h2)
{
  /// generic function to draw histograms versus absorber angular region
  TCanvas* c = new TCanvas(name, title);
  c->cd();
  h1->Draw();
  TH1D *proj1 = h2->ProjectionY(Form("%s_proj_0_2",h2->GetName()),1,2);
  proj1->Draw("sames");
  proj1->SetLineColor(2);
  TH1D *proj2 = h2->ProjectionY(Form("%s_proj_2_3",h2->GetName()),3,3);
  proj2->Draw("sames");
  proj2->SetLineColor(4);
  TH1D *proj3 = h2->ProjectionY(Form("%s__proj_3_10",h2->GetName()),4,10);
  proj3->Draw("sames");
  proj3->SetLineColor(3);
  return c;
}

//------------------------------------------------------------------------------------
TCanvas* DrawVsPos(const char* name, const char* title, TH2* h1, TH2* h2, TH2* h3)
{
  /// generic function to draw histograms versus position at absorber end
  TCanvas* c = new TCanvas(name, title);
  c->cd();
  h1->Draw();
  h1->SetMarkerColor(2);
  h2->Draw("sames");
  h2->SetMarkerColor(4);
  h3->Draw("sames");
  h3->SetMarkerColor(3);
  return c;
}

//------------------------------------------------------------------------------------
TCanvas* DrawResMomVsMom(const char* name, const char* title, TH2* h, Int_t nBins, TF1* f2, const char* fitting)
{
  /// generic function to draw and eventually fit momentum residuals versus momentum
  TLegend* l = new TLegend(0.15,0.25,0.3,0.85);
  TCanvas* c = new TCanvas(name, title);
  c->cd();
  TH1D* proj = 0x0;
  h->Sumw2();
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/nBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    if (f2) cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,nBins)<<flush;
    proj = h->ProjectionY(Form("%s_%d",h->GetName(),i/rebinFactorX),i-rebinFactorX+1,i);
    if (proj->GetEntries() > 0) proj->Scale(1./proj->GetEntries());
    proj->Draw((i==rebinFactorX)?"hist":"histsames");
    proj->SetLineColor(i/rebinFactorX);
    if (f2) {
      f2->SetParameters(0.2,0.,1.,1.);
      f2->SetLineColor(i/rebinFactorX);
      proj->Fit("f2","WWNQ","sames");
      Double_t fwhm = f2->GetParameter(0);
      Double_t sigma = f2->GetParameter(3);
      Double_t sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
      Int_t rebin = TMath::Max(Int_t(0.5*sigmaP/proj->GetBinWidth(1)),1);
      while (proj->GetNbinsX()%rebin!=0) rebin--;
      proj->Rebin(rebin);
      proj->Scale(1./rebin);
      proj->Fit("f2","Q","sames");
    } else proj->SetLineWidth(2);
    Double_t p = 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    l->AddEntry(proj,Form("%5.1f GeV",p));
  }
  if (f2) cout<<Form("\rFitting %s... %d/%d",fitting,nBins,nBins)<<endl;
  l->Draw("same");
  return c;
}

