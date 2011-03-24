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

/*
  //
  // FUNCTIONALITY:
  //
  // 1. The laser track is associated with the mirror
  //    see function FindMirror
  //
  // 2. The laser track is accepted for the analysis under certain condition
  //    (see function Accpet laser)
  //
  // 3. The drift velocity and jitter is calculated event by event
  //    (see function drift velocity)
  //
  // 4. The tracks are refitted at different sectors 
  //    Fit model 
  //      4.a) line
  //      4.b) parabola
  //      4.c) parabola with common P2 for inner and outer 
  //
  // To make laser scan the user interaction neccessary
  //
  .x ~/NimStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  TFile fcalib("CalibObjectsTrain2.root");
  AliTPCcalibLaser * laser = ( AliTPCcalibLaser *)fcalib->Get("laserTPC");
  laser->DumpMeanInfo(run)
  TFile fmean("laserMean.root")
  //
  //  laser track clasification;
  //
  TCut cutT("cutT","abs(Tr.fP[3])<0.06");
  TCut cutN("cutN","fTPCncls>70");
  TCut cutP("cutP","abs(atan2(x1,x0)-atan2(lx1,lx0))<0.03")
  TCut cutA = cutT+cutPt+cutP;
  TChain * chainTrL = tool.MakeChain("laser.txt","Track",0,10200);

  //
  //
  // Analyze  LASER scan 
  //
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  TChain * chain = tool.MakeChain("laserScan.txt","Mean",0,10200);
  chain->Lookup();
  AliTPCcalibLaser::DumpScanInfo(chain)
  TFile fscan("laserScan.root")
  TTree * treeT = (TTree*)fscan.Get("Mean")
  //
  // Analyze laser 
  //
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;
  AliXRDPROOFtoolkit::FilterList("laserDebug.list","* driftvN",1) 
  TChain * chainDrift = tool.MakeChainRandom("laser.txt.Good","driftv",0,50);
  chainDrift->Lookup();
  TChain * chainDriftN = tool.MakeChainRandom("laserDebug.list.Good","driftvN",0,300);
  chainDriftN->Lookup();


  TChain * chain = tool.MakeChain("laser.txt","Residuals",0,10200);
  chain->Lookup();
  TChain * chainFit = tool.MakeChain("laser.txt","FitModels",0,10200);
  chainFit->Lookup();
  TChain * chainTrack = tool.MakeChainRandom("laser.txt","Track",0,30);
  chainTrack->Lookup();

*/



#include "TLinearFitter.h"
#include "AliTPCcalibLaser.h"
#include "AliExternalTrackParam.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDtrack.h"
#include "AliTPCTracklet.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TVectorD.h"
#include "TTreeStream.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTracker.h"
#include "AliLog.h"
#include "TClonesArray.h"
#include "TPad.h"
#include "TSystem.h"
#include "TCut.h"
#include "TChain.h"
#include "TH2F.h"
#include "TStatToolkit.h"
#include "TROOT.h"
#include "TDatabasePDG.h"


#include "TTreeStream.h"
#include <iostream>
#include <sstream>
#include "AliTPCLaserTrack.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "TTimeStamp.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"
#include "AliGRPObject.h"
#include "AliTPCROC.h"

using namespace std;

ClassImp(AliTPCcalibLaser)

AliTPCcalibLaser::AliTPCcalibLaser():
  AliTPCcalibBase(),
  fESD(0),
  fESDfriend(0),
  fNtracks(0),
  fTracksMirror(336),
  fTracksEsd(336),
  fTracksEsdParam(336),
  fTracksTPC(336),
  fFullCalib(kTRUE),
  fDeltaZ(336),
  fDeltaP3(336),
  fDeltaP4(336),
  fDeltaPhi(336),
  fDeltaPhiP(336),
  fSignals(336),
  //
  fHisLaser(0),      //  N dim histogram of laser
  fHisLaserPad(0),      //  N dim histogram of laser
  fHisLaserTime(0),      //  N dim histogram of laser
  fHisNclIn(0),      //->Number of clusters inner
  fHisNclOut(0),     //->Number of clusters outer
  fHisNclIO(0),      //->Number of cluster inner outer
  fHisLclIn(0),      //->Level arm inner
  fHisLclOut(0),     //->Level arm outer
  fHisLclIO(0),      //->Number of cluster inner outer
  fHisdEdx(0),     //->dEdx histo
  fHisdZfit(0),     //->distance to the mirror after linear fit
  //
  //
  fHisChi2YIn1(0),      //->chi2 y inner - line
  fHisChi2YOut1(0),     //->chi2 y inner - line
  fHisChi2YIn2(0),      //->chi2 y inner - parabola
  fHisChi2YOut2(0),     //->chi2 y inner - parabola
  fHisChi2YIO1(0),      //->chi2 y IO    - common
  fHisChi2ZIn1(0),      //->chi2 z inner - line
  fHisChi2ZOut1(0),     //->chi2 z inner - line
  fHisChi2ZIn2(0),      //->chi2 z inner - parabola
  fHisChi2ZOut2(0),     //->chi2 z inner - parabola
  fHisChi2ZIO1(0),      //->chi2 z IO    - common
  //
  //
  fHisPy1vP0(0),     //-> delta y   P0outer-P0inner - line
  fHisPy2vP0(0),     //-> delta y   P0outer-P0inner - parabola
  fHisPy3vP0(0),     //-> delta y   P0outer-P0inner - common parabola
  fHisPy1vP1(0),     //-> delta ky  P1outer-P1inner - line
  fHisPy2vP1(0),     //-> delta ky  P1outer-P1inner - parabola
  fHisPy3vP1(0),     //-> delta ky  P1outer-P1inner - common parabola
  fHisPy2vP2In(0),   //-> Curv  P2inner - parabola
  fHisPy2vP2Out(0),  //-> Curv  P2outer - parabola
  fHisPy3vP2IO(0),   //-> Curv  P2outerinner - common parabola
  //
  //
  fHisPz1vP0(0),     //-> delta z   P0outer-P0inner - line
  fHisPz2vP0(0),     //-> delta z   P0outer-P0inner - parabola
  fHisPz3vP0(0),     //-> delta z   P0outer-P0inner - common parabola
  fHisPz1vP1(0),     //-> delta kz  P1outer-P1inner - line
  fHisPz2vP1(0),     //-> delta kz  P1outer-P1inner - parabola
  fHisPz3vP1(0),     //-> delta kz  P1outer-P1inner - common parabola
  fHisPz2vP2In(0),   //-> Curv  P2inner - parabola
  fHisPz2vP2Out(0),  //-> Curv  P2outer - parabola
  fHisPz3vP2IO(0),   //-> Curv  P2outerinner - common parabola
  //
  fDeltaYres(336),   //->2D histo of residuals
  fDeltaZres(336),   //->2D histo fo residuals
  fDeltaYres2(336),   //->2D histo of residuals
  fDeltaZres2(336),   //->2D histo fo residuals
  fDeltaYresAbs(336), //->2D histo of residuals
  fHisYAbsErrors(0),  //-> total errors per beam in the abs res y analysis
  fDeltaZresAbs(336), //->2D histo of residuals
  fHisZAbsErrors(0),  //-> total errors per beam in the abs res z analysis
  //fDeltaYres3(336),   //->2D histo of residuals
  //fDeltaZres3(336),   //->2D histo fo residuals
  fFitAside(new TVectorD(5)),
  fFitCside(new TVectorD(5)),      
  fFitACside(new TVectorD(6)),      
  fEdgeXcuts(3),    
  fEdgeYcuts(3),    
  fNClCuts(5),      
  fNcuts(0),
  fBeamSectorOuter(336),
  fBeamSectorInner(336),
  fBeamOffsetYOuter(336),
  fBeamSlopeYOuter(336),
  fBeamOffsetYInner(336),
  fBeamSlopeYInner(336),
  fBeamOffsetZOuter(336),
  fBeamSlopeZOuter(336),
  fBeamOffsetZInner(336),
  fBeamSlopeZInner(336),
  fInverseSlopeZ(kTRUE),
  fUseFixedDriftV(0),
  fFixedFitAside0(0.0),
  fFixedFitAside1(1.0),
  fFixedFitCside0(0.0),
  fFixedFitCside1(1.0)
{
  //
  // Constructor
  //
  fTracksEsdParam.SetOwner(kTRUE);
  for (Int_t i=0; i<336; i++) {
    fFitZ[i]=0;
    fCounter[i]=0;    //! counter of usage
    fClusterCounter[i]=0; //!couter of clusters in "sensitive are"
    fClusterSatur[i]=0;   //!couter of saturated clusters in "sensitive are"
  }
}

AliTPCcalibLaser::AliTPCcalibLaser(const Text_t *name, const Text_t *title, Bool_t full):
  AliTPCcalibBase(),
  fESD(0),
  fESDfriend(0),
  fNtracks(0),
  fTracksMirror(336),
  fTracksEsd(336),
  fTracksEsdParam(336),
  fTracksTPC(336),
  fFullCalib(full),
  //
  fDeltaZ(336),          // array of histograms of delta z for each track
  fDeltaP3(336),          // array of histograms of delta z for each track
  fDeltaP4(336),          // array of histograms of P3 for each track
  fDeltaPhi(336),          // array of histograms of P4 for each track
  fDeltaPhiP(336),          // array of histograms of delta z for each track
  fSignals(336),           // array of dedx signals
  //
  //
  fHisLaser(0),      //  N dim histogram of laser
  fHisLaserPad(0),      //  N dim histogram of laser
  fHisLaserTime(0),      //  N dim histogram of laser

  fHisNclIn(0),      //->Number of clusters inner
  fHisNclOut(0),     //->Number of clusters outer
  fHisNclIO(0),      //->Number of cluster inner outer
  fHisLclIn(0),      //->Level arm inner
  fHisLclOut(0),     //->Level arm outer
  fHisLclIO(0),      //->Number of cluster inner outer
  fHisdEdx(0), //->dEdx histo  
  fHisdZfit(0),     //->distance to the mirror after linear fit
  //
  //
  fHisChi2YIn1(0),      //->chi2 y inner - line
  fHisChi2YOut1(0),     //->chi2 y inner - line
  fHisChi2YIn2(0),      //->chi2 y inner - parabola
  fHisChi2YOut2(0),     //->chi2 y inner - parabola
  fHisChi2YIO1(0),      //->chi2 y IO    - common
  fHisChi2ZIn1(0),      //->chi2 z inner - line
  fHisChi2ZOut1(0),     //->chi2 z inner - line
  fHisChi2ZIn2(0),      //->chi2 z inner - parabola
  fHisChi2ZOut2(0),     //->chi2 z inner - parabola
  fHisChi2ZIO1(0),      //->chi2 z IO    - common
  //
  //
  fHisPy1vP0(0),     //-> delta y   P0outer-P0inner - line
  fHisPy2vP0(0),     //-> delta y   P0outer-P0inner - parabola
  fHisPy3vP0(0),     //-> delta y   P0outer-P0inner - common parabola
  fHisPy1vP1(0),     //-> delta ky  P1outer-P1inner - line
  fHisPy2vP1(0),     //-> delta ky  P1outer-P1inner - parabola
  fHisPy3vP1(0),     //-> delta ky  P1outer-P1inner - common parabola
  fHisPy2vP2In(0),   //-> Curv  P2inner - parabola
  fHisPy2vP2Out(0),  //-> Curv  P2outer - parabola
  fHisPy3vP2IO(0),   //-> Curv  P2outerinner - common parabola
  //
  //
  fHisPz1vP0(0),     //-> delta z   P0outer-P0inner - line
  fHisPz2vP0(0),     //-> delta z   P0outer-P0inner - parabola
  fHisPz3vP0(0),     //-> delta z   P0outer-P0inner - common parabola
  fHisPz1vP1(0),     //-> delta kz  P1outer-P1inner - line
  fHisPz2vP1(0),     //-> delta kz  P1outer-P1inner - parabola
  fHisPz3vP1(0),     //-> delta kz  P1outer-P1inner - common parabola
  fHisPz2vP2In(0),   //-> Curv  P2inner - parabola
  fHisPz2vP2Out(0),  //-> Curv  P2outer - parabola
  fHisPz3vP2IO(0),   //-> Curv  P2outerinner - common parabola
  //
  //
  //
  fDeltaYres(336),
  fDeltaZres(336),  
  fDeltaYres2(336),
  fDeltaZres2(336),  
  fDeltaYresAbs(336),
  fHisYAbsErrors(0),
  fDeltaZresAbs(336),
  fHisZAbsErrors(0),
  //  fDeltaYres3(336),
  //fDeltaZres3(336),  
  fFitAside(new TVectorD(5)),        // drift fit - A side
  fFitCside(new TVectorD(5)),        // drift fit - C- side
  fFitACside(new TVectorD(6)),        // drift fit - AC- side
  fEdgeXcuts(3),       // cuts in local x direction; used in the refit of the laser tracks
  fEdgeYcuts(3),       // cuts in local y direction; used in the refit of the laser tracks
  fNClCuts(5),         // cuts on the number of clusters per tracklet; used in the refit of the laser tracks
  fNcuts(0),           // number of cuts
  fBeamSectorOuter(336),
  fBeamSectorInner(336),
  fBeamOffsetYOuter(336),
  fBeamSlopeYOuter(336),
  fBeamOffsetYInner(336),
  fBeamSlopeYInner(336),
  fBeamOffsetZOuter(336),
  fBeamSlopeZOuter(336),
  fBeamOffsetZInner(336),
  fBeamSlopeZInner(336),
  fInverseSlopeZ(kTRUE),
  fUseFixedDriftV(0),
  fFixedFitAside0(0.0),
  fFixedFitAside1(1.0),
  fFixedFitCside0(0.0),
  fFixedFitCside1(1.0)
{
  SetName(name);
  SetTitle(title);
  //
  // Constructor
  //
  fTracksEsdParam.SetOwner(kTRUE);
  for (Int_t i=0; i<336; i++) {
    fFitZ[i]=0;
    fCounter[i]=0;    //! counter of usage
    fClusterCounter[i]=0; //!couter of clusters in "sensitive are"
    fClusterSatur[i]=0;   //!couter of saturated clusters in "sensitive are"
  }
}

AliTPCcalibLaser::AliTPCcalibLaser(const AliTPCcalibLaser& calibLaser):
  AliTPCcalibBase(calibLaser), 
  fESD(0),
  fESDfriend(0),
  fNtracks(0),
  fTracksMirror(336),
  fTracksEsd(336),
  fTracksEsdParam(336),
  fTracksTPC(336),
  fFullCalib(calibLaser.fFullCalib),
  //
  fDeltaZ(calibLaser.fDeltaZ),          // array of histograms of delta z for each track
  fDeltaP3(((calibLaser.fDeltaP3))),          // array of histograms of delta z for each track
  fDeltaP4(((calibLaser.fDeltaP4))),          // array of histograms of P3 for each track
  fDeltaPhi(((calibLaser.fDeltaPhi))),          // array of histograms of P4 for each track
  fDeltaPhiP(((calibLaser.fDeltaPhiP))),          // array of histograms of delta z for each track
  fSignals(((calibLaser.fSignals))),           // array of dedx signals
  //
  //
  fHisLaser(0),      //  N dim histogram of laser
  fHisLaserPad(0),      //  N dim histogram of laser
  fHisLaserTime(0),      //  N dim histogram of laser

  fHisNclIn(new TH2F(*(calibLaser.fHisNclIn))),      //->Number of clusters inner
  fHisNclOut(new TH2F(*(calibLaser.fHisNclOut))),     //->Number of clusters outer
  fHisNclIO(new TH2F(*(calibLaser.fHisNclIO))),      //->Number of cluster inner outer
  fHisLclIn(new TH2F(*(calibLaser.fHisLclIn))),      //->Level arm inner
  fHisLclOut(new TH2F(*(calibLaser.fHisLclOut))),     //->Level arm outer
  fHisLclIO(new TH2F(*(calibLaser.fHisLclIO))),      //->Number of cluster inner outer
  fHisdEdx(new TH2F(*(calibLaser.fHisdEdx))),      //
  fHisdZfit(new TH2F(*(calibLaser.fHisdZfit))),     //->distance to the mirror after linear fit
  //
  //
  fHisChi2YIn1(new TH2F(*(calibLaser.fHisChi2YIn1))),      //->chi2 y inner - line
  fHisChi2YOut1(new TH2F(*(calibLaser.fHisChi2YOut1))),     //->chi2 y inner - line
  fHisChi2YIn2(new TH2F(*(calibLaser.fHisChi2YIn2))),      //->chi2 y inner - parabola
  fHisChi2YOut2(new TH2F(*(calibLaser.fHisChi2YOut2))),     //->chi2 y inner - parabola
  fHisChi2YIO1(new TH2F(*(calibLaser.fHisChi2YIO1))),      //->chi2 y IO    - common
  fHisChi2ZIn1(new TH2F(*(calibLaser.fHisChi2ZIn1))),      //->chi2 z inner - line
  fHisChi2ZOut1(new TH2F(*(calibLaser.fHisChi2ZOut1))),     //->chi2 z inner - line
  fHisChi2ZIn2(new TH2F(*(calibLaser.fHisChi2ZIn2))),      //->chi2 z inner - parabola
  fHisChi2ZOut2(new TH2F(*(calibLaser.fHisChi2ZOut2))),     //->chi2 z inner - parabola
  fHisChi2ZIO1(new TH2F(*(calibLaser.fHisChi2ZIO1))),      //->chi2 z IO    - common
  //
  //
  fHisPy1vP0(new TH2F(*(calibLaser.fHisPy1vP0))),     //-> delta y   P0outer-P0inner - line
  fHisPy2vP0(new TH2F(*(calibLaser.fHisPy2vP0))),     //-> delta y   P0outer-P0inner - parabola
  fHisPy3vP0(new TH2F(*(calibLaser.fHisPy3vP0))),     //-> delta y   P0outer-P0inner - common parabola
  fHisPy1vP1(new TH2F(*(calibLaser.fHisPy1vP1))),     //-> delta ky  P1outer-P1inner - line
  fHisPy2vP1(new TH2F(*(calibLaser.fHisPy2vP1))),     //-> delta ky  P1outer-P1inner - parabola
  fHisPy3vP1(new TH2F(*(calibLaser.fHisPy3vP1))),     //-> delta ky  P1outer-P1inner - common parabola
  fHisPy2vP2In(new TH2F(*(calibLaser.fHisPy2vP2In))),   //-> Curv  P2inner - parabola
  fHisPy2vP2Out(new TH2F(*(calibLaser.fHisPy2vP2Out))),  //-> Curv  P2outer - parabola
  fHisPy3vP2IO(new TH2F(*(calibLaser.fHisPy3vP2IO))),   //-> Curv  P2outerinner - common parabola
  //
  //
  fHisPz1vP0(new TH2F(*(calibLaser.fHisPz1vP0))),     //-> delta z   P0outer-P0inner - line
  fHisPz2vP0(new TH2F(*(calibLaser.fHisPz2vP0))),     //-> delta z   P0outer-P0inner - parabola
  fHisPz3vP0(new TH2F(*(calibLaser.fHisPz3vP0))),     //-> delta z   P0outer-P0inner - common parabola
  fHisPz1vP1(new TH2F(*(calibLaser. fHisPz1vP1))),     //-> delta kz  P1outer-P1inner - line
  fHisPz2vP1(new TH2F(*(calibLaser.fHisPz2vP1))),     //-> delta kz  P1outer-P1inner - parabola
  fHisPz3vP1(new TH2F(*(calibLaser.fHisPz3vP1))),     //-> delta kz  P1outer-P1inner - common parabola
  fHisPz2vP2In(new TH2F(*(calibLaser.fHisPz2vP2In))),   //-> Curv  P2inner - parabola
  fHisPz2vP2Out(new TH2F(*(calibLaser. fHisPz2vP2Out))),  //-> Curv  P2outer - parabola
  fHisPz3vP2IO(new TH2F(*(calibLaser.fHisPz3vP2IO))),   //-> Curv  P2outerinner - common parabola
  //
  //
  fDeltaYres(((calibLaser.fDeltaYres))),
  fDeltaZres(((calibLaser.fDeltaZres))),  
  fDeltaYres2(((calibLaser.fDeltaYres))),
  fDeltaZres2(((calibLaser.fDeltaZres))),  
  fDeltaYresAbs(((calibLaser.fDeltaYresAbs))),
  fHisYAbsErrors(new TH1F(*(calibLaser.fHisYAbsErrors))),
  fDeltaZresAbs(((calibLaser.fDeltaZresAbs))),
  fHisZAbsErrors(new TH1F(*(calibLaser.fHisZAbsErrors))),
  //  fDeltaYres3(((calibLaser.fDeltaYres))),
  //fDeltaZres3(((calibLaser.fDeltaZres))),  
  fFitAside(new TVectorD(5)),        // drift fit - A side
  fFitCside(new TVectorD(5)),        // drift fit - C- side
  fFitACside(new TVectorD(6)),        // drift fit - C- side
  fEdgeXcuts(3),       // cuts in local x direction; used in the refit of the laser tracks
  fEdgeYcuts(3),       // cuts in local y direction; used in the refit of the laser tracks
  fNClCuts(5),         // cuts on the number of clusters per tracklet; used in the refit of the laser tracks
  fNcuts(0),           // number of cuts
  fBeamSectorOuter(336),
  fBeamSectorInner(336),
  fBeamOffsetYOuter(336),
  fBeamSlopeYOuter(336),
  fBeamOffsetYInner(336),
  fBeamSlopeYInner(336),
  fBeamOffsetZOuter(336),
  fBeamSlopeZOuter(336),
  fBeamOffsetZInner(336),
  fBeamSlopeZInner(336),
  fInverseSlopeZ(calibLaser.fInverseSlopeZ),
  fUseFixedDriftV(calibLaser.fUseFixedDriftV),
  fFixedFitAside0(calibLaser.fFixedFitAside0),
  fFixedFitAside1(calibLaser.fFixedFitAside1),
  fFixedFitCside0(calibLaser.fFixedFitCside0),
  fFixedFitCside1(calibLaser.fFixedFitCside1)
{
  //
  // copy constructor
  //
  for (Int_t i=0; i<336; i++) {
    fFitZ[i]=0;
    fCounter[i]=0;    //! counter of usage
    fClusterCounter[i]=0; //!couter of clusters in "sensitive are"
    fClusterSatur[i]=0;   //!couter of saturated clusters in "sensitive are"
  }
}



AliTPCcalibLaser & AliTPCcalibLaser::operator=(const AliTPCcalibLaser& calibLaser){
  //
  // assgnment operator
  //
  if (this != &calibLaser) {
    new (this) AliTPCcalibLaser(calibLaser);
  }
  return *this;

}




AliTPCcalibLaser::~AliTPCcalibLaser() {
  //
  // destructor
  //
  if ( fHisNclIn){
    delete fHisLaser;      //->
    delete fHisLaserPad;      //->
    delete fHisLaserTime;      //->

    delete fHisNclIn;      //->Number of clusters inner
    delete fHisNclOut;     //->Number of clusters outer
    delete fHisNclIO;      //->Number of cluster inner outer
    delete fHisLclIn;      //->Level arm inner
    delete fHisLclOut;     //->Level arm outer
    delete fHisLclIO;      //->Number of cluster inner outer
    delete fHisdEdx;
    delete fHisdZfit;     
    //
    //
    delete fHisChi2YIn1;      //->chi2 y inner - line
    delete fHisChi2YOut1;     //->chi2 y inner - line
    delete fHisChi2YIn2;      //->chi2 y inner - parabola
    delete fHisChi2YOut2;     //->chi2 y inner - parabola
    delete fHisChi2YIO1;      //->chi2 y IO    - common
    delete fHisChi2ZIn1;      //->chi2 z inner - line
    delete fHisChi2ZOut1;     //->chi2 z inner - line
    delete fHisChi2ZIn2;      //->chi2 z inner - parabola
    delete fHisChi2ZOut2;     //->chi2 z inner - parabola
    delete fHisChi2ZIO1;      //->chi2 z IO    - common
    //
    //
    delete fHisPy1vP0;     //-> delta y   P0outer-P0inner - line
    delete fHisPy2vP0;     //-> delta y   P0outer-P0inner - parabola
    delete fHisPy3vP0;     //-> delta y   P0outer-P0inner - common parabola
    delete fHisPy1vP1;     //-> delta ky  P1outer-P1inner - line
    delete fHisPy2vP1;     //-> delta ky  P1outer-P1inner - parabola
    delete fHisPy3vP1;     //-> delta ky  P1outer-P1inner - common parabola
    delete fHisPy2vP2In;   //-> Curv  P2inner - parabola
    delete fHisPy2vP2Out;  //-> Curv  P2outer - parabola
    delete fHisPy3vP2IO;   //-> Curv  P2outerinner - common parabola
    //
    //
    delete fHisPz1vP0;     //-> delta z   P0outer-P0inner - line
    delete fHisPz2vP0;     //-> delta z   P0outer-P0inner - parabola
    delete fHisPz3vP0;     //-> delta z   P0outer-P0inner - common parabola
    delete fHisPz1vP1;     //-> delta kz  P1outer-P1inner - line
    delete fHisPz2vP1;     //-> delta kz  P1outer-P1inner - parabola
    delete fHisPz3vP1;     //-> delta kz  P1outer-P1inner - common parabola
    delete fHisPz2vP2In;   //-> Curv  P2inner - parabola
    delete fHisPz2vP2Out;  //-> Curv  P2outer - parabola
    delete fHisPz3vP2IO;   //-> Curv  P2outerinner - common parabola
  }
  //
  //
  //
  fDeltaZ.SetOwner();          //-> array of histograms of delta z for each track
  fDeltaP3.SetOwner();         //-> array of histograms of P3      for each track
  fDeltaP4.SetOwner();         //-> array of histograms of P4      for each track
  fDeltaPhi.SetOwner();        //-> array of histograms of delta z for each track
  fDeltaPhiP.SetOwner();       //-> array of histograms of delta z for each track
  fSignals.SetOwner();         //->Array of dedx signals
  
  fDeltaZ.Delete();          //-> array of histograms of delta z for each track
  fDeltaP3.Delete();         //-> array of histograms of P3      for each track
  fDeltaP4.Delete();         //-> array of histograms of P4      for each track
  fDeltaPhi.Delete();        //-> array of histograms of delta z for each track
  fDeltaPhiP.Delete();       //-> array of histograms of delta z for each track
  fSignals.Delete();         //->Array of dedx signals

  fDeltaYres.SetOwner();
  fDeltaYres.Delete();
  delete fHisYAbsErrors;
  fDeltaZres.SetOwner();
  fDeltaZres.Delete();
  delete fHisZAbsErrors;
  fDeltaYres2.SetOwner();
  fDeltaYres2.Delete();
  fDeltaZres2.SetOwner();
  fDeltaZres2.Delete();
  
  fDeltaYresAbs.SetOwner();
  fDeltaYresAbs.Delete();
  fDeltaZresAbs.SetOwner();
  fDeltaZresAbs.Delete();
}



void AliTPCcalibLaser::Process(AliESDEvent * event) {
  //
  //
  // Loop over tracks and call  Process function
  //
  const Int_t  kMinTracks=20;
  const Int_t  kMinClusters=40;

  fESD = event;
  if (!fESD) {
    return;
  }
  fESDfriend=static_cast<AliESDfriend*>(fESD->FindListObject("AliESDfriend"));
  if (!fESDfriend) {
    return;
  }
  if (fESD->GetNumberOfTracks()<kMinTracks) return; //not enough tracks
  AliDebug(4,Form("Event number in current file: %d",event->GetEventNumberInFile()));
  //
  // find CE background if present
  //
  if (AliTPCLaserTrack::GetTracks()==0) AliTPCLaserTrack::LoadTracks();
  TH1D hisCE("hhisCE","hhisCE",100,-100,100);  
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    AliESDtrack *track=fESD->GetTrack(i);
    if (!track) continue;
    hisCE.Fill(track->GetZ());
    hisCE.Fill(track->GetZ()+2);
    hisCE.Fill(track->GetZ()-2);
  }
  //
  //


  fTracksTPC.Clear();
  fTracksEsd.Clear();
  fTracksEsdParam.Delete();
  for (Int_t id=0; id<336;id++) {
    fCounter[id]=0;
    fClusterCounter[id]=0;
    fClusterSatur[id]=0;
  }
  //
  Int_t n=fESD->GetNumberOfTracks();
  Int_t counter=0;
  for (Int_t i=0;i<n;++i) {
    AliESDfriendTrack *friendTrack=fESDfriend->GetTrack(i);
    if (!friendTrack) continue;
    AliESDtrack *track=fESD->GetTrack(i);
    if (!track) continue;
    Double_t binC = hisCE.GetBinContent(hisCE.FindBin(track->GetZ()));
    if (binC>336) continue; //remove CE background
    TObject *calibObject=0;
    AliTPCseed *seed=0;
    for (Int_t j=0;(calibObject=friendTrack->GetCalibObject(j));++j)
      if ((seed=dynamic_cast<AliTPCseed*>(calibObject)))
	break;
    if (track&&seed&&track->GetTPCNcls()>kMinClusters && seed->GetNumberOfClusters() >kMinClusters) {
      //filter CE tracks
      Int_t id = FindMirror(track,seed);
      if (id>=0) counter++;
    }
    //
  } 
  fNtracks=counter;
  if (counter<kMinTracks) return;

  //FitDriftV();
  FitDriftV(0.2);
  if (!fFullCalib) return;
  static Bool_t init=kFALSE;
  if (!init){
    init = kTRUE;  // way around for PROOF - to be investigated
    MakeFitHistos();
  }
  //
  for (Int_t id=0; id<336; id++){    
    //
    if (!fTracksEsdParam.At(id)) continue;
    if (fClusterSatur[id]>0.3)   continue;  // tracks in saturation
    DumpLaser(id);
    if ( AcceptLaser(id) && fFitZ[id]<0.5){
      RefitLaserJW(id);    
      MakeDistHisto(id);
    }
  }

}

void AliTPCcalibLaser::MakeDistHisto(Int_t id){
  //
  //
  //
  //  for (Int_t id=0; id<336; id++){
    //
    //
    if (!fTracksEsdParam.At(id)) return;
    if (!AcceptLaser(id)) return;
    Double_t xhis[12]={0,0,0,0,0,0,0,0,0,0,0,0};
    //
    //
    TH1F * hisdz = (TH1F*)fDeltaZ.At(id);
    if (!hisdz) MakeFitHistos();
    hisdz = (TH1F*)fDeltaZ.At(id);
    TH1F * hisP3 = (TH1F*)fDeltaP3.At(id);
    TH1F * hisP4 = (TH1F*)fDeltaP4.At(id);
    TH1F * hisdphi = (TH1F*)fDeltaPhi.At(id);
    TH1F * hisdphiP = (TH1F*)fDeltaPhiP.At(id);
    TH1F * hisSignal = (TH1F*)fSignals.At(id);
    //

    AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
    AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
    AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
    if (!param) return;
    if (!ltrp) return;
    if (!track) return;
    Double_t xyz[3];
    Double_t pxyz[3];
    Double_t lxyz[3];
    Double_t lpxyz[3];
    param->GetXYZ(xyz);
    param->GetPxPyPz(pxyz);
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);
    //
    Float_t dz   = param->GetZ()-ltrp->GetZ();
    Float_t dphi = (TMath::ATan2(xyz[1],xyz[0])- TMath::ATan2(lxyz[1],lxyz[0]))*254.;
    Float_t dphiP = param->GetParameter()[2]-ltrp->GetParameter()[2];
    if (hisdz) hisdz->Fill(dz);
    if (hisP3) hisP3->Fill(param->GetParameter()[3]);
    if (hisP4) hisP4->Fill(param->GetParameter()[4]);
    if (hisdphi) hisdphi->Fill(dphi);
    if (hisdphiP) hisdphiP->Fill(dphiP);
    if (hisSignal) hisSignal->Fill(TMath::Sqrt(TMath::Abs(track->GetTPCsignal())));
    // fill HisLaser
    xhis[0] = ltrp->GetId();
    xhis[1] = ltrp->GetSide();
    xhis[2] = ltrp->GetRod();
    xhis[3] = ltrp->GetBundle();
    xhis[4] = ltrp->GetBeam();
    xhis[5] = dphi;
    xhis[6] = fFitZ[id];
    xhis[7] = param->GetParameter()[2]-ltrp->GetParameter()[2]; //dp2
    xhis[8] = param->GetParameter()[3]-ltrp->GetParameter()[3]; //dp3
    xhis[9] = param->GetParameter()[4];
    xhis[10]= track->GetTPCNcls();
    xhis[11]= TMath::Sqrt(TMath::Abs(track->GetTPCsignal()));
    // }
    fHisLaser->Fill(xhis);
    //
    
}

void AliTPCcalibLaser::FitDriftV(){
  //
  // Fit corrections to the drift velocity - linear approximation in the z and global y
  //The transfromatiom from the drift time to the z position done in AliTPCTracnsform class
  // 
  /*
    Formulas:
    
    z  = s* (z0 - vd*(t-t0))
    
    s  - side -1 and +1 
    t0 - time 0
    vd - nominal drift velocity
    zs - miscalibrated position
    
    zs = s*(z0 - vd*(1+vr)*(t-(t0+dt))
    vr  - relative change of the drift velocity
    dzt - vd*dt
    dr  = zz0-s*z
    ..
    ==>
    zs ~ z - s*vr*(z0-s*z)+s*dzt
    --------------------------------
    1. Correction function vr constant:
    
    
    dz = zs-z = -s*vr *(z0-s*z)+s*dzt         
    dzs/dl = dz/dl +s*s*vr*dz/dl 
    d(dz/dl) = vr*dz/dl     
  */

  const Float_t kZCut        = 240;    // remove the closest laser beam
  const Float_t kSaturCut    = 0.05;   // remove saturated lasers - cut on fraction of saturated 
  const Float_t kDistCut     = 3;    // distance sigma cut
  const Float_t kDistCutAbs  = 0.25;  
  const Float_t kMinClusters = 60;  // minimal amount of the clusters
  const Float_t kMinSignal   = 16;  // minimal mean height of the signal
  const Float_t kChi2Cut     = 0.1; // chi2 cut to accept drift fit
  static TLinearFitter fdriftA(3,"hyp2");
  static TLinearFitter fdriftC(3,"hyp2");
  static TLinearFitter fdriftAC(4,"hyp3");
  TVectorD fitA(3),fitC(3),fitAC(4);
  
  AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();
  AliTPCParam  * tpcparam    = calib->GetParameters();    
  //
  for (Int_t id=0; id<336; id++) fFitZ[id]=0;

  Float_t chi2A = 10;
  Float_t chi2C = 10;
  Float_t chi2AC = 10;
  Int_t npointsA=0;
  Int_t npointsC=0;
  Int_t npointsAC=0;


  for (Int_t iter=0; iter<3; iter++){
    fdriftA.ClearPoints();
    fdriftC.ClearPoints();
    fdriftAC.ClearPoints();
    //
    for (Int_t id=0; id<336; id++){
      if (!fTracksEsdParam.At(id)) continue;
      if (!AcceptLaser(id)) continue;
      if ( fClusterSatur[id]>kSaturCut)  continue;
      if ( fClusterCounter[id]<kMinClusters)  continue;
      AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
      if (track->GetTPCsignal()<kMinSignal) continue;
      AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
      AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);

      Double_t xyz[3];
      Double_t pxyz[3];
      Double_t lxyz[3];
      Double_t lpxyz[3];
      param->GetXYZ(xyz);
      param->GetPxPyPz(pxyz);
      ltrp->GetXYZ(lxyz);
      ltrp->GetPxPyPz(lpxyz);
      if (TMath::Abs(lxyz[2])>kZCut) continue;
      Float_t sz = (ltrp->GetSide()==0) ? TMath::Sqrt(chi2A): TMath::Sqrt(chi2C);
      if (npointsAC>0) sz =TMath::Sqrt(chi2AC); 
      if (TMath::Abs(fFitZ[id])>sz*kDistCut) continue;
      if (iter>0 && TMath::Abs(fFitZ[id])>kDistCutAbs) continue;

      // drift distance
      Double_t zlength =  tpcparam->GetZLength(0);
      Double_t ldrift = zlength-TMath::Abs(lxyz[2]);
      Double_t mdrift = zlength-TMath::Abs(xyz[2]);
      Double_t xxx[2] = {ldrift,lxyz[1]*ldrift/(zlength*250.)};
      if (ltrp->GetSide()==0){
	fdriftA.AddPoint(xxx,mdrift,1);
      }else{
	fdriftC.AddPoint(xxx,mdrift,1);
      }
      Double_t xxx2[3] = {ltrp->GetSide(), ldrift,lxyz[1]*ldrift/(zlength*250.)};
      fdriftAC.AddPoint(xxx2,mdrift,1);
    }
    //
    if (fdriftA.GetNpoints()>10){
      //
      fdriftA.Eval();
      if (iter==0) fdriftA.EvalRobust(0.7);
      else fdriftA.EvalRobust(0.8);
      fdriftA.GetParameters(fitA);
      npointsA= fdriftA.GetNpoints();
      chi2A = fdriftA.GetChisquare()/fdriftA.GetNpoints();
      if (chi2A<kChi2Cut ||(*fFitAside)[0]==0 ) {
	if (fFitAside->GetNoElements()<5) fFitAside->ResizeTo(5);
	(*fFitAside)[0] = fitA[0];
	(*fFitAside)[1] = fitA[1];
	(*fFitAside)[2] = fitA[2];
	(*fFitAside)[3] = fdriftA.GetNpoints();
	(*fFitAside)[4] = chi2A;	
      }
    }
    if (fdriftC.GetNpoints()>10){
      fdriftC.Eval();
      if (iter==0) fdriftC.EvalRobust(0.7);
      else fdriftC.EvalRobust(0.8);
      //
      fdriftC.GetParameters(fitC);
      npointsC= fdriftC.GetNpoints();
      chi2C = fdriftC.GetChisquare()/fdriftC.GetNpoints();
      if (chi2C<kChi2Cut||(*fFitCside)[0]==0) {	
	if (fFitCside->GetNoElements()<5) fFitCside->ResizeTo(5);
	(*fFitCside)[0] = fitC[0];
	(*fFitCside)[1] = fitC[1];
	(*fFitCside)[2] = fitC[2];
	(*fFitCside)[3] = fdriftC.GetNpoints();
	(*fFitCside)[4] = chi2C;	
      }
    }

    if (fdriftAC.GetNpoints()>10&&fdriftC.GetNpoints()>10&&fdriftA.GetNpoints()>10){
      fdriftAC.Eval();
      if (iter==0) fdriftAC.EvalRobust(0.7);
      else fdriftAC.EvalRobust(0.8);
      //
      fdriftAC.GetParameters(fitAC);
      npointsAC= fdriftAC.GetNpoints();
      chi2AC = fdriftAC.GetChisquare()/fdriftAC.GetNpoints();
      if (chi2AC<kChi2Cut||(*fFitACside)[0]==0) (*fFitACside) = fitAC;
    }
    
    for (Int_t id=0; id<336; id++){
      if (!fTracksEsdParam.At(id)) continue;
      //
      AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
      AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
      Double_t xyz[3];
      Double_t pxyz[3];
      Double_t lxyz[3];
      Double_t lpxyz[3];
      param->GetXYZ(xyz);
      param->GetPxPyPz(pxyz);
      ltrp->GetXYZ(lxyz);
      ltrp->GetPxPyPz(lpxyz); 
      Double_t zlength =  tpcparam->GetZLength(0);
      Double_t ldrift = zlength-TMath::Abs(lxyz[2]);
      Double_t mdrift = zlength-TMath::Abs(xyz[2]);

      Float_t fz =0;
      if (ltrp->GetSide()==0){
	fz = (fitA)[0]+(fitA)[1]*ldrift+(fitA)[2]*lxyz[1]*ldrift/(zlength*250.);
      }else{
	fz = (fitC)[0]+(fitC)[1]*ldrift+(fitC)[2]*lxyz[1]*ldrift/(zlength*250.);	
      }
      if (npointsAC>10){
	fz = (fitAC)[0]+(fitAC)[1]*ltrp->GetSide()+(fitAC)[2]*ldrift+(fitAC)[3]*lxyz[1]*ldrift/(zlength*250.);
      }
      fFitZ[id]=mdrift-fz;
    }
    if (fStreamLevel>0){
      TTreeSRedirector *cstream = GetDebugStreamer();
      TTimeStamp tstamp(fTime);
      Float_t valuePressure0  = AliTPCcalibDB::GetPressure(tstamp,fRun,0);
      Float_t valuePressure1 = AliTPCcalibDB::GetPressure(tstamp,fRun,1);
      Double_t ptrelative0   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
      Double_t ptrelative1   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
      Double_t temp0         = AliTPCcalibDB::GetTemperature(tstamp,fRun,0);
      Double_t temp1         = AliTPCcalibDB::GetTemperature(tstamp,fRun,1);
      TVectorD vecGoofie(20);
      AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(fRun);
      if (goofieArray) 
	for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
	  AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
	  if (gsensor) vecGoofie[isensor]=gsensor->GetValue(tstamp);
	}

      if (cstream){
	(*cstream)<<"driftv"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "mag="<<fMagF<<             //  magnetic field
	  // Environment values
	  "press0="<<valuePressure0<<
	  "press1="<<valuePressure1<<
	  "pt0="<<ptrelative0<<
	  "pt1="<<ptrelative1<<
	  "temp0="<<temp0<<
	  "temp1="<<temp1<<
	  "vecGoofie.="<<&vecGoofie<<
	  //
	  //
	  "iter="<<iter<<
	  "driftA.="<<fFitAside<<
	  "driftC.="<<fFitCside<<
	  "driftAC.="<<fFitACside<<
	  "chi2A="<<chi2A<<
	  "chi2C="<<chi2C<<
	  "chi2AC="<<chi2AC<<
	  "nA="<<npointsA<<
	  "nC="<<npointsC<<
	  "nAC="<<npointsAC<<
	  "\n";
      }
    }
  }
}
Bool_t  AliTPCcalibLaser::FitDriftV(Float_t minFraction){
  //
  // Fit corrections to the drift velocity - linear approximation in the z and global y
  //The transfromatiom from the drift time to the z position done in AliTPCTracnsform class
  // 
  // Source of outlyers : 
  // 0. Track in the saturation - postpeak
  // 1. gating grid close the part of the signal for first bundle

  // The robust fit is performed in 2 itterations /robust fraction controlled by kFraction/
  // 1. Robust fit is used in the itteration number 0
  // only fraction of laser used
  // 2. Only the tracks close to the fit used in the second itteration
  /*
    Formulas:
    
    z  = s* (z0 - vd*(t-t0))
    
    s  - side -1 and +1 
    t0 - time 0
    vd - nominal drift velocity
    zs - miscalibrated position
    
    zs = s*(z0 - vd*(1+vr)*(t-(t0+dt))
    vr  - relative change of the drift velocity
    dzt - vd*dt
    dr  = zz0-s*z
    ..
    ==>
    zs ~ z - s*vr*(z0-s*z)+s*dzt
    --------------------------------
    1. Correction function vr constant:
    
    
    dz = zs-z = -s*vr *(z0-s*z)+s*dzt         
    dzs/dl = dz/dl +s*s*vr*dz/dl 
    d(dz/dl) = vr*dz/dl     
  */
  const Int_t   knLaser      = 336;     //n laser tracks
  const Float_t kFraction[2] = {0.70,0.95};    // robust fit fraction

  const Float_t kSaturCut    = 0.05;    // remove saturated lasers - cut on fraction of saturated 
  const Float_t kDistCut     = 3.;      // distance sigma cut - 3 sigma
  const Float_t kDistCutAbs  = 1.;      // absolute cut 1 cm
  const Float_t kMinClusters = 40.;      // minimal amount of the clusters
  const Float_t kMinSignal   = 2.5;      // minimal mean height of the signal
  const Float_t kChi2Cut     = 1.0;     // chi2 cut to accept drift fit
  //
  static TLinearFitter fdriftA(3,"hyp2");
  static TLinearFitter fdriftC(3,"hyp2");
  static TLinearFitter fdriftAC(4,"hyp3");
  TVectorD fitA(3),fitC(3),fitAC(4);
  
  AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();
  AliTPCParam  * tpcparam    = calib->GetParameters();    
  //
  // reset old data
  //
  for (Int_t id=0; id<336; id++) fFitZ[id]=0;
  if (fFitAside->GetNoElements()<5) fFitAside->ResizeTo(5);
  if (fFitCside->GetNoElements()<5) fFitCside->ResizeTo(5);
  for (Int_t i=0;i<5; i++){
    (*fFitCside)[i]=0;
    (*fFitAside)[i]=0;
  }
  //
  //
  Float_t chi2A = 10;
  Float_t chi2C = 10;
  Float_t chi2AC = 10;
  Int_t npointsA=0;
  Int_t npointsC=0;
  Int_t npointsAC=0;
  Int_t nbA[4]={0,0,0,0};
  Int_t nbC[4]={0,0,0,0};
  TVectorD vecZM(336);      // measured z potion of laser
  TVectorD vecA(336);       // accepted laser
  TVectorD vecZF(336);      // fitted position  
  TVectorD vecDz(336);      // deltaZ
  TVectorD vecZS(336);      // surveyed position of laser
  // additional variable to cut
  TVectorD vecdEdx(336);    // dEdx
  TVectorD vecSy(336);      // shape y
  TVectorD vecSz(336);      // shape z
  //
  //
  for (Int_t id=0; id<336; id++){
    Int_t reject=0;
    AliTPCLaserTrack *ltrp =
      (AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
    AliTPCseed * seed = (AliTPCseed*)fTracksTPC.At(id);
    vecZM(id)= (param==0) ? 0:param->GetZ();
    vecZS(id)= ltrp->GetZ();
    vecDz(id)= 0;
    vecA(id)=1;
    vecdEdx(id)=(seed)?seed->GetdEdx():0;
    vecSy(id) =(seed)?seed->CookShape(1):0;
    vecSz(id) =(seed)?seed->CookShape(2):0;
    //
    fFitZ[id]=0;
    if (!fTracksEsdParam.At(id)) reject|=1;
    if (!param) continue; 
    if (!AcceptLaser(id))        reject|=2;
    if ( fClusterSatur[id]>kSaturCut)  reject|=4;
    if ( fClusterCounter[id]<kMinClusters)  reject|=8;
    vecA(id)=reject;
    if (reject>0) continue;
    if (ltrp->GetSide()==0){
      npointsA++;
      nbA[ltrp->GetBundle()]++;
    }
    if (ltrp->GetSide()==1){
      npointsC++;
      nbC[ltrp->GetBundle()]++;
    }
  }
  //
  // reject  "bad events"
  //
  Bool_t isOK=kTRUE;
  Int_t  naA=0;
  Int_t  naC=0;
  if (npointsA<minFraction*0.5*knLaser && npointsC<minFraction*0.5*knLaser)
    isOK=kFALSE;
  for (Int_t i=0;i<4;i++){
    //count accepted for all layers
    if (nbA[i]>minFraction*0.5*0.25*knLaser) naA++;
    if (nbC[i]>minFraction*0.5*0.25*knLaser) naC++;
  }
  if (naA<3 &&naC<3) isOK=kFALSE;
  if (isOK==kFALSE) return kFALSE;
  
  //
  //
  //
  for (Int_t iter=0; iter<2; iter++){
    fdriftA.ClearPoints();
    fdriftC.ClearPoints();
    fdriftAC.ClearPoints(); 
    npointsA=0;  npointsC=0;  npointsAC=0;
    //
    for (Int_t id=0; id<336; id++){
      if (!fTracksEsdParam.At(id)) continue;
      if (!AcceptLaser(id)) continue;
      if ( fClusterSatur[id]>kSaturCut)  continue;
      if ( fClusterCounter[id]<kMinClusters)  continue;
      AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
      if (track->GetTPCsignal()<kMinSignal) continue;
      AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
      AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
      Double_t xyz[3];
      Double_t pxyz[3];
      Double_t lxyz[3];
      Double_t lpxyz[3];
      param->GetXYZ(xyz);
      param->GetPxPyPz(pxyz);
      ltrp->GetXYZ(lxyz);
      ltrp->GetPxPyPz(lpxyz);
      Float_t sz = (ltrp->GetSide()==0) ? TMath::Sqrt(chi2A): TMath::Sqrt(chi2C);
      //if (npointsAC>0) sz =TMath::Sqrt(chi2AC); 
      if (iter>0 && TMath::Abs(fFitZ[id])>sz*kDistCut) continue;
      if (iter>0 && TMath::Abs(fFitZ[id])>kDistCutAbs) continue;

//       // drift distance
//       Double_t zlength =  tpcparam->GetZLength(0);
//       Double_t ldrift  = zlength-TMath::Abs(lxyz[2]);
//       Double_t mdrift  = zlength-TMath::Abs(xyz[2]);
      //
      Double_t zlength = (ltrp->GetSide()==0)? tpcparam->GetZLength(36): tpcparam->GetZLength(71);
      Double_t ldrift  = (ltrp->GetSide()==0)?zlength-lxyz[2]:lxyz[2]+zlength;
      Double_t mdrift  = (ltrp->GetSide()==0)?zlength-xyz[2]:xyz[2]+zlength;


      Double_t xxx[2] = {ldrift,lxyz[1]*ldrift/(zlength*250.)};
      if (iter==0 &&ltrp->GetBundle()==0) continue;  
      // skip bundle 0 - can be wrong in the 0 iteration
      if (ltrp->GetSide()==0){
	fdriftA.AddPoint(xxx,mdrift,1);
      }else{
	fdriftC.AddPoint(xxx,mdrift,1);
      }
      Double_t xxx2[3] = {ltrp->GetSide(), ldrift,lxyz[1]*ldrift/(zlength*250.)};
      fdriftAC.AddPoint(xxx2,mdrift,1);
    }
    //
    if (fdriftA.GetNpoints()>minFraction*0.5*knLaser){
      //
      fdriftA.Eval();
      //if (iter==0) fdriftA.FixParameter(2,0); //fix global y gradient
      npointsA= fdriftA.GetNpoints();
      chi2A = fdriftA.GetChisquare()/fdriftA.GetNpoints();
      fdriftA.EvalRobust(kFraction[iter]);
      fdriftA.GetParameters(fitA);
      if (chi2A<kChi2Cut ||(*fFitAside)[0]==0 ) {
	if (fFitAside->GetNoElements()<5) fFitAside->ResizeTo(5);
	(*fFitAside)[0] = fitA[0];
	(*fFitAside)[1] = fitA[1];
	(*fFitAside)[2] = fitA[2];
	(*fFitAside)[3] = fdriftA.GetNpoints();
	(*fFitAside)[4] = chi2A;	
      }
    }
    if (fdriftC.GetNpoints()>minFraction*0.5*knLaser){
      fdriftC.Eval();
      //if (iter==0) fdriftC.FixParameter(2,0); //fix global y gradient
      npointsC= fdriftC.GetNpoints();
      chi2C = fdriftC.GetChisquare()/fdriftC.GetNpoints();
      fdriftC.EvalRobust(kFraction[iter]);
      fdriftC.GetParameters(fitC);
      if (chi2C<kChi2Cut||(*fFitCside)[0]==0) {	
	if (fFitCside->GetNoElements()<5) fFitCside->ResizeTo(5);
	(*fFitCside)[0] = fitC[0];
	(*fFitCside)[1] = fitC[1];
	(*fFitCside)[2] = fitC[2];
	(*fFitCside)[3] = fdriftC.GetNpoints();
	(*fFitCside)[4] = chi2C;	
      }
    }

    if (fdriftAC.GetNpoints()>minFraction*knLaser &&npointsA>0.5*minFraction*knLaser&&npointsC>0.5*minFraction*knLaser){
      fdriftAC.Eval();
      //if (iter==0) fdriftAC.FixParameter(2,0); //fix global y gradient
      npointsAC= fdriftAC.GetNpoints();
      chi2AC = fdriftAC.GetChisquare()/fdriftAC.GetNpoints();
      fdriftAC.EvalRobust(kFraction[iter]);
      fdriftAC.GetParameters(fitAC);
      if (chi2AC<kChi2Cut||(*fFitACside)[0]==0) (*fFitACside) = fitAC;
      (*fFitACside)[0] = fitAC[0];
      (*fFitACside)[1] = fitAC[1];
      (*fFitACside)[2] = fitAC[2];
      (*fFitACside)[3] = fdriftAC.GetNpoints();
      (*fFitACside)[4] = chi2AC;
    }
    
    for (Int_t id=0; id<336; id++){
      if (!fTracksEsdParam.At(id)) continue;
      //
      AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
      AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
      Double_t xyz[3];
      Double_t pxyz[3];
      Double_t lxyz[3];
      Double_t lpxyz[3];
      param->GetXYZ(xyz);
      param->GetPxPyPz(pxyz);
      ltrp->GetXYZ(lxyz);
      ltrp->GetPxPyPz(lpxyz); 
      //Double_t zlength =  tpcparam->GetZLength(0);
      //Double_t ldrift = zlength-TMath::Abs(lxyz[2]);
      //Double_t mdrift = zlength-TMath::Abs(xyz[2]);
      Double_t zlength = (ltrp->GetSide()==0)? tpcparam->GetZLength(36): tpcparam->GetZLength(71);
      Double_t ldrift  = (ltrp->GetSide()==0)?zlength-lxyz[2]:lxyz[2]+zlength;
      Double_t mdrift  = (ltrp->GetSide()==0)?zlength-xyz[2]:xyz[2]+zlength;


      Float_t fz =0;
      if (ltrp->GetSide()==0){
	fz = (fitA)[0]+(fitA)[1]*ldrift+(fitA)[2]*lxyz[1]*ldrift/(zlength*250.);
      }else{
	fz = (fitC)[0]+(fitC)[1]*ldrift+(fitC)[2]*lxyz[1]*ldrift/(zlength*250.);	
      }
      if (npointsAC>10){
	//fz = (fitAC)[0]+(fitAC)[1]*ltrp->GetSide()+(fitAC)[2]*ldrift+(fitAC)[3]*lxyz[1]*ldrift/(zlength*250.);
      }
      fFitZ[id]=mdrift-fz;
      vecZF[id]=fz;
      vecDz[id]=mdrift-fz;
    }
    if (fStreamLevel>0){
      TTreeSRedirector *cstream = GetDebugStreamer();
      TTimeStamp tstamp(fTime);
      Float_t valuePressure0  = AliTPCcalibDB::GetPressure(tstamp,fRun,0);
      Float_t valuePressure1 = AliTPCcalibDB::GetPressure(tstamp,fRun,1);
      Double_t ptrelative0   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,0);
      Double_t ptrelative1   = AliTPCcalibDB::GetPTRelative(tstamp,fRun,1);
      Double_t temp0         = AliTPCcalibDB::GetTemperature(tstamp,fRun,0);
      Double_t temp1         = AliTPCcalibDB::GetTemperature(tstamp,fRun,1);
      TVectorD vecGoofie(20);
      AliDCSSensorArray* goofieArray = AliTPCcalibDB::Instance()->GetGoofieSensors(fRun);
      if (goofieArray) 
	for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
	  AliDCSSensor *gsensor = goofieArray->GetSensor(isensor);
	  if (gsensor) vecGoofie[isensor]=gsensor->GetValue(tstamp);
	}

      if (cstream){
	(*cstream)<<"driftvN"<<
	  "run="<<fRun<<              //  run number
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<            //  time stamp of event
	  "trigger="<<fTrigger<<      //  trigger
	  "mag="<<fMagF<<             //  magnetic field
	  // Environment values
	  "press0="<<valuePressure0<<
	  "press1="<<valuePressure1<<
	  "pt0="<<ptrelative0<<
	  "pt1="<<ptrelative1<<
	  "temp0="<<temp0<<
	  "temp1="<<temp1<<
	  "vecGoofie.="<<&vecGoofie<<
	  //
	  //
	  "vecZM.="<<&vecZM<<   // measured z position
	  "vecZS.="<<&vecZS<<   // surveyed z position
	  "vecZF.="<<&vecZF<<   // fitted   z position
	  "vecDz.="<<&vecDz<<   // fitted   z position
	  "vecA.="<<&vecA<<     // accept laser flag
	  "vecdEdx.="<<&vecdEdx<<   // dEdx  - to cut on
	  "vecSy.="<<&vecSy<<       // shape y - to cut on
	  "vecSz.="<<&vecSz<<       // shape z - to cut on
	  //
	  "iter="<<iter<<
	  "driftA.="<<fFitAside<<
	  "driftC.="<<fFitCside<<
	  "driftAC.="<<fFitACside<<
	  "chi2A="<<chi2A<<
	  "chi2C="<<chi2C<<
	  "chi2AC="<<chi2AC<<
	  "nA="<<npointsA<<
	  "nC="<<npointsC<<
	  "nAC="<<npointsAC<<
	  "\n";
	/*
	  //
	  variables to check in debug mode:
	  //
	  chainDriftN->SetAlias("driftS","250-abs(vecZS.fElements)");
	  chainDriftN->SetAlias("driftM","250-abs(vecZM.fElements)");
	  chainDriftN->SetAlias("driftF","vecZF.fElements");	  
	  chainDriftN->SetAlias("deltaZ","driftF-driftM");  //deltaZ
	  TCut cutA="iter==1&&sqrt(chi2A)<0.1&&vecZS.fElements>0&&vecA.fElements==0";
	  TCut cutC="iter==1&&sqrt(chi2C)<0.1&&vecZS.fElements<0&&vecA.fElements==0";
	  
	 */
      }
    }
  }
  return kTRUE;
}

Float_t AliTPCcalibLaser::GetDistance(AliExternalTrackParam *param, AliTPCLaserTrack *ltrp){
  //
  // get distance between mirror and laser track
  //
  //
  Double_t xyz[3];
  Double_t lxyz[3];
  param->GetXYZ(xyz);
  ltrp->GetXYZ(lxyz);
  //
  //
  Double_t dist = 0;
  //radial distance
  dist+=TMath::Abs((TMath::ATan2(xyz[1],xyz[0])-TMath::ATan2(lxyz[1],lxyz[0]))*param->GetX());
  //
  // z distance
  // apply drift correction if already exist
  //
  Float_t dz = 0;
  if (ltrp->GetSide()==0){
    if ((*fFitAside)[1]>0.) dz = ((*fFitAside)[0]+(*fFitAside)[1]*lxyz[2]+(*fFitAside)[2]*lxyz[1])-xyz[2];
  }else{
    if ((*fFitCside)[1]>0.) dz = ((*fFitCside)[0]+(*fFitCside)[1]*lxyz[2]+(*fFitCside)[2]*lxyz[1])-xyz[2];
  }
  if (TMath::Abs(dz)>3*(TMath::Abs(lxyz[2]-xyz[2])+1)) dz= TMath::Abs(lxyz[2]-xyz[2]);
  dist+=TMath::Abs(dz);
  //
  // phi dist - divergence on 50 cm
  //
  dist = TMath::Abs((param->GetParameter()[2]-ltrp->GetParameter()[2])*50);
  return dist;
}


Bool_t  AliTPCcalibLaser::AcceptLaser(Int_t id){
  //
  //
  //
  /*
  TCut cutP0("cutP0","abs((atan2(x1,x0)-atan2(lx1,lx0))*254)<1.5");
  TCut cutP1("cutP1","abs(LTr.fP[1]-Tr.fP[1])<30");
  TCut cutP2("cutP2","abs(LTr.fP[2]-Tr.fP[2])<0.03");
  TCut cutP3("cutP3","abs(Tr.fP[3])<0.05");

  TCut cutA = cutP0+cutP1+cutP2+cutP3+cutP4;
  */
  AliExternalTrackParam *param =(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp  = ( AliTPCLaserTrack*)fTracksMirror.At(id);
  //AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
  Double_t xyz[3];
  Double_t lxyz[3];
  param->GetXYZ(xyz);
  ltrp->GetXYZ(lxyz);
  if (TMath::Abs((TMath::ATan2(xyz[1],xyz[0])-TMath::ATan2(lxyz[1],lxyz[0]))*254)>1.5) return kFALSE; //cut y- P0
  if (TMath::Abs(param->GetParameter()[1]-ltrp->GetParameter()[1])>30) return kFALSE;    // cutZ -P1
  if (TMath::Abs(param->GetParameter()[2]-ltrp->GetParameter()[2])>0.03) return kFALSE;  // cut -P2
  if (TMath::Abs(param->GetParameter()[3])>0.05) return kFALSE;   // cut Tl -P3
  //
  //

  return kTRUE;
}

Int_t  AliTPCcalibLaser::FindMirror(AliESDtrack *track, AliTPCseed *seed){
  //
  // Find corresponding mirror
  // add the corresponding tracks
 
  if (!track->GetOuterParam()) return -1;

  Float_t kRadius0  = 252;
  Float_t kRadius   = 254.2;
  Int_t countercl=0;
  Float_t counterSatur=0;
  Int_t csideA =0;
  Int_t csideC =0;
  for (Int_t irow=158;irow>-1;--irow) {
    AliTPCclusterMI *c=seed->GetClusterPointer(irow);
    if (!c) continue;
    Double_t pedgeY = c->GetX()*TMath::DegToRad()*(10)-TMath::Abs(c->GetY());
    Double_t pedgeX = TMath::Min((irow)*0.75, (159.-irow)*1.5);
    if (pedgeY<3) continue;
    if (pedgeX<3) continue;
    countercl++;
    if (c->GetDetector()%36<18) csideA++;
    if (c->GetDetector()%36>=18) csideC++;
    if (c->GetMax()>900) counterSatur++;
  }
  counterSatur/=(countercl+1);
  //
  //
  //
  if (csideA<0.9*seed->GetNumberOfClusters() && csideC<0.9*seed->GetNumberOfClusters()) return 0;  // cross laser track can not happen

  Int_t side= 0;
  if (csideC>0.5*seed->GetNumberOfClusters()) side=1;


  AliExternalTrackParam param(*(track->GetOuterParam()));
  AliTracker::PropagateTrackTo(&param,kRadius0,TDatabasePDG::Instance()->GetParticle("mu+")->Mass(),3,kTRUE);
  AliTracker::PropagateTrackTo(&param,kRadius,TDatabasePDG::Instance()->GetParticle("mu+")->Mass(),0.1,kTRUE);
  AliTPCLaserTrack ltr;
  AliTPCLaserTrack *ltrp=0x0;
  //  AliTPCLaserTrack *ltrpjw=0x0;
  //
  Int_t id   = AliTPCLaserTrack::IdentifyTrack(&param,side);
 // Int_t idjw = AliTPCLaserTrack::IdentifyTrackJW(&param);
  //AliDebug(4,Form("Identified Track: %03d (%03d)",id,idjw));

  if (id!=-1 && (AliTPCLaserTrack::GetTracks()->UncheckedAt(id)))
    ltrp=(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
  else
    ltrp=&ltr;

  if (id<0) return -1;
  if (ltrp->GetSide()!=side) return -1;
  fCounter[id]++;
  //
  //
  //
  //
  if (counterSatur>fClusterSatur[id]) fClusterSatur[id]=counterSatur;
  //
  //
  Float_t radius=TMath::Abs(ltrp->GetX());
  param.Rotate(ltrp->GetAlpha());
  AliTracker::PropagateTrackTo(&param,radius,0.10566,0.01,kFALSE);
  //
  if (!fTracksMirror.At(id)) fTracksMirror.AddAt(ltrp,id);
  Bool_t accept=kTRUE;  
  //
  // choose closer track
  //
  AliExternalTrackParam * param0  = (AliExternalTrackParam *)fTracksEsdParam.At(id);
  if (param0){
    Float_t dist0=GetDistance(param0,ltrp);
    Float_t dist1=GetDistance(&param,ltrp);
    if (dist0<dist1)    accept=kFALSE;       
  }
  
  if (accept){
    fClusterCounter[id]=countercl;
    fTracksEsdParam.AddAt(param.Clone(),id);
    fTracksEsd.AddAt(track,id);
    fTracksTPC.AddAt(seed,id);
  }
  return id;
}



void AliTPCcalibLaser::DumpLaser(Int_t id) {
  //
  //  Dump Laser info to the tree
  //
  AliESDtrack   *track    = (AliESDtrack*)fTracksEsd.At(id);
  AliExternalTrackParam *param=(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp = ( AliTPCLaserTrack*)fTracksMirror.At(id);
  //
  // Fast laser ID
  //
  Double_t xyz[3];
  Double_t pxyz[3];
  Double_t lxyz[3];
  Double_t lpxyz[3];
  param->GetXYZ(xyz);
  param->GetPxPyPz(pxyz);
  ltrp->GetXYZ(lxyz);
  ltrp->GetPxPyPz(lpxyz);
  Float_t dist3D   = GetDistance(param,ltrp);
  Float_t dist0    = (TMath::ATan2(xyz[1],xyz[0])-TMath::ATan2(lxyz[1],lxyz[0]))*param->GetX();
  Float_t distphi  = (param->GetParameter()[2]-ltrp->GetParameter()[2])*50;
  
  
  if (fStreamLevel>0){
    TTreeSRedirector *cstream = GetDebugStreamer();
    Int_t time = fESD->GetTimeStamp();
    Bool_t accept = AcceptLaser(id);
    if (cstream){
      (*cstream)<<"Track"<<
	//
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field

	"id="<<id<<
	"accept="<<accept<<
	"driftA.="<<fFitAside<<
	"driftC.="<<fFitCside<<
	"time="<<time<<
	"dist3D="<<dist3D<<
	"dist0="<<dist0<<
	"distphi="<<distphi<<
	//
	//
	"counter="<<fCounter[id]<<
	"clcounter="<<fClusterCounter[id]<<
	"satur="<<fClusterSatur[id]<<
	"fz="<<fFitZ[id]<<
	//
        "LTr.="<<ltrp<<
	"Esd.="<<track<<
	"Tr.="<<param<<
	"x0="<<xyz[0]<<
	"x1="<<xyz[1]<<
	"x2="<<xyz[2]<<
	"px0="<<pxyz[0]<<
	"px1="<<pxyz[1]<<
	"px2="<<pxyz[2]<<
	//
	"lx0="<<lxyz[0]<<
	"lx1="<<lxyz[1]<<
	"lx2="<<lxyz[2]<<
	"lpx0="<<lpxyz[0]<<
	"lpx1="<<lpxyz[1]<<
	"lpx2="<<lpxyz[2]<<
	"\n";
    }
  }
}

void AliTPCcalibLaser::RefitLaserJW(Int_t id){
  //
  // Refit the track with different tracklet models:
  // 1. Per ROC using the kalman filter, different edge cuts
  // 2. Per ROC linear in y and z
  // 3. Per ROC quadratic in y and z
  // 4. Per track offset for each sector, linear for each sector, common quadratic
  // store x, y, z information for all models and the cluster to calculate the residuals
  //
  
  // The clusters which do not fulfill given  criteria are skipped
  //
  // Cluters removed from fit
  //  1. Extended shape     > kShapeCut
  //  2. In saturation  Max > kMax
  //  3. Distance to edge   < cutEdge
  //
  // Clusters not used for the calibration:
  //
  //  1. Extended shape     > kShapeCut
  //  2. In saturation  Max > kMax


  AliTPCseed *track      = (AliTPCseed*)fTracksTPC.At(id);
  AliExternalTrackParam *extparam=(AliExternalTrackParam*)fTracksEsdParam.At(id);
  AliTPCLaserTrack *ltrp = (AliTPCLaserTrack*)fTracksMirror.At(id);

  AliTPCclusterMI dummyCl;

  //two tracklets
  Int_t kMaxTracklets=2;  
  Float_t kShapeCut = 1.3;
  Float_t kRatioCut = 2.;

  Float_t kMax      = 900.;


  //=============================================//
  // Linear Fitters for the Different Approaches //
  //=============================================//
  //linear fit model in y and z; inner - outer sector, combined with offset
  static TLinearFitter fy1I(2,"hyp1");
  static TLinearFitter fy1O(2,"hyp1");
  static TLinearFitter fz1I(2,"hyp1");
  static TLinearFitter fz1O(2,"hyp1");
  static TLinearFitter fy1IO(3,"hyp2");
  static TLinearFitter fz1IO(3,"hyp2");
  //quadratic fit model in y and z; inner - sector
  static TLinearFitter fy2I(3,"hyp2");
  static TLinearFitter fy2O(3,"hyp2");
  static TLinearFitter fz2I(3,"hyp2");
  static TLinearFitter fz2O(3,"hyp2");
  //common quadratic fit for IROC and OROC in y and z
  static TLinearFitter fy4(5,"hyp4");
  static TLinearFitter fz4(5,"hyp4");


  //set standard cuts
  if ( fNcuts==0 ){
      fNcuts=1;
      fEdgeXcuts[0]=4;
      fEdgeYcuts[0]=3;
      fNClCuts[0]=20;
  }
  //=============================//
  // Loop over all Tracklet Cuts //
  //=============================//
  for (Int_t icut=0; icut<fNcuts; icut++){
    Float_t xinMin = 2500, xinMax=-90;
    Float_t xoutMin=2500, xoutMax=-90;
    Float_t msigmaYIn=0;
    Float_t msigmaYOut=0;
    Float_t mqratioIn=0;
    Float_t mqratioOut=0;


      AliDebug(4,Form("Processing cut %d for track with ID %d",icut,id));
      //cut parameters
      Double_t edgeCutX = fEdgeXcuts[icut];
      Double_t edgeCutY = fEdgeYcuts[icut];
      Int_t    nclCut   = (Int_t)fNClCuts[icut];
      //=========================//
      // Parameters to calculate //
      //=========================//
      //fit parameter inner, outer tracklet and combined fit
      TVectorD vecy1resInner(2),vecz1resInner(2); //pol1 fit parameters inner
      TVectorD vecy2resInner(3),vecz2resInner(3); //pol2 fit parameters inner
      //
      TVectorD vecy1resOuter(2),vecz1resOuter(2); //pol1 fit parameters outer
      TVectorD vecy2resOuter(3),vecz2resOuter(3); //pol2 fit parameters outer
      TVectorD vecy4res(5),vecz4res(5);
      TVectorD vecy1resIO(3),vecz1resIO(3);
      // cluster and track positions for each row - used for residuals
      TVectorD vecgX(159);        // global X
      TVectorD vecgY(159);        // global Y
      TVectorD vecgZ(159);        // global Z

      TVectorD vecX(159);        // x is the same for all (row center)
      TVectorD vecYkalman(159);  // y from kalman fit
      TVectorD vecZkalman(159);  // z from kalman fit
      TVectorD vecY1(159);       // y from pol1 fit per ROC
      TVectorD vecZ1(159);       // z from pol1 fit per ROC
      TVectorD vecY1IO(159);     // y from pol1 fit per ROC
      TVectorD vecZ1IO(159);     // z from pol1 fit per ROC
      TVectorD vecY2(159);       // y from pol2 fit per ROC
      TVectorD vecZ2(159);       // z from pol2 fit per ROC
      TVectorD vecY4(159);       // y from sector fit
      TVectorD vecZ4(159);       // z from sector fit
      TVectorD vecClY(159);      // y cluster position
      TVectorD vecClZ(159);      // z cluster position
      TVectorD vecSec(159);      // sector for each row
      TVectorD  isReject(159);    // flag - cluster to be rejected 
      //chi2 of fits
      Double_t chi2I1z=0;       // chi2 of pol1 fit in z (inner)
      Double_t chi2I1y=0;       // chi2 of pol1 fit in y (inner)
      Double_t chi2O1z=0;       // chi2 of pol1 fit in z (outer)
      Double_t chi2O1y=0;       // chi2 of pol1 fit in y (outer)
      Double_t chi2IO1z=0;       // chi2 of pol1 fit in z (outer)
      Double_t chi2IO1y=0;       // chi2 of pol1 fit in y (outer)
      Double_t chi2I2z=0;       // chi2 of pol2 fit in z (inner)
      Double_t chi2I2y=0;       // chi2 of pol2 fit in y (inner)
      Double_t chi2O2z=0;       // chi2 of pol2 fit in z (outer)
      Double_t chi2O2y=0;       // chi2 of pol2 fit in y (outer)
      Double_t chi2IOz=0;       // chi2 of hyp4 fit in z (inner+outer)
      Double_t chi2IOy=0;       // chi2 of hyp4 fit in y (inner+outer)
      //more
      Int_t innerSector = -1;    // number of inner sector
      Int_t outerSector = -1;    // number of outer sector
      Int_t nclI=0;              // number of clusters (inner)
      Int_t nclO=0;              // number of clusters (outer)
      //=================================================//
      // Perform the Kalman Fit using the Tracklet Class //
      //=================================================//
      AliTPCTracklet::SetEdgeCut(edgeCutX,edgeCutY);
      TObjArray tracklets=
	  AliTPCTracklet::CreateTracklets(track,AliTPCTracklet::kKalman,
					  kFALSE,nclCut,kMaxTracklets);
      // tracklet pointers
      AliTPCTracklet *trInner = (AliTPCTracklet*)tracklets.At(0);
      AliTPCTracklet *trOuter = (AliTPCTracklet*)tracklets.At(1);
      AliTPCTracklet *tr=0x0;
      AliTPCTracklet dummy;
      //continue if we didn't find a tracklet
      if ( !trInner && !trOuter ) continue;
      //================================================================================//
      // Swap Inner and Outer if necessary (inner sector is definde by smaller local x) //
      //================================================================================//
      if ( trInner && trOuter ){
	  if ( !trInner->GetInner() || !trOuter->GetInner() ) continue;
	  if ( trInner->GetInner()->GetX() > trOuter->GetInner()->GetX() ){
	      tr = trInner;
	      trInner=trOuter;
	      trOuter=tr;
	      AliDebug(5,Form("Swapped Sectors: %02d (%f) <-> %02d (%f)", trOuter->GetSector(), trOuter->GetInner()->GetX(), trInner->GetSector(), trInner->GetInner()->GetX()));
	  }
      } else {
	  if ( trInner ){
              if ( !trInner->GetInner() ) continue;
	      trOuter=&dummy;
	      if ( trInner->GetSector()>35 ){
		  trOuter=trInner;
                  trInner=&dummy;
	      }
	  } else { //trOuter
              if ( !trOuter->GetInner() ) continue;
              trInner=&dummy;
	      if ( trOuter->GetSector()<36 ){
		  trInner=trOuter;
		  trOuter=&dummy;
	      }
	  }
      }
      innerSector = trInner->GetSector();
      if ( innerSector>=0 ) AliDebug(5,Form("Found inner Sector %02d at X %.2f", innerSector, trInner->GetInner()->GetX()));
      outerSector = trOuter->GetSector();
      if ( outerSector>=0 ) AliDebug(5,Form("Found outer Sector %02d at X %.2f", outerSector, trOuter->GetInner()->GetX()));

      // array of clusters
      TClonesArray arrCl("AliTPCclusterMI",159);
      arrCl.ExpandCreateFast(159);
      //=======================================//
      // fill fitters with cluster information //
      //=======================================//
      AliDebug(3,"Filing Cluster Information");

      //
      //

      for (Int_t irow=158;irow>-1;--irow) {
	  AliTPCclusterMI *c=track->GetClusterPointer(irow);
	  AliTPCclusterMI &cl = (AliTPCclusterMI&) (*arrCl[irow]);
	  cl=dummyCl;
	  vecX[irow]   = 0;
	  vecClY[irow] = 0;
	  vecClZ[irow] = 0;
	  Float_t meanY=0, sumY=0;
	  for (Int_t drow=-1;drow<=1;drow++) {
	    if (irow+drow<0) continue;
	    if (irow+drow>158) continue;
	    AliTPCclusterMI *ccurrent=track->GetClusterPointer(irow);
	    if (!ccurrent) continue;
	    Int_t roc = static_cast<Int_t>(ccurrent->GetDetector());
	    if ( roc!=innerSector && roc!=outerSector ) continue;
	    if (ccurrent->GetX()<10) continue;
	    if (ccurrent->GetY()==0) continue;
	    meanY+=ccurrent->GetY();
	    sumY++;
	  }
	  if (sumY>0)  meanY/=sumY; 
	 
          //
          vecSec[irow]=-1;
	  if (!c) continue;
	  Double_t pedgeY = c->GetX()*TMath::Tan(TMath::DegToRad()*(10))-TMath::Abs(meanY);
	  Double_t pedgeX = TMath::Min((irow)*0.75, (159.-irow)*1.5);
	  
          //
	  Int_t roc = static_cast<Int_t>(c->GetDetector());
	  if ( roc!=innerSector && roc!=outerSector ) continue;
          vecSec[irow]=roc;
	  //store clusters in clones array
	  cl=*c;	  
	  //
	  if (c->GetMax()<4)   continue;  // noise cluster?
	  if (TMath::Abs(c->GetY())<0.0001)   continue;  // noise cluster?
	  //cluster position
	  vecX[irow]   = c->GetX();
	  vecClY[irow] = c->GetY();
	  vecClZ[irow] = c->GetZ();
	  //
// 	  Float_t gxyz[3];
// 	  c->GetGlobalXYZ(gxyz);
// 	  vecgX[irow]   = gxyz[0];
// 	  vecgY[irow]   = gxyz[1];
// 	  vecgZ[irow]   = gxyz[2];
          //
          Double_t x = vecX[irow]-133.4; //reference is between IROC and OROC
          Double_t y = vecClY[irow];
	  Double_t z = vecClZ[irow];
          //
	  Double_t x2[2]={x,x*x};   //linear and parabolic parameters
	  Double_t x4[4]={0,0,0,0}; //hyp4 parameters
          Double_t xIO[2]={0,x};    //common linear + offset IROC-OROC
	  if ( roc == innerSector ) {
	      x4[0]=1; //offset inner - outer sector
	      x4[1]=x; //slope parameter inner sector
              xIO[0]=1;
	  } else {
	      x4[2]=x; //slope parameter outer sector
	  }
	  x4[3]=x*x;   //common parabolic shape
	  if (pedgeX < fEdgeXcuts[icut]) continue;
	  if (pedgeY < fEdgeYcuts[icut]) continue;
	  if (c->GetMax()>900) continue;  // cluster in saturation	  
          //
	  if ( roc==innerSector ){
	      fy1I.AddPoint(x2,y);
	      fz1I.AddPoint(x2,z);
	      fy2I.AddPoint(x2,y);
	      fz2I.AddPoint(x2,z);
              ++nclI;
	      if (c->GetX()<xinMin) xinMin=c->GetX();
	      if (c->GetX()>xinMax) xinMax=c->GetX();

	      msigmaYIn +=TMath::Sqrt(c->GetSigmaY2());
	      mqratioIn +=c->GetMax()/c->GetQ();

	  }
	  if ( roc==outerSector ){
	      fy1O.AddPoint(x2,y);
	      fz1O.AddPoint(x2,z);
	      fy2O.AddPoint(x2,y);
	      fz2O.AddPoint(x2,z);
              ++nclO;
	      if (c->GetX()<xoutMin) xoutMin=c->GetX();
	      if (c->GetX()>xoutMax) xoutMax=c->GetX();
	      msigmaYOut +=TMath::Sqrt(c->GetSigmaY2());
	      mqratioOut +=c->GetMax()/c->GetQ();

	  }
	  fy4.AddPoint(x4,y);
	  fz4.AddPoint(x4,z);
          fy1IO.AddPoint(xIO,y);
          fz1IO.AddPoint(xIO,z);
      }
      if (nclI>0)  {
	msigmaYIn/=nclI;
	mqratioIn/=nclI;
      }
      if (nclO>0)  {
	msigmaYOut/=nclO;
	mqratioOut/=nclO;
      }
      //======================================//
      // Evaluate and retrieve fit parameters //
      //======================================//
      AliDebug(5,Form("Evaluating fits with inner (outer) Sec: %02d (%02d)",innerSector,outerSector));
      //inner sector
      if (  (innerSector>-1) && (fy1I.GetNpoints()>0) ){
	  fy1I.Eval();
	  fz1I.Eval();
	  fy2I.Eval();
	  fz2I.Eval();
	  fy1I.GetParameters(vecy1resInner);
	  fz1I.GetParameters(vecz1resInner);
	  fy2I.GetParameters(vecy2resInner);
	  fz2I.GetParameters(vecz2resInner);
	  chi2I1y=fy1I.GetChisquare()/(fy1I.GetNpoints()-2);
	  chi2I1z=fz1I.GetChisquare()/(fz1I.GetNpoints()-2);
	  chi2I2y=fy2I.GetChisquare()/(fy2I.GetNpoints()-3);
	  chi2I2z=fz2I.GetChisquare()/(fz2I.GetNpoints()-3);
      }
      //outer sector
      if (  (outerSector>-1) && (fy1O.GetNpoints()>0) ){
	  fy1O.Eval();
	  fz1O.Eval();
	  fy2O.Eval();
	  fz2O.Eval();
	  fy1O.GetParameters(vecy1resOuter);
	  fz1O.GetParameters(vecz1resOuter);
	  fy2O.GetParameters(vecy2resOuter);
	  fz2O.GetParameters(vecz2resOuter);
	  chi2O1y=fy1O.GetChisquare()/(fy1O.GetNpoints()-2);
	  chi2O1z=fz1O.GetChisquare()/(fz1O.GetNpoints()-2);
	  chi2O2y=fy2O.GetChisquare()/(fy2O.GetNpoints()-3);
	  chi2O2z=fz2O.GetChisquare()/(fz2O.GetNpoints()-3);
      }
      //combined hyp4 fit
      if ( innerSector>0 && outerSector>0 ){
	  if (fy4.GetNpoints()>0) {
	      fy4.Eval();
	      fy4.GetParameters(vecy4res);
	      chi2IOy=fy4.GetChisquare()/(fy4.GetNpoints()-5);
	  }
	  if (fz4.GetNpoints()>0) {
	      fz4.Eval();
	      fz4.GetParameters(vecz4res);
	      chi2IOz=fz4.GetChisquare()/(fz4.GetNpoints()-5);
	  }
          if (fy1IO.GetNpoints()>0) {
            fy1IO.Eval();
            fy1IO.GetParameters(vecy1resIO);
            chi2IO1y=fy1IO.GetChisquare()/(fy1IO.GetNpoints()-3);
          }
          if (fz1IO.GetNpoints()>0) {
            fz1IO.Eval();
            fz1IO.GetParameters(vecz1resIO);
            chi2IO1z=fz1IO.GetChisquare()/(fz1IO.GetNpoints()-3);
          }
      }
      //clear points
      fy4.ClearPoints();  fz4.ClearPoints();
      fy1I.ClearPoints(); fy1O.ClearPoints();
      fz1I.ClearPoints(); fz1O.ClearPoints();
      fy2I.ClearPoints(); fy2O.ClearPoints();
      fz2I.ClearPoints(); fz2O.ClearPoints();
      fy1IO.ClearPoints(); fz1IO.ClearPoints();
      //==============================//
      // calculate tracklet positions //
      //==============================//
      AliDebug(4,"Calculate tracklet positions");
      for (Int_t irow=158;irow>-1;--irow) {
	isReject[irow]=0;
	 AliTPCclusterMI *c=track->GetClusterPointer(irow);
	 if ( vecSec[irow]!=innerSector && vecSec[irow]!=outerSector ) { // no cluster in given sectors
	   isReject[irow]+=1;
	 }

	 if (!c) { //no cluster	   
	   isReject[irow]+=2;
	 }else{
	   if (c->GetMax()>kMax){   //saturation
	     isReject[irow]+=4;
	   }
	   if ( vecSec[irow] == outerSector ) {   //extended shape
	     if (c->GetMax()/c->GetQ()> mqratioOut*kRatioCut)                isReject[irow]+=8;
	     if (TMath::Sqrt(c->GetSigmaY2())>msigmaYOut*kShapeCut)    isReject[irow]+=16;
	   }else{
	     if (c->GetMax()/c->GetQ()> mqratioIn*kRatioCut)                isReject[irow]+=8;
	     if (TMath::Sqrt(c->GetSigmaY2())>msigmaYIn*kShapeCut)    isReject[irow]+=16;
	   }	   
	 }



	  if ( vecSec[irow]==-1 ) continue;  //no cluster info
          if ( vecSec[irow]!=innerSector && vecSec[irow]!=outerSector ) continue;
	  tr=&dummy;
	  Double_t x=vecX[irow];
	  Double_t xref=x-133.4;
          //
          Double_t yoffInner=0;
          Double_t zoffInner=0;
          Double_t yoffInner1=0;
          Double_t zoffInner1=0;
          Double_t yslopeInner=0;
          Double_t yslopeOuter=0;
          Double_t zslopeInner=0;
	  Double_t zslopeOuter=0;
          //positions of hyperplane fits
	  if ( vecSec[irow] == outerSector ) {
	      tr=trOuter;
              vecY1[irow]=vecy1resOuter[0]+vecy1resOuter[1]*xref;
              vecZ1[irow]=vecz1resOuter[0]+vecz1resOuter[1]*xref;
              vecY2[irow]=vecy2resOuter[0]+vecy2resOuter[1]*xref+vecy2resOuter[2]*xref*xref;
	      vecZ2[irow]=vecz2resOuter[0]+vecz2resOuter[1]*xref+vecz2resOuter[2]*xref*xref;
              yslopeOuter=vecy4res[3];
	      zslopeOuter=vecz4res[3];
	  } else {
	      tr=trInner;
              vecY1[irow]=vecy1resInner[0]+vecy1resInner[1]*xref;
              vecZ1[irow]=vecz1resInner[0]+vecz1resInner[1]*xref;
              vecY2[irow]=vecy2resInner[0]+vecy2resInner[1]*xref+vecy2resInner[2]*xref*xref;
              vecZ2[irow]=vecz2resInner[0]+vecz2resInner[1]*xref+vecz2resInner[2]*xref*xref;
              yoffInner=vecy4res[1];
              zoffInner=vecz4res[1];
              yoffInner1=vecy1resIO[1];
              zoffInner1=vecz1resIO[1];
              yslopeInner=vecy4res[2];
	      zslopeInner=vecz4res[2];
	  }
          vecY1IO[irow]=vecy1resIO[0]+yoffInner1+vecy1resIO[2]*xref;
          vecZ1IO[irow]=vecz1resIO[0]+zoffInner1+vecz1resIO[2]*xref;
	  vecY4[irow]=vecy4res[0]+yoffInner+yslopeInner*xref+yslopeOuter*xref+vecy4res[4]*xref*xref;
	  vecZ4[irow]=vecz4res[0]+zoffInner+zslopeInner*xref+zslopeOuter*xref+vecz4res[4]*xref*xref;
          //positions of kalman fits
	  Double_t gxyz[3],xyz[3];
	  AliExternalTrackParam *param = 0x0;
          //
	  param=tr->GetInner();
	  if (param){
	      param->GetXYZ(gxyz);
	      Float_t bz = AliTracker::GetBz(gxyz);
	      param->GetYAt(x, bz, xyz[1]);
	      param->GetZAt(x, bz, xyz[2]);
	      vecYkalman[irow]=xyz[1];
	      vecZkalman[irow]=xyz[2];
	  }
	  //
	  //
	  //
	 
      }
      //=====================================================================//
      // write results from the different tracklet fits with debug streamers //
      //=====================================================================//
      if (fStreamLevel>4){
	  TTreeSRedirector *cstream = GetDebugStreamer();
	  if (cstream){
	    Float_t dedx = track->GetdEdx();
	    (*cstream)<<"FitModels"<<
	      "run="<<fRun<<              //  run number
	      "event="<<fEvent<<          //  event number
	      "time="<<fTime<<            //  time stamp of event
	      "trigger="<<fTrigger<<      //  trigger
	      "mag="<<fMagF<<             //  magnetic field	      
	      //
	      "cutNr="      << icut <<
	      "edgeCutX="   << edgeCutX <<
	      "edgeCutY="   << edgeCutY <<
	      "nclCut="     << nclCut <<
	      "innerSector="<< innerSector <<
	      "outerSector="<< outerSector <<
	      "dEdx="       << dedx <<
	      "LTr.="       << ltrp <<
	      "Tr.="        << extparam <<
              "yPol1In.="   << &vecy1resInner <<
              "zPol1In.="   << &vecz1resInner <<
              "yPol1InOut.="<< &vecy1resIO <<
              "zPol1InOut.="<< &vecz1resIO <<
              "yPol2In.="   << &vecy2resInner <<
	      "zPol2In.="   << &vecz2resInner <<
	      "yPol1Out.="  << &vecy1resOuter <<
	      "zPol1Out.="  << &vecz1resOuter <<
	      "yPol2Out.="  << &vecy2resOuter <<
	      "zPol2Out.="  << &vecz2resOuter <<
	      "yInOut.="    << &vecy4res <<
	      "zInOut.="    << &vecz4res <<
              "chi2y1In="   << chi2I1y <<
              "chi2z1In="   << chi2I1z <<
              "chi2y1InOut="<< chi2IO1y <<
              "chi2z1InOut="<< chi2IO1z <<
              "chi2y1Out="  << chi2O1y <<
	      "chi2z1Out="  << chi2O1z <<
	      "chi2y2In="   << chi2I2y <<
	      "chi2z2In="   << chi2I2z <<
	      "chi2y2Out="  << chi2O2y <<
	      "chi2z2Out="  << chi2O2z <<
	      "chi2yInOut=" << chi2IOy <<
	      "chi2zInOut=" << chi2IOz <<
	      "trletIn.="   << trInner <<
	      "trletOut.="  << trOuter <<
	      "nclI="       << nclI <<
	      "nclO="       << nclO <<
 	      "xinMin="     << xinMin<<
	      "xinMax="     << xinMax<<
 	      "xoutMin="    << xoutMin<<
	      "xoutMax="     << xoutMax<<
	      "msigmaYIn=" <<msigmaYIn<<
	      "msigmaYOut=" <<msigmaYOut<<
	      "mqratioIn=" <<mqratioIn<<
	      "mqratioOut=" <<	mqratioOut      <<
	      "\n";
	  }
      }

      // wirte residuals information
      if (fStreamLevel>5){
	  TTreeSRedirector *cstream = GetDebugStreamer();
	  if (cstream){
	    Float_t dedx = track->GetdEdx();
	    (*cstream)<<"Residuals"<<
	      "run="<<fRun<<              //  run number
	      "event="<<fEvent<<          //  event number
	      "time="<<fTime<<            //  time stamp of event
	      "trigger="<<fTrigger<<      //  trigger
	      "mag="<<fMagF<<             //  magnetic field	      
	      //
	      "cutNr="      << icut <<
	      "edgeCutX="   << edgeCutX <<
	      "edgeCutY="   << edgeCutY <<
	      "nclCut="     << nclCut   <<
	      "LTr.="       << ltrp <<
	      "Tr.="        << extparam<<
	      "dEdx="       << dedx <<
	      "Cl.="        << &arrCl <<
	      "vX.="        << &vecgX<<   // global x
	      "vY.="        << &vecgY<<   // global y
	      "vZ.="        << &vecgZ<<   // global z
	      "TrX.="       << &vecX <<
	      "TrYpol1.="   << &vecY1 <<
	      "TrZpol1.="   << &vecZ1 <<
	      "TrYpol2.="   << &vecY2 <<
	      "TrZpol2.="   << &vecZ2 <<
              "TrYpol1InOut.="<< &vecY1IO <<
              "TrZpol1InOut.="<< &vecZ1IO <<
              "TrYInOut.="  << &vecY4 <<
              "TrZInOut.="  << &vecZ4 <<
              "ClY.="       << &vecClY <<
	      "ClZ.="       << &vecClZ <<
	      "isReject.="  << &isReject<<
	      "sec.="       << &vecSec <<
	      "nclI="       << nclI <<
	      "nclO="       << nclO <<
 	      "xinMin="     << xinMin<<
	      "xinMax="     << xinMax<<
 	      "xoutMin="    << xoutMin<<
	      "xoutMax="     << xoutMax<<
	      "msigmaYIn=" <<msigmaYIn<<
	      "msigmaYOut=" <<msigmaYOut<<
	      "mqratioIn=" <<mqratioIn<<
	      "mqratioOut=" <<	mqratioOut      <<
	      "yInOut.="    << &vecy4res <<
	      "zInOut.="    << &vecz4res <<
	      //chi2s
	      "chi2y1In="   << chi2I1y <<    //
	      "chi2z1In="   << chi2I1z <<
	      "chi2y1Out="  << chi2O1y <<
	      "chi2z1Out="  << chi2O1z <<
              "chi2y1InOut="<< chi2IO1y <<
              "chi2z1InOut="<< chi2IO1z <<
              "chi2y2In="   << chi2I2y <<
	      "chi2z2In="   << chi2I2z <<
	      "chi2y2Out="  << chi2O2y <<
	      "chi2z2Out="  << chi2O2z <<
	      "chi2yInOut=" << chi2IOy <<
	      "chi2zInOut=" << chi2IOz <<
	      // fit parameters
	      "yPol1In.="   << &vecy1resInner <<
	      "zPol1In.="   << &vecz1resInner <<
	      "yPol2In.="   << &vecy2resInner <<
	      "zPol2In.="   << &vecz2resInner <<
	      "yPol1Out.="  << &vecy1resOuter <<
	      "zPol1Out.="  << &vecz1resOuter <<
              "yPol1InOut.="<< &vecy1resIO <<
              "zPol1InOut.="<< &vecz1resIO <<
              "yPol2Out.="  << &vecy2resOuter <<
	      "zPol2Out.="  << &vecz2resOuter <<

	      "\n";

	  }
      }
      //==========================//
      // Fill Residual Histograms //
      //==========================//
      if (!fHisNclIn) MakeFitHistos(); 

      TH2F *profy = (TH2F*)fDeltaYres.UncheckedAt(id);
      TH2F *profz = (TH2F*)fDeltaZres.UncheckedAt(id);
      TH2F *profy2 = (TH2F*)fDeltaYres2.UncheckedAt(id);
      TH2F *profz2 = (TH2F*)fDeltaZres2.UncheckedAt(id);
      TH2F *profyabs = (TH2F*)fDeltaYresAbs.UncheckedAt(id);
      TH2F *profzabs = (TH2F*)fDeltaZresAbs.UncheckedAt(id);
      //      TH2F *profy3 = (TH2F*)fDeltaYres3.UncheckedAt(id);
      //TH2F *profz3 = (TH2F*)fDeltaZres3.UncheckedAt(id);
      //
      for (Int_t irow=158;irow>-1;--irow) {
	if (vecSec[irow]==-1)continue; // no cluster info
	if (isReject[irow]>0.5) continue;  // 
	Double_t ycl  = vecClY[irow];
	Double_t yfit = vecY1[irow];
	Double_t yfit2 = vecY2[irow];
	Double_t x    = vecX[irow];
	Double_t yabsbeam = -1000;
	if(vecSec[irow]==outerSector && outerSector==fBeamSectorOuter[id])
	  yabsbeam = fBeamSlopeYOuter[id]*x + fBeamOffsetYOuter[id];
	else if(innerSector==fBeamSectorInner[id])
	  yabsbeam = fBeamSlopeYInner[id]*x + fBeamOffsetYInner[id];

	//	Double_t yfit3 = vecY2[irow];
	Double_t zcl  = vecClZ[irow];
	Double_t zfit = vecZ1[irow];
 	Double_t zfit2 = vecZ2[irow];
	//Double_t zfit3 = vecZ2[irow];

	// dz abs 
	// The expressions for zcorrected has been obtained by
	// inverting the fits in the FitDriftV() method (ignoring the
	// global y dependence for now):
	// A side:
	// 250 - zmeasured = [0] + [1]*(250 - zreal) + .... (yglobal)
	// =>
	// zreal = (zmeasured + [0] - (1-[1])*250.0)/[1]
	// 
	// C side:
	// 250 + zmeasured  = [0] + [1]*(250+zreal) + .... (yglobal)
	// =>
	// zreal = (zmeasured - [0] + (1 - [1])*250.0)/[1]

	Double_t dzabs = -1000;
	Double_t zcorrected = -1000;
	if (ltrp->GetSide()==0){
	  if ((*fFitAside)[1]>0. || fUseFixedDriftV) { 
	    // ignore global y dependence for now
	    zcorrected = 0;	    
	    if(!fUseFixedDriftV) 
	      zcorrected = (zcl + (*fFitAside)[0] - 
			    (1.0-(*fFitAside)[1])*250.0)/(*fFitAside)[1];
	    else
	      zcorrected = (zcl + fFixedFitAside0 - 
			    (1.0-fFixedFitAside1)*250.0)/fFixedFitAside1;
	    //	    zcorrected = zcl;
	    if(vecSec[irow]==outerSector && outerSector==fBeamSectorOuter[id])
	      dzabs = zcorrected -fBeamSlopeZOuter[id]*x -fBeamOffsetZOuter[id];
	    else if(innerSector==fBeamSectorInner[id])
	      dzabs = zcorrected -fBeamSlopeZInner[id]*x -fBeamOffsetZInner[id];
	  }
	} else {
	  if ((*fFitCside)[1]>0. || fUseFixedDriftV) {
	    
	    if(!fUseFixedDriftV) 
	      zcorrected = (zcl - (*fFitCside)[0] + 
			    (1.0-(*fFitCside)[1])*250.0)/(*fFitCside)[1];
	    else
	      zcorrected = (zcl - fFixedFitCside0 + 
			    (1.0-fFixedFitCside1)*250.0)/fFixedFitCside1;
	    
	    //	    zcorrected = zcl;
	    if(vecSec[irow]==outerSector && outerSector==fBeamSectorOuter[id])
	      dzabs = zcorrected -fBeamSlopeZOuter[id]*x -fBeamOffsetZOuter[id];
	    else if(innerSector==fBeamSectorInner[id])
	      dzabs = zcorrected -fBeamSlopeZInner[id]*x -fBeamOffsetZInner[id];
	  }
	}
	
	if (TMath::Abs(yfit-ycl)<2&&TMath::Abs(zfit-zcl)<2){
	  if (profy){	    
	      profy->Fill(irow,ycl-yfit);
	      profy2->Fill(irow,ycl-yfit2);
	      if(yabsbeam<-100) {
		fHisYAbsErrors->Fill(id);
		//		profyabs->Fill(irow,-0.99);
	      }	else
		profyabs->Fill(irow,ycl-yabsbeam);

	      //	      profy3->Fill(irow,ycl-yfit3);
	  }
	  if (profz) {
	      profz->Fill(irow,zcl-zfit);
	      profz2->Fill(irow,zcl-zfit2);
	      //profz3->Fill(irow,zcl-zfit3);
	      if(dzabs<-100) {

		fHisZAbsErrors->Fill(id);
	      }else
		profzabs->Fill(irow,dzabs);
	  }
	}
      }
      //
      //
      // Fill laser fit histograms
      //
      Float_t dedx = track->GetdEdx();
      if (nclI>20&&nclO>20){
	fHisNclIn->Fill(id,nclI);      //->Number of clusters inner
	fHisNclOut->Fill(id,nclO);     //->Number of clusters outer
	fHisNclIO->Fill(id,nclI+nclO);      //->Number of cluster inner outer
	//
	fHisLclIn->Fill(id,xinMax-xinMin);        //->Level arm inner
	fHisLclOut->Fill(id,xoutMax-xoutMin);     //->Level arm outer
	fHisLclIO->Fill(id,xoutMax-xinMin);       //->Number of cluster inner outer
	//
	fHisdEdx->Fill(id,TMath::Sqrt(TMath::Abs(dedx)));
	fHisdZfit->Fill(id,fFitZ[id]);
	// 
	//
	fHisChi2YIn1->Fill(id,TMath::Sqrt(chi2I1y));      //->chi2 y inner - line
	fHisChi2YOut1->Fill(id,TMath::Sqrt(chi2O1y));     //->chi2 y inner - line
	fHisChi2YIn2->Fill(id,TMath::Sqrt(chi2I2y));      //->chi2 y inner - parabola
	fHisChi2YOut2->Fill(id,TMath::Sqrt(chi2O2y));     //->chi2 y inner - parabola
	fHisChi2YIO1->Fill(id,TMath::Sqrt(chi2IOy));      //->chi2 y IO    - common
	
	
	fHisChi2ZIn1->Fill(id,TMath::Sqrt(chi2I1z));      //->chi2 z inner - line
	fHisChi2ZOut1->Fill(id,TMath::Sqrt(chi2O1z));     //->chi2 z inner - line
	fHisChi2ZIn2->Fill(id,TMath::Sqrt(chi2O2y));      //->chi2 z inner - parabola
	fHisChi2ZOut2->Fill(id,TMath::Sqrt(chi2O2z));     //->chi2 z inner - parabola
	fHisChi2ZIO1->Fill(id,TMath::Sqrt(chi2IOz));      //->chi2 z IO    - common
	//
	//
	fHisPy1vP0->Fill(id,vecy1resOuter[0]-vecy1resInner[0]);     //-> delta y   P0outer-P0inner - line
	fHisPy2vP0->Fill(id,vecy2resOuter[0]-vecy2resInner[0]);     //-> delta y   P0outer-P0inner - parabola
	fHisPy3vP0->Fill(id,vecy4res[1]);                           //-> delta y   P0outer-P0inner - common parabola
	//
	fHisPy1vP1->Fill(id,vecy1resOuter[1]-vecy1resInner[1]);     //-> delta ky  P1outer-P1inner - line
	fHisPy2vP1->Fill(id,vecy2resOuter[1]-vecy2resInner[1]);     //-> delta ky  P1outer-P1inner - parabola
	fHisPy3vP1->Fill(id,vecy4res[3]-vecy4res[2]);               //-> delta ky  P1outer-P1inner - common parabola
	//
	fHisPy3vP2IO->Fill(id,vecy4res[4]);        //-> Curv  P2outerinner - common parabola
	fHisPz1vP0->Fill(id,vecz1resOuter[0]-vecz1resInner[0]);     //-> delta z   P0outer-P0inner - line
	fHisPz2vP0->Fill(id,vecz2resOuter[0]-vecz2resInner[0]);     //-> delta z   P0outer-P0inner - parabola
	fHisPz3vP0->Fill(id,vecz4res[1]);                           //-> delta z   P0outer-P0inner - common parabola
	//
	fHisPz1vP1->Fill(id,vecz1resOuter[1]-vecz1resInner[1]);     //-> delta kz  P1outer-P1inner - line
	fHisPz2vP1->Fill(id,vecz2resOuter[1]-vecz2resInner[1]);     //-> delta kz  P1outer-P1inner - parabola
	fHisPz3vP1->Fill(id,vecz4res[3]-vecz4res[2]);               //-> delta kz  P1outer-P1inner - common parabola
	fHisPz3vP2IO->Fill(id,vecz4res[4]);        //-> Curv  P2outerinner - common parabola
      }
      if(nclI>20){
	fHisPy2vP2In->Fill(id,vecy2resInner[2]);   //-> Curv  P2inner - parabola
	fHisPz2vP2In->Fill(id,vecz2resInner[2]);   //-> Curv  P2inner - parabola
      }
      //
      if (nclO>20){
	fHisPz2vP2Out->Fill(id,vecz2resOuter[2]);  //-> Curv  P2outer - parabola
	fHisPy2vP2Out->Fill(id,vecy2resOuter[2]);  //-> Curv  P2outer - parabola
      }
      
  }
  //
  // Fill raw THnSparses
  //
  for (Int_t irow=0;irow<159;irow++) {
    AliTPCclusterMI *c=track->GetClusterPointer(irow);
    if (!c) continue;
    if (c->GetMax()>800) continue; // saturation cut
    //if (TMath::Sqrt(TMath::Abs(c->GetSigmaY2()))>1) continue; // saturation cut
    //
    Double_t deltaY=c->GetY()-(*ltrp->GetVecLY())[irow];
    Double_t deltaZ=c->GetZ()-(*ltrp->GetVecLZ())[irow];
    //TString axisName[6]={"Delta","bin", "rms shape", "Q", "row","trackID"}
    Double_t xyz[6]={0, 0, 0,TMath::Sqrt(c->GetMax()),irow,id};
    xyz[0]=deltaY;
    xyz[1]=c->GetPad();
    xyz[2]=TMath::Sqrt(TMath::Abs(c->GetSigmaY2()));
    fHisLaserPad->Fill(xyz);
    xyz[0]=deltaZ;
    xyz[1]=c->GetTimeBin();
    xyz[2]=TMath::Sqrt(TMath::Abs(c->GetSigmaZ2()));
    fHisLaserTime->Fill(xyz);
  }
}



void AliTPCcalibLaser::DumpMeanInfo(Int_t run){
  //
  //  Dump information about laser beams
  //  isOK variable indicates usability of the beam  
  //  Beam is not usable if:
  //     a.  No entries in range (krmsCut0)
  //     b.  Big sperad          (krmscut1)
  //     c.  RMSto fit sigma bigger then (kmultiCut)
  //     d.  Too big angular spread 
  //  

  const Float_t krmsCut0=0.001;
  const Float_t krmsCut1=0.16;
  const Float_t kmultiCut=2;
  const Float_t kcutP0=0.002;
  AliMagF* magF=  dynamic_cast<AliMagF*> (TGeoGlobalMagField::Instance()->GetField());
  Double_t xyz[3]={90,0,10};         // tmp. global position 
  Double_t bxyz[3]={90,0,10};        // tmp. mag field  integral - cylindrical
  Double_t bgxyz[3]={90,0,10};       // tmp. mag field  integral - local
  //
  AliTPCcalibLaser *laser = this;
  TTreeSRedirector *pcstream = new TTreeSRedirector("laserMean.root");
  TF1 fg("fg","gaus");
  AliTPCParam  * tpcparam    = 0;   
  // start set up for absolute residuals analysis
  //
  AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();
  tpcparam    = calib->GetParameters(); 
  if (!tpcparam) tpcparam    = new AliTPCParamSR;
  tpcparam->Update();
  AliGRPObject *grp = AliTPCcalibDB::GetGRP(run);
  Float_t current=0;
  Float_t bfield      = 0, bz=0;

  if (grp){
    Float_t polarity = (grp->GetL3Polarity()>0) ? -1.:1;
    current = grp->GetL3Current((AliGRPObject::Stats)0);
    bfield = polarity*5*current/30000.;
    bz = polarity*5*current/30000.;
    printf("Run%d\tL3 current%f\tBz\t%f\n",run,current,bz);
  }

  SetBeamParameters(fBeamOffsetYOuter, fBeamSlopeYOuter, fBeamSectorOuter,0);
  SetBeamParameters(fBeamOffsetYInner, fBeamSlopeYInner, fBeamSectorInner,1);
  TLinearFitter lfabsyInner(2);
  lfabsyInner.SetFormula("1 ++ x");
  TLinearFitter lfabszInner(2);
  lfabszInner.SetFormula("1 ++ x");

  TLinearFitter lfabsyOuter(2);
  lfabsyOuter.SetFormula("1 ++ x");
  TLinearFitter lfabszOuter(2);
  lfabszOuter.SetFormula("1 ++ x");
  // end set up for absolute residuals analysis

  //
  //
  for (Int_t id=0; id<336; id++){
    Bool_t isOK=kTRUE;
    TH1F * hisphi  = (TH1F*)laser->fDeltaPhi.At(id);
    TH1F * hisphiP = (TH1F*)laser->fDeltaPhiP.At(id);
    TH1F * hisZ    = (TH1F*)laser->fDeltaZ.At(id);
    TH1F * hisP3    = (TH1F*)laser->fDeltaP3.At(id);
    TH1F * hisP4    = (TH1F*)laser->fDeltaP4.At(id);
    TH1F * hisS    = (TH1F*)laser->fSignals.At(id);
    //if (!hisphi) continue;
    Double_t entries = (hisphi==0)? 0: hisphi->GetEntries();
    //if (entries<minEntries) continue;
    //
    AliTPCLaserTrack *ltrp = (AliTPCLaserTrack*)fTracksMirror.At(id);
    if (!ltrp) {
     AliTPCLaserTrack::LoadTracks();
      ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    }
    ltrp->UpdatePoints();
    pcstream->GetFile()->cd();
    if (hisphi)  hisphi->Write();
    if (hisphiP) hisphiP->Write();
    if (hisZ)    hisZ->Write();
    if (hisP3)    hisP3->Write();
    if (hisP4)    hisP4->Write();
    
    Float_t meanphi = hisphi->GetMean();
    Float_t rmsphi = hisphi->GetRMS();
    //
    Float_t meanphiP = hisphiP->GetMean();
    Float_t rmsphiP = hisphiP->GetRMS();
    Float_t meanZ = hisZ->GetMean();
    Float_t rmsZ = hisZ->GetRMS();
    if (hisphi->GetRMS()>hisphi->GetBinWidth(1))
      hisphi->Fit(&fg,"","",hisphi->GetMean()-4*hisphi->GetRMS(),hisphi->GetMean()+4*hisphi->GetRMS());
    Double_t gphi1 = fg.GetParameter(1);
    Double_t gphi2 = fg.GetParameter(2);
    if (hisphiP->GetRMS()>0)
      hisphiP->Fit(&fg,"","",hisphiP->GetMean()-4*hisphiP->GetRMS(),hisphiP->GetMean()+4*hisphiP->GetRMS());
    Double_t gphiP1 = fg.GetParameter(1);
    Double_t gphiP2 = fg.GetParameter(2);
    //
    if (hisZ->GetRMS()>hisZ->GetBinWidth(1))
      hisZ->Fit(&fg,"","");
    Double_t gz1 = fg.GetParameter(1);
    Double_t gz2 = fg.GetParameter(2);
    //
    if (hisP3->GetRMS()>hisP3->GetBinWidth(1))
      hisP3->Fit(&fg,"","",hisP3->GetMean()-4*hisP3->GetRMS(),hisP3->GetMean()+4*hisP3->GetRMS());
    Double_t gp31 = fg.GetParameter(1);
    Double_t gp32 = fg.GetParameter(2);
    Double_t meanp3 = hisP3->GetMean();
    Double_t rmsp3  = hisP3->GetRMS();
    //
    if (hisP4->GetRMS()>hisP4->GetBinWidth(1))
      hisP4->Fit(&fg,"","",hisP4->GetMean()-4*hisP4->GetRMS(),hisP4->GetMean()+4*hisP4->GetRMS());
    Double_t gp41 = fg.GetParameter(1);
    Double_t gp42 = fg.GetParameter(2);
    Double_t meanp4 = hisP4->GetMean();
    Double_t rmsp4  = hisP4->GetRMS();
    //
    Float_t meanS=hisS->GetMean();
    //
    Double_t lxyz[3];
    Double_t lpxyz[3];
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);

    if (rmsphi<krmsCut0) isOK=kFALSE; // empty in range - not entries inside
    if (rmsphi>krmsCut1) isOK=kFALSE; // empty in range - not entries inside
    if (rmsphi>krmsCut0) if (gphi2/rmsphi>kmultiCut) isOK=kFALSE;   // multi peak structure
    if (gphiP2>kcutP0) isOK=kFALSE;
    //
    //
    //
    TH1 *his =0;
    //
    his = fHisNclIn->ProjectionY("aaa",id+1,id+1);
    Float_t mnclIn =  his->GetMean();
    delete his;
    his = fHisNclOut->ProjectionY("aaa",id+1,id+1);
    Float_t mnclOut =  his->GetMean();
    delete his;
    his = fHisNclIO->ProjectionY("aaa",id+1,id+1);
    Float_t mnclIO =  his->GetMean();
    delete his;
    his = fHisLclIn->ProjectionY("aaa",id+1,id+1);
    Float_t mLclIn =  his->GetMean();
    delete his;
    his = fHisLclOut->ProjectionY("aaa",id+1,id+1);
    Float_t mLclOut =  his->GetMean();
    delete his;
    his = fHisLclIO->ProjectionY("aaa",id+1,id+1);
    Float_t mLclIO =  his->GetMean();
    delete his;
    //
    his = fHisdEdx->ProjectionY("aaa",id+1,id+1);
    Float_t mdEdx =  his->GetMean();
    delete his;
    //
    //
    //
    //
    his = fHisChi2YIn1->ProjectionY("aaa",id+1,id+1);      //->chi2 y inner - line
    Float_t  mChi2YIn1=  his->GetMean();
    delete his;
    his = fHisChi2YOut1->ProjectionY("aaa",id+1,id+1);     //->chi2 y inner - line
    Float_t mChi2YOut1 =  his->GetMean();
    delete his;
    his = fHisChi2YIn2->ProjectionY("aaa",id+1,id+1);      //->chi2 y inner - parabola
    Float_t mChi2YIn2 =  his->GetMean();
    delete his;
    his = fHisChi2YOut2->ProjectionY("aaa",id+1,id+1);     //->chi2 y inner - parabola
    Float_t mChi2YOut2 =  his->GetMean();
    delete his;
    his = fHisChi2YIO1->ProjectionY("aaa",id+1,id+1);      //->chi2 y IO    - common
    Float_t mChi2YIO1 =  his->GetMean();
    delete his;
    his = fHisChi2ZIn1->ProjectionY("aaa",id+1,id+1);      //->chi2 z inner - line
    Float_t mChi2ZIn1 =  his->GetMean();
    delete his;
    his = fHisChi2ZOut1->ProjectionY("aaa",id+1,id+1);     //->chi2 z inner - line
    Float_t mChi2ZOut1 =  his->GetMean();
    delete his;
    his = fHisChi2ZIn2->ProjectionY("aaa",id+1,id+1);      //->chi2 z inner - parabola
    Float_t mChi2ZIn2 =  his->GetMean();
    delete his;
    his = fHisChi2ZOut2->ProjectionY("aaa",id+1,id+1);     //->chi2 z inner - parabola
    Float_t mChi2ZOut2 =  his->GetMean();
    delete his;
    his = fHisChi2ZIO1->ProjectionY("aaa",id+1,id+1);      //->chi2 z IO    - common
    Float_t mChi2ZIO1 =  his->GetMean();
    delete his;
    //
    // fit res. histos
    // 
    his = fHisdZfit->ProjectionY("aaa",id+1,id+1);
    Float_t  edZfit  = his->GetEntries();
    Float_t  mdZfit =  his->GetMean();
    Float_t  rdZfit  = his->GetRMS();
    delete his;

    his = fHisPy1vP0->ProjectionY("aaa",id+1,id+1);     //-> delta y   P0outer-P0inner - line
    Float_t  ePy1vP0  = his->GetEntries();
    Float_t  mPy1vP0  = his->GetMean();
    Float_t  rPy1vP0  = his->GetRMS();
    delete his;
    
    his = fHisPy2vP0->ProjectionY("aaa",id+1,id+1);     //-> delta y   P0outer-P0inner - parabola
    Float_t  ePy2vP0  = his->GetEntries();
    Float_t  mPy2vP0  = his->GetMean();
    Float_t  rPy2vP0  = his->GetRMS();
    delete his;
    
    his = fHisPy3vP0->ProjectionY("aaa",id+1,id+1);     //-> delta y   P0outer-P0inner - common parabola
    Float_t  ePy3vP0  = his->GetEntries();
    Float_t  mPy3vP0  = his->GetMean();
    Float_t  rPy3vP0  = his->GetRMS();
    delete his;
    
    his = fHisPy1vP1->ProjectionY("aaa",id+1,id+1);     //-> delta ky  P1outer-P1inner - line
    Float_t  ePy1vP1  = his->GetEntries();
    Float_t  mPy1vP1  = his->GetMean();
    Float_t  rPy1vP1  = his->GetRMS();
    delete his;
    
    his = fHisPy2vP1->ProjectionY("aaa",id+1,id+1);     //-> delta ky  P1outer-P1inner - parabola
    Float_t  ePy2vP1  = his->GetEntries();
    Float_t  mPy2vP1  = his->GetMean();
    Float_t  rPy2vP1  = his->GetRMS();
    delete his;
    
    his = fHisPy3vP1->ProjectionY("aaa",id+1,id+1);     //-> delta ky  P1outer-P1inner - common parabola
    Float_t  ePy3vP1  = his->GetEntries();
    Float_t  mPy3vP1  = his->GetMean();
    Float_t  rPy3vP1  = his->GetRMS();
    delete his;
    
    his = fHisPy2vP2In->ProjectionY("aaa",id+1,id+1);   //-> Curv  P2inner - parabola
    Float_t  ePy2vP2In  = his->GetEntries();
    Float_t  mPy2vP2In  = his->GetMean();
    Float_t  rPy2vP2In  = his->GetRMS();
    delete his;
    
    his = fHisPy2vP2Out->ProjectionY("aaa",id+1,id+1);  //-> Curv  P2outer - parabola
    Float_t  ePy2vP2Out  = his->GetEntries();
    Float_t  mPy2vP2Out  = his->GetMean();
    Float_t  rPy2vP2Out  = his->GetRMS();
    delete his;
    
    his = fHisPy3vP2IO->ProjectionY("aaa",id+1,id+1);   //-> Curv  P2outerinner - common parabola
    Float_t  ePy3vP2IO  = his->GetEntries();
    Float_t  mPy3vP2IO  = his->GetMean();
    Float_t  rPy3vP2IO  = his->GetRMS();
    delete his;
    
    //
    //
    his = fHisPz1vP0->ProjectionY("aaa",id+1,id+1);     //-> delta z   P0outer-P0inner - line
    Float_t  ePz1vP0  = his->GetEntries();
    Float_t  mPz1vP0  = his->GetMean();
    Float_t  rPz1vP0  = his->GetRMS();
    delete his;
    
    his = fHisPz2vP0->ProjectionY("aaa",id+1,id+1);     //-> delta z   P0outer-P0inner - parabola
    Float_t  ePz2vP0  = his->GetEntries();
    Float_t  mPz2vP0  = his->GetMean();
    Float_t  rPz2vP0  = his->GetRMS();
    delete his;
    
    his = fHisPz3vP0->ProjectionY("aaa",id+1,id+1);     //-> delta z   P0outer-P0inner - common parabola
    Float_t  ePz3vP0  = his->GetEntries();
    Float_t  mPz3vP0  = his->GetMean();
    Float_t  rPz3vP0  = his->GetRMS();
    delete his;
    
    his = fHisPz1vP1->ProjectionY("aaa",id+1,id+1);     //-> delta kz  P1outer-P1inner - line
    Float_t  ePz1vP1  = his->GetEntries();
    Float_t  mPz1vP1  = his->GetMean();
    Float_t  rPz1vP1  = his->GetRMS();
    delete his;
    
    his = fHisPz2vP1->ProjectionY("aaa",id+1,id+1);     //-> delta kz  P1outer-P1inner - parabola
    Float_t  ePz2vP1  = his->GetEntries();
    Float_t  mPz2vP1  = his->GetMean();
    Float_t  rPz2vP1  = his->GetRMS();
    delete his;
    
    his = fHisPz3vP1->ProjectionY("aaa",id+1,id+1);     //-> delta kz  P1outer-P1inner - common parabola
    Float_t  ePz3vP1  = his->GetEntries();
    Float_t  mPz3vP1  = his->GetMean();
    Float_t  rPz3vP1  = his->GetRMS();
    delete his;
    
    his = fHisPz2vP2In->ProjectionY("aaa",id+1,id+1);   //-> Curv  P2inner - parabola
    Float_t  ePz2vP2In  = his->GetEntries();
    Float_t  mPz2vP2In  = his->GetMean();
    Float_t  rPz2vP2In  = his->GetRMS();
    delete his;
    
    his = fHisPz2vP2Out->ProjectionY("aaa",id+1,id+1);  //-> Curv  P2outer - parabola
    Float_t  ePz2vP2Out  = his->GetEntries();
    Float_t  mPz2vP2Out  = his->GetMean();
    Float_t  rPz2vP2Out  = his->GetRMS();
    delete his;
    
    his = fHisPz3vP2IO->ProjectionY("aaa",id+1,id+1);   //-> Curv  P2outerinner - common parabola
    Float_t  ePz3vP2IO  = his->GetEntries();
    Float_t  mPz3vP2IO  = his->GetMean();
    Float_t  rPz3vP2IO  = his->GetRMS();
    delete his;

    // Fit absolute laser residuals
    TH2F* histAbsY = (TH2F*)laser->fDeltaYresAbs[id];
    TH2F* histAbsZ = (TH2F*)laser->fDeltaZresAbs[id];

    Int_t secInner = TMath::Nint(fBeamSectorInner[id]); 
    Int_t secOuter = TMath::Nint(fBeamSectorOuter[id]); 

    TVectorD vecX(159);        // X
    TVectorD vecY(159);        // Y
    TVectorD vecR(159);        // R
    TVectorD vecDY(159);       // absolute residuals in Y
    TVectorD vecDZ(159);       // absolute residuals in Z
    TVectorD vecN(159);        // number of clusters
    TVectorD vecEy(159);       //error y
    TVectorD vecEz(159);       //error z
    TVectorD vecPhi(159);      // local tangent
    TVectorD vecPhiR(159);     // local tangent
    // magnetic field integrals
    TVectorD vecIBR(159);        // radial
    TVectorD vecIBRPhi(159);     // r-phi
    TVectorD vecIBLX(159);       // local x
    TVectorD vecIBLY(159);       // local y
    TVectorD vecIBGX(159);       // local x
    TVectorD vecIBGY(159);       // local y
    TVectorD vecIBZ(159);        // z
    //
    for (Int_t irow=0;irow<159;irow++){
      vecIBR[irow]=0;
      vecIBRPhi[irow]=0;
      vecIBLX[irow]=0;
      vecIBLY[irow]=0;
      vecIBGX[irow]=0;
      vecIBGY[irow]=0;
      vecIBZ[irow]=0;
      Double_t gx    =(*(ltrp->fVecGX))[irow];
      Double_t gy    =(*(ltrp->fVecGY))[irow];
      Int_t    lsec  =TMath::Nint((*(ltrp->fVecSec))[irow]);
      Double_t   ca  =TMath::Cos(TMath::Pi()*(lsec+0.5)/9.);
      Double_t   sa  =TMath::Sin(TMath::Pi()*(lsec+0.5)/9.);
      xyz[2]=(*(ltrp->fVecGZ))[irow];
      xyz[0]=TMath::Sqrt(gx*gx+gy*gy);
      xyz[1]=TMath::ATan2(gy,gx);
      Double_t gxyz[3]={gx,gy,(*(ltrp->fVecGZ))[irow]};
      if (magF){
	magF->GetTPCIntCyl(xyz,bxyz);
	magF->GetTPCInt(gxyz,bgxyz);
	vecIBR[irow]=bxyz[0];
	vecIBRPhi[irow]=bxyz[1];
	//
	vecIBGX[irow]=bgxyz[0];
	vecIBGY[irow]=bgxyz[1];
	//
	vecIBLX[irow]=  bgxyz[0]*ca+bgxyz[1]*sa;
	vecIBLY[irow]= -bgxyz[0]*sa+bgxyz[1]*ca;
	//

	vecIBZ[irow]=bxyz[2];
      }
    }


    lfabsyInner.ClearPoints();    
    lfabszInner.ClearPoints();    
    lfabsyOuter.ClearPoints();    
    lfabszOuter.ClearPoints();    
    // dummy fit values
    Int_t    nClInY         = 0;
    Float_t  yAbsInOffset  = -100;
    Float_t  yAbsInSlope   = -100;
    Float_t  yAbsInDeltay  = -100;
    Int_t    nClInZ         = 0;
    Float_t  zAbsInOffset  = -100;
    Float_t  zAbsInSlope   = -100;
    Float_t  zAbsInDeltaz  = -100;
    Int_t    nClOutY         = 0;
    Float_t  yAbsOutOffset  = -100;
    Float_t  yAbsOutSlope   = -100;
    Float_t  yAbsOutDeltay  = -100;
    Int_t    nClOutZ         = 0;
    Float_t  zAbsOutOffset  = -100;
    Float_t  zAbsOutSlope   = -100;
    Float_t  zAbsOutDeltaz  = -100;

    Float_t  lasTanPhiLocIn = -100;
    Float_t  lasTanPhiLocOut = -100;

    if(histAbsY && histAbsY->GetEntries()>0) {
      
      Double_t rotAngOut = 10;
      Double_t rotAngIn = 10;
      if((secInner%36)!=(secOuter%36))
	rotAngIn += 20; // 30 degrees
      
      // Calculate laser mirror X position in local frame
      Double_t laserposOut = 
	TMath::Abs(ltrp->GetX()*TMath::Cos(rotAngOut*TMath::DegToRad()));
      Double_t laserposIn = 
	TMath::Abs(ltrp->GetX()*TMath::Cos(rotAngIn*TMath::DegToRad()));
      
      // Calculate laser tan phi in local frame
      lasTanPhiLocOut  = TMath::ASin(ltrp->GetSnp());
      if(lasTanPhiLocOut<0) {
	lasTanPhiLocIn  = lasTanPhiLocOut - rotAngIn*TMath::DegToRad();
	lasTanPhiLocOut -= rotAngOut*TMath::DegToRad();
      } else {
	
	lasTanPhiLocIn  = lasTanPhiLocOut + rotAngIn*TMath::DegToRad();
	lasTanPhiLocOut += rotAngOut*TMath::DegToRad();
      }
      
      lasTanPhiLocIn = TMath::Tan(lasTanPhiLocIn);
      lasTanPhiLocOut = TMath::Tan(lasTanPhiLocOut);
      
      TProfile* yprof  = histAbsY->ProfileX("yprof");
      TProfile* zprof  = histAbsZ->ProfileX("zprof");
      
      for(Int_t bin = 1; bin<=159; bin++) {
	
	if(yprof->GetBinEntries(bin)<5&&
	   zprof->GetBinEntries(bin)<5) {
	  continue;	
	}
	
	// There is a problem in defining inner and outer sectors for
	// the outer beams (0 and 6) where both sectors are OROCs. To
	// make sure there is no overlap row 94 to 99 are cutted.
	if(((ltrp->GetBeam()==0)||(ltrp->GetBeam()==6))&&bin>=95&&bin<=100)
	  continue;
	
	Int_t row = (bin-1);
	if(row>62)
	  row -= 63;

	Bool_t isOuter = kTRUE;
	Int_t sector = TMath::Nint(fBeamSectorOuter[id]);
	
	if(bin<=62 ||                                        
	   (((ltrp->GetBeam()==0)||(ltrp->GetBeam()==6))&&bin<=95)) {

	  isOuter = kFALSE;
	  sector = TMath::Nint(fBeamSectorInner[id]);
	}


	Double_t x = tpcparam->GetPadRowRadii(sector, row); // slope
	vecN[bin-1] =yprof->GetBinEntries(bin);
	vecEy[bin-1]=yprof->GetBinError(bin);
	vecEz[bin-1]=zprof->GetBinError(bin);
	vecX[bin-1] = x;
	vecDY[bin-1] = yprof->GetBinContent(bin);
	vecDZ[bin-1] = zprof->GetBinContent(bin);
	if (bin>0&&bin<159){
	  //
	  //truncated mean - skip first and the last  pad row
	  //
	  Int_t firstBin=TMath::Max(bin-5,0);
	  Int_t lastBin =TMath::Min(bin+5,158);
	  histAbsY->GetXaxis()->SetRangeUser(firstBin,lastBin);
	  histAbsY->GetYaxis()->SetRangeUser(-2,2);
	  vecEy[bin-1]=histAbsY->GetRMS(2);
	  vecDY[bin-1]=histAbsY->GetMean(2);
	  histAbsY->GetXaxis()->SetRangeUser(firstBin+2,lastBin-2);//use+-2 bins
	  histAbsY->GetYaxis()->SetRangeUser(vecDY[bin-1]-4*vecEy[bin-1],vecDY[bin-1]+4*vecEy[bin-1]);
	  if (yprof->GetBinEntries(bin-1)>0) vecEy[bin-1]=histAbsY->GetRMS(2)/TMath::Sqrt(yprof->GetBinEntries(bin-1));
	  vecDY[bin-1]=histAbsY->GetMean(2);
	}

	if(!isOuter) { // inner	  
	  vecPhi[bin-1]=lasTanPhiLocIn;
	  // Calculate local y from residual and database
	  Double_t y  = fBeamSlopeYInner[id]*x + fBeamOffsetYInner[id]
	    + vecDY[bin-1];
	  vecY[bin-1] = y;
	  Double_t r  = TMath::Sqrt(x*x + y*y);
	  vecR[bin-1] = r;
	  // Find angle between laser vector and R vector
	  // cos angle = (x,y).(1,fBeamSlopeYInner)/R/sqrt(1+fBeamSlopeYInner**2)
	  Double_t cosPhi = x + y*fBeamSlopeYInner[id];
	  cosPhi /= r*TMath::Sqrt(1+fBeamSlopeYInner[id]*fBeamSlopeYInner[id]);
	  vecPhiR[bin-1] = TMath::Tan(TMath::ACos(cosPhi));
	  if(lasTanPhiLocIn<0)
	    vecPhiR[bin-1]*=-1; // to have the same sign
	  
	  if(yprof->GetBinEntries(bin)>=10) {
	    lfabsyInner.AddPoint(&x, yprof->GetBinContent(bin), 
				  TMath::Max(yprof->GetBinError(bin), 0.001));
	  } 
	  if(zprof->GetBinEntries(bin)>=10) {
	    lfabszInner.AddPoint(&x, zprof->GetBinContent(bin), 
				  TMath::Max(zprof->GetBinError(bin), 0.001));
	  }
	} else { // outer
	  vecPhi[bin-1]=lasTanPhiLocOut;	  
	  // Calculate local y from residual and database
	  Double_t y  = fBeamSlopeYOuter[id]*x + fBeamOffsetYOuter[id]
	    + vecDY[bin-1];
	  vecY[bin-1] = y;
	  Double_t r  = TMath::Sqrt(x*x + y*y);
	  vecR[bin-1] = r;

	  Double_t cosPhi = x + y*fBeamSlopeYOuter[id];
	  cosPhi /= r*TMath::Sqrt(1+fBeamSlopeYOuter[id]*fBeamSlopeYOuter[id]);
	  vecPhiR[bin-1] = TMath::Tan(TMath::ACos(cosPhi));
	  if(lasTanPhiLocOut<0)
	    vecPhiR[bin-1]*=-1; // to have the same sign

	  if(yprof->GetBinEntries(bin)>=10) {
	    lfabsyOuter.AddPoint(&x, yprof->GetBinContent(bin), 
				  TMath::Max(yprof->GetBinError(bin), 0.001));
	  }
	  if(zprof->GetBinEntries(bin)>=10) {
	    lfabszOuter.AddPoint(&x, zprof->GetBinContent(bin), 
				  TMath::Max(zprof->GetBinError(bin), 0.001));
	  }
	}
	// global position
	
      }
	
      delete yprof; delete zprof;

      
      // Fit laser abs residuals with linear robust fit (exclude 5% outliers)
      nClInY = lfabsyInner.GetNpoints();
      if(lfabsyInner.GetNpoints()>10) {
	lfabsyInner.EvalRobust(0.95);
	
	TVectorD result(2);
	lfabsyInner.GetParameters(result);
	yAbsInOffset = result[0];
	yAbsInSlope  = result[1];
	yAbsInDeltay = yAbsInSlope*laserposIn + yAbsInOffset;
      }
      nClInZ = lfabszInner.GetNpoints();
      if(lfabszInner.GetNpoints()>10) {
	lfabszInner.EvalRobust(0.95);
	
	TVectorD result(2);
	lfabszInner.GetParameters(result);
	zAbsInOffset = result[0];
	zAbsInSlope  = result[1];
	zAbsInDeltaz = zAbsInSlope*laserposIn + zAbsInOffset;
      } 
      nClOutY = lfabsyOuter.GetNpoints();
      if(lfabsyOuter.GetNpoints()>10) {
	lfabsyOuter.EvalRobust(0.95);
	
	TVectorD result(2);
	lfabsyOuter.GetParameters(result);
	yAbsOutOffset = result[0];
	yAbsOutSlope  = result[1];
	yAbsOutDeltay = yAbsOutSlope*laserposOut + yAbsOutOffset;
      }
      nClOutZ = lfabszOuter.GetNpoints();
      if(lfabszOuter.GetNpoints()>10) {
	lfabszOuter.EvalRobust(0.95);
	
	TVectorD result(2);
	lfabszOuter.GetParameters(result);
	zAbsOutOffset = result[0];
	zAbsOutSlope  = result[1];
	zAbsOutDeltaz = zAbsOutSlope*laserposOut + zAbsOutOffset;
      }
    }

    
    Int_t itime=-1;
    Float_t  coverIA=AliTPCcalibDB::GetCoverVoltage(run,0,itime);
    Float_t  coverIC=AliTPCcalibDB::GetCoverVoltage(run,18,itime);
    Float_t  coverOA=AliTPCcalibDB::GetCoverVoltage(run,36,itime);
    Float_t  coverOC=AliTPCcalibDB::GetCoverVoltage(run,54,itime);
    //
    Float_t  skirtA=AliTPCcalibDB::GetSkirtVoltage(run,0,itime);
    Float_t  skirtC=AliTPCcalibDB::GetSkirtVoltage(run,18,itime);
    //
    Float_t  ggOffA=AliTPCcalibDB::GetGGoffsetVoltage(run,0,itime);
    Float_t  ggOffC=AliTPCcalibDB::GetGGoffsetVoltage(run,18,itime);

    //
    (*pcstream)<<"Mean"<<
      "run="<<run<<               //
      //voltages
      "UcIA="<<coverIA<<
      "UcIC="<<coverIC<<
      "UcOA="<<coverOA<<
      "UcOC="<<coverOC<<
      "UsA="<<skirtA<<
      "UsC="<<skirtC<<
      "UggA="<<ggOffA<<
      "UggC="<<ggOffC<<
      //
      "isOK="<<isOK<<             //
      "id="<<id<<                 // track id
      "entries="<<entries<<       // number of entries
      "bz="<<bfield<<             // bfield
      "LTr.="<<ltrp<<             // refernece track
      //
      "mnclIn="<<mnclIn<<         // mean number of clusters in inner
      "mnclOut="<<mnclOut<<       // mean number of clusters in outer
      "mnclIO="<<mnclIO<<         // mean number of clusters in inner+outer
      "mLclIn="<<mLclIn<<         // mean number of clusters in inner
      "mLclOut="<<mLclOut<<       // mean number of clusters in outer
      "mLclIO="<<mLclIO<<         // mean number of clusters in inner+outer
      "mdEdx="<<mdEdx<<           // mean dEdx
      "edZfit="<<edZfit<<           // entries z fit
      "mdZfit="<<mdZfit<<           // mean z fit
      "rdZfit="<<rdZfit<<           // RMS z fit
      //
      //
      "mChi2YIn1="<<mChi2YIn1<<       //->chi2 y inner - line
      "mChi2YOut1="<<mChi2YOut1<<     //->chi2 y inner - line
      "mChi2YIn2="<<mChi2YIn2<<       //->chi2 y inner - parabola
      "mChi2YOut2="<<mChi2YOut2<<     //->chi2 y inner - parabola
      "mChi2YIO1="<<mChi2YIO1<<       //->chi2 y IO    - common
      "mChi2ZIn1="<<mChi2ZIn1<<       //->chi2 z inner - line
      "mChi2ZOut1="<<mChi2ZOut1<<     //->chi2 z inner - line
      "mChi2ZIn2="<<mChi2ZIn2<<       //->chi2 z inner - parabola
      "mChi2ZOut2="<<mChi2ZOut2<<     //->chi2 z inner - parabola
      "mChi2ZIO1="<<mChi2ZIO1<<       //->chi2 z IO    - common
      //
      //
      "ePy1vP0="<<ePy1vP0<<
      "mPy1vP0="<<mPy1vP0<<
      "rPy1vP0="<<rPy1vP0<<
      "ePy2vP0="<<ePy2vP0<<
      "mPy2vP0="<<mPy2vP0<<
      "rPy2vP0="<<rPy2vP0<<      
      "ePy3vP0="<<ePy3vP0<<
      "mPy3vP0="<<mPy3vP0<<
      "rPy3vP0="<<rPy3vP0<<
      "ePy1vP1="<<ePy1vP1<<
      "mPy1vP1="<<mPy1vP1<<
      "rPy1vP1="<<rPy1vP1<<
      "ePy2vP1="<<ePy2vP1<<
      "mPy2vP1="<<mPy2vP1<<
      "rPy2vP1="<<rPy2vP1<<
      "ePy3vP1="<<ePy3vP1<<
      "mPy3vP1="<<mPy3vP1<<
      "rPy3vP1="<<rPy3vP1<<
      "ePy2vP2In="<<ePy2vP2In<<
      "mPy2vP2In="<<mPy2vP2In<<
      "rPy2vP2In="<<rPy2vP2In<<
      "ePy2vP2Out="<<ePy2vP2Out<<
      "mPy2vP2Out="<<mPy2vP2Out<<
      "rPy2vP2Out="<<rPy2vP2Out<<
      "ePy3vP2IO="<<ePy3vP2IO<<
      "mPy3vP2IO="<<mPy3vP2IO<<
      "rPy3vP2IO="<<rPy3vP2IO<<
      //
      //
      "ePz1vP0="<<ePz1vP0<<
      "mPz1vP0="<<mPz1vP0<<
      "rPz1vP0="<<rPz1vP0<<
      "ePz2vP0="<<ePz2vP0<<
      "mPz2vP0="<<mPz2vP0<<
      "rPz2vP0="<<rPz2vP0<<
      "ePz3vP0="<<ePz3vP0<<
      "mPz3vP0="<<mPz3vP0<<
      "rPz3vP0="<<rPz3vP0<<
      "ePz1vP1="<<ePz1vP1<<
      "mPz1vP1="<<mPz1vP1<<
      "rPz1vP1="<<rPz1vP1<<
      "ePz2vP1="<<ePz2vP1<<
      "mPz2vP1="<<mPz2vP1<<
      "rPz2vP1="<<rPz2vP1<<
      "ePz3vP1="<<ePz3vP1<<
      "mPz3vP1="<<mPz3vP1<<
      "rPz3vP1="<<rPz3vP1<<
      "ePz2vP2In="<<ePz2vP2In<<
      "mPz2vP2In="<<mPz2vP2In<<
      "rPz2vP2In="<<rPz2vP2In<<
      "ePz2vP2Out="<<ePz2vP2Out<<
      "mPz2vP2Out="<<mPz2vP2Out<<
      "rPz2vP2Out="<<rPz2vP2Out<<
      "ePz3vP2IO="<<ePz3vP2IO<<  
      "mPz3vP2IO="<<mPz3vP2IO<<
      "rPz3vP2IO="<<rPz3vP2IO<<     
      //
      //
      //
      "lx0="<<lxyz[0]<<            // reference x
      "lx1="<<lxyz[1]<<            // reference y
      "lx2="<<lxyz[2]<<            // refernece z
      "lpx0="<<lpxyz[0]<<           // reference x
      "lpx1="<<lpxyz[1]<<          // reference y
      "lpx2="<<lpxyz[2]<<          // refernece z
      //
      "msig="<<meanS<<
      //
      "mphi="<<meanphi<<         //
      "rmsphi="<<rmsphi<<        //
      "gphi1="<<gphi1<<
      "gphi2="<<gphi2<<
      //
      "mphiP="<<meanphiP<<       //
      "rmsphiP="<<rmsphiP<<      //
      "gphiP1="<<gphiP1<<
      "gphiP2="<<gphiP2<<
      //
      "meanZ="<<meanZ<<
      "rmsZ="<<rmsZ<<
      "gz1="<<gz1<<
      "gz2="<<gz2<<
      //
      "gp31="<<gp31<<            //gaus mean - tgl
      "gp32="<<gp32<<            //gaus rms  -tgl
      "meanp3="<<meanp3<<            //mean - tgl
      "rmsp3="<<rmsp3<<            //rms  -tgl
      "gp41="<<gp41<<            //gaus mean - P4
      "gp42="<<gp42<<            //gaus rms  - P4
      "meanp4="<<meanp4<<            //mean - P4
      "rmsp4="<<rmsp4<<            //rms  - P4
      // Parameters from abs res analysis
      "SecIn="<<secInner<<              // inner sector
      "SecOut="<<secOuter<<             // outer sector
      "lasTanPhiLocIn="<<lasTanPhiLocIn<< // laser tan phi in local frame (inner)
      "lasTanPhiLocOut="<<lasTanPhiLocOut<<// laser tan phi in local frame (outer)
      "ibr.="<<&vecIBR<<   // radial filed integral
      "ibrphi.="<<&vecIBRPhi<<   // r=phifiled integral
      "ibr.="<<&vecIBR<<   // radial filed integral
      "ibz.="<<&vecIBZ<<   // z filed integral
      //
      "iblx.="<<&vecIBLX<<   // local bx  integral
      "ibly.="<<&vecIBLY<<   // local by integral
      "ibgx.="<<&vecIBGX<<   // global bx  integral
      "ibgy.="<<&vecIBGY<<   // global by integral
      //
      "X.="<<&vecX<<       // local x 
      "Y.="<<&vecY<<       // local y 
      "R.="<<&vecR<<       // radius 
      "dY.="<<&vecDY<<     // abs y residuals
      "dZ.="<<&vecDZ<<     // abs z residuals
      "eY.="<<&vecEy<<     // error of y residuals
      "eZ.="<<&vecEz<<     // error of z residuals
      "kY.="<<&vecPhi<<    // local tangent y (fixed for sector)
      "kYR.="<<&vecPhiR<<  // tangent between laser and R vector (varies inside sector) 
      "nCl.="<<&vecN<<     // number of clusters
      //
      "nClInY="<<nClInY<<               // Number of clusters for inner y
      "yAbsInOffset="<<yAbsInOffset<< // fitted offset absres (inner y)
      "yAbsInSlope="<<yAbsInSlope <<  // fitted slope absres (inner y)
      "yAbsInDeltay="<<yAbsInDeltay<< // fitted laser offset absres (inner y)
      "nClInZ="<<nClInZ<<               // Number of clusters for inner z
      "zAbsInOffset="<<zAbsInOffset<< // fitted offset absres (inner z)
      "zAbsInSlope="<<zAbsInSlope <<  // fitted slope absres (inner z)
      "zAbsInDeltaz="<<zAbsInDeltaz<< // fitted laser offset absres (inner z)
      //
      "nClOutY="<<nClOutY<<               // Number of clusters for outer y
      "yAbsOutOffset="<<yAbsOutOffset<< // fitted offset absres (outer y)
      "yAbsOutSlope="<<yAbsOutSlope <<  // fitted slope absres (outer y)
      "yAbsOutDeltay="<<yAbsOutDeltay<< // fitted laser offset absres (outer y)
      "nClOutZ="<<nClOutZ<<               // Number of clusters for outer z
      "zAbsOutOffset="<<zAbsOutOffset<< // fitted offset absres (outer z)
      "zAbsOutSlope="<<zAbsOutSlope <<  // fitted slope absres (outer z)
      "zAbsOutDeltaz="<<zAbsOutDeltaz<< // fitted laser offset absres (outer z)
      //
      "\n";
  }
  delete pcstream;
  /*
    Browse the content
    TFile fmean("laserMean.root");
    

   */


}



void AliTPCcalibLaser::DumpScanInfo(TTree * chain, const char * cutUser){
  //
  //
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector("laserScan.root");
  TFile * f = pcstream->GetFile();
  f->mkdir("dirphi");
  f->mkdir("dirphiP");
  f->mkdir("dirZ");
  TF1 fp("p1","pol1");
  //
  //
  char    cut[1000];
  char grname[1000];
  char grnamefull[1000];

  Double_t mphi[100];
  Double_t mphiP[100];
  Double_t smphi[100];
  Double_t smphiP[100];
  Double_t mZ[100];
  Double_t smZ[100];
  Double_t bz[100];
  Double_t sbz[100];
  // fit parameters
  Double_t pphi[3];
  Double_t pphiP[3];
  Double_t pmZ[3];
  
  //
  for (Int_t id=0; id<336; id++){
    // id =205;
    snprintf(cut,1000,"fId==%d&&%s",id,cutUser);
    Int_t entries = chain->Draw("bz",cut,"goff");
    if (entries<3) continue;
    AliTPCLaserTrack *ltrp = 0;
    if (!AliTPCLaserTrack::GetTracks()) AliTPCLaserTrack::LoadTracks();
    ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    Double_t lxyz[3];
    Double_t lpxyz[3];
    ltrp->GetXYZ(lxyz);
    ltrp->GetPxPyPz(lpxyz);

    chain->Draw("bz",cut,"goff");
    memcpy(bz, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.01*abs(bz)+0.02",cut,"goff");
    memcpy(sbz, chain->GetV1(), entries*sizeof(Double_t));
    //
    chain->Draw("gphi1",cut,"goff");
    memcpy(mphi, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.05*abs(mphi)+abs(gphi2)*0.5+0.05",cut,"goff");
    memcpy(smphi, chain->GetV1(), entries*sizeof(Double_t));
    //
    chain->Draw("gphiP1",cut,"goff");
    memcpy(mphiP, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.05*abs(mphiP)+abs(gphiP2)*0.5+0.001",cut,"goff");
    memcpy(smphiP, chain->GetV1(), entries*sizeof(Double_t));
    //
    chain->Draw("gz1",cut,"goff");
    memcpy(mZ, chain->GetV1(), entries*sizeof(Double_t));
    chain->Draw("0.01*abs(meanZ)+abs(gz2)*0.5+0.1",cut,"goff");
    memcpy(smZ, chain->GetV1(), entries*sizeof(Double_t));
    //
    //
    snprintf(grnamefull,1000,"Side_%d_Bundle_%d_Rod_%d_Beam_%d",
	    ltrp->GetSide(),  ltrp->GetBundle(), ltrp->GetRod(), ltrp->GetBeam());
    // store data  
    // phi
    f->cd("dirphi");
    Float_t phiB0 =0;
    Float_t phiPB0=0;
    Float_t zB0=0;
    for (Int_t ientry=0; ientry<entries; ientry++){
      if (TMath::Abs(bz[ientry])<0.05){
	phiB0  = mphi[ientry];
	phiPB0 = mphiP[ientry];
	zB0    = mZ[ientry];
      }
    }
    TGraphErrors *grphi = new TGraphErrors(entries,bz,mphi,sbz,smphi);
    grphi->Draw("a*");
    grphi->Fit(&fp);
    pphi[0] = fp.GetParameter(0);                          // offset
    pphi[1] = fp.GetParameter(1);                          // slope
    pphi[2] = TMath::Sqrt(fp.GetChisquare()/(entries-2.));  // normalized chi2
    snprintf(grname,1000,"phi_id%d",id);
    grphi->SetName(grname);  grphi->SetTitle(grnamefull);
    grphi->GetXaxis()->SetTitle("b_{z} (T)");
    grphi->GetYaxis()->SetTitle("#Delta r#phi (cm)");
    grphi->SetMaximum(1.2);
    grphi->SetMinimum(-1.2);
    grphi->Draw("a*");

    grphi->Write();
    gPad->SaveAs(Form("pic/phi/phi_%s.gif",grnamefull));
    // phiP
    f->cd("dirphiP)");
    TGraphErrors *grphiP = new TGraphErrors(entries,bz,mphiP,sbz,smphiP);
    grphiP->Draw("a*");
    grphiP->Fit(&fp);
    pphiP[0] = fp.GetParameter(0);                          // offset
    pphiP[1] = fp.GetParameter(1);                          // slope
    pphiP[2] = TMath::Sqrt(fp.GetChisquare()/(entries-2.));  // normalized chi2
    snprintf(grname,1000,"phiP_id%d",id);
    grphiP->SetName(grname);  grphiP->SetTitle(grnamefull);
    grphiP->GetXaxis()->SetTitle("b_{z} (T)");
    grphiP->GetYaxis()->SetTitle("#Delta #phi (rad)");
    grphiP->SetMaximum(pphiP[0]+0.005);
    grphiP->SetMinimum(pphiP[0]-0.005);

    gPad->SaveAs(Form("pic/phiP/phiP_%s.gif",grnamefull));
    grphiP->Write();
    //
    //Z
    f->cd("dirZ");
    TGraphErrors *grmZ = new TGraphErrors(entries,bz,mZ,sbz,smZ);
    grmZ->Draw("a*");
    grmZ->Fit(&fp);
    pmZ[0] = fp.GetParameter(0);                          // offset
    pmZ[1] = fp.GetParameter(1);                          // slope
    pmZ[2] = TMath::Sqrt(fp.GetChisquare()/(entries-2.));  // normalized chi2
    snprintf(grname,1000,"mZ_id%d",id);
    grmZ->SetName(grname);  grmZ->SetTitle(grnamefull);
    grmZ->GetXaxis()->SetTitle("b_{z} (T)");
    grmZ->GetYaxis()->SetTitle("#Delta z (cm)");

    gPad->SaveAs(Form("pic/z/z_%s.gif",grnamefull));
    grmZ->Write();
    //
    // P4
    //

    for (Int_t ientry=0; ientry<entries; ientry++){
      (*pcstream)<<"Mean"<<
	"id="<<id<<
	"LTr.="<<ltrp<<
	"entries="<<entries<<
	"bz="<<bz[ientry]<<
	"lx0="<<lxyz[0]<<          // reference x
	"lx1="<<lxyz[1]<<          // reference y
	"lx2="<<lxyz[2]<<          // refernece z      
	"lpx0="<<lpxyz[0]<<          // reference x
	"lpx1="<<lpxyz[1]<<          // reference y
	"lpx2="<<lpxyz[2]<<          // refernece z            
	//values
	"phiB0="<<phiB0<<          // position shift at 0 field
	"phiPB0="<<phiPB0<<        // angular  shift at 0 field
	"zB0="<<zB0<<              // z shift for 0 field
	//
	"gphi1="<<mphi[ientry]<< // mean - from gaus fit
	"pphi0="<<pphi[0]<<   // offset
	"pphi1="<<pphi[1]<<   // slope
	"pphi2="<<pphi[2]<<   // norm chi2
	//
	"gphiP1="<<mphiP[ientry]<< // mean - from gaus fit
	"pphiP0="<<pphiP[0]<< // offset
	"pphiP1="<<pphiP[1]<< // slope
	"pphiP2="<<pphiP[2]<< // norm chi2
	//
	"gz1="<<mZ[ientry]<<
	"pmZ0="<<pmZ[0]<<     // offset
	"pmZ1="<<pmZ[1]<<     // slope
	"pmZ2="<<pmZ[2]<<     // norm chi2
	"\n";
    }
  }
  
  delete pcstream;
  
}


void AliTPCcalibLaser::Analyze(){
  //
  //
  //
}


Long64_t AliTPCcalibLaser::Merge(TCollection *li) {

  TIterator* iter = li->MakeIterator();
  AliTPCcalibLaser* cal = 0;
  static Int_t counter0=0;
  while ((cal = (AliTPCcalibLaser*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibLaser::Class())) {
      //      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    printf("Marging number %d\n", counter0);
    counter0++;
    //
    MergeFitHistos(cal);
    TH1F *h=0x0;
    TH1F *hm=0x0;
    TH2F *h2=0x0;
    TH2F *h2m=0x0;
    //    TProfile *hp=0x0;
    //TProfile *hpm=0x0;

    for (Int_t id=0; id<336; id++){
      // merge fDeltaZ histograms
      hm = (TH1F*)cal->fDeltaZ.At(id);
      h  = (TH1F*)fDeltaZ.At(id);
      if (!h &&hm &&hm->GetEntries()>0) {
	h=(TH1F*)hm->Clone();
	h->SetDirectory(0);
	fDeltaZ.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      // merge fP3 histograms
      hm = (TH1F*)cal->fDeltaP3.At(id);
      h  = (TH1F*)fDeltaP3.At(id);
      if (!h&&hm &&hm->GetEntries()>0) {
	h=(TH1F*)hm->Clone();
	h->SetDirectory(0);
	fDeltaP3.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      // merge fP4 histograms
      hm = (TH1F*)cal->fDeltaP4.At(id);
      h  = (TH1F*)fDeltaP4.At(id);
      if (!h &&hm &&hm->GetEntries()>0) {
	h=(TH1F*)hm->Clone();
	h->SetDirectory(0);
	fDeltaP4.AddAt(h,id);
      }
      if (hm) h->Add(hm);

      //
      // merge fDeltaPhi histograms
      hm = (TH1F*)cal->fDeltaPhi.At(id);
      h  = (TH1F*)fDeltaPhi.At(id);
      if (!h &&hm &&hm->GetEntries()>0) {
	h= (TH1F*)hm->Clone();
	h->SetDirectory(0);
	fDeltaPhi.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      // merge fDeltaPhiP histograms
      hm = (TH1F*)cal->fDeltaPhiP.At(id);
      h  = (TH1F*)fDeltaPhiP.At(id);
      if (!h&&hm &&hm->GetEntries()>0) {
	h=(TH1F*)hm->Clone();
	h->SetDirectory(0);
	fDeltaPhiP.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      // merge fSignals histograms
      hm = (TH1F*)cal->fSignals.At(id);
      h  = (TH1F*)fSignals.At(id);
      if (!h&&hm &&hm->GetEntries()>0) {
	h=(TH1F*)hm->Clone();
	h->SetDirectory(0);
	fSignals.AddAt(h,id);
      }
      if (hm) h->Add(hm);
      //
      //
      // merge ProfileY histograms -0
      h2m = (TH2F*)cal->fDeltaYres.At(id);
      h2  = (TH2F*)fDeltaYres.At(id);
      if (h2m&&h2) h2->Add(h2m);
      //
      h2m = (TH2F*)cal->fDeltaZres.At(id);
      h2  = (TH2F*)fDeltaZres.At(id);
      if (h2m&&h2) h2->Add(h2m);
      // merge ProfileY histograms - 2
      h2m = (TH2F*)cal->fDeltaYres2.At(id);
      h2  = (TH2F*)fDeltaYres2.At(id);
      if (h2m&&h2) h2->Add(h2m);
      //
      h2m = (TH2F*)cal->fDeltaZres2.At(id);
      h2  = (TH2F*)fDeltaZres2.At(id);
      if (h2m&&h2) h2->Add(h2m);

      // merge ProfileY histograms - abs
      h2m = (TH2F*)cal->fDeltaYresAbs.At(id);
      h2  = (TH2F*)fDeltaYresAbs.At(id);
      if (h2m&&h2) h2->Add(h2m);
      if (h2m&&!h2) { h2=(TH2F*)h2m->Clone(); h2->SetDirectory(0); fDeltaYresAbs.AddAt(h2,id);}
      h2m = (TH2F*)cal->fDeltaZresAbs.At(id);
      h2  = (TH2F*)fDeltaZresAbs.At(id);
      if (h2m&&h2) h2->Add(h2m);
      if (h2m&&!h2) { h2=(TH2F*)h2m->Clone(); h2->SetDirectory(0); fDeltaZresAbs.AddAt(h2,id);}
      // merge ProfileY histograms - 3
      //h2m = (TH2F*)cal->fDeltaYres3.At(id);
      //h2  = (TH2F*)fDeltaYres3.At(id);
      //if (h2m) h2->Add(h2m);
      //
      //h2m = (TH2F*)cal->fDeltaZres3.At(id);
      //h2  = (TH2F*)fDeltaZres3.At(id);
      //if (h2m) h->Add(h2m);
      //
      //
    }
  }
  return 0;
}

void   AliTPCcalibLaser::MakeFitHistos(){
  //
  // Make a fit histograms
  // 
  // Number of clusters
  //
  //TH2F           *fHisNclIn;      //->Number of clusters inner
  //TH2F           *fHisNclOut;     //->Number of clusters outer
  //TH2F           *fHisNclIO;      //->Number of cluster inner outer
  // TH2F           *fHisdEdx;       //->dEdx histo
  fHisNclIn = new TH2F("HisNclIn","HisNclIn",336,0,336,64,10,64);
  fHisNclOut = new TH2F("HisNclOut","HisNclOut",336,0,336,100,10,100);
  fHisNclIO = new TH2F("HisNclIO","HisNclIO",336,0,336,160,10,160);
  //
  fHisLclIn  = new TH2F("HisLclIn","HisLclIn",336,0,336,64,10,64);
  fHisLclOut = new TH2F("HisLclOut","HisLclOut",336,0,336,100,10,150);
  fHisLclIO  = new TH2F("HisLclIO","HisLclIO",336,0,336,160,10,160);
  //
  fHisdEdx = new TH2F("HisdEdx","HisdEdx",336,0,336,160,1,50);
  fHisdZfit = new TH2F("HisdZfit","HisdZfit",336,0,336,300,-0.3,0.3);
  
  //
  //  Chi2
  //
  //   TH2F           *fHisChi2YIn1;      //->chi2 y inner - line
  //   TH2F           *fHisChi2YOut1;     //->chi2 y inner - line
  //   TH2F           *fHisChi2YIn2;      //->chi2 y inner - parabola
  //   TH2F           *fHisChi2YOut2;     //->chi2 y inner - parabola
  //   TH2F           *fHisChi2YIO1;      //->chi2 y IO    - common
  fHisChi2YIn1   = new TH2F("Chi2YIn1","Chi2YIn1",336,0,336,500,0.001,0.5);
  fHisChi2YOut1  = new TH2F("Chi2YOut1","Chi2YOut1",336,0,336,500,0.001,0.5);
  fHisChi2YIn2   = new TH2F("Chi2YIn2","Chi2YIn2",336,0,336,500,0.001,0.5);
  fHisChi2YOut2  = new TH2F("Chi2YOut2","Chi2YOut2",336,0,336,500,0.001,0.5);
  fHisChi2YIO1   = new TH2F("Chi2YIO1","Chi2YIO1",336,0,336,500,0.001,0.5);
  //   TH2F           *fHisChi2ZIn1;      //->chi2 z inner - line
  //   TH2F           *fHisChi2ZOut1;     //->chi2 z inner - line
  //   TH2F           *fHisChi2ZIn2;      //->chi2 z inner - parabola
  //   TH2F           *fHisChi2ZOut2;     //->chi2 z inner - parabola
  //   TH2F           *fHisChi2ZIO1;      //->chi2 z IO    - common
  fHisChi2ZIn1   = new TH2F("Chi2ZIn1","Chi2ZIn1",336,0,336,500,0.001,0.5);
  fHisChi2ZOut1  = new TH2F("Chi2ZOut1","Chi2ZOut1",336,0,336,500,0.001,0.5);
  fHisChi2ZIn2   = new TH2F("Chi2ZIn2","Chi2ZIn2",336,0,336,500,0.001,0.5);
  fHisChi2ZOut2  = new TH2F("Chi2ZOut2","Chi2ZOut2",336,0,336,500,0.001,0.5);
  fHisChi2ZIO1   = new TH2F("Chi2ZIO1","Chi2ZIO1",336,0,336,500,0.001,0.5);
  //
  // Fit
  //
  //
  //   TH2F           *fHisPy1vP0;     //-> delta y   P0outer-P0inner - line
  //   TH2F           *fHisPy2vP0;     //-> delta y   P0outer-P0inner - parabola
  //   TH2F           *fHisPy3vP0;     //-> delta y   P0outer-P0inner - common parabola
  //   TH2F           *fHisPy1vP1;     //-> delta ky  P1outer-P1inner - line
  //   TH2F           *fHisPy2vP1;     //-> delta ky  P1outer-P1inner - parabola
  //   TH2F           *fHisPy3vP1;     //-> delta ky  P1outer-P1inner - common parabola
  //   TH2F           *fHisPy2vP2In;   //-> Curv  P2inner - parabola
  //   TH2F           *fHisPy2vP2Out;  //-> Curv  P2outer - parabola
  //   TH2F           *fHisPy3vP2IO;   //-> Curv  P2outerinner - common parabola
  fHisPy1vP0   = new TH2F("HisPy1vP0",   "HisPy1vP0",336,0,336,500,-0.3,0.3);
  fHisPy2vP0   = new TH2F("HisPy2vP0",   "HisPy2vP0",336,0,336,500,-0.3,0.3);
  fHisPy3vP0   = new TH2F("HisPy3vP0",   "HisPy3vP0",336,0,336,500,-0.3,0.3);
  fHisPy1vP1   = new TH2F("HisPy1vP1",   "HisPy1vP1",336,0,336,500,-0.01,0.01);
  fHisPy2vP1   = new TH2F("HisPy2vP1",   "HisPy2vP1",336,0,336,500,-0.01,0.01);
  fHisPy3vP1   = new TH2F("HisPy3vP1",   "HisPy3vP1",336,0,336,500,-0.01,0.01);
  fHisPy2vP2In = new TH2F("HisPy2vP2In", "HisPy2vP2In",336,0,336,500,-0.0003,0.0003);
  fHisPy2vP2Out= new TH2F("HisPy2vP2Out","HisPy2vP2Out",336,0,336,500,-0.0003,0.0003);
  fHisPy3vP2IO = new TH2F("HisPy3vP2IO", "HisPy3vP2IO",336,0,336,500,-0.0003,0.0003);
  //
  //
  //   TH2F           *fHisPz1vP0;     //-> delta z   P0outer-P0inner - line
  //   TH2F           *fHisPz2vP0;     //-> delta z   P0outer-P0inner - parabola
  //   TH2F           *fHisPz3vP0;     //-> delta z   P0outer-P0inner - common parabola
  //   TH2F           *fHisPz1vP1;     //-> delta kz  P1outer-P1inner - line
  //   TH2F           *fHisPz2vP1;     //-> delta kz  P1outer-P1inner - parabola
  //   TH2F           *fHisPz3vP1;     //-> delta kz  P1outer-P1inner - common parabola
  //   TH2F           *fHisPz2vP2In;   //-> Curv  P2inner - parabola
  //   TH2F           *fHisPz2vP2Out;  //-> Curv  P2outer - parabola
  //   TH2F           *fHisPz3vP2IO;   //-> Curv  P2outerinner - common parabola
  fHisPz1vP0   = new TH2F("HisPz1vP0",   "HisPz1vP0",336,0,336,500,-0.3,0.3);
  fHisPz2vP0   = new TH2F("HisPz2vP0",   "HisPz2vP0",336,0,336,500,-0.3,0.3);
  fHisPz3vP0   = new TH2F("HisPz3vP0",   "HisPz3vP0",336,0,336,500,-0.3,0.3);
  fHisPz1vP1   = new TH2F("HisPz1vP1",   "HisPz1vP1",336,0,336,500,-0.01,0.01);
  fHisPz2vP1   = new TH2F("HisPz2vP1",   "HisPz2vP1",336,0,336,500,-0.01,0.01);
  fHisPz3vP1   = new TH2F("HisPz3vP1",   "HisPz3vP1",336,0,336,500,-0.01,0.01);
  fHisPz2vP2In = new TH2F("HisPz2vP2In", "HisPz2vP2In",336,0,336,500,-0.0003,0.0003);
  fHisPz2vP2Out= new TH2F("HisPz2vP2Out","HisPz2vP2Out",336,0,336,500,-0.0002,0.0002);
  fHisPz3vP2IO = new TH2F("HisPz3vP2IO", "HisPz3vP2IO",336,0,336,500,-0.0002,0.0002);
  
  fHisYAbsErrors = new TH1F("HisYabsErrors", "Errors per beam (y)", 336, 0, 336);
  fHisZAbsErrors = new TH1F("HisZabsErrors", "Errors per beam (z)", 336, 0, 336);

  fHisNclIn->SetDirectory(0);      //->Number of clusters inner
  fHisNclOut->SetDirectory(0);     //->Number of clusters outer
  fHisNclIO->SetDirectory(0);      //->Number of cluster inner outer
  fHisLclIn->SetDirectory(0);      //->Level arm inner
  fHisLclOut->SetDirectory(0);     //->Level arm outer
  fHisLclIO->SetDirectory(0);      //->Number of cluster inner outer
  fHisdEdx->SetDirectory(0);      //->Number of cluster inner outer
  fHisdZfit->SetDirectory(0);      //->Number of cluster inner outer
  //
  //
  fHisChi2YIn1->SetDirectory(0);      //->chi2 y inner - line
  fHisChi2YOut1->SetDirectory(0);     //->chi2 y inner - line
  fHisChi2YIn2->SetDirectory(0);      //->chi2 y inner - parabola
  fHisChi2YOut2->SetDirectory(0);     //->chi2 y inner - parabola
  fHisChi2YIO1->SetDirectory(0);      //->chi2 y IO    - common
  fHisChi2ZIn1->SetDirectory(0);      //->chi2 z inner - line
  fHisChi2ZOut1->SetDirectory(0);     //->chi2 z inner - line
  fHisChi2ZIn2->SetDirectory(0);      //->chi2 z inner - parabola
  fHisChi2ZOut2->SetDirectory(0);     //->chi2 z inner - parabola
  fHisChi2ZIO1->SetDirectory(0);      //->chi2 z IO    - common
  //
  //
  fHisPy1vP0->SetDirectory(0);     //-> delta y   P0outer-P0inner - line
  fHisPy2vP0->SetDirectory(0);     //-> delta y   P0outer-P0inner - parabola
  fHisPy3vP0->SetDirectory(0);     //-> delta y   P0outer-P0inner - common parabola
  fHisPy1vP1->SetDirectory(0);     //-> delta ky  P1outer-P1inner - line
  fHisPy2vP1->SetDirectory(0);     //-> delta ky  P1outer-P1inner - parabola
  fHisPy3vP1->SetDirectory(0);     //-> delta ky  P1outer-P1inner - common parabola
  fHisPy2vP2In->SetDirectory(0);   //-> Curv  P2inner - parabola
  fHisPy2vP2Out->SetDirectory(0);  //-> Curv  P2outer - parabola
  fHisPy3vP2IO->SetDirectory(0);   //-> Curv  P2outerinner - common parabola
  //
  //
  fHisPz1vP0->SetDirectory(0);     //-> delta z   P0outer-P0inner - line
  fHisPz2vP0->SetDirectory(0);     //-> delta z   P0outer-P0inner - parabola
  fHisPz3vP0->SetDirectory(0);     //-> delta z   P0outer-P0inner - common parabola
  fHisPz1vP1->SetDirectory(0);     //-> delta kz  P1outer-P1inner - line
  fHisPz2vP1->SetDirectory(0);     //-> delta kz  P1outer-P1inner - parabola
  fHisPz3vP1->SetDirectory(0);     //-> delta kz  P1outer-P1inner - common parabola
  fHisPz2vP2In->SetDirectory(0);   //-> Curv  P2inner - parabola
  fHisPz2vP2Out->SetDirectory(0);  //-> Curv  P2outer - parabola
  fHisPz3vP2IO->SetDirectory(0);   //-> Curv  P2outerinner - common parabola

  fHisYAbsErrors->SetDirectory(0); //-> total errors per beam in the abs res y analysis
  fHisZAbsErrors->SetDirectory(0); //-> total errors per beam in the abs res z analysis

  

  //
  //
  //
  for (Int_t id=0; id<336;id++){
    TH2F *profy = (TH2F*)fDeltaYres.UncheckedAt(id);
    TH2F *profz = (TH2F*)fDeltaZres.UncheckedAt(id);
    //TH2F *profy2 = (TH2F*)fDeltaYres2.UncheckedAt(id);
    TH2F *profy2 = 0;
    TH2F *profz2 = (TH2F*)fDeltaZres2.UncheckedAt(id);
    TH2F *profyabs = (TH2F*)fDeltaYresAbs.UncheckedAt(id);
    TH2F *profzabs = (TH2F*)fDeltaYresAbs.UncheckedAt(id);
    //    TH2F *profy3 = (TH2F*)fDeltaYres3.UncheckedAt(id);
    //TH2F *profz3 = (TH2F*)fDeltaZres3.UncheckedAt(id);
    if (!profy){
      profy=new TH2F(Form("pry%03d",id),Form("Y Residuals for Laser Beam %03d -Linear",id),160,0,160,50,-0.5,0.5);
      profy->SetDirectory(0);
      fDeltaYres.AddAt(profy,id);
      profy2=new TH2F(Form("pry%03d",id),Form("Y Residuals for Laser Beam %03d -Parabolic",id),160,0,160,50,-0.5,0.5);
      profy2->SetDirectory(0);
      fDeltaYres2.AddAt(profy2,id);
      if(!fUseFixedDriftV)
	profyabs=new TH2F(Form("pryabs%03d",id),Form("Y Residuals for Laser Beam %03d -Absolute",id),160,0,160,100,-1.0,1.0); // has to be bigger based on earlier studies
      else
	profyabs=new TH2F(Form("pryabs%03d",id),Form("Y Residuals for Laser Beam %03d -Absolute",id),160,0,160,200,-2.0,2.0); // has to be bigger based on earlier studies
      profyabs->SetDirectory(0);
      fDeltaYresAbs.AddAt(profyabs,id);
      //profy3=new TH2F(Form("pry%03d",id),Form("Y Residuals for Laser Beam %03d- Parabolic2",id),160,0,160,100,-0.5,0.5);
      //profy3->SetDirectory(0);
      //fDeltaYres3.AddAt(profy3,id);
    }
    if (!profz){
      profz=new TH2F(Form("prz%03d",id),Form("Z Residuals for Laser Beam %03d Linear",id),160,0,160,50,-0.5,0.5);
      profz->SetDirectory(0);
      fDeltaZres.AddAt(profz,id);
      profz2=new TH2F(Form("prz%03d",id),Form("Z Residuals for Laser Beam %03d - Parabolic",id),160,0,160,50,-0.5,0.5);
      profz2->SetDirectory(0);
      fDeltaZres2.AddAt(profz2,id);
      if(!fUseFixedDriftV)
	profzabs=new TH2F(Form("przabs%03d",id),Form("Z Residuals for Laser Beam %03d -Absolute",id),160,0,160,100,-1.0,1.0); // has to be bigger based on earlier studies
      else
	profzabs=new TH2F(Form("przabs%03d",id),Form("Z Residuals for Laser Beam %03d -Absolute",id),160,0,160,200,-2.0,2.0); // has to be bigger based on earlier studies
      profzabs->SetDirectory(0);
      fDeltaZresAbs.AddAt(profzabs,id);
      //profz3=new TH2F(Form("prz%03d",id),Form("Z Residuals for Laser Beam %03d- Parabolic2",id),160,0,160,100,-0.5,0.5);
      //profz3->SetDirectory(0);
      //fDeltaZres3.AddAt(profz3,id);
    }
  }
  //
  //
  for (Int_t id=0; id<336;id++){
    TH1F * hisdz = (TH1F*)fDeltaZ.At(id);
    //TH1F * hisP3 = (TH1F*)fDeltaP3.At(id);
    TH1F * hisP3 = 0;
    TH1F * hisP4 = (TH1F*)fDeltaP4.At(id);
    
    TH1F * hisdphi = (TH1F*)fDeltaPhi.At(id);
    TH1F * hisdphiP = (TH1F*)fDeltaPhiP.At(id);
    TH1F * hisSignal = (TH1F*)fSignals.At(id);

    if (!hisdz){
      hisdz = new TH1F(Form("hisdz%d",id),Form("hisdz%d",id),1000,-10,10);
      hisdz->SetDirectory(0);
      fDeltaZ.AddAt(hisdz,id);

      hisP3 = new TH1F(Form("hisPar3v%d",id),Form("hisPar3v%d",id),400,-0.06,0.06);
      hisP3->SetDirectory(0);
      fDeltaP3.AddAt(hisP3,id);
      //
      hisP4 = new TH1F(Form("hisPar4v%d",id),Form("hisPar4v%d",id),200,-0.06,0.06);
      hisP4->SetDirectory(0);
      fDeltaP4.AddAt(hisP4,id);

      //
      hisdphi = new TH1F(Form("hisdphi%d",id),Form("hisdphi%d",id),1000,-1,1);
      hisdphi->SetDirectory(0);
      fDeltaPhi.AddAt(hisdphi,id);
      //
      hisdphiP = new TH1F(Form("hisdphiP%d",id),Form("hisdphiP%d",id),1000,-0.01,0.01);
      hisdphiP->SetDirectory(0);
      fDeltaPhiP.AddAt(hisdphiP,id);
      hisSignal = new TH1F(Form("hisSignal%d",id),Form("hisSignal%d",id),100,0,300);
      hisSignal->SetDirectory(0);
      fSignals.AddAt(hisSignal,id);
    }
  }

  SetBeamParameters(fBeamOffsetZOuter, fBeamSlopeZOuter, fBeamSectorOuter,2);
  SetBeamParameters(fBeamOffsetZInner, fBeamSlopeZInner, fBeamSectorInner,3);
  SetBeamParameters(fBeamOffsetYOuter, fBeamSlopeYOuter, fBeamSectorOuter,0);
  SetBeamParameters(fBeamOffsetYInner, fBeamSlopeYInner, fBeamSectorInner,1);
  //
  // Make THnSparse
  //
  //                   id   side   rod  bundle  beam  dP0  dP1  dP2  dP3  dP4 ncl dEdx 
  Int_t binsLaser[12]= {336,  //id
			2,    //side
			6,    //rod
			4,    //bundle
			7,    //beam
			300,  //dP0
			300,  //dP1
			300,  //dP2
			300,  //dP3
			300,  //dP4
			80,   //ncl
			50};   //dEdx
  Double_t xminLaser[12]= {0,  //id
			0,     //side
			0,     //rod
			0,     //bundle
			0,     //beam
			-1,    //dP0
			-1,    //dP1
			-0.01,   //dP2
			-0.01,  //dP3
			-0.1,  //dP4
			0,   //ncl
			0};   //sqrt dEdx
  Double_t xmaxLaser[12]= {336,  //id
			2,    //side
			6,    //rod
			4,    //bundle
			7,    //beam
			1,  //dP0
			1,  //dP1
			0.01,  //dP2
			0.01,  //dP3
			0.1,  //dP4
			160,   //ncl
			40};   //sqrt dEdx
  
  TString nameLaser[12]= {"id",
			  "side",
			  "rod",
			  "bundle",
			  "beam",
			  "dP0",
			  "dP1",
			  "dP2",
			  "dP3",
			  "dP4",
			  "ncl",
			  "sqrt dEdx"};
  TString titleLaser[12]= {"id",
			  "side",
			  "rod",
			  "bundle",
			  "beam",
			  "#Delta_{P0}",
			  "#Delta_{P1}",
			  "#Delta_{P2}",
			  "#Delta_{P3}",
			  "#Delta_{P4}",
			  "N_{cl}",
			  "#sqrt{dEdx}"};
  fHisLaser = new THnSparseS("dLaser","#Delta_{Laser}", 12, binsLaser,xminLaser, xmaxLaser);
  for (Int_t iaxis=1; iaxis<12; iaxis++){
    fHisLaser->GetAxis(iaxis)->SetName(nameLaser[iaxis]);
    fHisLaser->GetAxis(iaxis)->SetTitle(titleLaser[iaxis]);
  }
  //
  //                  Delta       Time bin
  //                              Pad        SigmaShape      Q charge  pad row  trackID
  Int_t   binsRow[6]={200,        10000,           20,            30,     159,  336};
  Double_t axisMin[6]={-1,             0,           0,            1,     0  ,    0};
  Double_t axisMax[6]={ 1,          1000,           1,           30,     159,  336};
  TString axisName[6]={"Delta","bin", "rms shape", "sqrt(Q)", "row","trackID"};

  binsRow[1]=2000;
  axisMin[1]=0;
  axisMax[1]=200;
  fHisLaserPad = new THnSparseS("laserPad","#Delta_{Laser}", 6, binsRow,axisMin, axisMax);  
  //
  binsRow[0]=1000;
  axisMin[0]=-20;
  axisMax[0]=20;
  binsRow[1]=10000;
  axisMin[1]=0;
  axisMax[1]=1000;
  //
  fHisLaserTime= new THnSparseS("laserTime","#Delta_{Laser}", 6, binsRow,axisMin, axisMax);
  //
  for (Int_t iaxis=0; iaxis<6; iaxis++){
    fHisLaserPad->GetAxis(iaxis)->SetName(axisName[iaxis]);
    fHisLaserTime->GetAxis(iaxis)->SetTitle(axisName[iaxis]);
  }
}

void AliTPCcalibLaser::MergeFitHistos(AliTPCcalibLaser * laser){
  //
  // Merge content of histograms 
  //
  // Only first histogram is checked - all other should be the same
  if (fHisLaser &&laser->fHisLaser) fHisLaser->Add(laser->fHisLaser);
  if (fHisLaserPad &&laser->fHisLaserPad) fHisLaserPad->Add(laser->fHisLaserPad);
  if (!fHisLaserPad &&laser->fHisLaserPad) fHisLaserPad=(THnSparseS*)laser->fHisLaserPad->Clone();
  if (fHisLaserTime &&laser->fHisLaserTime) fHisLaserTime->Add(laser->fHisLaserTime);
  if (!fHisLaserTime &&laser->fHisLaserTime) fHisLaserTime=(THnSparseS*)laser->fHisLaserTime->Clone();
  
  if (!laser->fHisNclIn) laser->MakeFitHistos();  // empty histograms
  if (!fHisNclIn) MakeFitHistos();
  if (fHisNclIn->GetEntries()<1) MakeFitHistos();
  //
  fHisNclIn->Add(laser->fHisNclIn  );      //->Number of clusters inner
  fHisNclOut->Add(laser->fHisNclOut  );     //->Number of clusters outer
  fHisNclIO->Add(laser->fHisNclIO  );      //->Number of cluster inner outer
  fHisLclIn->Add(laser->fHisLclIn  );      //->Level arm inner
  fHisLclOut->Add(laser->fHisLclOut  );     //->Level arm outer
  fHisLclIO->Add(laser->fHisLclIO  );      //->Number of cluster inner outer
  fHisdEdx->Add(laser->fHisdEdx  );      //->dedx
  fHisdZfit->Add(laser->fHisdZfit  );      //->dedx
  //
  //
  fHisChi2YIn1->Add(laser->fHisChi2YIn1  );      //->chi2 y inner - line
  fHisChi2YOut1->Add(laser->fHisChi2YOut1  );     //->chi2 y inner - line
  fHisChi2YIn2->Add(laser->fHisChi2YIn2  );      //->chi2 y inner - parabola
  fHisChi2YOut2->Add(laser->fHisChi2YOut2  );     //->chi2 y inner - parabola
  fHisChi2YIO1->Add(laser->fHisChi2YIO1  );      //->chi2 y IO    - common
  fHisChi2ZIn1->Add(laser->fHisChi2ZIn1  );      //->chi2 z inner - line
  fHisChi2ZOut1->Add(laser->fHisChi2ZOut1  );     //->chi2 z inner - line
  fHisChi2ZIn2->Add(laser->fHisChi2ZIn2  );      //->chi2 z inner - parabola
  fHisChi2ZOut2->Add(laser->fHisChi2ZOut2  );     //->chi2 z inner - parabola
  fHisChi2ZIO1->Add(laser->fHisChi2ZIO1  );      //->chi2 z IO    - common
  //
  //
  fHisPy1vP0->Add(laser->fHisPy1vP0  );     //-> delta y   P0outer-P0inner - line
  fHisPy2vP0->Add(laser->fHisPy2vP0  );     //-> delta y   P0outer-P0inner - parabola
  fHisPy3vP0->Add(laser->fHisPy3vP0  );     //-> delta y   P0outer-P0inner - common parabola
  fHisPy1vP1->Add(laser->fHisPy1vP1  );     //-> delta ky  P1outer-P1inner - line
  fHisPy2vP1->Add(laser->fHisPy2vP1  );     //-> delta ky  P1outer-P1inner - parabola
  fHisPy3vP1->Add(laser->fHisPy3vP1  );     //-> delta ky  P1outer-P1inner - common parabola
  fHisPy2vP2In->Add(laser-> fHisPy2vP2In  );   //-> Curv  P2inner - parabola
  fHisPy2vP2Out->Add(laser->fHisPy2vP2Out  );  //-> Curv  P2outer - parabola
  fHisPy3vP2IO->Add(laser->fHisPy3vP2IO  );   //-> Curv  P2outerinner - common parabola
  //
  //
  fHisPz1vP0->Add(laser->fHisPz1vP0  );     //-> delta z   P0outer-P0inner - line
  fHisPz2vP0->Add(laser->fHisPz2vP0  );     //-> delta z   P0outer-P0inner - parabola
  fHisPz3vP0->Add(laser->fHisPz3vP0  );     //-> delta z   P0outer-P0inner - common parabola
  fHisPz1vP1->Add(laser->fHisPz1vP1  );     //-> delta kz  P1outer-P1inner - line
  fHisPz2vP1->Add(laser->fHisPz2vP1  );     //-> delta kz  P1outer-P1inner - parabola
  fHisPz3vP1->Add(laser->fHisPz3vP1  );     //-> delta kz  P1outer-P1inner - common parabola
  fHisPz2vP2In->Add(laser->fHisPz2vP2In  );   //-> Curv  P2inner - parabola
  fHisPz2vP2Out->Add(laser->fHisPz2vP2Out  );  //-> Curv  P2outer - parabola
  fHisPz3vP2IO->Add(laser->fHisPz3vP2IO  );   //-> Curv  P2outerinner - common parabola
  fHisYAbsErrors->Add(laser->fHisYAbsErrors); //-> total errors per beam in the abs res y analysis
  fHisZAbsErrors->Add(laser->fHisZAbsErrors); //-> total errors per beam in the abs res z analysis
  //
  //
  //
  



}




void AliTPCcalibLaser::DumpFitInfo(TTree * chainFit,Int_t id){
  //
  // Dump fit information - collect information from the streamers 
  //
  /*
    TChain * chainFit=0;
    TChain * chainTrack=0;
    TChain * chain=0;
    //
    gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
    gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
    AliXRDPROOFtoolkit tool;
    chainTrack = tool.MakeChain("laser.txt","Track",0,10200);
    chainTrack->Lookup();
    chainTrack->SetProof(kTRUE);
    chainDrift = tool.MakeChain("laser.txt","driftv",0,10200);
    chainDrift->Lookup();
    chainDrift->SetProof(kTRUE);
 
    chain = tool.MakeChain("laser.txt","Residuals",0,10200);
    chain->Lookup();
    chainFit = tool.MakeChain("laser.txt","FitModels",0,10200);
    chainFit->Lookup();
    chainFit->SetProof(kTRUE);
    chain->SetProof(kTRUE);
    AliTPCLaserTrack::LoadTracks();  
    //AliTPCcalibLaser::DumpFitInfo(chainFit,0);

  */
  //
  // Fit cuts
  //
  TCut cutP3("abs(Tr.fP[3])<0.1");
  TCut cutP4("abs(Tr.fP[4])<0.5");
  TCut cutPx = cutP3+cutP4;
  TCut cutChi2YOut("sqrt(chi2y2Out*dEdx)<5&&chi2y2Out>0");
  TCut cutChi2ZOut("sqrt(chi2z2Out*dEdx)<5&&chi2z2Out>0");
  TCut cutChi2YIn("sqrt(chi2y2In*dEdx)<5&&chi2y2In>0");
  TCut cutChi2ZIn("sqrt(chi2z2In*dEdx)<5&&chi2z2In>0");
  //
  TCut cutdEdx("sqrt(dEdx)>3");
  TCut cutDY("abs(yPol2In.fElements[2]*nclI*nclI/4.)<3");
  TCut cutN("nclO>20&&nclI>20");
  TCut cutA = cutChi2YOut+cutChi2ZOut+cutChi2YIn+cutChi2ZIn+cutN+cutdEdx+cutPx+"accept";
  //
  // Cluster cuts
  //
  TCut cutClY("abs(Cl[].fY-TrYpol2.fElements)<0.15");
  TCut cutClZ("abs(Cl[].fZ-TrZpol2.fElements)<0.15");
  TCut cutClX("abs(Cl[].fX)>10");
  TCut cutE("abs(Cl[].fY/Cl[].fX)<0.14");
  TCut cutSY("sqrt(Cl[].fSigmaY2)>0.05");
  TCut cutSZ("sqrt(Cl[].fSigmaZ2)>0.05");
  TCut cutQ("sqrt(Cl[].fMax)>4");
  TCut cutCl=cutClY+cutClZ+cutClX+cutE+cutSY+cutSZ+cutQ;


  TH1F * phisAl     = 0;
  TH1F * phisAccept = 0;
  TH1F * phisOut    = 0;
  TProfile * pdEdx  = 0;

  TProfile * pP0    = 0;
  TProfile * pP1    = 0;
  TProfile * pP2    = 0;
  TProfile * pP3    = 0;
  TProfile * pP4    = 0;
  //
  TProfile * pNclI  = 0;
  TProfile * pNclO  = 0;
  //
  TProfile * pchi2YIn   =0;
  TProfile * pchi2ZIn   =0;
  TProfile * pchi2YOut  =0;
  TProfile * pchi2ZOut  =0;
  TProfile * pchi2YInOut =0;
  TProfile * pchi2ZInOut =0;;
  // laser counters
  chainFit->Draw("LTr.fId>>hisAl(350,0,350)","LTr.fId<350");
  phisAl = (TH1F*)gROOT->FindObject("hisAl");
  chainFit->Draw("LTr.fId>>hisAccept(350,0,350)","LTr.fId<350"+cutA);
  phisAccept = (TH1F*)gROOT->FindObject("hisAccept");
  chainFit->Draw("LTr.fId>>hisOut(350,0,350)","LTr.fId<350"+!cutA);
  phisOut = (TH1F*)gROOT->FindObject("hisOut");
  //
  chainFit->Draw("sqrt(dEdx):LTr.fId>>hdEdx(350,0,350)","","prof");
  pdEdx   = (TProfile*)gROOT->FindObject("hdEdx");
  // track param
  //
  chainFit->Draw("Tr.fP[0]:LTr.fId>>hP0(350,0,350)","Tr.fP[4]/sqrt(Tr.fC[14])<20"+cutA,"prof");
  pP0   = (TProfile*)gROOT->FindObject("hP0");
  chainFit->Draw("Tr.fP[1]:LTr.fId>>hP1(350,0,350)","Tr.fP[4]/sqrt(Tr.fC[14])<20"+cutA,"prof");
  pP1   = (TProfile*)gROOT->FindObject("hP1");
  chainFit->Draw("Tr.fP[2]:LTr.fId>>hP2(350,0,350)","Tr.fP[4]/sqrt(Tr.fC[14])<20"+cutA,"prof");
  pP2   = (TProfile*)gROOT->FindObject("hP2");
  chainFit->Draw("Tr.fP[3]:LTr.fId>>hP3(350,0,350)","Tr.fP[4]/sqrt(Tr.fC[14])<20"+cutA,"prof");
  pP3   = (TProfile*)gROOT->FindObject("hP3");
  chainFit->Draw("Tr.fP[4]:LTr.fId>>hP4(350,0,350)","Tr.fP[4]/sqrt(Tr.fC[14])<20"+cutA,"prof");
  pP4   = (TProfile*)gROOT->FindObject("hP4");
  //
  chainFit->Draw("nclI:LTr.fId>>hNclI(350,0,350)","","prof");
  pNclI   = (TProfile*)gROOT->FindObject("hNclI");
  chainFit->Draw("nclO:LTr.fId>>hNclO(350,0,350)","","prof");
  pNclO   = (TProfile*)gROOT->FindObject("hNclO");
  //
  //
  chainFit->Draw("sqrt(chi2y2In):LTr.fId>>hChi2YIn(350,0,350)",cutA+"","prof");
  pchi2YIn   = (TProfile*)gROOT->FindObject("hChi2YIn");
  chainFit->Draw("sqrt(chi2y2Out):LTr.fId>>hChi2YOut(350,0,350)",cutA+"","prof");
  pchi2YOut   = (TProfile*)gROOT->FindObject("hChi2YOut");
  chainFit->Draw("sqrt(chi2yInOut):LTr.fId>>hChi2YInOut(350,0,350)",cutA+"","prof");
  pchi2YInOut   = (TProfile*)gROOT->FindObject("hChi2YInOut");
  chainFit->Draw("sqrt(chi2z2In):LTr.fId>>hChi2ZIn(350,0,350)",cutA+"","prof");
  pchi2ZIn   = (TProfile*)gROOT->FindObject("hChi2ZIn");
  chainFit->Draw("sqrt(chi2z2Out):LTr.fId>>hChi2ZOut(350,0,350)",cutA+"","prof");
  pchi2ZOut   = (TProfile*)gROOT->FindObject("hChi2ZOut");
  chainFit->Draw("sqrt(chi2zInOut):LTr.fId>>hChi2ZInOut(350,0,350)",cutA+"","prof");
  pchi2ZInOut   = (TProfile*)gROOT->FindObject("hChi2ZInOut");
  //
  // second derivatives
  //
  TH2F * phisPy2In = new TH2F("Py2Inner","Py2Inner",350,0,350,100,-0.001,0.001);
  chainFit->Draw("yPol2In.fElements[2]:LTr.fId>>Py2Inner",cutA,"");
  TH2F * phisPy2Out = new TH2F("Py2Outer","Py2Outer",350,0,350,200,-0.0005,0.0005);
  chainFit->Draw("yPol2Out.fElements[2]:LTr.fId>>Py2Outer",cutA,"");
  TH2F * phisPy2InOut = new TH2F("Py2InOut","Py2InOut",350,0,350,200,-0.0005,0.0005);
  chainFit->Draw("yInOut.fElements[4]:LTr.fId>>Py2InOut",cutA,"");
  //
  phisPy2In->FitSlicesY();
  TH1D * phisPy2InEntries = (TH1D*)gROOT->FindObject("Py2Inner_0");
  TH1D * phisPy2InMean = (TH1D*)gROOT->FindObject("Py2Inner_1");
  TH1D * phisPy2InSigma = (TH1D*)gROOT->FindObject("Py2Inner_2");
  //
  phisPy2Out->FitSlicesY();
  TH1D * phisPy2OutEntries = (TH1D*)gROOT->FindObject("Py2Outer_0");
  TH1D * phisPy2OutMean = (TH1D*)gROOT->FindObject("Py2Outer_1");
  TH1D * phisPy2OutSigma = (TH1D*)gROOT->FindObject("Py2Outer_2");
  //
  phisPy2InOut->FitSlicesY();
  TH1D * phisPy2InOutEntries = (TH1D*)gROOT->FindObject("Py2InOut_0");
  TH1D * phisPy2InOutMean = (TH1D*)gROOT->FindObject("Py2InOut_1");
  TH1D * phisPy2InOutSigma = (TH1D*)gROOT->FindObject("Py2InOut_2");
  //
  TH2F * phisPz2In = new TH2F("Pz2Inner","Pz2Inner",350,0,350,100,-0.001,0.001);
  chainFit->Draw("zPol2In.fElements[2]:LTr.fId>>Pz2Inner",cutA,"");
  TH2F * phisPz2Out = new TH2F("Pz2Outer","Pz2Outer",350,0,350,200,-0.0005,0.0005);
  chainFit->Draw("zPol2Out.fElements[2]:LTr.fId>>Pz2Outer",cutA,"");
  TH2F * phisPz2InOut = new TH2F("Pz2InOut","Pz2InOut",350,0,350,200,-0.0005,0.0005);
  chainFit->Draw("zInOut.fElements[4]:LTr.fId>>Pz2InOut",cutA,"");
  //
  phisPz2In->FitSlicesY();
  TH1D * phisPz2InEntries = (TH1D*)gROOT->FindObject("Pz2Inner_0");
  TH1D * phisPz2InMean = (TH1D*)gROOT->FindObject("Pz2Inner_1");
  TH1D * phisPz2InSigma = (TH1D*)gROOT->FindObject("Pz2Inner_2");
  //
  phisPz2Out->FitSlicesY();
  TH1D * phisPz2OutEntries = (TH1D*)gROOT->FindObject("Pz2Outer_0");
  TH1D * phisPz2OutMean = (TH1D*)gROOT->FindObject("Pz2Outer_1");
  TH1D * phisPz2OutSigma = (TH1D*)gROOT->FindObject("Pz2Outer_2");
  //
  phisPz2InOut->FitSlicesY();
  TH1D * phisPz2InOutEntries = (TH1D*)gROOT->FindObject("Pz2InOut_0");
  TH1D * phisPz2InOutMean = (TH1D*)gROOT->FindObject("Pz2InOut_1");
  TH1D * phisPz2InOutSigma = (TH1D*)gROOT->FindObject("Pz2InOut_2");
  //
  //
  //


  {
    TTreeSRedirector *pcstream = new TTreeSRedirector("vscan.root");
    for (Int_t ilaser=0; ilaser<336; ilaser++){
      Float_t all=phisAl->GetBinContent(ilaser+1);
      Float_t accept=phisAccept->GetBinContent(ilaser+1);
      Float_t out=phisOut->GetBinContent(ilaser+1);
      Float_t sdedx = pdEdx->GetBinContent(ilaser+1);
      Float_t mP0 = pP0->GetBinContent(ilaser+1);
      Float_t mP1 = pP1->GetBinContent(ilaser+1);
      Float_t mP2 = pP2->GetBinContent(ilaser+1);
      Float_t mP3 = pP3->GetBinContent(ilaser+1);
      Float_t mP4 = pP4->GetBinContent(ilaser+1);


      Float_t nclI  = pNclI->GetBinContent(ilaser+1); 
      Float_t nclO  = pNclO->GetBinContent(ilaser+1); 
      //
      Float_t chi2YIn=pchi2YIn->GetBinContent(ilaser+1); 
      Float_t chi2YOut=pchi2YOut->GetBinContent(ilaser+1); 
      Float_t chi2YInOut=pchi2YInOut->GetBinContent(ilaser+1); 
      Float_t chi2ZIn=pchi2ZIn->GetBinContent(ilaser+1); 
      Float_t chi2ZOut=pchi2ZOut->GetBinContent(ilaser+1); 
      Float_t chi2ZInOut=pchi2ZInOut->GetBinContent(ilaser+1); 
      //
      Float_t entriesPy2In  = phisPy2InEntries->GetBinContent(ilaser+1);
      Float_t meanPy2In     = phisPy2InMean->GetBinContent(ilaser+1);
      Float_t sigmaPy2In    = phisPy2InSigma->GetBinContent(ilaser+1);
      //
      Float_t entriesPy2Out  = phisPy2OutEntries->GetBinContent(ilaser+1);
      Float_t meanPy2Out     = phisPy2OutMean->GetBinContent(ilaser+1);
      Float_t sigmaPy2Out    = phisPy2OutSigma->GetBinContent(ilaser+1);
      //
      Float_t entriesPy2InOut  = phisPy2InOutEntries->GetBinContent(ilaser+1);
      Float_t meanPy2InOut     = phisPy2InOutMean->GetBinContent(ilaser+1);
      Float_t sigmaPy2InOut    = phisPy2InOutSigma->GetBinContent(ilaser+1);
      //
      Float_t entriesPz2In  = phisPz2InEntries->GetBinContent(ilaser+1);
      Float_t meanPz2In     = phisPz2InMean->GetBinContent(ilaser+1);
      Float_t sigmaPz2In    = phisPz2InSigma->GetBinContent(ilaser+1);
      //
      Float_t entriesPz2Out  = phisPz2OutEntries->GetBinContent(ilaser+1);
      Float_t meanPz2Out     = phisPz2OutMean->GetBinContent(ilaser+1);
      Float_t sigmaPz2Out    = phisPz2OutSigma->GetBinContent(ilaser+1);
      //
      Float_t entriesPz2InOut  = phisPz2InOutEntries->GetBinContent(ilaser+1);
      Float_t meanPz2InOut     = phisPz2InOutMean->GetBinContent(ilaser+1);
      Float_t sigmaPz2InOut    = phisPz2InOutSigma->GetBinContent(ilaser+1);
      
      AliTPCLaserTrack* ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(ilaser);
      (*pcstream)<<"Scan"<<
	"Run="<<id<<
	"LTr.="<<ltrp<<
	"all="<<all<<
	"accept="<<accept<<
	"out="<<out<<
	"sdedx="<<sdedx<<
	"mP0="<<mP0<<
	"mP1="<<mP1<<
	"mP2="<<mP2<<
	"mP3="<<mP3<<
	"mP4="<<mP4<<
	"nclI="<<nclI<<
	"nclO="<<nclO<<
	"chi2YIn="<<chi2YIn<<
	"chi2YOut="<<chi2YOut<<
	"chi2YInOut="<<chi2YInOut<<
	"chi2ZIn="<<chi2ZIn<<
	"chi2ZOut="<<chi2ZOut<<
	"chi2ZInOut="<<chi2ZInOut<<
	//
	"nPy2In="<<entriesPy2In<<
	"mPy2In="<<meanPy2In<<
	"sPy2In="<<sigmaPy2In<<
	//
	"nPy2Out="<<entriesPy2Out<<
	"mPy2Out="<<meanPy2Out<<
	"sPy2Out="<<sigmaPy2Out<<
	//
	"nPy2InOut="<<entriesPy2InOut<<
	"mPy2InOut="<<meanPy2InOut<<
	"sPy2InOut="<<sigmaPy2InOut<<
	//
	"nPz2In="<<entriesPz2In<<
	"mPz2In="<<meanPz2In<<
	"sPz2In="<<sigmaPz2In<<
	//
	"nPz2Out="<<entriesPz2Out<<
	"mPz2Out="<<meanPz2Out<<
	"sPz2Out="<<sigmaPz2Out<<
	//
	"nPz2InOut="<<entriesPz2InOut<<
	"mPz2InOut="<<meanPz2InOut<<
	"sPz2InOut="<<sigmaPz2InOut<<
	"\n";
    }
    
    delete pcstream;
  }
}

void AliTPCcalibLaser::SetBeamParameters(TVectorD& meanOffset,
					 TVectorD& meanSlope,
					 TVectorD& sectorArray,
					 Int_t option)
{
  // This method should ideally go in AliTPCLaser
  // option == 0 (pads outer - closest to beam)
  // option == 1 (pads inner)
  // option == 2 (time outer)
  // option == 3 (time inner)
  Int_t nFailures = 0;

  for(Int_t id = 0; id < 336; id++) {
    
    if (!AliTPCLaserTrack::GetTracks()) 
      AliTPCLaserTrack::LoadTracks();
    AliTPCLaserTrack *ltrp =
      (AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    
    AliExternalTrackParam trackParam(*ltrp);
    
    Double_t deltaangle = 10.0; // sector is 1/2 a sector away from mirror
    if((option==1 || option==3)&& (ltrp->GetBeam()<=1 || ltrp->GetBeam()>=5))
      deltaangle = 30.0; // inner sector is 1 sector further away from mirror
    
    Double_t angle = trackParam.GetAlpha();
    if(angle<0)
      angle += 2*TMath::Pi();
    if(trackParam.GetSnp()>0) // track points to sector "before"
      angle -= deltaangle*TMath::DegToRad();
    else // track points to sector "after"
      angle += deltaangle*TMath::DegToRad();
    
    Bool_t success = trackParam.Rotate(angle);
    
    if(!success) {
      //      cout << "WARNING: Rotate failed for ID: " << id << endl;
      nFailures++;
    }
    
    angle *= TMath::RadToDeg();
    
    Int_t sector = TMath::Nint((angle-10.0)/20.0);
    if(sector<0)
      sector += 18;
    else if(sector>=18)
      sector -= 18;
    if(ltrp->GetSide()==1) // C side
      sector += 18;
    if(option==0 || option==2)
      sector += 36;
    if((option==1||option==3)&&(ltrp->GetBeam()==0||ltrp->GetBeam()==6)) 
      sector += 36;

    sectorArray[id] = sector;

    const Double_t x0 = 0;
    
    Double_t slopey  = TMath::Tan(TMath::ASin(trackParam.GetSnp()));
    Double_t slopez  = trackParam.GetTgl();
    //  One needs a factor sqrt(1+slopey**2) to take into account the
    //  longer path length
    slopez *= TMath::Sqrt(1.0 + slopey*slopey);    
    if(fInverseSlopeZ) // wrong sign in database, should be fixed there
      slopez *= -1;
    //    Double_t offsetz = trackParam.GetZ();
    Double_t offsety = trackParam.GetY() + slopey*(x0-trackParam.GetX());
    Double_t offsetz = trackParam.GetZ() + slopez*(x0-trackParam.GetX());
    if(option==2 || option==3) {
      meanOffset[id] = offsetz; meanSlope[id] = slopez;
    } else {
      meanOffset[id] = offsety; meanSlope[id] = slopey;
    }
  }

  if(nFailures>0)
    AliWarning(Form("Rotate method failed %d times", nFailures));
}



void AliTPCcalibLaser::DumpLaser(const char *finput, Int_t run){
  //
  //
  //input="TPCLaserObjects.root"
  //
  // 0. OBJ: TAxis     Delta
  // 1. OBJ: TAxis     bin
  // 2. OBJ: TAxis     rms shape
  // 3. OBJ: TAxis     sqrt(Q)
  // 4. OBJ: TAxis     row
  // 5. OBJ: TAxis     trackID

  const Double_t kSigma=4.;
  TFile f(finput);
  AliTPCcalibLaser *laserTPC = (AliTPCcalibLaser*) f.Get("laserTPC");
  THnSparse * hisPadInput   = laserTPC->fHisLaserPad;
  THnSparse * hisTimeInput = laserTPC->fHisLaserTime;
  TTreeSRedirector *pcstream= new TTreeSRedirector("hisLasers.root");
  TVectorD meanY(159), sigmaY(159);
  TVectorD meanZ(159), sigmaZ(159);
  TVectorD meanPad(159), sigmaPad(159);
  TVectorD meanTime(159), sigmaTime(159);
  TVectorD meanDPad(159), sigmaDPad(159);
  TVectorD meanDTime(159), sigmaDTime(159);
  TVectorD meandEdx(159), sigmadEdx(159);
  TVectorD meanSTime(159), sigmaSTime(159);
  TVectorD meanSPad(159), sigmaSPad(159);
  TVectorD entries(159);
  //
  Int_t indexes[10]={0,1,2,3,4,5,6};
  TH1 *his=0;
  AliTPCLaserTrack::LoadTracks();
  //
  for (Int_t id=0; id<336; id++){ // llop over laser beams 
    printf("id=\t%d\n",id);
    //
    AliTPCLaserTrack *ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    //
    hisPadInput->GetAxis(5)->SetRange(id+1,id+1);
    hisTimeInput->GetAxis(5)->SetRange(id+1,id+1);
    //
    his=hisTimeInput->Projection(3);
    Int_t firstBindEdx=his->FindFirstBinAbove(0);
    Int_t lastBindEdx=his->FindLastBinAbove(0);
    hisPadInput->GetAxis(3)->SetRange(firstBindEdx, lastBindEdx);
    hisTimeInput->GetAxis(3)->SetRange(firstBindEdx, lastBindEdx);
    delete his;
    //
    his=hisTimeInput->Projection(1);
    //    Int_t firstBinTime=his->FindFirstBinAbove(0);
    //Int_t lastBinTime=his->FindLastBinAbove(0);
    //hisTimeInput->GetAxis(1)->SetRange(firstBinTime, lastBinTime);
    delete his;
    //
    //
    his=hisTimeInput->Projection(2);
    //Int_t firstBinZ=his->FindFirstBinAbove(0);
    //Int_t lastBinZ=his->FindLastBinAbove(0);
    //hisTimeInput->GetAxis(2)->SetRange(firstBinZ, lastBinZ);
    delete his;
    //
    his=hisPadInput->Projection(2);
    //    Int_t firstBinY=his->FindFirstBinAbove(0);
    //Int_t lastBinY=his->FindLastBinAbove(0);
    //hisPadInput->GetAxis(2)->SetRange(firstBinY, lastBinY);
    delete his;
    //
    //
    //
    THnSparse *hisPad0  = hisPadInput->Projection(5,indexes);
    THnSparse *hisTime0 = hisTimeInput->Projection(5,indexes);
    //
    //    
    for (Int_t irow=0; irow<159; irow++){
      entries[irow]=0;
      if ((*(ltrp->GetVecSec()))[irow] <0) continue;
      if ((*(ltrp->GetVecLX()))[irow] <80) continue;

      hisPad0->GetAxis(4)->SetRange(irow+1,irow+1);
      hisTime0->GetAxis(4)->SetRange(irow+1,irow+1);
      //THnSparse *hisPad  = hisPad0->Projection(4,indexes);
      //THnSparse *hisTime = hisTime0->Projection(4,indexes);
      THnSparse *hisPad  = hisPad0;
      THnSparse *hisTime = hisTime0;
      //
      // Get mean value of QA variables
      //
      // dEdx
      his=hisTime->Projection(3);
      his->GetXaxis()->SetRangeUser(his->GetMean()-kSigma*his->GetRMS(), his->GetMean()+kSigma*his->GetRMS());
      meandEdx[irow] =his->GetMean();
      sigmadEdx[irow]=his->GetRMS();
      //      Int_t bindedx0= his->FindBin(meandEdx[irow]-kSigma*sigmadEdx[irow]);
      //Int_t bindedx1= his->FindBin(meandEdx[irow]+kSigma*sigmadEdx[irow]);
      //      hisPad->GetAxis(3)->SetRange(bindedx0,bindedx1);
      //hisTime->GetAxis(3)->SetRange(bindedx0,bindedx1 );
      delete his;
      //
      // sigma Time
      //
      his=hisTime->Projection(2);
      his->GetXaxis()->SetRangeUser(his->GetMean()-kSigma*his->GetRMS(), his->GetMean()-kSigma*his->GetRMS());
      meanSTime[irow] =his->GetMean();
      sigmaSTime[irow]=his->GetRMS();
      //Int_t binSTime0= his->FindBin(his->GetMean()-kSigma*his->GetRMS());
      //Int_t binSTime1= his->FindBin(his->GetMean()+kSigma*his->GetRMS());
      //      hisTime->GetAxis(2)->SetRange(binSTime0, binSTime1);
      delete his;
      //
      // sigma Pad
      his=hisPad->Projection(2);
      his->GetXaxis()->SetRangeUser(his->GetMean()-kSigma*his->GetRMS(), his->GetMean()+kSigma*his->GetRMS());
      meanSPad[irow] =his->GetMean();
      sigmaSPad[irow]=his->GetRMS();   
      //      Int_t binSPad0= his->FindBin(his->GetMean()-kSigma*his->GetRMS());
      //Int_t binSPad1= his->FindBin(his->GetMean()+kSigma*his->GetRMS());
      //      hisPad->GetAxis(2)->SetRange(binSPad0, binSPad1);
      delete his;
      //
      // apply selection on QA variables
      //
      //
      //
      // Y
      his=hisPad->Projection(0);
      entries[irow]=his->GetEntries();
      his->GetXaxis()->SetRangeUser(his->GetMean()-kSigma*his->GetRMS(), his->GetMean()+kSigma*his->GetRMS());
      meanY[irow] =his->GetMean();
      sigmaY[irow]=his->GetRMS();
      delete his;
      // Z
      his=hisTime->Projection(0);
      his->GetXaxis()->SetRangeUser(his->GetMean()-kSigma*his->GetRMS(), his->GetMean()+kSigma*his->GetRMS());
      meanZ[irow] =his->GetMean();
      sigmaZ[irow]=his->GetRMS();
      delete his;
      // Pad
      his=hisPad->Projection(1);
      his->GetXaxis()->SetRangeUser(his->GetMean()-kSigma*his->GetRMS(), his->GetMean()+kSigma*his->GetRMS());
      meanPad[irow] =his->GetMean();
      meanDPad[irow] =his->GetMean()-Int_t(his->GetMean());
      sigmaPad[irow]=his->GetRMS();
      delete his;
      // Time
      his=hisTime->Projection(1);
      his->GetXaxis()->SetRangeUser(his->GetMean()-kSigma*his->GetRMS(), his->GetMean()+kSigma*his->GetRMS());
      meanTime[irow]  = his->GetMean();
      meanDTime[irow] = his->GetMean()-Int_t(his->GetMean());
      sigmaTime[irow]=his->GetRMS();
      delete his;
      //
      //delete hisTime;
      //delete hisPad;
    }
    //
    //
    //
    (*pcstream)<<"laserClusters"<<
      "id="<<id<<      
      "run="<<run<<
      "LTr.="<<ltrp<<
      //
      "entries.="<<&entries<<
      "my.="<<&meanY<<           //mean delta y
      "rmsy.="<<&sigmaY<<        //rms deltay
      "mz.="<<&meanZ<<           //mean deltaz
      "rmsz.="<<&sigmaZ<<        //rms z
      //
      "mPad.="<<&meanPad<<       // mean pad
      "mDPad.="<<&meanDPad<<     // mead dpad
      "rmsPad.="<<&sigmaPad<<    // rms pad
      "mTime.="<<&meanTime<<     
      "mDTime.="<<&meanTime<<
      "rmsTime.="<<&sigmaTime<<
      //
      "mdEdx.="<<&meandEdx<<      //mean dedx
      "rmsdEdx.="<<&sigmadEdx<<   //rms dedx
      "mSPad.="<<&meanSPad<<      //mean sigma pad
      "rmsSPad.="<<&sigmaSPad<<   //rms sigma pad
      "mSTime.="<<&meanSTime<<    //mean sigma time
      "rmsSTime.="<<&sigmaSTime<<
      "\n";
    //
    delete hisPad0;
    delete hisTime0;
  }
  delete pcstream;

  /*
    
  */
}

void AliTPCcalibLaser::FitLaserClusters(Int_t run){
  //
  //
  //input="TPCLaserObjects.root"
  //Algorithm:
  //   1. Select cluster candidates, remove outlyers
  //             edge clusters
  //             clusters with atypical spread (e.g due track overlaps)
  //             small amount of entries clusters (absolute minimal cut + raltive -to mean cut) 
  //   2. Fit the tracklets -per sector - in pad and time coordinate frame
  //             Remove outlyers
  //             Store info distance of track to pad, time center
  //             Fit the correction for distance to the center (sin,cos)
  //   3. Do local fit
  const Double_t kEpsilon=0.000001;
  const Int_t     kMinClusters=20;
  const Double_t kEdgeCut=3;
  const Double_t kDistCut=1.5;                 // cut distance to the ideal track
  const Double_t kDistCutFit=0.5;              
  const Double_t kDistCutFitPad=0.25;
  const Double_t kDistCutFitTime=0.25;
  const Int_t kSmoothRow=5;
  TFile f("hisLasers.root");  // Input file
  TTree * treeInput=(TTree*)f.Get("laserClusters"); 
  TTreeSRedirector *pcstream=new TTreeSRedirector("fitLasers.root");
  TVectorD *vecN=0;
  TVectorD *vecMY=0;
  TVectorD *vecMZ=0;
  TVectorD *vecPad=0;
  TVectorD *vecTime=0;
  TVectorD *vecSY=0;
  TVectorD *vecSZ=0;
  TVectorD *meandEdx=0;
  TVectorD  isOK(159);
  TVectorD  fitPad(159);
  TVectorD  fitTime(159);
  TVectorD  fitPadLocal(159);
  TVectorD  fitTimeLocal(159);
  TVectorD  fitDPad(159);
  TVectorD  fitDTime(159);
  TVectorD  fitIPad(159);
  TVectorD  fitITime(159);
  Double_t  chi2PadIROC=0;
  Double_t  chi2PadOROC=0;
  //
  treeInput->SetBranchAddress("my.",&vecMY);
  treeInput->SetBranchAddress("mz.",&vecMZ);
  treeInput->SetBranchAddress("mPad.",&vecPad);
  treeInput->SetBranchAddress("mTime.",&vecTime);
  treeInput->SetBranchAddress("rmsy.",&vecSY);
  treeInput->SetBranchAddress("rmsz.",&vecSZ);
  treeInput->SetBranchAddress("entries.",&vecN);
  treeInput->SetBranchAddress("mdEdx.",&meandEdx);

  AliTPCLaserTrack::LoadTracks();
  //
  //
  TVectorD fitPadIROC(3),     fitPadOROC(3);
  TVectorD fitPadIROCSin(3),  fitPadOROCSin(3);
  TVectorD fitTimeIROC(3),    fitTimeOROC(3);
  TVectorD fitTimeIROCSin(3), fitTimeOROCSin(3);
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  Double_t refX=roc->GetPadRowRadii(0,roc->GetNRows(0)-1);
  
  //
  for (Int_t id=0; id<336; id++){
    //
    treeInput->GetEntry(id);
    AliTPCLaserTrack *ltrp =(AliTPCLaserTrack*)AliTPCLaserTrack::GetTracks()->UncheckedAt(id);
    Int_t medianEntries = TMath::Nint(TMath::Median(159,vecN->GetMatrixArray()));
    Double_t medianRMSY = TMath::Median(159,vecSY->GetMatrixArray());
    Double_t rmsRMSY    = TMath::RMS(159,vecSY->GetMatrixArray());
    Double_t medianRMSZ = TMath::Median(159,vecSZ->GetMatrixArray());
    Double_t rmsRMSZ    = TMath::RMS(159,vecSZ->GetMatrixArray());
    Double_t mdEdx      = TMath::Median(159,meandEdx->GetMatrixArray());
    Int_t sectorInner= TMath::Nint(ltrp->GetVecSec()->GetMatrixArray()[63/2]);
    Int_t sectorOuter= TMath::Nint(ltrp->GetVecSec()->GetMatrixArray()[64+96/2]);
    TLinearFitter fitterY(2,"pol1");
    TLinearFitter fitterZ(2,"pol1");
    TLinearFitter fitterPad(2,"pol1");
    TLinearFitter fitterTime(2,"pol1");
    TLinearFitter fitterPadSin(2,"hyp1");
    TLinearFitter fitterTimeSin(3,"hyp2");    
    //
    //
    for (UInt_t irow=0; irow<159; irow++){
      fitPad[irow]=0; fitIPad[irow]=0; fitDPad[irow]=0;
      fitTime[irow]=0; fitITime[irow]=0; fitDTime[irow]=0;
      Double_t sign=(ltrp->GetZ()>0) ? 1.:-1.;
      isOK[irow]=kFALSE;    
      fitPad[irow]=0;
      fitTime[irow]=0;
      Int_t sector=(irow<roc->GetNRows(0))? sectorInner:sectorOuter;      
      Int_t npads=(irow<roc->GetNRows(0))? roc->GetNPads(sector,irow):roc->GetNPads(sector,irow-roc->GetNRows(0));
      (*vecPad)[irow]-=npads*0.5;	
      //
      if ((irow<roc->GetNRows(0)) &&TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[irow]-sectorInner)>0.1) continue;
      if ((irow>=roc->GetNRows(0)) &&TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[irow]-sectorOuter)>0.1) continue;
      //
      if (TMath::Abs((*vecMY)[irow])<kEpsilon) continue;   //not determined position
      if (TMath::Abs((*vecMZ)[irow])<kEpsilon) continue;   //not determined position
      if (TMath::Abs((*vecPad)[irow])<kEpsilon) continue;   //not determined position
      if (TMath::Abs((*vecTime)[irow])<kEpsilon) continue;   //not determined position
      if ((*vecN)[irow]<0.5*medianEntries) continue;       //small amount of clusters
      if ((*vecSY)[irow]>medianRMSY+3*rmsRMSY) continue;   //big sigma
      if ((*vecSZ)[irow]>medianRMSZ+3*rmsRMSZ) continue;   //big sigma
      Double_t dEdge= TMath::Abs((*(ltrp->GetVecLY()))[irow])-(*(ltrp->GetVecLX()))[irow]*TMath::Tan(TMath::Pi()/18.); //edge cut
      if (TMath::Abs(dEdge)<kEdgeCut) continue;
      if (irow<roc->GetNRows(0)){
	if (TMath::Abs(((*ltrp->GetVecLY())[irow])-sign*(*vecPad)[irow]*0.4)>kDistCut) continue;
      }
      if (irow>roc->GetNRows(0)){
	if (TMath::Abs(((*ltrp->GetVecLY())[irow])-sign*(*vecPad)[irow]*0.6)>kDistCut) continue;
      }
      
      isOK[irow]=kTRUE;     
    }
    //
    //fit OROC - get delta pad and delta time 
    //
    fitterPad.ClearPoints();
    fitterTime.ClearPoints();    
    fitterPadSin.ClearPoints();
    fitterTimeSin.ClearPoints();    
    {for (Int_t irow=2; irow<157; irow++){
	if (isOK[irow]<0.5) continue;	
	if (TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[irow]-sectorOuter)>0.1) continue;
	if (TMath::Abs(ltrp->GetVecLX()->GetMatrixArray()[irow])<80) continue;
	Double_t y=(*vecPad)[irow];
	Double_t z=(*vecTime)[irow];
	Double_t x=ltrp->GetVecLX()->GetMatrixArray()[irow]-refX;
	fitterPad.AddPoint(&x,y);
	fitterTime.AddPoint(&x,z);
      }}
    chi2PadOROC=0;
    if (fitterPad.GetNpoints()>kMinClusters&&fitterTime.GetNpoints()>kMinClusters){
      fitterPad.Eval();
      fitterTime.Eval();
      chi2PadOROC=TMath::Sqrt(fitterPad.GetChisquare()/fitterPad.GetNpoints());
      for (Int_t irow=2; irow<157; irow++){
	if (isOK[irow]<0.5) continue;	
	if (TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[irow]-sectorOuter)>0.1) continue;
	if (TMath::Abs(ltrp->GetVecLX()->GetMatrixArray()[irow])<80) continue;
	Double_t y=(*vecPad)[irow];
	Double_t z=(*vecTime)[irow];
	Double_t x=ltrp->GetVecLX()->GetMatrixArray()[irow]-refX;
	Double_t fitP=fitterPad.GetParameter(0)+fitterPad.GetParameter(1)*x;
	Double_t fitT=fitterTime.GetParameter(0)+fitterTime.GetParameter(1)*x;
	fitPad[irow]=fitterPad.GetParameter(0)+fitterPad.GetParameter(1)*x;
	fitTime[irow]=fitterTime.GetParameter(0)+fitterTime.GetParameter(1)*x;
	fitDPad[irow]=y-(fitterPad.GetParameter(0)+fitterPad.GetParameter(1)*x);
	fitDTime[irow]=z-(fitterTime.GetParameter(0)+fitterTime.GetParameter(1)*x);
	fitIPad[irow]=fitP-TMath::Nint(fitP-0.5);
	fitITime[irow]=fitT-TMath::Nint(fitT-0.5);
	if (fitDPad[irow]>kDistCutFit) isOK[irow]=kFALSE;
	if (fitDTime[irow]>kDistCutFit) isOK[irow]=kFALSE;
	if (isOK[irow]>0){
	  Double_t xxxPad[2]={TMath::Sin(2*TMath::Pi()*fitIPad[irow])};
	  Double_t xxxTime[3]={TMath::Sin(2*TMath::Pi()*fitITime[irow]),
			       TMath::Cos(2*TMath::Pi()*fitITime[irow])};
	  fitterPadSin.AddPoint(xxxPad,fitDPad[irow]);
	  fitterTimeSin.AddPoint(xxxTime,fitDTime[irow]);
	}
      }
      fitterPadSin.Eval();
      fitterTimeSin.Eval();
      fitterPadSin.FixParameter(0,0);
      fitterTimeSin.FixParameter(0,0);
      fitterPadSin.Eval();
      fitterTimeSin.Eval();
      //
      fitterPad.GetParameters(fitPadOROC);
      fitterTime.GetParameters(fitTimeOROC);
      fitterPadSin.GetParameters(fitPadOROCSin);
      fitterTimeSin.GetParameters(fitTimeOROCSin);
    }
    //
    //
    //fit IROC
    //
    fitterPad.ClearPoints();
    fitterTime.ClearPoints();    
    fitterPadSin.ClearPoints();
    fitterTimeSin.ClearPoints();    
    for (Int_t irow=2; irow<157; irow++){
      if (isOK[irow]<0.5) continue;	
      if (TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[irow]-sectorInner)>0.1) continue;
      if (TMath::Abs(ltrp->GetVecLX()->GetMatrixArray()[irow])<80) continue;
      Double_t y=(*vecPad)[irow];
      Double_t z=(*vecTime)[irow];
      Double_t x=ltrp->GetVecLX()->GetMatrixArray()[irow]-refX;
      fitterPad.AddPoint(&x,y);
      fitterTime.AddPoint(&x,z);
    }
    chi2PadIROC=0;
    if (fitterPad.GetNpoints()>kMinClusters&&fitterTime.GetNpoints()>kMinClusters){
      fitterPad.Eval();
      fitterTime.Eval();
      chi2PadIROC=TMath::Sqrt(fitterPad.GetChisquare()/fitterPad.GetNpoints());
      for (Int_t irow=2; irow<157; irow++){
	if (isOK[irow]<0.5) continue;	
	if (TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[irow]-sectorInner)>0.1) continue;
	if (TMath::Abs(ltrp->GetVecLX()->GetMatrixArray()[irow])<80) continue;
	Double_t y=(*vecPad)[irow];
	Double_t z=(*vecTime)[irow];
	Double_t x=ltrp->GetVecLX()->GetMatrixArray()[irow]-refX;
	Double_t fitP=fitterPad.GetParameter(0)+fitterPad.GetParameter(1)*x;
	Double_t fitT=fitterTime.GetParameter(0)+fitterTime.GetParameter(1)*x;
	fitPad[irow]=fitterPad.GetParameter(0)+fitterPad.GetParameter(1)*x;
	fitTime[irow]=fitterTime.GetParameter(0)+fitterTime.GetParameter(1)*x;
	fitDPad[irow]=y-(fitterPad.GetParameter(0)+fitterPad.GetParameter(1)*x);
	fitDTime[irow]=z-(fitterTime.GetParameter(0)+fitterTime.GetParameter(1)*x);
	fitIPad[irow]=fitP-TMath::Nint(fitP-0.5);
	fitITime[irow]=fitT-TMath::Nint(fitT-0.5);
	if (fitDPad[irow]>kDistCutFit) isOK[irow]=kFALSE;
	if (fitDTime[irow]>kDistCutFit) isOK[irow]=kFALSE;
	if (isOK[irow]>0.5){
	  Double_t xxxPad[3]={TMath::Sin(2*TMath::Pi()*fitIPad[irow]),
			 TMath::Cos(2*TMath::Pi()*fitIPad[irow])};
	  Double_t xxxTime[3]={TMath::Sin(2*TMath::Pi()*fitITime[irow]),
			 TMath::Cos(2*TMath::Pi()*fitITime[irow])};
	  fitterPadSin.AddPoint(xxxPad,fitDPad[irow]);
	  fitterTimeSin.AddPoint(xxxTime,fitDTime[irow]);
	}
      }
      fitterPadSin.Eval();
      fitterTimeSin.Eval();
      fitterPadSin.FixParameter(0,0);
      fitterTimeSin.FixParameter(0,0);
      fitterPadSin.Eval();
      fitterTimeSin.Eval();
      fitterPad.GetParameters(fitPadIROC);
      fitterTime.GetParameters(fitTimeIROC);
      fitterPadSin.GetParameters(fitPadIROCSin);
      fitterTimeSin.GetParameters(fitTimeIROCSin);
    }
    for (Int_t irow=0; irow<159; irow++){
      if (TMath::Abs(fitDPad[irow])<kEpsilon)  isOK[irow]=kFALSE;
      if (TMath::Abs(fitDTime[irow])<kEpsilon) isOK[irow]=kFALSE;
      if (TMath::Abs(fitDPad[irow])>kDistCutFitPad)  isOK[irow]=kFALSE;
      if (TMath::Abs(fitDTime[irow])>kDistCutFitTime) isOK[irow]=kFALSE;
    }
    for (Int_t irow=kSmoothRow/2; irow<159-kSmoothRow/2; irow++){
      fitPadLocal[irow]=0;
      fitTimeLocal[irow]=0;
      if (isOK[irow]<0.5) continue;     
      Int_t sector=(irow<Int_t(roc->GetNRows(0)))? sectorInner:sectorOuter;
      if (TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[irow]-sector)>0.1) continue;
      //
      TLinearFitter fitterPadLocal(2,"pol1");
      TLinearFitter fitterTimeLocal(2,"pol1");
      Double_t xref=ltrp->GetVecLX()->GetMatrixArray()[irow];
      for (Int_t delta=-kSmoothRow; delta<=kSmoothRow; delta++){
	Int_t jrow=irow+delta;
	if (jrow<0) jrow=0;
	if (jrow>159) jrow=159;	
	if (isOK[jrow]<0.5) continue;	
	if (TMath::Abs(ltrp->GetVecSec()->GetMatrixArray()[jrow]-sector)>0.1) continue;
      	Double_t y=(*vecPad)[jrow];
	Double_t z=(*vecTime)[jrow];
	Double_t x=ltrp->GetVecLX()->GetMatrixArray()[jrow]-xref;
	fitterPadLocal.AddPoint(&x,y);
	fitterTimeLocal.AddPoint(&x,z);
      }      
      if (fitterPadLocal.GetNpoints()<kSmoothRow) continue;
      fitterPadLocal.Eval();
      fitterTimeLocal.Eval();
      fitPadLocal[irow]=fitterPadLocal.GetParameter(0);
      fitTimeLocal[irow]=fitterTimeLocal.GetParameter(0);
    }
    //
    //
    (*pcstream)<<"fit"<<
      "run="<<run<<
      "id="<<id<<
      "chi2PadIROC="<<chi2PadIROC<<
      "chi2PadOROC="<<chi2PadOROC<<
      "mdEdx="<<mdEdx<<
      "LTr.="<<ltrp<<
      "isOK.="<<&isOK<<
      // mean measured-ideal values
      "mY.="<<vecMY<<
      "mZ.="<<vecMZ<<
      // local coordinate fit
      "mPad.="<<vecPad<<
      "mTime.="<<vecTime<<
      "fitPad.="<<&fitPad<<
      "fitTime.="<<&fitTime<<
      "fitPadLocal.="<<&fitPadLocal<<
      "fitTimeLocal.="<<&fitTimeLocal<<
      "fitDPad.="<<&fitDPad<<
      "fitDTime.="<<&fitDTime<<
      "fitIPad.="<<&fitIPad<<
      "fitITime.="<<&fitITime<<         
      //
      "fitPadIROC.="<<&fitPadIROC<<           // pad fit linear IROC
      "fitPadIROCSin.="<<&fitPadIROCSin<<     // pad fit linear+ pad correction
      "fitPadOROC.="<<&fitPadOROC<<
      "fitPadOROCSin.="<<&fitPadOROCSin<<
      //
      "fitTimeIROC.="<<&fitTimeIROC<<
      "fitTimeIROCSin.="<<&fitTimeIROCSin<<
      "fitTimeOROC.="<<&fitTimeOROC<<
      "fitTimeOROCSin.="<<&fitTimeOROCSin<<
      "\n";    
  }  
  delete pcstream;
}
