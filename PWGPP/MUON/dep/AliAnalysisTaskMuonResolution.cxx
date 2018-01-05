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

/* $Id$ */

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <Riostream.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TList.h>
#include <TObjString.h>
#include <TRegexp.h>

// STEER includes
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliGeomManager.h"
#include "AliVVertex.h"

// ANALYSIS includes
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisTaskMuonResolution.h"
#include "AliMultSelection.h"

// MUON includes
#include "AliMUONCDB.h"
#include "AliMUONConstants.h"
#include "AliMUONRecoParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONVCluster.h"
#include "AliMUONVDigit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMpDEIterator.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpPad.h"
#include "AliMpDetElement.h"
#include "AliMpCathodType.h"

#ifndef SafeDelete
#define SafeDelete(x) if (x != NULL) { delete x; x = NULL; }
#endif

using std::cout;
using std::endl;
using std::flush;

ClassImp(AliAnalysisTaskMuonResolution)

const Int_t AliAnalysisTaskMuonResolution::fgkMinEntries = 10;

//________________________________________________________________________
AliAnalysisTaskMuonResolution::AliAnalysisTaskMuonResolution() :
  AliAnalysisTaskSE(), 
  fResiduals(NULL),
  fResidualsVsP(NULL),
  fResidualsVsCent(NULL),
  fResidualsVsAngle(NULL),
  fLocalChi2(NULL),
  fChamberRes(NULL),
  fTrackRes(NULL),
  fCanvases(NULL),
  fTmpHists(NULL),
  fDefaultStorage(""),
  fNEvents(0),
  fShowProgressBar(kFALSE),
  fPrintClResPerCh(kFALSE),
  fPrintClResPerDE(kFALSE),
  fGaus(NULL),
  fMinMomentum(0.),
  fMinPt(0.),
  fSign(0),
  fUseMCLabel(kFALSE),
  fExtrapMode(1),
  fCorrectForSystematics(kTRUE),
  fRemoveMonoCathCl(kFALSE),
  fCheckAllPads(kFALSE),
  fImproveTracks(kFALSE),
  fShiftHalfCh(kFALSE),
  fPrintHalfChShift(kFALSE),
  fShiftDE(kFALSE),
  fPrintDEShift(kFALSE),
  fOCDBLoaded(kFALSE),
  fNDE(0),
  fReAlign(kFALSE),
  fOldAlignStorage(""),
  fOldAlignVersion(-1),
  fOldAlignSubVersion(-1),
  fNewAlignStorage(""),
  fNewAlignVersion(-1),
  fNewAlignSubVersion(-1),
  fOldGeoTransformer(NULL),
  fNewGeoTransformer(NULL),
  fMuonEventCuts(0x0),
  fMuonTrackCuts(0x0)
{
  /// Default constructor
  
  for (Int_t i = 0; i < 10; i++) SetStartingResolution(i, -1., -1.);
  for (Int_t i = 0; i < 1100; i++) fDEIndices[i] = 0;
  for (Int_t i = 0; i < 200; i++) fDEIds[i] = 0;
  for (Int_t i = 0; i < 20; i++) SetHalfChShift(i, 0., 0.);
  for (Int_t i = 0; i < 200; i++) SetDEShift(i, 0., 0.);
  
}

//________________________________________________________________________
AliAnalysisTaskMuonResolution::AliAnalysisTaskMuonResolution(const char *name) :
  AliAnalysisTaskSE(name), 
  fResiduals(NULL),
  fResidualsVsP(NULL),
  fResidualsVsCent(NULL),
  fResidualsVsAngle(NULL),
  fLocalChi2(NULL),
  fChamberRes(NULL),
  fTrackRes(NULL),
  fCanvases(NULL),
  fTmpHists(NULL),
  fDefaultStorage("raw://"),
  fNEvents(0),
  fShowProgressBar(kFALSE),
  fPrintClResPerCh(kFALSE),
  fPrintClResPerDE(kFALSE),
  fGaus(NULL),
  fMinMomentum(0.),
  fMinPt(0.),
  fSign(0),
  fUseMCLabel(kFALSE),
  fExtrapMode(1),
  fCorrectForSystematics(kTRUE),
  fRemoveMonoCathCl(kFALSE),
  fCheckAllPads(kFALSE),
  fImproveTracks(kFALSE),
  fShiftHalfCh(kFALSE),
  fPrintHalfChShift(kFALSE),
  fShiftDE(kFALSE),
  fPrintDEShift(kFALSE),
  fOCDBLoaded(kFALSE),
  fNDE(0),
  fReAlign(kFALSE),
  fOldAlignStorage(""),
  fOldAlignVersion(-1),
  fOldAlignSubVersion(-1),
  fNewAlignStorage(""),
  fNewAlignVersion(-1),
  fNewAlignSubVersion(-1),
  fOldGeoTransformer(NULL),
  fNewGeoTransformer(NULL),
  fMuonEventCuts(0x0),
  fMuonTrackCuts(0x0)
{
  /// Constructor
  
  for (Int_t i = 0; i < 10; i++) SetStartingResolution(i, -1., -1.);
  for (Int_t i = 0; i < 1100; i++) fDEIndices[i] = 0;
  for (Int_t i = 0; i < 200; i++) fDEIds[i] = 0;
  for (Int_t i = 0; i < 20; i++) SetHalfChShift(i, 0., 0.);
  for (Int_t i = 0; i < 200; i++) SetDEShift(i, 0., 0.);
  
  FitResiduals();
  
  // Output slot #1 writes into a TObjArray container
  DefineOutput(1,TObjArray::Class());
  // Output slot #2 writes into a TObjArray container
  DefineOutput(2,TObjArray::Class());
  // Output slot #3 writes into a TObjArray container
  DefineOutput(3,TObjArray::Class());
  // Output slot #4 writes into a TObjArray container
  DefineOutput(4,TObjArray::Class());
  // Output slot #5 writes into a TObjArray container
  DefineOutput(5,TObjArray::Class());
  // Output slot #6 writes into a TObjArray container
  DefineOutput(6,TObjArray::Class());
  // Output slot #7 writes into a TObjArray container
  DefineOutput(7,TObjArray::Class());
}

//________________________________________________________________________
AliAnalysisTaskMuonResolution::~AliAnalysisTaskMuonResolution()
{
  /// Destructor
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    SafeDelete(fResiduals);
    SafeDelete(fResidualsVsP);
    SafeDelete(fResidualsVsCent);
    SafeDelete(fResidualsVsAngle);
    SafeDelete(fTrackRes);
  }
  SafeDelete(fLocalChi2);
  SafeDelete(fChamberRes);
  SafeDelete(fCanvases);
  SafeDelete(fTmpHists);
  SafeDelete(fGaus);
  SafeDelete(fOldGeoTransformer);
  SafeDelete(fNewGeoTransformer);
  SafeDelete(fMuonEventCuts);
  SafeDelete(fMuonTrackCuts);
}

//___________________________________________________________________________
void AliAnalysisTaskMuonResolution::UserCreateOutputObjects()
{
  /// Create histograms
  
  // do it once the OCDB has been loaded (i.e. from NotifyRun())
  if (!fOCDBLoaded) return;
  
  fResiduals = new TObjArray(1000);
  fResiduals->SetOwner();
  fResidualsVsP = new TObjArray(1000);
  fResidualsVsP->SetOwner();
  fResidualsVsCent = new TObjArray(1000);
  fResidualsVsCent->SetOwner();
  fResidualsVsAngle = new TObjArray(1000);
  fResidualsVsAngle->SetOwner();
  fTrackRes = new TObjArray(1000);
  fTrackRes->SetOwner();
  TH2F* h2;
  
  // find the highest chamber resolution and set histogram bins
  const AliMUONRecoParam* recoParam = AliMUONESDInterface::GetTracker()->GetRecoParam();
  Double_t maxSigma[2] = {-1., -1.};
  for (Int_t i = 0; i < 10; i++) {
    if (recoParam->GetDefaultNonBendingReso(i) > maxSigma[0]) maxSigma[0] = recoParam->GetDefaultNonBendingReso(i);
    if (recoParam->GetDefaultBendingReso(i) > maxSigma[1]) maxSigma[1] = recoParam->GetDefaultBendingReso(i);
  }
  const char* axes[2] = {"X", "Y"};
  const char* side[2] = {"I", "O"};
  const Int_t nBins = 5000;
  const Int_t nSigma = 10;
  const Int_t pNBins = 20;
  const Double_t pEdges[2] = {0., 50.};
  Int_t nCentBins = 12;
  Double_t centRange[2] = {-10., 110.};
  Int_t nAngleBins = 20;
  Double_t angleRange[2][2] = {{-15., 15.}, {-40., 40.}};
  TString name, title;
  
  for (Int_t ia = 0; ia < 2; ia++) {
    
    Double_t maxRes = nSigma*maxSigma[ia];
    
    // List of residual histos
    name = Form("hResidual%sPerCh_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s distribution per chamber (cluster attached to the track);chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 10, 0.5, 10.5, nBins, -maxRes, maxRes);
    fResiduals->AddAtAndExpand(h2, kResidualPerChClusterIn+ia);
    name = Form("hResidual%sPerCh_ClusterOut",axes[ia]); title = Form("cluster-track residual-%s distribution per chamber (cluster not attached to the track);chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 10, 0.5, 10.5, nBins, -2.*maxRes, 2.*maxRes);
    fResiduals->AddAtAndExpand(h2, kResidualPerChClusterOut+ia);
    
    name = Form("hResidual%sPerHalfCh_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s distribution per half chamber (cluster attached to the track);half chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 20, 0.5, 20.5, nBins, -maxRes, maxRes);
    for (Int_t i = 0; i < 10; i++) { h2->GetXaxis()->SetBinLabel(2*i+1, Form("%d-I",i+1)); h2->GetXaxis()->SetBinLabel(2*i+2, Form("%d-O",i+1)); }
    fResiduals->AddAtAndExpand(h2, kResidualPerHalfChClusterIn+ia);
    name = Form("hResidual%sPerHalfCh_ClusterOut",axes[ia]); title = Form("cluster-track residual-%s distribution per half chamber (cluster not attached to the track);half chamber ID;#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 20, 0.5, 20.5, nBins, -2.*maxRes, 2.*maxRes);
    for (Int_t i = 0; i < 10; i++) { h2->GetXaxis()->SetBinLabel(2*i+1, Form("%d-I",i+1)); h2->GetXaxis()->SetBinLabel(2*i+2, Form("%d-O",i+1)); }
    fResiduals->AddAtAndExpand(h2, kResidualPerHalfChClusterOut+ia);
    
    name = Form("hResidual%sPerDE_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s distribution per DE (cluster attached to the track);DE ID;#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), fNDE, 0.5, fNDE+0.5, nBins, -maxRes, maxRes);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kResidualPerDEClusterIn+ia);
    name = Form("hResidual%sPerDE_ClusterOut",axes[ia]); title = Form("cluster-track residual-%s distribution per DE (cluster not attached to the track);DE ID;#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), fNDE, 0.5, fNDE+0.5, nBins, -2.*maxRes, 2.*maxRes);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kResidualPerDEClusterOut+ia);
    
    name = Form("hTrackRes%sPerCh",axes[ia]); title = Form("track #sigma_{%s} per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 10, 0.5, 10.5, nBins, 0., maxRes);
    fResiduals->AddAtAndExpand(h2, kTrackResPerCh+ia);
    name = Form("hTrackRes%sPerHalfCh",axes[ia]); title = Form("track #sigma_{%s} per half Ch;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 20, 0.5, 20.5, nBins, 0., maxRes);
    for (Int_t i = 0; i < 10; i++) { h2->GetXaxis()->SetBinLabel(2*i+1, Form("%d-I",i+1)); h2->GetXaxis()->SetBinLabel(2*i+2, Form("%d-O",i+1)); }
    fResiduals->AddAtAndExpand(h2, kTrackResPerHalfCh+ia);
    name = Form("hTrackRes%sPerDE",axes[ia]); title = Form("track #sigma_{%s} per DE;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), fNDE, 0.5, fNDE+0.5, nBins, 0., maxRes);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kTrackResPerDE+ia);
    
    name = Form("hMCS%sPerCh",axes[ia]); title = Form("MCS %s-dispersion of extrapolated track per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 10, 0.5, 10.5, nBins, 0., 0.2);
    fResiduals->AddAtAndExpand(h2, kMCSPerCh+ia);
    name = Form("hMCS%sPerHalfCh",axes[ia]); title = Form("MCS %s-dispersion of extrapolated track per half Ch;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 20, 0.5, 20.5, nBins, 0., 0.2);
    for (Int_t i = 0; i < 10; i++) { h2->GetXaxis()->SetBinLabel(2*i+1, Form("%d-I",i+1)); h2->GetXaxis()->SetBinLabel(2*i+2, Form("%d-O",i+1)); }
    fResiduals->AddAtAndExpand(h2, kMCSPerHalfCh+ia);
    name = Form("hMCS%sPerDE",axes[ia]); title = Form("MCS %s-dispersion of extrapolated track per DE;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), fNDE, 0.5, fNDE+0.5, nBins, 0., 0.2);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kMCSPerDE+ia);
    
    name = Form("hClusterRes2%sPerCh",axes[ia]); title = Form("cluster #sigma_{%s}^{2} per Ch;chamber ID;#sigma_{%s}^{2} (cm^{2})",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 10, 0.5, 10.5, nSigma*nBins, -0.1*maxRes*maxRes, maxRes*maxRes);
    fResiduals->AddAtAndExpand(h2, kClusterRes2PerCh+ia);
    
    // List of residual histos vs. p or vs. centrality or vs. angle
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      name = Form("hResidual%sInCh%dVsP_ClusterIn",axes[ia],i+1); title = Form("cluster-track residual-%s distribution in chamber %d versus momentum (cluster attached to the track);p (GeV/c^{2});#Delta_{%s} (cm)",axes[ia],i+1,axes[ia]);
      h2 = new TH2F(name.Data(), title.Data(), pNBins, pEdges[0], pEdges[1], nBins, -maxRes, maxRes);
      fResidualsVsP->AddAtAndExpand(h2, kResidualInChVsPClusterIn+10*ia+i);
      name = Form("hResidual%sInCh%dVsP_ClusterOut",axes[ia],i+1); title = Form("cluster-track residual-%s distribution in chamber %d versus momentum (cluster not attached to the track);p (GeV/c^{2});#Delta_{%s} (cm)",axes[ia],i+1,axes[ia]);
      h2 = new TH2F(name.Data(), title.Data(), pNBins, pEdges[0], pEdges[1], nBins, -2.*maxRes, 2.*maxRes);
      fResidualsVsP->AddAtAndExpand(h2, kResidualInChVsPClusterOut+10*ia+i);
      
      name = Form("hResidual%sInCh%dVsCent_ClusterIn",axes[ia],i+1); title = Form("cluster-track residual-%s distribution in chamber %d versus centrality (cluster attached to the track);centrality (%%);#Delta_{%s} (cm)",axes[ia],i+1,axes[ia]);
      h2 = new TH2F(name.Data(), title.Data(), nCentBins, centRange[0], centRange[1], nBins, -maxRes, maxRes);
      fResidualsVsCent->AddAtAndExpand(h2, kResidualInChVsCentClusterIn+10*ia+i);
      name = Form("hResidual%sInCh%dVsCent_ClusterOut",axes[ia],i+1); title = Form("cluster-track residual-%s distribution in chamber %d versus centrality (cluster not attached to the track);centrality (%%);#Delta_{%s} (cm)",axes[ia],i+1,axes[ia]);
      h2 = new TH2F(name.Data(), title.Data(), nCentBins, centRange[0], centRange[1], nBins, -2.*maxRes, 2.*maxRes);
      fResidualsVsCent->AddAtAndExpand(h2, kResidualInChVsCentClusterOut+10*ia+i);
      
      name = Form("hResidual%sInCh%dVsAngle_ClusterIn",axes[ia],i+1); title = Form("cluster-track residual-%s distribution in chamber %d versus track angle (cluster attached to the track);%s-angle (deg);#Delta_{%s} (cm)",axes[ia],i+1,axes[ia],axes[ia]);
      h2 = new TH2F(name.Data(), title.Data(), nAngleBins, angleRange[ia][0], angleRange[ia][1], nBins, -maxRes, maxRes);
      fResidualsVsAngle->AddAtAndExpand(h2, kResidualInChVsAngleClusterIn+10*ia+i);
      name = Form("hResidual%sInCh%dVsAngle_ClusterOut",axes[ia],i+1); title = Form("cluster-track residual-%s distribution in chamber %d versus track angle (cluster not attached to the track);%s-angle (deg);#Delta_{%s} (cm)",axes[ia],i+1,axes[ia],axes[ia]);
      h2 = new TH2F(name.Data(), title.Data(), nAngleBins, angleRange[ia][0], angleRange[ia][1], nBins, -2.*maxRes, 2.*maxRes);
      fResidualsVsAngle->AddAtAndExpand(h2, kResidualInChVsAngleClusterOut+10*ia+i);
      
      for (Int_t j = 0; j < 2; j++) {
	name = Form("hResidual%sInHalfCh%d%sVsAngle_ClusterIn",axes[ia],i+1,side[j]); title = Form("cluster-track residual-%s distribution in half-chamber %d-%s versus track angle (cluster attached to the track);%s-angle (deg);#Delta_{%s} (cm)",axes[ia],i+1,side[j],axes[ia],axes[ia]);
	h2 = new TH2F(name.Data(), title.Data(), nAngleBins, angleRange[ia][0], angleRange[ia][1], nBins, -maxRes, maxRes);
	fResidualsVsAngle->AddAtAndExpand(h2, kResidualInHalfChVsAngleClusterIn+20*ia+2*i+j);
      }
    }
    
    name = Form("hResidual%sVsP_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s distribution integrated over chambers versus momentum (cluster attached to the track);p (GeV/c^{2});#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), pNBins, pEdges[0], pEdges[1], nBins, -maxRes, maxRes);
    fResidualsVsP->AddAtAndExpand(h2, kResidualVsPClusterIn+ia);
    name = Form("hResidual%sVsP_ClusterOut",axes[ia]); title = Form("cluster-track residual-%s distribution integrated over chambers versus momentum (cluster not attached to the track);p (GeV/c^{2});#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), pNBins, pEdges[0], pEdges[1], nBins, -2.*maxRes, 2.*maxRes);
    fResidualsVsP->AddAtAndExpand(h2, kResidualVsPClusterOut+ia);
    
    name = Form("hResidual%sVsCent_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s distribution integrated over chambers versus centrality (cluster attached to the track);centrality (%%);#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), nCentBins, centRange[0], centRange[1], nBins, -maxRes, maxRes);
    fResidualsVsCent->AddAtAndExpand(h2, kResidualVsCentClusterIn+ia);
    name = Form("hResidual%sVsCent_ClusterOut",axes[ia]); title = Form("cluster-track residual-%s distribution integrated over chambers versus centrality (cluster not attached to the track);centrality (%%);#Delta_{%s} (cm)",axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), nCentBins, centRange[0], centRange[1], nBins, -2.*maxRes, 2.*maxRes);
    fResidualsVsCent->AddAtAndExpand(h2, kResidualVsCentClusterOut+ia);
    
    name = Form("hResidual%sVsAngle_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s distribution integrated over chambers versus track angle (cluster attached to the track);%s-angle (deg);#Delta_{%s} (cm)",axes[ia],axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), nAngleBins, angleRange[ia][0], angleRange[ia][1], nBins, -maxRes, maxRes);
    fResidualsVsAngle->AddAtAndExpand(h2, kResidualVsAngleClusterIn+ia);
    name = Form("hResidual%sVsAngle_ClusterOut",axes[ia]); title = Form("cluster-track residual-%s distribution integrated over chambers versus track angle (cluster not attached to the track);%s-angle (deg);#Delta_{%s} (cm)",axes[ia],axes[ia],axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), nAngleBins, angleRange[ia][0], angleRange[ia][1], nBins, -2.*maxRes, 2.*maxRes);
    fResidualsVsAngle->AddAtAndExpand(h2, kResidualVsAngleClusterOut+ia);
    
    // local chi2
    name = Form("hLocalChi2%sPerCh",axes[ia]); title = Form("local chi2-%s distribution per chamber;chamber ID;local #chi^{2}_{%s}", axes[ia], axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 10, 0.5, 10.5, 1000, 0., 25.);
    fResiduals->AddAtAndExpand(h2, kLocalChi2PerCh+ia);
    name = Form("hLocalChi2%sPerDE",axes[ia]); title = Form("local chi2-%s distribution per DE;DE ID;local #chi^{2}_{%s}", axes[ia], axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), fNDE, 0.5, fNDE+0.5, 1000, 0., 25.);
    for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
    fResiduals->AddAtAndExpand(h2, kLocalChi2PerDE+ia);
    
    // track resolution
    name = Form("hUncorrSlope%sRes",axes[ia]); title = Form("muon slope_{%s} reconstructed resolution at first cluster vs p;p (GeV/c); #sigma_{slope_{%s}}", axes[ia], axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 300, 0., 300., 1000, 0., 0.003);
    fTrackRes->AddAtAndExpand(h2, kUncorrSlopeRes+ia);
    name = Form("hSlope%sRes",axes[ia]); title = Form("muon slope_{%s} reconstructed resolution at vertex vs p;p (GeV/c); #sigma_{slope_{%s}}", axes[ia], axes[ia]);
    h2 = new TH2F(name.Data(), title.Data(), 300, 0., 300., 1000, 0., 0.02);
    fTrackRes->AddAtAndExpand(h2, kSlopeRes+ia);
  }
  
  // local chi2 X+Y
  h2 = new TH2F("hLocalChi2PerCh", "local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) distribution per chamber;chamber ID;local #chi^{2}", 10, 0.5, 10.5, 1000, 0., 25.);
  fResiduals->AddAtAndExpand(h2, kLocalChi2PerCh+2);
  h2 = new TH2F("hLocalChi2PerDE", "local chi2 (~0.5*(#chi^{2}_{X}+#chi^{2}_{Y})) distribution per chamber;DE ID;local #chi^{2}", fNDE, 0.5, fNDE+0.5, 1000, 0., 25.);
  for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
  fResiduals->AddAtAndExpand(h2, kLocalChi2PerDE+2);
  
  // track resolution
  h2 = new TH2F("hUncorrPRes", "muon momentum reconstructed resolution at first cluster vs p;p (GeV/c); #sigma_{p}/p (%)", 300, 0., 300., 1000, 0., 10.);
  fTrackRes->AddAtAndExpand(h2, kUncorrPRes);
  h2 = new TH2F("hPRes", "muon momentum reconstructed resolution at vertex vs p;p (GeV/c); #sigma_{p}/p (%)", 300, 0., 300., 1000, 0., 10.);
  fTrackRes->AddAtAndExpand(h2, kPRes);
  h2 = new TH2F("hUncorrPtRes", "muon transverse momentum reconstructed resolution at first cluster vs p_{t};p_{t} (GeV/c); #sigma_{p_{t}}/p_{t} (%)", 300, 0., 30., 300, 0., 30.);
  fTrackRes->AddAtAndExpand(h2, kUncorrPtRes);
  h2 = new TH2F("hPtRes", "muon transverse momentum reconstructed resolution at vertex vs p_{t};p_{t} (GeV/c); #sigma_{p_{t}}/p_{t} (%)", 300, 0., 30., 300, 0., 30.);
  fTrackRes->AddAtAndExpand(h2, kPtRes);
  
  // Post data at least once per task to ensure data synchronisation (required for merging)
  PostData(1, fResiduals);
  PostData(2, fResidualsVsP);
  PostData(3, fResidualsVsCent);
  PostData(4, fResidualsVsAngle);
  PostData(5, fTrackRes);
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::UserExec(Option_t *)
{
  /// Main event loop
  
  // check if OCDB properly loaded
  if (!fOCDBLoaded) AliFatal("Problem occur while loading OCDB objects");
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!esd) return;
  
  if (fShowProgressBar && (++fNEvents)%100 == 0) cout<<"\rEvent processing... "<<fNEvents<<"\r"<<flush;
  
  // skip events that do not pass required selections
  if (fMuonEventCuts && !fMuonEventCuts->IsSelected(fInputHandler)) return;
  
  // get the centrality
  AliMultSelection *multSelection = static_cast<AliMultSelection*>(esd->FindListObject("MultSelection"));
  Float_t centrality = multSelection ? multSelection->GetMultiplicityPercentile("V0M") : -1.;
  
  // get tracker to refit
  AliMUONVTrackReconstructor* tracker = AliMUONESDInterface::GetTracker();
  
  // loop over tracks
  Int_t nTracks = (Int_t) esd->GetNumberOfMuonTracks(); 
  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {
    
    // get the ESD track
    AliESDMuonTrack* esdTrack = esd->GetMuonTrack(iTrack);
    
    // skip ghost tracks
    if (!esdTrack->ContainTrackerData()) continue;
    
    // apply standard track cuts if any
    if (fMuonTrackCuts && !fMuonTrackCuts->IsSelected(esdTrack)) continue;
    
    // select positive and/or negative tracks
    if (fSign*esdTrack->Charge() < 0) continue;
    
    // select only tracks with MC label or not
    if (fUseMCLabel && esdTrack->GetLabel() < 0) continue;
    
    // skip tracks with not enough clusters
    //if (esdTrack->GetNClusters() < 8) continue;
    
    // get the corresponding MUON track
    AliMUONTrack track;
    AliMUONESDInterface::ESDToMUON(*esdTrack, track, kFALSE);
    
    // change the cluster resolution (and position)
    ModifyClusters(track);
    
    // refit the track
    if (!tracker->RefitTrack(track, fImproveTracks) || (fImproveTracks && !track.IsImproved())) continue;
    
    // save track unchanged
    AliMUONTrack referenceTrack(track);
    
    // get track param at first cluster and add MCS in first chamber
    AliMUONTrackParam trackParamAtFirstCluster(*(static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->First())));
    Int_t firstCh = trackParamAtFirstCluster.GetClusterPtr()->GetChamberId();
    AliMUONTrackExtrap::AddMCSEffect(&trackParamAtFirstCluster, AliMUONConstants::ChamberThicknessInX0(firstCh)/2., -1.);
    
    // fill momentum error at first cluster
    Double_t pXUncorr = trackParamAtFirstCluster.Px();
    Double_t pYUncorr = trackParamAtFirstCluster.Py();
    Double_t pZUncorr = trackParamAtFirstCluster.Pz();
    Double_t pUncorr = trackParamAtFirstCluster.P();
    TMatrixD covUncorr(5,5);
    Cov2CovP(trackParamAtFirstCluster,covUncorr);
    Double_t sigmaPUncorr = TMath::Sqrt(pXUncorr * (pXUncorr*covUncorr(2,2) + pYUncorr*covUncorr(2,3) + pZUncorr*covUncorr(2,4)) +
					pYUncorr * (pXUncorr*covUncorr(3,2) + pYUncorr*covUncorr(3,3) + pZUncorr*covUncorr(3,4)) +
					pZUncorr * (pXUncorr*covUncorr(4,2) + pYUncorr*covUncorr(4,3) + pZUncorr*covUncorr(4,4))) / pUncorr;
    ((TH2F*)fTrackRes->UncheckedAt(kUncorrPRes))->Fill(pUncorr,100.*sigmaPUncorr/pUncorr);
    
    // fill transverse momentum error at first cluster
    Double_t ptUncorr = TMath::Sqrt(pXUncorr*pXUncorr + pYUncorr*pYUncorr);
    Double_t sigmaPtUncorr = TMath::Sqrt(pXUncorr * (pXUncorr*covUncorr(2,2)  + pYUncorr*covUncorr(2,3)) + pYUncorr * (pXUncorr*covUncorr(3,2) + pYUncorr*covUncorr(3,3))) / ptUncorr;
    ((TH2F*)fTrackRes->UncheckedAt(kUncorrPtRes))->Fill(ptUncorr,100.*sigmaPtUncorr/ptUncorr);
    
    // fill slopeX-Y error at first cluster
    ((TH2F*)fTrackRes->UncheckedAt(kUncorrSlopeRes))->Fill(pUncorr,TMath::Sqrt(trackParamAtFirstCluster.GetCovariances()(1,1)));
    ((TH2F*)fTrackRes->UncheckedAt(kUncorrSlopeRes+1))->Fill(pUncorr,TMath::Sqrt(trackParamAtFirstCluster.GetCovariances()(3,3)));
    
    // fill momentum error at vertex
    AliMUONTrackParam trackParamAtVtx(trackParamAtFirstCluster);
    AliMUONTrackExtrap::ExtrapToVertex(&trackParamAtVtx, esdTrack->GetNonBendingCoor(), esdTrack->GetBendingCoor(), esdTrack->GetZ(), 0., 0.);
    Double_t pXVtx = trackParamAtVtx.Px();
    Double_t pYVtx = trackParamAtVtx.Py();
    Double_t pZVtx = trackParamAtVtx.Pz();
    Double_t pVtx = trackParamAtVtx.P();
    TMatrixD covVtx(5,5);
    Cov2CovP(trackParamAtVtx,covVtx);
    Double_t sigmaPVtx = TMath::Sqrt(pXVtx * (pXVtx*covVtx(2,2) + pYVtx*covVtx(2,3) + pZVtx*covVtx(2,4)) +
				     pYVtx * (pXVtx*covVtx(3,2) + pYVtx*covVtx(3,3) + pZVtx*covVtx(3,4)) +
				     pZVtx * (pXVtx*covVtx(4,2) + pYVtx*covVtx(4,3) + pZVtx*covVtx(4,4))) / pVtx;
    ((TH2F*)fTrackRes->UncheckedAt(kPRes))->Fill(pVtx,100.*sigmaPVtx/pVtx);
    
    // fill transverse momentum error at vertex
    Double_t ptVtx = TMath::Sqrt(pXVtx*pXVtx + pYVtx*pYVtx);
    Double_t sigmaPtVtx = TMath::Sqrt(pXVtx * (pXVtx*covVtx(2,2)  + pYVtx*covVtx(2,3)) + pYVtx * (pXVtx*covVtx(3,2) + pYVtx*covVtx(3,3))) / ptVtx;
    ((TH2F*)fTrackRes->UncheckedAt(kPtRes))->Fill(ptVtx,100.*sigmaPtVtx/ptVtx);
    
    // fill slopeX-Y error at vertex
    ((TH2F*)fTrackRes->UncheckedAt(kSlopeRes))->Fill(pVtx,TMath::Sqrt(trackParamAtVtx.GetCovariances()(1,1)));
    ((TH2F*)fTrackRes->UncheckedAt(kSlopeRes+1))->Fill(pVtx,TMath::Sqrt(trackParamAtVtx.GetCovariances()(3,3)));
    
    // skip low momentum tracks in the computation of resolution (cut on initial momentum to avoid correlations)
    if (esdTrack->PUncorrected() < fMinMomentum) continue;
    
    // skip low pT tracks in the computation of resolution (cut on initial pT to avoid correlations)
    Double_t pXUncorr0 = esdTrack->PxUncorrected();
    Double_t pYUncorr0 = esdTrack->PyUncorrected();
    if (TMath::Sqrt(pXUncorr0*pXUncorr0+pYUncorr0*pYUncorr0) < fMinPt) continue;
    
    // loop over clusters
    Int_t nClusters = track.GetNClusters();
    for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
      
      // Get current, previous and next trackParams
      AliMUONTrackParam* trackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster));
      AliMUONTrackParam* previousTrackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->Before(trackParam));
      AliMUONTrackParam* nextTrackParam = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->After(trackParam));
      if (!previousTrackParam && !nextTrackParam) {
	AliError(Form("unable to find a cluster neither before nor after the one at position %d !?!", iCluster));
	track.RecursiveDump();
	break;
      }
      
      // remove mono-cathod clusters on stations 3-4-5 if required
      if (fRemoveMonoCathCl && trackParam->GetClusterPtr()->GetChamberId() > 3) {
	Bool_t hasBending, hasNonBending;
	if (fCheckAllPads) CheckPads(trackParam->GetClusterPtr(), hasBending, hasNonBending);
	else CheckPadsBelow(trackParam->GetClusterPtr(), hasBending, hasNonBending);
	if (!hasBending || !hasNonBending) continue;
      }
      
      // save current trackParam and remove it from the track
      AliMUONTrackParam currentTrackParam(*trackParam);
      track.RemoveTrackParamAtCluster(trackParam);
      
      // get cluster info
      AliMUONVCluster* cluster = currentTrackParam.GetClusterPtr();
      Int_t chId = cluster->GetChamberId();
      Int_t halfChId = (cluster->GetX() > 0) ? 2*chId : 2*chId+1;
      Int_t deId = cluster->GetDetElemId();
      
      // compute residuals with cluster still attached to the track
      AliMUONTrackParam* referenceTrackParam = static_cast<AliMUONTrackParam*>(referenceTrack.GetTrackParamAtCluster()->UncheckedAt(iCluster));
      Double_t deltaX = cluster->GetX() - referenceTrackParam->GetNonBendingCoor();
      Double_t deltaY = cluster->GetY() - referenceTrackParam->GetBendingCoor();
      
      // compute local track angles
      Double_t angleX = TMath::ATan(referenceTrackParam->GetNonBendingSlope())*TMath::RadToDeg();
      Double_t angleY = TMath::ATan(referenceTrackParam->GetBendingSlope())*TMath::RadToDeg();
      
      // compute local chi2
      Double_t sigmaDeltaX2 = cluster->GetErrX2() - referenceTrackParam->GetCovariances()(0,0);
      Double_t sigmaDeltaY2 = cluster->GetErrY2() - referenceTrackParam->GetCovariances()(2,2);
      Double_t localChi2X = (sigmaDeltaX2 > 0.) ? deltaX*deltaX/sigmaDeltaX2 : 0.;
      Double_t localChi2Y = (sigmaDeltaY2 > 0.) ? deltaY*deltaY/sigmaDeltaY2 : 0.;
      Double_t localChi2 = 0.5 * referenceTrackParam->GetLocalChi2();
      
      // fill local chi2 info at every clusters
      ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerCh))->Fill(chId+1, localChi2X);
      ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerCh+1))->Fill(chId+1, localChi2Y);
      ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerCh+2))->Fill(chId+1, localChi2);
      ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerDE))->Fill(fDEIndices[deId], localChi2X);
      ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerDE+1))->Fill(fDEIndices[deId], localChi2Y);
      ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerDE+2))->Fill(fDEIndices[deId], localChi2);
      
      // make sure the track has another cluster in the same station and can still be refitted
      Bool_t refit = track.IsValid( 1 << (chId/2) );
      if (refit) {
	
	// refit the track and proceed if everything goes fine
	if (tracker->RefitTrack(track, kFALSE)) {
	  
	  // fill histograms of residuals with cluster still attached to the track
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerChClusterIn))->Fill(chId+1, deltaX);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerChClusterIn+1))->Fill(chId+1, deltaY);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterIn))->Fill(halfChId+1, deltaX);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterIn+1))->Fill(halfChId+1, deltaY);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterIn))->Fill(fDEIndices[deId], deltaX);
	  ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterIn+1))->Fill(fDEIndices[deId], deltaY);
	  ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterIn+chId))->Fill(pUncorr, deltaX);
	  ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterIn+10+chId))->Fill(pUncorr, deltaY);
	  ((TH2F*)fResidualsVsP->UncheckedAt(kResidualVsPClusterIn))->Fill(pUncorr, deltaX);
	  ((TH2F*)fResidualsVsP->UncheckedAt(kResidualVsPClusterIn+1))->Fill(pUncorr, deltaY);
	  ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterIn+chId))->Fill(centrality, deltaX);
	  ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterIn+10+chId))->Fill(centrality, deltaY);
	  ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualVsCentClusterIn))->Fill(centrality, deltaX);
	  ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualVsCentClusterIn+1))->Fill(centrality, deltaY);
	  ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterIn+chId))->Fill(angleX, deltaX);
	  ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterIn+10+chId))->Fill(angleY, deltaY);
	  ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualVsAngleClusterIn))->Fill(angleX, deltaX);
	  ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualVsAngleClusterIn+1))->Fill(angleY, deltaY);
	  ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInHalfChVsAngleClusterIn+halfChId))->Fill(angleX, deltaX);
	  ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInHalfChVsAngleClusterIn+20+halfChId))->Fill(angleY, deltaY);
	  
	  // find the track parameters closest to the current cluster position
	  Double_t dZWithPrevious = (previousTrackParam) ? TMath::Abs(previousTrackParam->GetClusterPtr()->GetZ() - cluster->GetZ()) : FLT_MAX;
	  Int_t previousChId = (previousTrackParam) ? previousTrackParam->GetClusterPtr()->GetChamberId() : -1;
	  Double_t dZWithNext = (nextTrackParam) ? TMath::Abs(nextTrackParam->GetClusterPtr()->GetZ() - cluster->GetZ()) : FLT_MAX;
	  AliMUONTrackParam* startingTrackParam = (nextTrackParam) ? nextTrackParam : previousTrackParam;
	  if ((fExtrapMode == 0 && previousTrackParam && dZWithPrevious < dZWithNext) ||
	      (fExtrapMode == 1 && previousTrackParam && !(chId/2 == 2 && previousChId/2 == 1) &&
	       !(chId/2 == 3 && previousChId/2 == 2))) startingTrackParam = previousTrackParam;
	  
	  // reset current parameters
	  currentTrackParam.SetParameters(startingTrackParam->GetParameters());
	  currentTrackParam.SetZ(startingTrackParam->GetZ());
	  currentTrackParam.SetCovariances(startingTrackParam->GetCovariances());
	  currentTrackParam.ResetPropagator();
	  
	  // extrapolate to the current cluster position and fill histograms of residuals if everything goes fine
	  if (AliMUONTrackExtrap::ExtrapToZCov(&currentTrackParam, currentTrackParam.GetClusterPtr()->GetZ(), kTRUE)) {
	    
	    // compute MCS dispersion on the first chamber
	    TMatrixD mcsCov(5,5);
	    if (startingTrackParam == nextTrackParam && chId == 0) {
	      AliMUONTrackParam trackParamForMCS;
	      trackParamForMCS.SetParameters(nextTrackParam->GetParameters());
	      AliMUONTrackExtrap::AddMCSEffect(&trackParamForMCS,AliMUONConstants::ChamberThicknessInX0(nextTrackParam->GetClusterPtr()->GetChamberId()),-1.);
	      const TMatrixD &propagator = currentTrackParam.GetPropagator();
	      TMatrixD tmp(trackParamForMCS.GetCovariances(),TMatrixD::kMultTranspose,propagator);
	      mcsCov.Mult(propagator,tmp);
	    } else mcsCov.Zero();
	    
	    // compute residuals
	    Double_t trackResX2 = currentTrackParam.GetCovariances()(0,0) + mcsCov(0,0);
	    Double_t trackResY2 = currentTrackParam.GetCovariances()(2,2) + mcsCov(2,2);
	    deltaX = cluster->GetX() - currentTrackParam.GetNonBendingCoor();
	    deltaY = cluster->GetY() - currentTrackParam.GetBendingCoor();
	    
	    // fill histograms with cluster not attached to the track
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerChClusterOut))->Fill(chId+1, deltaX);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerChClusterOut+1))->Fill(chId+1, deltaY);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterOut))->Fill(halfChId+1, deltaX);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterOut+1))->Fill(halfChId+1, deltaY);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterOut))->Fill(fDEIndices[deId], deltaX);
	    ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterOut+1))->Fill(fDEIndices[deId], deltaY);
	    ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterOut+chId))->Fill(pUncorr, deltaX);
	    ((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterOut+10+chId))->Fill(pUncorr, deltaY);
	    ((TH2F*)fResidualsVsP->UncheckedAt(kResidualVsPClusterOut))->Fill(pUncorr, deltaX);
	    ((TH2F*)fResidualsVsP->UncheckedAt(kResidualVsPClusterOut+1))->Fill(pUncorr, deltaY);
	    ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterOut+chId))->Fill(centrality, deltaX);
	    ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterOut+10+chId))->Fill(centrality, deltaY);
	    ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualVsCentClusterOut))->Fill(centrality, deltaX);
	    ((TH2F*)fResidualsVsCent->UncheckedAt(kResidualVsCentClusterOut+1))->Fill(centrality, deltaY);
	    ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterOut+chId))->Fill(angleX, deltaX);
	    ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterOut+10+chId))->Fill(angleY, deltaY);
	    ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualVsAngleClusterOut))->Fill(angleX, deltaX);
	    ((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualVsAngleClusterOut+1))->Fill(angleY, deltaY);
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerCh))->Fill(chId+1, TMath::Sqrt(trackResX2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerCh+1))->Fill(chId+1, TMath::Sqrt(trackResY2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerHalfCh))->Fill(halfChId+1, TMath::Sqrt(trackResX2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerHalfCh+1))->Fill(halfChId+1, TMath::Sqrt(trackResY2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerDE))->Fill(fDEIndices[deId], TMath::Sqrt(trackResX2));
	    ((TH2F*)fResiduals->UncheckedAt(kTrackResPerDE+1))->Fill(fDEIndices[deId], TMath::Sqrt(trackResY2));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerCh))->Fill(chId+1, TMath::Sqrt(mcsCov(0,0)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerCh+1))->Fill(chId+1, TMath::Sqrt(mcsCov(2,2)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerHalfCh))->Fill(halfChId+1, TMath::Sqrt(mcsCov(0,0)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerHalfCh+1))->Fill(halfChId+1, TMath::Sqrt(mcsCov(2,2)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerDE))->Fill(fDEIndices[deId], TMath::Sqrt(mcsCov(0,0)));
	    ((TH2F*)fResiduals->UncheckedAt(kMCSPerDE+1))->Fill(fDEIndices[deId], TMath::Sqrt(mcsCov(2,2)));
	    ((TH2F*)fResiduals->UncheckedAt(kClusterRes2PerCh))->Fill(chId+1, deltaX*deltaX - trackResX2);
	    ((TH2F*)fResiduals->UncheckedAt(kClusterRes2PerCh+1))->Fill(chId+1, deltaY*deltaY - trackResY2);
	  }
	  
	}
	
      }
      
      // restore the track
      track.AddTrackParamAtCluster(currentTrackParam, *(currentTrackParam.GetClusterPtr()), kTRUE);
      
    }
    
  }
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fResiduals);
  PostData(2, fResidualsVsP);
  PostData(3, fResidualsVsCent);
  PostData(4, fResidualsVsAngle);
  PostData(5, fTrackRes);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMuonResolution::UserNotify()
{
  /// setup OCDB default storage
  
  // do it only once
  if (fOCDBLoaded) return kTRUE;
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage(fDefaultStorage.Data());
  
  return kTRUE;
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::NotifyRun()
{
  /// load necessary data from OCDB corresponding to the first run number and initialize analysis
  
  if (fOCDBLoaded) return;
  
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetRun(fCurrentRunNumber);
  
  if (!AliMUONCDB::LoadField()) return;
  
  if (!AliMUONCDB::LoadMapping()) return;
  
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) return;
  
  if (fImproveTracks) recoParam->ImproveTracks(kTRUE);
  else recoParam->ImproveTracks(kFALSE);
  
  AliMUONESDInterface::ResetTracker(recoParam);
  
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    
    // set the cluster resolution to default if not already set and create output objects
    if (fClusterResNB[i] < 0.) fClusterResNB[i] = recoParam->GetDefaultNonBendingReso(i);
    if (fClusterResB[i] < 0.) fClusterResB[i] = recoParam->GetDefaultBendingReso(i);
    
    // fill correspondence between DEId and indices for histo (starting from 1)
    AliMpDEIterator it;
    it.First(i);
    while (!it.IsDone()) {
      fNDE++;
      fDEIndices[it.CurrentDEId()] = fNDE;
      fDEIds[fNDE] = it.CurrentDEId();
      it.Next();
    }
    
  }
  
  // recover default storage full name (raw:// cannot be used to set specific storage)
  TString defaultStorage(cdbm->GetDefaultStorage()->GetType());
  if (defaultStorage == "alien") defaultStorage += Form("://folder=%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
  else defaultStorage += Form("://%s", cdbm->GetDefaultStorage()->GetBaseFolder().Data());
  
  if (fReAlign) {
    
    // reset existing geometry/alignment if any
    if (cdbm->GetEntryCache()->Contains("GRP/Geometry/Data")) cdbm->UnloadFromCache("GRP/Geometry/Data");
    if (cdbm->GetEntryCache()->Contains("MUON/Align/Data")) cdbm->UnloadFromCache("MUON/Align/Data");
    if (AliGeomManager::GetGeometry()) AliGeomManager::GetGeometry()->UnlockGeometry();
    
    // get original geometry transformer
    AliGeomManager::LoadGeometry();
    if (!AliGeomManager::GetGeometry()) return;
    if (fOldAlignStorage != "none") {
      cdbm->SetSpecificStorage("MUON/Align/Data", fOldAlignStorage.IsNull() ? defaultStorage.Data() : fOldAlignStorage.Data(), fOldAlignVersion, fOldAlignSubVersion);
      AliGeomManager::ApplyAlignObjsFromCDB("MUON");
    }
    fOldGeoTransformer = new AliMUONGeometryTransformer();
    fOldGeoTransformer->LoadGeometryData();
    
  }
  
  // reset existing geometry/alignment if any
  if (cdbm->GetEntryCache()->Contains("GRP/Geometry/Data")) cdbm->UnloadFromCache("GRP/Geometry/Data");
  if (cdbm->GetEntryCache()->Contains("MUON/Align/Data")) cdbm->UnloadFromCache("MUON/Align/Data");
  if (AliGeomManager::GetGeometry()) AliGeomManager::GetGeometry()->UnlockGeometry();
  
  // load original or realigned geometry (also used to check pads below clusters and for track extrapolation to vertex)
  AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry()) return;
  if (fNewAlignStorage != "none") {
    cdbm->SetSpecificStorage("MUON/Align/Data", fNewAlignStorage.IsNull() ? defaultStorage.Data() : fNewAlignStorage.Data(), fNewAlignVersion, fNewAlignSubVersion);
    AliGeomManager::ApplyAlignObjsFromCDB("MUON");
  }
  fNewGeoTransformer = new AliMUONGeometryTransformer();
  fNewGeoTransformer->LoadGeometryData();
  
  // print starting chamber resolution if required
  if (fPrintClResPerCh) {
    printf("\nstarting chamber resolution:\n");
    printf(" - non-bending:");
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %5.3f":", %5.3f",fClusterResNB[i]);
    printf("\n -     bending:");
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %6.4f":", %6.4f",fClusterResB[i]);
    printf("\n\n");
  }
  
  if (fMuonEventCuts) fMuonEventCuts->Print();
  
  // get the trackCuts for this run
  if (fMuonTrackCuts) {
    fMuonTrackCuts->SetRun(fInputHandler);
    fMuonTrackCuts->Print();
  }
  
  fOCDBLoaded = kTRUE;
  
  UserCreateOutputObjects();
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::Terminate(Option_t *)
{
  /// compute final results
  
  // recover output objects
  fResiduals = static_cast<TObjArray*> (GetOutputData(1));
  fResidualsVsP = static_cast<TObjArray*> (GetOutputData(2));
  fResidualsVsCent = static_cast<TObjArray*> (GetOutputData(3));
  fResidualsVsAngle = static_cast<TObjArray*> (GetOutputData(4));
  fTrackRes = static_cast<TObjArray*> (GetOutputData(5));
  if (!fResiduals || !fResidualsVsP || !fResidualsVsCent || !fResidualsVsAngle || !fTrackRes) return;
  
  // summary graphs
  fLocalChi2 = new TObjArray(1000);
  fLocalChi2->SetOwner();
  fChamberRes = new TObjArray(1000);
  fChamberRes->SetOwner();
  fCanvases = new TObjArray(1000);
  fCanvases->SetOwner();
  fTmpHists = new TObjArray(1000);
  fTmpHists->SetOwner();
  TGraphErrors* g;
  TMultiGraph* mg;
  TCanvas* c;
  
  const char* axes[2] = {"X", "Y"};
  const char* side[2] = {"I", "O"};
  Double_t newClusterRes[2][10], newClusterResErr[2][10];
  fNDE = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterIn))->GetXaxis()->GetNbins();
  TString name, title;
  
  for (Int_t ia = 0; ia < 2; ia++) {
    
    // shifts per chamber
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChMean_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: mean (cluster in);chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerChMeanClusterIn+ia);
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChMean_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: mean (cluster out);chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerChMeanClusterOut+ia);
    
    // shifts per half-chamber
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerHalfChMean_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per half Ch: mean (cluster in);half chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerHalfChMeanClusterIn+ia);
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerHalfChMean_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per half Ch: mean (cluster out);half chamber ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerHalfChMeanClusterOut+ia);
    
    // shifts per DE
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gResidual%sPerDEMean_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per DE: mean (cluster in);DE ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerDEMeanClusterIn+ia);
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gResidual%sPerDEMean_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per DE: mean (cluster out);DE ID;<#Delta_{%s}> (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerDEMeanClusterOut+ia);
    
    // resolution per chamber
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChSigma_ClusterIn",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: sigma (cluster in);chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerChSigmaClusterIn+ia);
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChSigma_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: sigma (cluster out);chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerChSigmaClusterOut+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gResidual%sPerChDispersion_ClusterOut",axes[ia]));
    g->SetTitle(Form("cluster-track residual-%s per Ch: dispersion (cluster out);chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kResidualPerChDispersionClusterOut+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gCombinedResidual%sPerChSigma",axes[ia]));
    g->SetTitle(Form("combined cluster-track residual-%s per Ch: sigma;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kCombinedResidualPerChSigma+ia);
    
    // resolution per half-chamber
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gCombinedResidual%sPerHalfChSigma",axes[ia]));
    g->SetTitle(Form("combined cluster-track residual-%s per half Ch: sigma;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kCombinedResidualPerHalfChSigma+ia);
    
    // resolution per DE
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gCombinedResidual%sPerDESigma",axes[ia]));
    g->SetTitle(Form("combined cluster-track residual-%s per DE: sigma;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kCombinedResidualPerDESigma+ia);
    
    // shifts versus p
    name = Form("mgResidual%sMeanVsP_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s per chamber versus momentum: mean (cluster in);p (GeV/c^{2});<#Delta_{%s}> (cm)",axes[ia],axes[ia]);
    mg = new TMultiGraph(name.Data(), title.Data());
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      g = new TGraphErrors(((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterIn+10*ia+i))->GetNbinsX());
      g->SetName(Form("gShift%sVsP_ch%d",axes[ia],i+1));
      g->SetMarkerStyle(kFullDotMedium);
      g->SetMarkerColor(i+1+i/9);
      mg->Add(g,"p");
    }
    fChamberRes->AddAtAndExpand(mg, kResidualMeanClusterInVsP+ia);
    
    // resolutions versus p
    name = Form("mgCombinedResidual%sSigmaVsP",axes[ia]); title = Form("cluster %s-resolution per chamber versus momentum;p (GeV/c^{2});#sigma_{%s} (cm)",axes[ia],axes[ia]);
    mg = new TMultiGraph(name.Data(), title.Data());
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      g = new TGraphErrors(((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterIn+10*ia+i))->GetNbinsX());
      g->SetName(Form("gRes%sVsP_ch%d",axes[ia],i+1));
      g->SetMarkerStyle(kFullDotMedium);
      g->SetMarkerColor(i+1+i/9);
      mg->Add(g,"p");
    }
    fChamberRes->AddAtAndExpand(mg, kCombinedResidualSigmaVsP+ia);
    
    g = new TGraphErrors(((TH2F*)fResidualsVsP->UncheckedAt(kResidualVsPClusterIn+ia))->GetNbinsX());
    g->SetName(Form("gCombinedResidual%sSigmaVsP",axes[ia]));
    g->SetTitle(Form("cluster %s-resolution integrated over chambers versus momentum: sigma;p (GeV/c^{2});#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kCombinedResidualAllChSigmaVsP+ia);
    
    // shifts versus centrality
    name = Form("mgResidual%sMeanVsCent_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s per chamber versus centrality: mean (cluster in);centrality (%%);<#Delta_{%s}> (cm)",axes[ia],axes[ia]);
    mg = new TMultiGraph(name.Data(), title.Data());
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      g = new TGraphErrors(((TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterIn+10*ia+i))->GetNbinsX());
      g->SetName(Form("gShift%sVsCent_ch%d",axes[ia],i+1));
      g->SetMarkerStyle(kFullDotMedium);
      g->SetMarkerColor(i+1+i/9);
      mg->Add(g,"p");
    }
    fChamberRes->AddAtAndExpand(mg, kResidualMeanClusterInVsCent+ia);
    
    // resolutions versus centrality
    name = Form("mgCombinedResidual%sSigmaVsCent",axes[ia]); title = Form("cluster %s-resolution per chamber versus centrality;centrality (%%);#sigma_{%s} (cm)",axes[ia],axes[ia]);
    mg = new TMultiGraph(name.Data(), title.Data());
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      g = new TGraphErrors(((TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterIn+10*ia+i))->GetNbinsX());
      g->SetName(Form("gRes%sVsCent_ch%d",axes[ia],i+1));
      g->SetMarkerStyle(kFullDotMedium);
      g->SetMarkerColor(i+1+i/9);
      mg->Add(g,"p");
    }
    fChamberRes->AddAtAndExpand(mg, kCombinedResidualSigmaVsCent+ia);
    
    g = new TGraphErrors(((TH2F*)fResidualsVsCent->UncheckedAt(kResidualVsCentClusterIn+ia))->GetNbinsX());
    g->SetName(Form("gCombinedResidual%sSigmaVsCent",axes[ia]));
    g->SetTitle(Form("cluster %s-resolution integrated over chambers versus centrality: sigma;centrality (%%);#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kCombinedResidualAllChSigmaVsCent+ia);
    
    // shifts versus track angle
    name = Form("mgResidual%sMeanVsAngle_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s per chamber versus track angle: mean (cluster in);%s-angle (deg);<#Delta_{%s}> (cm)",axes[ia],axes[ia],axes[ia]);
    mg = new TMultiGraph(name.Data(), title.Data());
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      g = new TGraphErrors(((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterIn+10*ia+i))->GetNbinsX());
      g->SetName(Form("gShift%sVsAngle_ch%d",axes[ia],i+1));
      g->SetMarkerStyle(kFullDotMedium);
      g->SetMarkerColor(i+1+i/9);
      mg->Add(g,"p");
    }
    fChamberRes->AddAtAndExpand(mg, kResidualMeanClusterInVsAngle+ia);
    
    name = Form("mgHChResidual%sMeanVsAngle_ClusterIn",axes[ia]); title = Form("cluster-track residual-%s per half-chamber versus track angle: mean (cluster in);%s-angle (deg);<#Delta_{%s}> (cm)",axes[ia],axes[ia],axes[ia]);
    mg = new TMultiGraph(name.Data(), title.Data());
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      for (Int_t j = 0; j < 2; j++) {
	g = new TGraphErrors(((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInHalfChVsAngleClusterIn+20*ia+2*i+j))->GetNbinsX());
	g->SetName(Form("gShift%sVsAngle_halfCh%d%s",axes[ia],i+1,side[j]));
	g->SetMarkerStyle(kFullDotMedium);
	g->SetMarkerColor(2*i+j+1+(2*i+j)/9);
	mg->Add(g,"p");
      }
    }
    fChamberRes->AddAtAndExpand(mg, kHChResidualMeanClusterInVsAngle+ia);
    
    // resolutions versus track angle
    name = Form("mgCombinedResidual%sSigmaVsAngle",axes[ia]); title = Form("cluster %s-resolution per chamber versus track angle;%s-angle (deg);#sigma_{%s} (cm)",axes[ia],axes[ia],axes[ia]);
    mg = new TMultiGraph(name.Data(), title.Data());
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      g = new TGraphErrors(((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterIn+10*ia+i))->GetNbinsX());
      g->SetName(Form("gRes%sVsAngle_ch%d",axes[ia],i+1));
      g->SetMarkerStyle(kFullDotMedium);
      g->SetMarkerColor(i+1+i/9);
      mg->Add(g,"p");
    }
    fChamberRes->AddAtAndExpand(mg, kCombinedResidualSigmaVsAngle+ia);
    
    g = new TGraphErrors(((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualVsAngleClusterIn+ia))->GetNbinsX());
    g->SetName(Form("gCombinedResidual%sSigmaVsAngle",axes[ia]));
    g->SetTitle(Form("cluster %s-resolution integrated over chambers versus track angle: sigma;%s-angle (deg);#sigma_{%s} (cm)",axes[ia],axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kCombinedResidualAllChSigmaVsAngle+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gTrackRes%sPerCh",axes[ia]));
    g->SetTitle(Form("track <#sigma_{%s}> per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kTrackResPerChMean+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gMCS%sPerCh",axes[ia]));
    g->SetTitle(Form("MCS %s-dispersion of extrapolated track per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kMCSPerChMean+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gClusterRes%sPerCh",axes[ia]));
    g->SetTitle(Form("cluster <#sigma_{%s}> per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kClusterResPerCh+ia);
    
    g = new TGraphErrors(2*AliMUONConstants::NTrackingCh());
    g->SetName(Form("gClusterRes%sPerHalfCh",axes[ia]));
    g->SetTitle(Form("cluster <#sigma_{%s}> per half Ch;half chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kClusterResPerHalfCh+ia);
    
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gClusterRes%sPerDE",axes[ia]));
    g->SetTitle(Form("cluster <#sigma_{%s}> per DE;DE ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kClusterResPerDE+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gCalcClusterRes%sPerCh",axes[ia]));
    g->SetTitle(Form("calculated cluster <#sigma_{%s}> per Ch;chamber ID;#sigma_{%s} (cm)",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fChamberRes->AddAtAndExpand(g, kCalcClusterResPerCh+ia);
    
    g = new TGraphErrors(AliMUONConstants::NTrackingCh());
    g->SetName(Form("gLocalChi2%sPerChMean",axes[ia]));
    g->SetTitle(Form("local chi2-%s per Ch: mean;chamber ID;<local #chi^{2}_{%s}>",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fLocalChi2->AddAtAndExpand(g, kLocalChi2PerChMean+ia);
    
    g = new TGraphErrors(fNDE);
    g->SetName(Form("gLocalChi2%sPerDEMean",axes[ia]));
    g->SetTitle(Form("local chi2-%s per DE: mean;DE ID;<local #chi^{2}_{%s}>",axes[ia],axes[ia]));
    g->SetMarkerStyle(kFullDotLarge);
    fLocalChi2->AddAtAndExpand(g, kLocalChi2PerDEMean+ia);
    
    // canvases
    name = Form("cDetailRes%sPerChClIn",axes[ia]); title = Form("cDetailRes%sPerChClIn",axes[ia]);
    c = new TCanvas(name.Data(), title.Data(), 1200, 500);
    c->Divide(5,2);
    fCanvases->AddAtAndExpand(c, kDetailResPerCh+2*ia);
    
    name = Form("cDetailRes%sPerChClOut",axes[ia]); title = Form("cDetailRes%sPerChClOut",axes[ia]);
    c = new TCanvas(name.Data(), title.Data(), 1200, 500);
    c->Divide(5,2);
    fCanvases->AddAtAndExpand(c, kDetailResPerCh+1+2*ia);
    
    name = Form("cDetailRes%sPerHalfChClIn",axes[ia]); title = Form("cDetailRes%sPerHalfChClIn",axes[ia]);
    c = new TCanvas(name.Data(), title.Data(), 1200, 800);
    c->Divide(5,4);
    fCanvases->AddAtAndExpand(c, kDetailResPerHalfCh+2*ia);
    
    name = Form("cDetailRes%sPerHalfChClOut",axes[ia]); title = Form("cDetailRes%sPerHalfChClOut",axes[ia]);
    c = new TCanvas(name.Data(), title.Data(), 1200, 800);
    c->Divide(5,4);
    fCanvases->AddAtAndExpand(c, kDetailResPerHalfCh+1+2*ia);
    
    name = Form("cDetailRes%sPerDEClIn",axes[ia]); title = Form("cDetailRes%sPerDEClIn",axes[ia]);
    c = new TCanvas(name.Data(), title.Data(), 1200, 800);
    c->Divide(13,12);
    fCanvases->AddAtAndExpand(c, kDetailResPerDE+2*ia);
    
    name = Form("cDetailRes%sPerDEClOut",axes[ia]); title = Form("cDetailRes%sPerDEClOut",axes[ia]);
    c = new TCanvas(name.Data(), title.Data(), 1200, 800);
    c->Divide(13,12);
    fCanvases->AddAtAndExpand(c, kDetailResPerDE+1+2*ia);
    
    // compute residual mean and dispersion integrated over chambers versus p
    FillMeanSigmaClusterVsX((TH2F*)fResidualsVsP->UncheckedAt(kResidualVsPClusterIn+ia),
			    (TH2F*)fResidualsVsP->UncheckedAt(kResidualVsPClusterOut+ia),
			    0x0, (TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualAllChSigmaVsP+ia));
    
    // compute residual mean and dispersion integrated over chambers versus centrality
    FillMeanSigmaClusterVsX((TH2F*)fResidualsVsCent->UncheckedAt(kResidualVsCentClusterIn+ia),
			    (TH2F*)fResidualsVsCent->UncheckedAt(kResidualVsCentClusterOut+ia),
			    0x0, (TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualAllChSigmaVsCent+ia));
    
    // compute residual mean and dispersion integrated over chambers versus track angle
    FillMeanSigmaClusterVsX((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualVsAngleClusterIn+ia),
			    (TH2F*)fResidualsVsAngle->UncheckedAt(kResidualVsAngleClusterOut+ia),
			    0x0, (TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualAllChSigmaVsAngle+ia));
    
    // compute residual mean and dispersion and averaged local chi2 per chamber and half chamber
    Double_t meanIn, meanInErr, meanOut, meanOutErr, sigma, sigmaIn, sigmaInErr, sigmaOut, sigmaOutErr;
    Double_t sigmaTrack, sigmaTrackErr, sigmaMCS, sigmaMCSErr, clusterRes, clusterResErr, sigmaCluster, sigmaClusterErr;
    Double_t dumy1, dumy2;
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
      
      // method 1
      TH1D *tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerChClusterIn+ia))->ProjectionY(Form("hRes%sCh%dClIn",axes[ia],i+1),i+1,i+1,"e");
      tmp->SetTitle(Form("chamber %d",i+1));
      GetMeanRMS(tmp, meanIn, meanInErr, sigmaIn, sigmaInErr, (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChMeanClusterIn+ia),
		 (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChSigmaClusterIn+ia), i, i+1);
      fTmpHists->AddLast(tmp);
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerChClusterOut+ia))->ProjectionY(Form("hRes%sCh%dClOut",axes[ia],i+1),i+1,i+1,"e");
      tmp->SetTitle(Form("chamber %d",i+1));
      GetMeanRMS(tmp, meanOut, meanOutErr, sigmaOut, sigmaOutErr, (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChMeanClusterOut+ia),
		 (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChSigmaClusterOut+ia), i, i+1);
      fTmpHists->AddLast(tmp);
      
      if (fCorrectForSystematics) {
	sigma = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
	sigmaInErr = (sigma>0) ? TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigma : 0.;
	sigmaIn = sigma;
	sigma = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
	sigmaOutErr = (sigma>0) ? TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigma : 0.;
	sigmaOut = sigma;
      }
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChDispersionClusterOut+ia))->SetPoint(i, i+1, sigmaOut);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChDispersionClusterOut+ia))->SetPointError(i, 0., sigmaOutErr);
      
      clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
      //      clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
      clusterResErr = TMath::Sqrt(sigmaInErr*sigmaOutErr);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerChSigma+ia))->SetPoint(i, i+1, clusterRes);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerChSigma+ia))->SetPointError(i, 0., clusterResErr);
      newClusterRes[ia][i] = clusterRes;
      newClusterResErr[ia][i] = clusterResErr;
      
      // method 2
      tmp = ((TH2F*)fResiduals->UncheckedAt(kTrackResPerCh+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMeanRMS(tmp, sigmaTrack, sigmaTrackErr, dumy1, dumy2, (TGraphErrors*)fChamberRes->UncheckedAt(kTrackResPerChMean+ia), 0x0, i, i+1, kFALSE, kFALSE);
      delete tmp;
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kMCSPerCh+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMeanRMS(tmp, sigmaMCS, sigmaMCSErr, dumy1, dumy2, (TGraphErrors*)fChamberRes->UncheckedAt(kMCSPerChMean+ia), 0x0, i, i+1, kFALSE, kFALSE);
      delete tmp;
      
      sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
      if (sigmaCluster > 0.) {
	sigmaCluster = TMath::Sqrt(sigmaCluster);
	sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
      } else {
	sigmaCluster = 0.;
	sigmaClusterErr = 0.;
      }
      ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerCh+ia))->SetPoint(i, i+1, sigmaCluster);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerCh+ia))->SetPointError(i, 0., sigmaClusterErr);
      
      // method 3
      tmp = ((TH2F*)fResiduals->UncheckedAt(kClusterRes2PerCh+ia))->ProjectionY("tmp",i+1,i+1,"e");
      ZoomRight(tmp);
      clusterRes = tmp->GetMean();
      if (clusterRes > 0.) {
	((TGraphErrors*)fChamberRes->UncheckedAt(kCalcClusterResPerCh+ia))->SetPoint(i, i+1, TMath::Sqrt(clusterRes));
	((TGraphErrors*)fChamberRes->UncheckedAt(kCalcClusterResPerCh+ia))->SetPointError(i, 0., 0.5 * tmp->GetMeanError() / TMath::Sqrt(clusterRes));
      } else {
	((TGraphErrors*)fChamberRes->UncheckedAt(kCalcClusterResPerCh+ia))->SetPoint(i, i+1, 0.);
	((TGraphErrors*)fChamberRes->UncheckedAt(kCalcClusterResPerCh+ia))->SetPointError(i, 0., 0.);
      }
      delete tmp;
      
      // method 1 versus p
      FillMeanSigmaClusterVsX((TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterIn+10*ia+i),
			      (TH2F*)fResidualsVsP->UncheckedAt(kResidualInChVsPClusterOut+10*ia+i),
			      (TGraphErrors*)((TMultiGraph*)fChamberRes->UncheckedAt(kResidualMeanClusterInVsP+ia))->GetListOfGraphs()->FindObject(Form("gShift%sVsP_ch%d",axes[ia],i+1)),
			      (TGraphErrors*)((TMultiGraph*)fChamberRes->UncheckedAt(kCombinedResidualSigmaVsP+ia))->GetListOfGraphs()->FindObject(Form("gRes%sVsP_ch%d",axes[ia],i+1)));
      
      // method 1 versus centrality
      FillMeanSigmaClusterVsX((TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterIn+10*ia+i),
			      (TH2F*)fResidualsVsCent->UncheckedAt(kResidualInChVsCentClusterOut+10*ia+i),
			      (TGraphErrors*)((TMultiGraph*)fChamberRes->UncheckedAt(kResidualMeanClusterInVsCent+ia))->GetListOfGraphs()->FindObject(Form("gShift%sVsCent_ch%d",axes[ia],i+1)),
			      (TGraphErrors*)((TMultiGraph*)fChamberRes->UncheckedAt(kCombinedResidualSigmaVsCent+ia))->GetListOfGraphs()->FindObject(Form("gRes%sVsCent_ch%d",axes[ia],i+1)));
      
      // method 1 versus track angle
      FillMeanSigmaClusterVsX((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterIn+10*ia+i),
			      (TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInChVsAngleClusterOut+10*ia+i),
			      (TGraphErrors*)((TMultiGraph*)fChamberRes->UncheckedAt(kResidualMeanClusterInVsAngle+ia))->GetListOfGraphs()->FindObject(Form("gShift%sVsAngle_ch%d",axes[ia],i+1)),
			      (TGraphErrors*)((TMultiGraph*)fChamberRes->UncheckedAt(kCombinedResidualSigmaVsAngle+ia))->GetListOfGraphs()->FindObject(Form("gRes%sVsAngle_ch%d",axes[ia],i+1)));
      
      // compute residual mean and dispersion per half chamber
      for (Int_t j = 0; j < 2; j++) {
	Int_t k = 2*i+j;
	
	// method 1
	tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterIn+ia))->ProjectionY(Form("hRes%sHalfCh%dClIn",axes[ia],k+1),k+1,k+1,"e");
	tmp->SetTitle(Form("half chamber %s",((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterIn+ia))->GetXaxis()->GetBinLabel(k+1)));
	GetMeanRMS(tmp, meanIn, meanInErr, sigmaIn, sigmaInErr, (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterIn+ia), 0x0, k, k+1);
	fTmpHists->AddLast(tmp);
	
	tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterOut+ia))->ProjectionY(Form("hRes%sHalfCh%dClOut",axes[ia],k+1),k+1,k+1,"e");
	tmp->SetTitle(Form("half chamber %s",((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterOut+ia))->GetXaxis()->GetBinLabel(k+1)));
	GetMeanRMS(tmp, meanOut, meanOutErr, sigmaOut, sigmaOutErr, (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterOut+ia), 0x0, k, k+1);
	fTmpHists->AddLast(tmp);
	
	if (fCorrectForSystematics) {
	  sigma = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
	  sigmaInErr = (sigma>0) ? TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigma : 0.;
	  sigmaIn = sigma;
	  sigma = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
	  sigmaOutErr = (sigma>0) ? TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigma : 0.;
	  sigmaOut = sigma;
	}
	
	clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
	//	clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
	clusterResErr = TMath::Sqrt(sigmaInErr*sigmaOutErr);
	((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerHalfChSigma+ia))->SetPoint(k, k+1, clusterRes);
	((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerHalfChSigma+ia))->SetPointError(k, 0., clusterResErr);
	
	// method 2
	tmp = ((TH2F*)fResiduals->UncheckedAt(kTrackResPerHalfCh+ia))->ProjectionY("tmp",k+1,k+1,"e");
	GetMeanRMS(tmp, sigmaTrack, sigmaTrackErr, dumy1, dumy2, 0x0, 0x0, 0, 0, kFALSE, kFALSE);
	delete tmp;
	
	tmp = ((TH2F*)fResiduals->UncheckedAt(kMCSPerHalfCh+ia))->ProjectionY("tmp",k+1,k+1,"e");
	GetMeanRMS(tmp, sigmaMCS, sigmaMCSErr, dumy1, dumy2, 0x0, 0x0, 0, 0, kFALSE, kFALSE);
	delete tmp;
	
	sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
	if (sigmaCluster > 0.) {
	  sigmaCluster = TMath::Sqrt(sigmaCluster);
	  sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
	} else {
	  sigmaCluster = 0.;
	  sigmaClusterErr = 0.;
	}
	((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerHalfCh+ia))->SetPoint(k, k+1, sigmaCluster);
	((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerHalfCh+ia))->SetPointError(k, 0., sigmaClusterErr);
	
	// method 1 versus track angle
	FillMeanSigmaClusterVsX((TH2F*)fResidualsVsAngle->UncheckedAt(kResidualInHalfChVsAngleClusterIn+20*ia+2*i+j), 0x0,
				(TGraphErrors*)((TMultiGraph*)fChamberRes->UncheckedAt(kHChResidualMeanClusterInVsAngle+ia))->GetListOfGraphs()->FindObject(Form("gShift%sVsAngle_halfCh%d%s",axes[ia],i+1,side[j])), 0x0);
	
      }
      
      // compute averaged local chi2
      tmp = ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerCh+ia))->ProjectionY("tmp",i+1,i+1,"e");
      ((TGraphErrors*)fLocalChi2->UncheckedAt(kLocalChi2PerChMean+ia))->SetPoint(i, i+1, tmp->GetMean());
      ((TGraphErrors*)fLocalChi2->UncheckedAt(kLocalChi2PerChMean+ia))->SetPointError(i, 0., tmp->GetMeanError());
      delete tmp;
      
    }
    
    // compute residual mean and dispersion per DE
    for (Int_t i = 0; i < fNDE; i++) {
      
      // method 1
      TH1D *tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterIn+ia))->ProjectionY(Form("hRes%sDE%dClIn",axes[ia],i+1),i+1,i+1,"e");
      tmp->SetTitle(Form("DE %s",((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterIn+ia))->GetXaxis()->GetBinLabel(i+1)));
      GetMeanRMS(tmp, meanIn, meanInErr, sigmaIn, sigmaInErr, (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterIn+ia), 0x0, i, i+1);
      fTmpHists->AddLast(tmp);
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterOut+ia))->ProjectionY(Form("hRes%sDE%dClOut",axes[ia],i+1),i+1,i+1,"e");
      tmp->SetTitle(Form("DE %s",((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterOut+ia))->GetXaxis()->GetBinLabel(i+1)));
      GetMeanRMS(tmp, meanOut, meanOutErr, sigmaOut, sigmaOutErr, (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterOut+ia), 0x0, i, i+1);
      fTmpHists->AddLast(tmp);
      
      if (fCorrectForSystematics) {
	sigma = TMath::Sqrt(sigmaIn*sigmaIn + meanIn*meanIn);
	sigmaInErr = (sigma>0) ? TMath::Sqrt(sigmaIn*sigmaIn*sigmaInErr*sigmaInErr + meanIn*meanIn*meanInErr*meanInErr) / sigma : 0.;
	sigmaIn = sigma;
	sigma = TMath::Sqrt(sigmaOut*sigmaOut + meanOut*meanOut);
	sigmaOutErr = (sigma>0) ? TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + meanOut*meanOut*meanOutErr*meanOutErr) / sigma : 0.;
	sigmaOut = sigma;
      }
      
      clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
      //      clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
      clusterResErr = TMath::Sqrt(sigmaInErr*sigmaOutErr);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerDESigma+ia))->SetPoint(i, i+1, clusterRes);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerDESigma+ia))->SetPointError(i, 0., clusterResErr);
      
      // method 2
      tmp = ((TH2F*)fResiduals->UncheckedAt(kTrackResPerDE+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMeanRMS(tmp, sigmaTrack, sigmaTrackErr, dumy1, dumy2, 0x0, 0x0, 0, 0, kFALSE, kFALSE);
      delete tmp;
      
      tmp = ((TH2F*)fResiduals->UncheckedAt(kMCSPerDE+ia))->ProjectionY("tmp",i+1,i+1,"e");
      GetMeanRMS(tmp, sigmaMCS, sigmaMCSErr, dumy1, dumy2, 0x0, 0x0, 0, 0, kFALSE, kFALSE);
      delete tmp;
      
      sigmaCluster = sigmaOut*sigmaOut - sigmaTrack*sigmaTrack;
      if (sigmaCluster > 0.) {
	sigmaCluster = TMath::Sqrt(sigmaCluster);
	sigmaClusterErr = TMath::Sqrt(sigmaOut*sigmaOut*sigmaOutErr*sigmaOutErr + sigmaTrack*sigmaTrack*sigmaTrackErr*sigmaTrackErr) / sigmaCluster;
      } else {
	sigmaCluster = 0.;
	sigmaClusterErr = 0.;
      }
      ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerDE+ia))->SetPoint(i, i+1, sigmaCluster);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerDE+ia))->SetPointError(i, 0., sigmaClusterErr);
      
      // compute averaged local chi2
      tmp = ((TH2F*)fResiduals->UncheckedAt(kLocalChi2PerDE+ia))->ProjectionY("tmp",i+1,i+1,"e");
      ((TGraphErrors*)fLocalChi2->UncheckedAt(kLocalChi2PerDEMean+ia))->SetPoint(i, i+1, tmp->GetMean());
      ((TGraphErrors*)fLocalChi2->UncheckedAt(kLocalChi2PerDEMean+ia))->SetPointError(i, 0., tmp->GetMeanError());
      delete tmp;
      
    }
    
    // set half-chamber graph labels
    TAxis* xAxis = ((TH2F*)fResiduals->UncheckedAt(kResidualPerHalfChClusterIn+ia))->GetXaxis();
    ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterIn+ia))->GetXaxis()->Set(20, 0.5, 20.5);
    ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterOut+ia))->GetXaxis()->Set(20, 0.5, 20.5);
    ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerHalfChSigma+ia))->GetXaxis()->Set(20, 0.5, 20.5);
    ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerHalfCh+ia))->GetXaxis()->Set(20, 0.5, 20.5);
    for (Int_t i = 1; i <= 20; i++) {
      const char* label = xAxis->GetBinLabel(i);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterIn+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterOut+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerHalfChSigma+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerHalfCh+ia))->GetXaxis()->SetBinLabel(i, label);
    }
    
    // set DE graph labels
    xAxis = ((TH2F*)fResiduals->UncheckedAt(kResidualPerDEClusterOut+ia))->GetXaxis();
    ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterIn+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterOut+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerDESigma+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerDE+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    ((TGraphErrors*)fLocalChi2->UncheckedAt(kLocalChi2PerDEMean+ia))->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
    for (Int_t i = 1; i <= fNDE; i++) {
      const char* label = xAxis->GetBinLabel(i);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterIn+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterOut+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerDESigma+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerDE+ia))->GetXaxis()->SetBinLabel(i, label);
      ((TGraphErrors*)fLocalChi2->UncheckedAt(kLocalChi2PerDEMean+ia))->GetXaxis()->SetBinLabel(i, label);
    }
    
  }
  
  // compute averaged local chi2 per chamber (X+Y)
  TH2F* h2 = (TH2F*)fResiduals->UncheckedAt(kLocalChi2PerCh+2);
  g = new TGraphErrors(AliMUONConstants::NTrackingCh());
  g->SetName("gLocalChi2PerChMean");
  g->SetTitle("local chi2 per Ch: mean;chamber ID;<local #chi^{2}>");
  g->SetMarkerStyle(kFullDotLarge);
  fLocalChi2->AddAtAndExpand(g, kLocalChi2PerChMean+2);
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    TH1D* tmp = h2->ProjectionY("tmp",i+1,i+1,"e");
    g->SetPoint(i, i+1, tmp->GetMean());
    g->SetPointError(i, 0., tmp->GetMeanError());
    delete tmp;
  }
  
  // compute averaged local chi2 per DE (X+Y)
  h2 = (TH2F*)fResiduals->UncheckedAt(kLocalChi2PerDE+2);
  g = new TGraphErrors(fNDE);
  g->SetName("gLocalChi2PerDEMean");
  g->SetTitle("local chi2 per DE: mean;DE ID;<local #chi^{2}>");
  g->SetMarkerStyle(kFullDotLarge);
  fLocalChi2->AddAtAndExpand(g, kLocalChi2PerDEMean+2);
  for (Int_t i = 0; i < fNDE; i++) {
    TH1D* tmp = h2->ProjectionY("tmp",i+1,i+1,"e");
    g->SetPoint(i, i+1, tmp->GetMean());
    g->SetPointError(i, 0., tmp->GetMeanError());
    delete tmp;
  }
  
  // set graph labels
  g->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
  for (Int_t i = 1; i <= fNDE; i++) {
    const char* label = h2->GetXaxis()->GetBinLabel(i);
    g->GetXaxis()->SetBinLabel(i, label);
  }
  
  // display
  TLegend *lResPerChMean = new TLegend(0.75,0.85,0.99,0.99);
  TLegend *lResPerChSigma1 = new TLegend(0.75,0.85,0.99,0.99);
  TLegend *lResPerChSigma2 = new TLegend(0.75,0.85,0.99,0.99);
  TLegend *lResPerChSigma3 = new TLegend(0.75,0.85,0.99,0.99);
  
  TCanvas* cResPerCh = new TCanvas("cResPerCh","cResPerCh",1200,500);
  cResPerCh->Divide(4,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cResPerCh->cd(1+4*ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChMeanClusterOut+ia);
    g->Draw("ap");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    if (ia == 0) lResPerChMean->AddEntry(g,"cluster out","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChMeanClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    if (ia == 0) lResPerChMean->AddEntry(g,"cluster in","PL");
    if (ia == 0) lResPerChMean->Draw();
    else lResPerChMean->DrawClone();
    gROOT->SetSelectedPad(cResPerCh->cd(2+4*ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChSigmaClusterOut+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    if (ia == 0) lResPerChSigma1->AddEntry(g,"cluster out","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChSigmaClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    if (ia == 0) lResPerChSigma1->AddEntry(g,"cluster in","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kMCSPerChMean+ia);
    g->Draw("p");
    g->SetMarkerColor(5);
    g->SetLineColor(5);
    if (ia == 0) lResPerChSigma1->AddEntry(g,"MCS","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerChSigma+ia);
    g->Draw("p");
    g->SetMarkerColor(3);
    g->SetLineColor(3);
    if (ia == 0) lResPerChSigma1->AddEntry(g,"combined 1","PL");
    if (ia == 0) lResPerChSigma1->Draw();
    else lResPerChSigma1->DrawClone();
    gROOT->SetSelectedPad(cResPerCh->cd(3+4*ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerChDispersionClusterOut+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    if (ia == 0) lResPerChSigma2->AddEntry(g,"cluster out","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kMCSPerChMean+ia);
    g->Draw("p");
    if (ia == 0) lResPerChSigma2->AddEntry(g,"MCS","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kTrackResPerChMean+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    if (ia == 0) lResPerChSigma2->AddEntry(g,"track res.","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerCh+ia);
    g->Draw("p");
    if (ia == 0) lResPerChSigma2->AddEntry(g,"combined 2","PL");
    if (ia == 0) lResPerChSigma2->Draw();
    else lResPerChSigma2->DrawClone();
    gROOT->SetSelectedPad(cResPerCh->cd(4+4*ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerChSigma+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    if (ia == 0) lResPerChSigma3->AddEntry(g,"combined 1","PL");
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerCh+ia);
    g->Draw("p");
    if (ia == 0) lResPerChSigma3->AddEntry(g,"combined 2","PL");
    if (ia == 0) lResPerChSigma3->Draw();
    else lResPerChSigma3->DrawClone();
  }
  fCanvases->AddAtAndExpand(cResPerCh, kResPerCh);
  
  TCanvas* cResPerHalfCh = new TCanvas("cResPerHalfCh","cResPerHalfCh",1200,500);
  cResPerHalfCh->Divide(2,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cResPerHalfCh->cd(1+2*ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterOut+ia);
    g->Draw("ap");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    lResPerChMean->DrawClone();
    gROOT->SetSelectedPad(cResPerHalfCh->cd(2+2*ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerHalfChSigma+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(3);
    g->SetLineColor(3);
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerHalfCh+ia);
    g->Draw("p");
    lResPerChSigma3->DrawClone();
  }
  fCanvases->AddAtAndExpand(cResPerHalfCh, kResPerHalfCh);
  
  TCanvas* cResPerDE = new TCanvas("cResPerDE","cResPerDE",1200,800);
  cResPerDE->Divide(1,4);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cResPerDE->cd(1+ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterOut+ia);
    g->Draw("ap");
    g->SetMarkerColor(2);
    g->SetLineColor(2);
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterIn+ia);
    g->Draw("p");
    g->SetMarkerColor(4);
    g->SetLineColor(4);
    lResPerChMean->DrawClone();
    gROOT->SetSelectedPad(cResPerDE->cd(3+ia));
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerDESigma+ia);
    g->Draw("ap");
    g->SetMinimum(0.);
    g->SetMarkerColor(3);
    g->SetLineColor(3);
    g = (TGraphErrors*)fChamberRes->UncheckedAt(kClusterResPerDE+ia);
    g->Draw("p");
    lResPerChSigma3->DrawClone();
  }
  fCanvases->AddAtAndExpand(cResPerDE, kResPerDE);
  
  TCanvas* cResPerChVsP = new TCanvas("cResPerChVsP","cResPerChVsP");
  cResPerChVsP->Divide(1,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cResPerChVsP->cd(1+ia));
    mg = (TMultiGraph*)fChamberRes->UncheckedAt(kCombinedResidualSigmaVsP+ia);
    mg->Draw("ap");
  }
  fCanvases->AddAtAndExpand(cResPerChVsP, kResPerChVsP);
  
  TCanvas* cResPerChVsCent = new TCanvas("cResPerChVsCent","cResPerChVsCent");
  cResPerChVsCent->Divide(1,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cResPerChVsCent->cd(1+ia));
    mg = (TMultiGraph*)fChamberRes->UncheckedAt(kCombinedResidualSigmaVsCent+ia);
    mg->Draw("ap");
  }
  fCanvases->AddAtAndExpand(cResPerChVsCent, kResPerChVsCent);
  
  TCanvas* cResPerChVsAngle = new TCanvas("cResPerChVsAngle","cResPerChVsAngle");
  cResPerChVsAngle->Divide(1,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cResPerChVsAngle->cd(1+ia));
    mg = (TMultiGraph*)fChamberRes->UncheckedAt(kCombinedResidualSigmaVsAngle+ia);
    mg->Draw("ap");
  }
  fCanvases->AddAtAndExpand(cResPerChVsAngle, kResPerChVsAngle);
  
  TCanvas* cShiftPerChVsP = new TCanvas("cShiftPerChVsP","cShiftPerChVsP");
  cShiftPerChVsP->Divide(1,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cShiftPerChVsP->cd(1+ia));
    mg = (TMultiGraph*)fChamberRes->UncheckedAt(kResidualMeanClusterInVsP+ia);
    mg->Draw("ap");
  }
  fCanvases->AddAtAndExpand(cShiftPerChVsP, kShiftPerChVsP);
  
  TCanvas* cShiftPerChVsCent = new TCanvas("cShiftPerChVsCent","cShiftPerChVsCent");
  cShiftPerChVsCent->Divide(1,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cShiftPerChVsCent->cd(1+ia));
    mg = (TMultiGraph*)fChamberRes->UncheckedAt(kResidualMeanClusterInVsCent+ia);
    mg->Draw("ap");
  }
  fCanvases->AddAtAndExpand(cShiftPerChVsCent, kShiftPerChVsCent);
  
  TCanvas* cShiftPerChVsAngle = new TCanvas("cShiftPerChVsAngle","cShiftPerChVsAngle");
  cShiftPerChVsAngle->Divide(1,2);
  for (Int_t ia = 0; ia < 2; ia++) {
    gROOT->SetSelectedPad(cShiftPerChVsAngle->cd(1+ia));
    mg = (TMultiGraph*)fChamberRes->UncheckedAt(kResidualMeanClusterInVsAngle+ia);
    mg->Draw("ap");
  }
  fCanvases->AddAtAndExpand(cShiftPerChVsAngle, kShiftPerChVsAngle);
  
  // print results
  if (fPrintClResPerCh) {
    printf("\nchamber resolution:\n");
    printf(" - non-bending:");
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %5.3f":", %5.3f",newClusterRes[0][i]);
    printf("\n -     bending:");
    for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) printf((i==0)?" %6.4f":", %6.4f",newClusterRes[1][i]);
    printf("\n\n");
  }
  
  if (fPrintClResPerDE) {
    Double_t iDE, clRes;
    printf("\nDE resolution:\n");
    printf(" - non-bending:");
    for (Int_t i = 0; i < fNDE; i++) {
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerDESigma))->GetPoint(i, iDE, clRes);
      printf((i==0)?" %5.3f":", %5.3f", clRes);
    }
    printf("\n -     bending:");
    for (Int_t i = 0; i < fNDE; i++) {
      ((TGraphErrors*)fChamberRes->UncheckedAt(kCombinedResidualPerDESigma+1))->GetPoint(i, iDE, clRes);
      printf((i==0)?" %6.4f":", %6.4f", clRes);
    }
    printf("\n\n");
  }
  
  if (fPrintHalfChShift) {
    Double_t iHCh, hChShift;
    printf("\nhalf-chamber residual displacements:\n");
    printf(" - non-bending:");
    for (Int_t i = 0; i < 2*AliMUONConstants::NTrackingCh(); i++) {
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterIn))->GetPoint(i, iHCh, hChShift);
      printf((i==0)?" %6.4f":", %6.4f", hChShift);
    }
    printf("\n -     bending:");
    for (Int_t i = 0; i < 2*AliMUONConstants::NTrackingCh(); i++) {
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerHalfChMeanClusterIn+1))->GetPoint(i, iHCh, hChShift);
      printf((i==0)?" %6.4f":", %6.4f", hChShift);
    }
    printf("\n\n");
  }
  
  if (fPrintDEShift) {
    Double_t iDE, deShift;
    printf("\nDE residual displacements:\n");
    printf(" - non-bending:");
    for (Int_t i = 0; i < fNDE; i++) {
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterIn))->GetPoint(i, iDE, deShift);
      printf((i==0)?" %6.4f":", %6.4f", deShift);
    }
    printf("\n -     bending:");
    for (Int_t i = 0; i < fNDE; i++) {
      ((TGraphErrors*)fChamberRes->UncheckedAt(kResidualPerDEMeanClusterIn+1))->GetPoint(i, iDE, deShift);
      printf((i==0)?" %6.4f":", %6.4f", deShift);
    }
    printf("\n\n");
  }
  
  // Post final data.
  PostData(6, fLocalChi2);
  PostData(7, fChamberRes);
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::ModifyClusters(AliMUONTrack& track)
{
  /// Reset the clusters resolution from the ones given to the task and change
  /// the cluster position according to the new alignment parameters if required
  
  Double_t gX,gY,gZ,lX,lY,lZ;
  
  // loop over clusters
  Int_t nClusters = track.GetNClusters();
  for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
    
    AliMUONVCluster* cl = static_cast<AliMUONTrackParam*>(track.GetTrackParamAtCluster()->UncheckedAt(iCluster))->GetClusterPtr();
    Int_t chId = cl->GetChamberId();
    Int_t halfChId = (cl->GetX() > 0) ? 2*chId : 2*chId+1;
    Int_t deId = cl->GetDetElemId();
    
    // change their resolution
    cl->SetErrXY(fClusterResNB[chId], fClusterResB[chId]);
    
    // change their position
    gX = cl->GetX();
    gY = cl->GetY();
    gZ = cl->GetZ();
    if (fReAlign) { // change the alignement
      fOldGeoTransformer->Global2Local(deId,gX,gY,gZ,lX,lY,lZ);
      fNewGeoTransformer->Local2Global(deId,lX,lY,lZ,gX,gY,gZ);
    }
    if (fShiftHalfCh) { // correct for half-chamber displacement
      gX -= fHalfChShiftNB[halfChId];
      gY -= fHalfChShiftB[halfChId];
    }
    if (fShiftDE) { // correct for DE displacement
      gX -= fDEShiftNB[fDEIndices[deId]-1];
      gY -= fDEShiftB[fDEIndices[deId]-1];
    }
    cl->SetXYZ(gX,gY,gZ);
    
    // "remove" mono-cathod clusters on stations 3-4-5 if required
    // (to be done after moving clusters to the new position)
    if (fRemoveMonoCathCl && chId > 3) {
      Bool_t hasBending, hasNonBending;
      if (fCheckAllPads) CheckPads(cl, hasBending, hasNonBending);
      else CheckPadsBelow(cl, hasBending, hasNonBending);
      if (!hasNonBending) cl->SetErrXY(10., cl->GetErrY());
      if (!hasBending) cl->SetErrXY(cl->GetErrX(), 10.);
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::Zoom(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic at each edge
  ZoomLeft(h, fractionCut);
  ZoomRight(h, fractionCut);
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::ZoomLeft(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic on the left side
  Int_t maxEventsCut = (Int_t) (fractionCut * h->GetEntries());
  Int_t nBins = h->GetNbinsX();
  
  // set low edge  
  Int_t minBin;
  Int_t eventsCut = 0;
  for (minBin = 1; minBin <= nBins; minBin++) {
    eventsCut += (Int_t) h->GetBinContent(minBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set new axis range
  h->GetXaxis()->SetRange(--minBin, h->GetXaxis()->GetLast());
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::ZoomRight(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic on the right side
  Int_t maxEventsCut = (Int_t) (fractionCut * h->GetEntries());
  Int_t nBins = h->GetNbinsX();
  
  // set high edge
  Int_t maxBin;
  Int_t eventsCut = 0;
  for (maxBin = nBins; maxBin >= 1; maxBin--) {
    eventsCut += (Int_t) h->GetBinContent(maxBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set new axis range
  h->GetXaxis()->SetRange(h->GetXaxis()->GetFirst(), ++maxBin);
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::GetMeanRMS(TH1* h, Double_t& mean, Double_t& meanErr,
					       Double_t& rms, Double_t& rmsErr,
					       TGraphErrors* gMean, TGraphErrors* gRMS,
					       Int_t i, Double_t x, Bool_t zoom, Bool_t enableFit)
{
  /// Fill graphs with the mean and rms values and the corresponding error (zooming if required)
  
  if (h->GetEntries() < fgkMinEntries) { // not enough entries
    
    mean = 0.;
    meanErr = 0.;
    rms = 0.;
    rmsErr = 0.;
    
  } else if (enableFit && fGaus) { // take the mean of a gaussian fit
    
    // first fit
    Double_t xMin = h->GetXaxis()->GetXmin();
    Double_t xMax = h->GetXaxis()->GetXmax();
    fGaus->SetRange(xMin, xMax);
    fGaus->SetParameters(h->GetEntries(), 0., 0.1);
    fGaus->SetParLimits(1, xMin, xMax);
    h->Fit("fGaus", "WWNQ");
    
    // rebin histo
    Int_t rebin = static_cast<Int_t>(TMath::Min(0.1*h->GetNbinsX(),TMath::Max(0.3*fGaus->GetParameter(2)/h->GetBinWidth(1),1.)));
    while (h->GetNbinsX()%rebin!=0) rebin--;
    h->Rebin(rebin);
    
    // second fit
    xMin = TMath::Max(fGaus->GetParameter(1)-10.*fGaus->GetParameter(2), h->GetXaxis()->GetXmin());
    xMax = TMath::Min(fGaus->GetParameter(1)+10.*fGaus->GetParameter(2), h->GetXaxis()->GetXmax());
    fGaus->SetRange(xMin, xMax);
    fGaus->SetParLimits(1, xMin, xMax);
    h->Fit("fGaus","NQR");
    
    mean = fGaus->GetParameter(1);
    meanErr = fGaus->GetParError(1);
    rms = fGaus->GetParameter(2);
    rmsErr = fGaus->GetParError(2);
    
    // display the detail of the fit
    if (!strstr(h->GetName(),"tmp")) {
      Int_t ia = (strstr(h->GetName(),"ResX")) ? 0 : 1;
      Int_t ib = (strstr(h->GetName(),"ClIn")) ? 0 : 1;
      if (strstr(h->GetName(),"Half"))
        gROOT->SetSelectedPad(((TCanvas*)fCanvases->UncheckedAt(kDetailResPerHalfCh+ib+2*ia))->cd(i+1));
      else if (strstr(h->GetName(),"Ch"))
        gROOT->SetSelectedPad(((TCanvas*)fCanvases->UncheckedAt(kDetailResPerCh+ib+2*ia))->cd(i+1));
      else
        gROOT->SetSelectedPad(((TCanvas*)fCanvases->UncheckedAt(kDetailResPerDE+ib+2*ia))->cd(i+1));
      gPad->SetLogy();
      h->Draw("hist");
      TF1* f = (TF1*)fGaus->DrawClone("same");
      f->SetLineWidth(1);
      f->SetLineColor(2);
      fTmpHists->AddLast(f);
    }
    
  } else { // take the mean of the distribution
    
    if (zoom) Zoom(h);
    
    mean = h->GetMean();
    meanErr = h->GetMeanError();
    rms = h->GetRMS();
    rmsErr = h->GetRMSError();
    
    if (zoom) h->GetXaxis()->SetRange(0,0);
    
  }
  
  // fill graph if required
  if (gMean) {
    gMean->SetPoint(i, x, mean);
    gMean->SetPointError(i, 0., meanErr);
  }
  if (gRMS) {
    gRMS->SetPoint(i, x, rms);
    gRMS->SetPointError(i, 0., rmsErr);
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::FillMeanSigmaClusterVsX(const TH2* hIn, const TH2* hOut,
							    TGraphErrors* gMean, TGraphErrors* gSigma)
{
  /// Fill graph with cluster shift (cluster in) if gMean != OxO
  /// and resolution from combined residuals with cluster in/out (zooming if required)
  Double_t meanIn, meanInErr, sigmaIn, sigmaInErr, sigmaOut, sigmaOutErr, clusterRes, clusterResErr, dumy1, dumy2;
  for (Int_t j = 1; j <= hIn->GetNbinsX(); j++) {
    TH1D* tmp = hIn->ProjectionY("tmp",j,j,"e");
    GetMeanRMS(tmp, meanIn, meanInErr, sigmaIn, sigmaInErr, 0x0, 0x0, 0, 0.);
    delete tmp;
    if (hOut) {
      tmp = hOut->ProjectionY("tmp",j,j,"e");
      GetMeanRMS(tmp, dumy1, dumy2, sigmaOut, sigmaOutErr, 0x0, 0x0, 0, 0.);
      delete tmp;
    }
    Double_t x = 0.5 * (hIn->GetXaxis()->GetBinLowEdge(j) + hIn->GetXaxis()->GetBinLowEdge(j+1));
    Double_t xErr = x - hIn->GetXaxis()->GetBinLowEdge(j);
    if (gMean) {
      gMean->SetPoint(j-1, x, meanIn);
      gMean->SetPointError(j-1, xErr, meanInErr);
    }
    if (hOut && gSigma) {
      clusterRes = TMath::Sqrt(sigmaIn*sigmaOut);
      //clusterResErr = (clusterRes > 0.) ? 0.5 * TMath::Sqrt(sigmaInErr*sigmaInErr*sigmaOut*sigmaOut + sigmaIn*sigmaIn*sigmaOutErr*sigmaOutErr) / clusterRes : 0.;
      clusterResErr = TMath::Sqrt(sigmaInErr*sigmaOutErr);
      gSigma->SetPoint(j-1, x, clusterRes);
      gSigma->SetPointError(j-1, xErr, clusterResErr);
    }
  }
}

//__________________________________________________________________________
void AliAnalysisTaskMuonResolution::Cov2CovP(const AliMUONTrackParam &param, TMatrixD &covP)
{
  /// change coordinate system: (X, SlopeX, Y, SlopeY, q/Pyz) -> (X, Y, pX, pY, pZ)
  /// parameters (param) are given in the (X, SlopeX, Y, SlopeY, q/Pyz) coordinate system
  
  // Get useful parameters
  Double_t slopeX = param.GetNonBendingSlope();
  Double_t slopeY = param.GetBendingSlope();
  Double_t qOverPYZ = param.GetInverseBendingMomentum();
  Double_t pZ = param.Pz();
  
  // compute Jacobian to change the coordinate system from (X,SlopeX,Y,SlopeY,c/pYZ) to (X,Y,pX,pY,pZ)
  Double_t dpZdSlopeY = - qOverPYZ * qOverPYZ * pZ * pZ * pZ * slopeY;
  Double_t dpZdQOverPYZ = (qOverPYZ != 0.) ? - pZ / qOverPYZ : - FLT_MAX;
  TMatrixD jacob(5,5);
  jacob.Zero();
  jacob(0,0) = 1.;
  jacob(1,2) = 1.;
  jacob(2,1) = pZ;
  jacob(2,3) = slopeX * dpZdSlopeY;
  jacob(2,4) = slopeX * dpZdQOverPYZ;
  jacob(3,3) = pZ + slopeY * dpZdSlopeY;
  jacob(3,4) = slopeY * dpZdQOverPYZ;
  jacob(4,3) = dpZdSlopeY;
  jacob(4,4) = dpZdQOverPYZ;
  
  // compute covariances in new coordinate system
  TMatrixD tmp(param.GetCovariances(),TMatrixD::kMultTranspose,jacob);
  covP.Mult(jacob,tmp);
}

//__________________________________________________________________________
void AliAnalysisTaskMuonResolution::CheckPads(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const
{
  /// Check that this cluster contains pads on both cathods
  
  // reset
  hasBending = kFALSE;
  hasNonBending = kFALSE;
  
  // loop over digits contained in the cluster
  for (Int_t iDigit = 0; iDigit < cl->GetNDigits(); iDigit++) {
    
    Int_t manuId = AliMUONVDigit::ManuId(cl->GetDigitId(iDigit));
    
    // check the location of the manu the digit belongs to
    if (manuId > 0) {
      if (manuId & AliMpConstants::ManuMask(AliMp::kNonBendingPlane)) hasNonBending = kTRUE;
      else hasBending = kTRUE;
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonResolution::CheckPadsBelow(AliMUONVCluster *cl, Bool_t &hasBending, Bool_t &hasNonBending) const
{
  /// Check that this cluster contains pads on both cathods just under its position
  
  // reset
  hasBending = kFALSE;
  hasNonBending = kFALSE;
  
  // get the cathod corresponding to the bending/non-bending plane
  Int_t deId = cl->GetDetElemId();
  AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(deId, kFALSE);
  if (!de) return;
  AliMp::CathodType cath1 = de->GetCathodType(AliMp::kBendingPlane); 
  AliMp::CathodType cath2 = de->GetCathodType(AliMp::kNonBendingPlane); 
  
  // get the corresponding segmentation
  const AliMpVSegmentation* seg1 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, cath1);
  const AliMpVSegmentation* seg2 = AliMpSegmentation::Instance()->GetMpSegmentation(deId, cath2);
  if (!seg1 || !seg2) return;
  
  // get global coordinate of the cluster
  Double_t gX = cl->GetX();
  Double_t gY = cl->GetY();
  Double_t gZ = cl->GetZ();
  
  // revert half-chamber or DE displacement if any
  Int_t chId = cl->GetChamberId();
  Int_t halfChId = (cl->GetX() > 0) ? 2*chId : 2*chId+1;
  if (fShiftHalfCh) {
    gX += fHalfChShiftNB[halfChId];
    gY += fHalfChShiftB[halfChId];
  }
  if (fShiftDE) {
    gX += fDEShiftNB[fDEIndices[deId]-1];
    gY += fDEShiftB[fDEIndices[deId]-1];
  }
  
  // get local coordinate of the cluster
  Double_t lX,lY,lZ;
  fNewGeoTransformer->Global2Local(deId,gX,gY,gZ,lX,lY,lZ);
  
  // find pads below the cluster
  AliMpPad pad1 = seg1->PadByPosition(lX, lY, kFALSE);
  AliMpPad pad2 = seg2->PadByPosition(lX, lY, kFALSE);
  
  // build their ID if pads are valid
  UInt_t padId1 = (pad1.IsValid()) ? AliMUONVDigit::BuildUniqueID(deId, pad1.GetManuId(), pad1.GetManuChannel(), cath1) : 0;
  UInt_t padId2 = (pad2.IsValid()) ? AliMUONVDigit::BuildUniqueID(deId, pad2.GetManuId(), pad2.GetManuChannel(), cath2) : 0;
  
  // check if the cluster contains these pads 
  for (Int_t iDigit = 0; iDigit < cl->GetNDigits(); iDigit++) {
    
    UInt_t digitId = cl->GetDigitId(iDigit);
    
    if (digitId == padId1) {
      
      hasBending = kTRUE;
      if (hasNonBending) break;
      
    } else if (digitId == padId2) {
      
      hasNonBending = kTRUE;
      if (hasBending) break;
      
    }
    
  }
  
}

