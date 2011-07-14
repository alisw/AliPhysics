/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskMuonPerformance
/// Analysis task to chek the tracker and trigger reconstruction
/// The output is a list of histograms.
///
/// \author Diego Stocco and Philippe Pillot
//-----------------------------------------------------------------------------

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TMCProcess.h"
#include "TGeoManager.h"

// STEER includes
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

// CORRFW includes
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "AliCFEffGrid.h"

// MUON includes
#include "AliMUONConstants.h"
#include "AliMUONTrack.h"
#include "AliMUONTriggerTrack.h"
#include "AliMUONLocalTrigger.h"
#include "AliMUONVTrackStore.h"
#include "AliMUONVTriggerTrackStore.h"
#include "AliMUONESDInterface.h"
#include "AliMUONRecoParam.h"
#include "AliMUONCDB.h"
#include "AliMUONTrackExtrap.h"
#include "AliMUONTrackParam.h"
#include "AliMUONRecoCheck.h"
#include "AliMUONVCluster.h"
#include "AliMpDEIterator.h"

#include "AliAnalysisTaskMuonPerformance.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskMuonPerformance) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskMuonPerformance::AliAnalysisTaskMuonPerformance() :
AliAnalysisTaskSE(),
fDefaultStorage(""),
fNPBins(30),
fCorrectForSystematics(kFALSE),
fFitResiduals(kFALSE),
fRequestedStationMask(0),
fRequest2ChInSameSt45(0),
fSigmaCut(-1.),
fSigmaCutTrig(-1.),
fNDE(0),
fCFContainer(0x0),
fEfficiencyList(0x0),
fTriggerList(0x0),
fTrackerList(0x0),
fPAtVtxList(0x0),
fSlopeAtVtxList(0x0),
fEtaAtVtxList(0x0),
fPhiAtVtxList(0x0),
fPAt1stClList(0x0),
fSlopeAt1stClList(0x0),
fDCAList(0x0),
fClusterList(0x0)
{
  //
  /// Default Constructor.
  //
}


//________________________________________________________________________
AliAnalysisTaskMuonPerformance::AliAnalysisTaskMuonPerformance(const char *name) :
AliAnalysisTaskSE(name),
fDefaultStorage("raw://"),
fNPBins(30),
fCorrectForSystematics(kFALSE),
fFitResiduals(kFALSE),
fRequestedStationMask(0),
fRequest2ChInSameSt45(0),
fSigmaCut(-1.),
fSigmaCutTrig(-1.),
fNDE(0),
fCFContainer(0x0),
fEfficiencyList(0x0),
fTriggerList(0x0),
fTrackerList(0x0),
fPAtVtxList(0x0),
fSlopeAtVtxList(0x0),
fEtaAtVtxList(0x0),
fPhiAtVtxList(0x0),
fPAt1stClList(0x0),
fSlopeAt1stClList(0x0),
fDCAList(0x0),
fClusterList(0x0)
{
  //
  /// Constructor.
  //
  fPRange[0] = 0.;
  fPRange[1] = 300.;
  fClusterMaxRes[0] = 0.;
  fClusterMaxRes[1] = 0.;
  
  DefineOutput(1, AliCFContainer::Class());
  DefineOutput(2, TObjArray::Class());
  DefineOutput(3, TObjArray::Class());
  DefineOutput(4, TObjArray::Class());
  DefineOutput(5, TObjArray::Class());
  DefineOutput(6, TObjArray::Class());
  DefineOutput(7, TObjArray::Class());
  DefineOutput(8, TObjArray::Class());
  DefineOutput(9, TObjArray::Class());
  DefineOutput(10, TObjArray::Class());
  DefineOutput(11, TObjArray::Class());
  DefineOutput(12, TObjArray::Class());
}


//________________________________________________________________________
AliAnalysisTaskMuonPerformance::~AliAnalysisTaskMuonPerformance()
{
  //
  /// Destructor
  //
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fCFContainer;
    delete fEfficiencyList;
    delete fTriggerList;
    delete fTrackerList;
    delete fPAtVtxList;
    delete fSlopeAtVtxList;
    delete fEtaAtVtxList;
    delete fPhiAtVtxList;
    delete fPAt1stClList;
    delete fSlopeAt1stClList;
    delete fDCAList;
    delete fClusterList;
  }
}


//___________________________________________________________________________
void AliAnalysisTaskMuonPerformance::NotifyRun()
{
  //
  /// Use the event handler information to correctly fill the analysis flags:
  /// - check if Monte Carlo information is present
  //
  
  // load OCDB objects only once
  if (fSigmaCut > 0) return;
  
  // set OCDB location
  AliCDBManager* cdbm = AliCDBManager::Instance();
  cdbm->SetDefaultStorage(fDefaultStorage.Data());
  cdbm->SetRun(fCurrentRunNumber);
  
  // load magnetic field for track extrapolation
  if (!AliMUONCDB::LoadField()) return;
  
  // load mapping
  if (!AliMUONCDB::LoadMapping()) return;

  // load geometry for track extrapolation to vertex
  if (!AliGeomManager::GetGeometry()) AliGeomManager::LoadGeometry();
  if (!AliGeomManager::GetGeometry()) return;
  
  // load recoParam
  AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();
  if (!recoParam) {
    fRequestedStationMask = 0;
    fRequest2ChInSameSt45 = kFALSE;
    fSigmaCut = -1.;
    fSigmaCutTrig = -1.;
    AliError("--> skip this run");
    return;
  }
  
  // compute the mask of requested stations from recoParam
  fRequestedStationMask = 0;
  for (Int_t i = 0; i < 5; i++) if (recoParam->RequestStation(i)) fRequestedStationMask |= ( 1 << i );
  
  // get from recoParam whether a track need 2 chambers hit in the same station (4 or 5) or not to be reconstructible
  fRequest2ChInSameSt45 = !recoParam->MakeMoreTrackCandidates();
  
  // get sigma cut from recoParam to associate clusters with TrackRefs in case the labels are not used
  fSigmaCut = (recoParam->ImproveTracks()) ? recoParam->GetSigmaCutForImprovement() : recoParam->GetSigmaCutForTracking();
  
  // get sigma cut from recoParam to associate trigger track to triggerable track
  fSigmaCutTrig = recoParam->GetSigmaCutForTrigger();

  AliMUONESDInterface::ResetTracker(recoParam);
  
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    
    // find the highest chamber resolution and set histogram bins
    if (recoParam->GetDefaultNonBendingReso(i) > fClusterMaxRes[0]) fClusterMaxRes[0] = recoParam->GetDefaultNonBendingReso(i);
    if (recoParam->GetDefaultBendingReso(i) > fClusterMaxRes[1]) fClusterMaxRes[1] = recoParam->GetDefaultBendingReso(i);
    
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
  
  UserCreateOutputObjects();
}


//___________________________________________________________________________
void AliAnalysisTaskMuonPerformance::UserCreateOutputObjects() 
{
  //
  /// Create output objects
  //

  // do it once the OCDB has been loaded (i.e. from NotifyRun())
  if (fSigmaCut < 0) return;
  
  // ------ CFContainer ------
  
  // define axes of particle container
  Int_t nPtBins = 60;
  Double_t ptMin = 0., ptMax = 30.;
  TString ptName("Pt"), ptTitle("p_{t}"), ptUnits("GeV/c");

  Int_t nEtaBins = 25;
  Double_t etaMin = -4.5, etaMax = -2.;
  TString etaName("Eta"), etaTitle("#eta"), etaUnits("a.u.");

  Int_t nPhiBins = 36;
  Double_t phiMin = 0.; Double_t phiMax = 2*TMath::Pi();
  TString phiName("Phi"), phiTitle("#phi"), phiUnits("rad");

  Int_t nThetaAbsEndBins = 4;
  Double_t thetaAbsEndMin = -0.5, thetaAbsEndMax = 3.5;
  TString thetaAbsEndName("ThetaAbsEnd"), thetaAbsEndTitle("#theta_{abs}"), thetaAbsEndUnits("a.u.");  

  Int_t nChargeBins = 2;
  Double_t chargeMin = -2, chargeMax = 2.;
  TString chargeName("Charge"), chargeTitle("charge"), chargeUnits("e");

  Int_t nHasTrackerBins = 2;
  Double_t hasTrackerMin = -0.5, hasTrackerMax = (Double_t)nHasTrackerBins - 0.5;
  TString hasTrackerName("HasTracker"), hasTrackerTitle("Has tracker"), hasTrackerUnits("");
  
  Int_t nTriggerBins = kNtrigCuts;
  Double_t triggerMin = -0.5, triggerMax = (Double_t)nTriggerBins - 0.5;
  TString triggerName("MatchTrig"), triggerTitle("Trigger match"), triggerUnits("");

  Int_t nMotherTypeBins = kNtrackSources;
  Double_t motherTypeMin = -0.5, motherTypeMax = (Double_t)kNtrackSources - 0.5;
  TString motherTypeName("MotherType"), motherTypeTitle("motherType"), motherTypeUnits("");

  Int_t nMatchMCBins = kNMatchMC;
  Double_t matchMCMin = -0.5, matchMCMax = (Double_t)kNMatchMC - 0.5;
  TString matchMCName("MatchMC"), matchMCTitle("MatchMC"), matchMCUnits("");
  
  Int_t nbins[kNvars] = {nPtBins, nEtaBins, nPhiBins, nThetaAbsEndBins, nChargeBins, nHasTrackerBins, nTriggerBins, nMotherTypeBins, nMatchMCBins};
  Double_t xmin[kNvars] = {ptMin, etaMin, phiMin, thetaAbsEndMin, chargeMin, hasTrackerMin, triggerMin, motherTypeMin, matchMCMin};
  Double_t xmax[kNvars] = {ptMax, etaMax, phiMax, thetaAbsEndMax, chargeMax, hasTrackerMax, triggerMax, motherTypeMax, matchMCMax};
  TString axisTitle[kNvars] = {ptTitle, etaTitle, phiTitle, thetaAbsEndTitle, chargeTitle, hasTrackerTitle, triggerTitle, motherTypeTitle, matchMCTitle};
  TString axisUnits[kNvars] = {ptUnits, etaUnits, phiUnits, thetaAbsEndUnits, chargeUnits, hasTrackerUnits, triggerUnits, motherTypeUnits, matchMCUnits};
  
  // create particle container
  fCFContainer = new AliCFContainer(GetOutputSlot(1)->GetContainer()->GetName(),"container for tracks",kNsteps,kNvars,nbins);
  
  // set axes title and limits
  for ( Int_t idim = 0; idim<kNvars; idim++) {
    TString histoTitle = Form("%s (%s)", axisTitle[idim].Data(), axisUnits[idim].Data());
    histoTitle.ReplaceAll("()","");
    fCFContainer->SetVarTitle(idim, histoTitle.Data());
    fCFContainer->SetBinLimits(idim, xmin[idim], xmax[idim]);
  }
  
  // define histo names
  TString stepTitle[kNsteps] = {"reconstructed", "generated"};
  
  // define axes labels if any
  TString trigName[kNtrigCuts];
  trigName[kNoMatchTrig] = "NoMatch";
  trigName[kAllPtTrig]   = "AllPt";
  trigName[kLowPtTrig]   = "LowPt";
  trigName[kHighPtTrig]  = "HighPt";
  
  TString srcName[kNtrackSources];
  srcName[kCharmMu]     = "Charm";
  srcName[kBeautyMu]    = "Beauty";
  srcName[kPrimaryMu]   = "Decay";
  srcName[kSecondaryMu] = "Secondary";
  srcName[kRecoHadron]  = "Hadrons";
  srcName[kUnknownPart] = "Fakes";
  
  TString mMCName[kNMatchMC];
  mMCName[kNoMatch]     = "NoMatch";
  mMCName[kTrackerOnly] = "TrackerOnly";
  mMCName[kMatchedSame] = "MatchedSame";
  mMCName[kMatchedDiff] = "MatchedDiff";
  mMCName[kTriggerOnly] = "TriggerOnly";
  
  // set histo name and axis labels if any
  for (Int_t istep=0; istep<kNsteps; istep++) {
    
    // histo name
    fCFContainer->SetStepTitle(istep, stepTitle[istep].Data());
    AliCFGridSparse* gridSparse = fCFContainer->GetGrid(istep);
    
    // trigger info
    TAxis* triggerAxis =  gridSparse->GetAxis(kVarTrigger);
    for ( Int_t ibin=0; ibin<kNtrigCuts; ibin++ ) {
      triggerAxis->SetBinLabel(ibin+1,trigName[ibin]);
    }
    
    // mother type
    TAxis* motherTypeAxis = gridSparse->GetAxis(kVarMotherType);
    for ( Int_t ibin=0; ibin<kNtrackSources; ibin++ ) {
      motherTypeAxis->SetBinLabel(ibin+1,srcName[ibin]);
    }
    
    // MC matching flag
    TAxis* matchMCAxis = gridSparse->GetAxis(kVarMatchMC);
    for ( Int_t ibin=0; ibin<kNMatchMC; ibin++ ) {
      matchMCAxis->SetBinLabel(ibin+1,mMCName[ibin]);
    }
    
  }
  
  // ------ trigger histograms ------
  
  fTriggerList = new TObjArray(100);
  fTriggerList->SetOwner();
  
  TH1F* h1 = new TH1F("hResTrigX11", "Residual X11;X11_{reco} - X11_{MC} (cm)", 100, -10., 10.);
  fTriggerList->AddAt(h1, kResTrigX11);
  h1 = new TH1F("hResTrigY11", "Residual Y11;Y11_{reco} - Y11_{MC} (cm)", 100, -10., 10.);
  fTriggerList->AddAt(h1, kResTrigY11);
  h1 = new TH1F("hResTrigSlopeY", "Residual slope y;ySlope_{reco} - ySlope_{MC} (rad)", 100, -0.1, 0.1);
  fTriggerList->AddAt(h1, kResTrigSlopeY);
  
  // ------ tracker histograms ------
  
  fTrackerList = new TObjArray(100);
  fTrackerList->SetOwner();
  
  // momentum resolution at vertex
  const Int_t deltaPAtVtxNBins = 250;
  Double_t deltaPAtVtxEdges[2];
  deltaPAtVtxEdges[0] = -20. - 0.05 * fPRange[1];
  deltaPAtVtxEdges[1] =   5. + 0.05 * fPRange[1];
  
  h1 = new TH1F("hResPAtVtx"," delta P at vertex;#Delta_{p} (GeV/c)",deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h1, kResPAtVtx);
  TH2F *h2 = new TH2F("hResPAtVtxVsP","#Delta_{p} at vertex versus p;p (GeV/c);#Delta_{p} (GeV/c)",2*fNPBins,fPRange[0],fPRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsP);
  h2 = new TH2F("hResPAtVtxVsPIn23deg","#Delta_{p} at vertex versus p for tracks between 2 and 3 degrees at absorber end;p (GeV/c);#Delta_{p} (GeV/c)",2*fNPBins,fPRange[0],fPRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsPIn23deg);
  h2 = new TH2F("hResPAtVtxVsPIn310deg","#Delta_{p} at vertex versus p for tracks between 3 and 10 degrees at absorber end;p (GeV/c);#Delta_{p} (GeV/c)",2*fNPBins,fPRange[0],fPRange[1],deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsPIn310deg);
  h2 = new TH2F("hResPAtVtxVsPIn02degMC","#Delta_{p} at vertex versus p for tracks with MC angle below 2 degrees;p (GeV/c);#Delta_{p} (GeV/c)",2*fNPBins,fPRange[0],fPRange[1],deltaPAtVtxNBins/10,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsPIn02degMC);
  h2 = new TH2F("hResPAtVtxVsPosAbsEndIn02degMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{p} (GeV/c)",1000,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsPosAbsEndIn02degMC);
  h2 = new TH2F("hResPAtVtxVsPosAbsEndIn23degMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{p} (GeV/c)",1000,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsPosAbsEndIn23degMC);
  h2 = new TH2F("hResPAtVtxVsPosAbsEndIn310degMC","#Delta_{p} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{p} (GeV/c)",1000,0.,100.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsPosAbsEndIn310degMC);
  h2 = new TH2F("hResPAtVtxVsAngleAtAbsEnd","#Delta_{p} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{p} (GeV/c)",10,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsAngleAtAbsEnd);
  h2 = new TH2F("hResPAtVtxVsMCAngle","#Delta_{p} at vertex versus MC angle;MC angle (Deg);#Delta_{p} (GeV/c)",10,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPAtVtxVsMCAngle);
  TH3F *h3 = new TH3F("hResPAtVtxVsAngleAtAbsEndVsP","#Delta_{p} at vertex versus track position at absorber end converted to degrees versus momentum;p (GeV/c);angle (Deg);#Delta_{p} (GeV/c)",2*fNPBins,fPRange[0],fPRange[1],100,0.,10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0],deltaPAtVtxEdges[1]);
  fTrackerList->AddAt(h3, kResPAtVtxVsAngleAtAbsEndVsP);
  
  // transverse momentum resolution at vertex
  h2 = new TH2F("hResPtAtVtxVsPt","#Delta_{p_{t}} at vertex versus p_{t};p_{t} (GeV/c);#Delta_{p_{t}} (GeV/c)",2*fNPBins,fPRange[0]/10.,fPRange[1]/10.,deltaPAtVtxNBins,deltaPAtVtxEdges[0]/10.,deltaPAtVtxEdges[1]/10.);
  fTrackerList->AddAt(h2, kResPtAtVtxVsPt);
  
  // momentum resolution at first cluster
  const Int_t deltaPAtFirstClNBins = 500;
  Double_t deltaPAtFirstClEdges[2];
  deltaPAtFirstClEdges[0] = -5. - 0.05 * fPRange[1];
  deltaPAtFirstClEdges[1] =  5. + 0.05 * fPRange[1];
  
  h1 = new TH1F("hResPAt1stCl"," delta P at first cluster;#Delta_{p} (GeV/c)",deltaPAtFirstClNBins,deltaPAtFirstClEdges[0],deltaPAtFirstClEdges[1]);
  fTrackerList->AddAt(h1, kResPAt1stCl);
  h2 = new TH2F("hResPAt1stClVsP","#Delta_{p} at first cluster versus p;p (GeV/c);#Delta_{p} (GeV/c)",2*fNPBins,fPRange[0],fPRange[1],deltaPAtFirstClNBins,deltaPAtFirstClEdges[0],deltaPAtFirstClEdges[1]);
  fTrackerList->AddAt(h2, kResPAt1stClVsP);
  
  // transverse momentum resolution at vertex
  h2 = new TH2F("hResPtAt1stClVsPt","#Delta_{p_{t}} at first cluster versus p_{t};p_{t} (GeV/c);#Delta_{p_{t}} (GeV/c)",2*fNPBins,fPRange[0]/10.,fPRange[1]/10.,deltaPAtFirstClNBins,deltaPAtFirstClEdges[0]/10.,deltaPAtFirstClEdges[1]/10.);
  fTrackerList->AddAt(h2, kResPtAt1stClVsPt);
  
  // angular resolution at vertex
  const Int_t deltaSlopeAtVtxNBins = 500;
  const Double_t deltaSlopeAtVtxEdges[2] = {-0.05, 0.05};
  
  h1 = new TH1F("hResSlopeXAtVtx","#Delta_{slope_{X}} at vertex;#Delta_{slope_{X}}", deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h1, kResSlopeXAtVtx);
  h1 = new TH1F("hResSlopeYAtVtx","#Delta_{slope_{Y}} at vertex;#Delta_{slope_{Y}}", deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h1, kResSlopeYAtVtx);
  h2 = new TH2F("hResSlopeXAtVtxVsP","#Delta_{slope_{X}} at vertex versus p;p (GeV/c);#Delta_{slope_{X}}",2*fNPBins,fPRange[0],fPRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeXAtVtxVsP);
  h2 = new TH2F("hResSlopeYAtVtxVsP","#Delta_{slope_{Y}} at vertex versus p;p (GeV/c);#Delta_{slope_{Y}}",2*fNPBins,fPRange[0],fPRange[1], deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeYAtVtxVsP);
  h2 = new TH2F("hResSlopeXAtVtxVsPosAbsEndIn02degMC","#Delta_{slope_{X}} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{slope_{X}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeXAtVtxVsPosAbsEndIn02degMC);
  h2 = new TH2F("hResSlopeYAtVtxVsPosAbsEndIn02degMC","#Delta_{slope_{Y}} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{slope_{Y}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeYAtVtxVsPosAbsEndIn02degMC);
  h2 = new TH2F("hResSlopeXAtVtxVsPosAbsEndIn23degMC","#Delta_{slope_{X}} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{slope_{X}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeXAtVtxVsPosAbsEndIn23degMC);
  h2 = new TH2F("hResSlopeYAtVtxVsPosAbsEndIn23degMC","#Delta_{slope_{Y}} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{slope_{Y}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeYAtVtxVsPosAbsEndIn23degMC);
  h2 = new TH2F("hResSlopeXAtVtxVsPosAbsEndIn310degMC","#Delta_{slope_{X}} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{slope_{X}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeXAtVtxVsPosAbsEndIn310degMC);
  h2 = new TH2F("hResSlopeYAtVtxVsPosAbsEndIn310degMC","#Delta_{slope_{Y}} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{slope_{Y}}",1000,0.,100.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeYAtVtxVsPosAbsEndIn310degMC);
  h2 = new TH2F("hResSlopeXAtVtxVsAngleAtAbsEnd","#Delta_{slope_{X}} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{slope_{X}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeXAtVtxVsAngleAtAbsEnd);
  h2 = new TH2F("hResSlopeYAtVtxVsAngleAtAbsEnd","#Delta_{slope_{Y}} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{slope_{Y}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeYAtVtxVsAngleAtAbsEnd);
  h2 = new TH2F("hResSlopeXAtVtxVsMCAngle","#Delta_{slope_{X}} at vertex versus MC angle;MC angle (Deg);#Delta_{slope_{X}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeXAtVtxVsMCAngle);
  h2 = new TH2F("hResSlopeYAtVtxVsMCAngle","#Delta_{slope_{Y}} at vertex versus MC angle;MC angle (Deg);#Delta_{slope_{Y}}",10,0.,10.,deltaSlopeAtVtxNBins, deltaSlopeAtVtxEdges[0], deltaSlopeAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeYAtVtxVsMCAngle);
  
  // angular resolution at first cluster
  const Int_t deltaSlopeAtFirstClNBins = 500;
  const Double_t deltaSlopeAtFirstClEdges[2] = {-0.01, 0.01};
  
  h1 = new TH1F("hResSlopeXAt1stCl","#Delta_{slope_{X}} at first cluster;#Delta_{slope_{X}}", deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  fTrackerList->AddAt(h1, kResSlopeXAt1stCl);
  h1 = new TH1F("hResSlopeYAt1stCl","#Delta_{slope_{Y}} at first cluster;#Delta_{slope_{Y}}", deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  fTrackerList->AddAt(h1, kResSlopeYAt1stCl);
  h2 = new TH2F("hResSlopeXAt1stClVsP","#Delta_{slope_{X}} at first cluster versus p;p (GeV/c);#Delta_{slope_{X}}",2*fNPBins,fPRange[0],fPRange[1], deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeXAt1stClVsP);
  h2 = new TH2F("hResSlopeYAt1stClVsP","#Delta_{slope_{Y}} at first cluster versus p;p (GeV/c);#Delta_{slope_{Y}}",2*fNPBins,fPRange[0],fPRange[1], deltaSlopeAtFirstClNBins, deltaSlopeAtFirstClEdges[0], deltaSlopeAtFirstClEdges[1]);
  fTrackerList->AddAt(h2, kResSlopeYAt1stClVsP);
  
  // eta resolution at vertex
  const Int_t deltaEtaAtVtxNBins = 500;
  const Double_t deltaEtaAtVtxEdges[2] = {-0.5, 0.5};
  
  h1 = new TH1F("hResEtaAtVtx","#Delta_{eta} at vertex;#Delta_{eta}", deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  fTrackerList->AddAt(h1, kResEtaAtVtx);
  h2 = new TH2F("hResEtaAtVtxVsP","#Delta_{eta} at vertex versus p;p (GeV/c);#Delta_{eta}",2*fNPBins,fPRange[0],fPRange[1], deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResEtaAtVtxVsP);
  h2 = new TH2F("hResEtaAtVtxVsPosAbsEndIn02degMC","#Delta_{eta} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{eta}",1000,0.,100.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResEtaAtVtxVsPosAbsEndIn02degMC);
  h2 = new TH2F("hResEtaAtVtxVsPosAbsEndIn23degMC","#Delta_{eta} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{eta}",1000,0.,100.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResEtaAtVtxVsPosAbsEndIn23degMC);
  h2 = new TH2F("hResEtaAtVtxVsPosAbsEndIn310degMC","#Delta_{eta} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{eta}",1000,0.,100.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResEtaAtVtxVsPosAbsEndIn310degMC);
  h2 = new TH2F("hResEtaAtVtxVsAngleAtAbsEnd","#Delta_{eta} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{eta}",10,0.,10.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResEtaAtVtxVsAngleAtAbsEnd);
  h2 = new TH2F("hResEtaAtVtxVsMCAngle","#Delta_{eta} at vertex versus MC angle;MC angle (Deg);#Delta_{eta}",10,0.,10.,deltaEtaAtVtxNBins, deltaEtaAtVtxEdges[0], deltaEtaAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResEtaAtVtxVsMCAngle);
  
  // phi resolution at vertex
  const Int_t deltaPhiAtVtxNBins = 500;
  const Double_t deltaPhiAtVtxEdges[2] = {-0.5, 0.5};
  
  h1 = new TH1F("hResPhiAtVtx","#Delta_{phi} at vertex;#Delta_{phi}", deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  fTrackerList->AddAt(h1, kResPhiAtVtx);
  h2 = new TH2F("hResPhiAtVtxVsP","#Delta_{phi} at vertex versus p;p (GeV/c);#Delta_{phi}",2*fNPBins,fPRange[0],fPRange[1], deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPhiAtVtxVsP);
  h2 = new TH2F("hResPhiAtVtxVsPosAbsEndIn02degMC","#Delta_{phi} at vertex versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);#Delta_{phi}",1000,0.,100.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPhiAtVtxVsPosAbsEndIn02degMC);
  h2 = new TH2F("hResPhiAtVtxVsPosAbsEndIn23degMC","#Delta_{phi} at vertex versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);#Delta_{phi}",1000,0.,100.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPhiAtVtxVsPosAbsEndIn23degMC);
  h2 = new TH2F("hResPhiAtVtxVsPosAbsEndIn310degMC","#Delta_{phi} at vertex versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);#Delta_{phi}",1000,0.,100.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPhiAtVtxVsPosAbsEndIn310degMC);
  h2 = new TH2F("hResPhiAtVtxVsAngleAtAbsEnd","#Delta_{phi} at vertex versus track position at absorber end converted to degrees;angle (Deg);#Delta_{phi}",10,0.,10.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPhiAtVtxVsAngleAtAbsEnd);
  h2 = new TH2F("hResPhiAtVtxVsMCAngle","#Delta_{phi} at vertex versus MC angle;MC angle (Deg);#Delta_{phi}",10,0.,10.,deltaPhiAtVtxNBins, deltaPhiAtVtxEdges[0], deltaPhiAtVtxEdges[1]);
  fTrackerList->AddAt(h2, kResPhiAtVtxVsMCAngle);
  
  // DCA resolution
  const Int_t deltaPDCANBins = 500;
  const Double_t deltaPDCAEdges[2] = {0., 1000.};
  const Double_t deltaPMCSAngEdges[2] = {-1.5, 1.5};
  
  h1 = new TH1F("hPDCA","p #times DCA at vertex;p #times DCA (GeV #times cm)", deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h1, kPDCA);
  h2 = new TH2F("hPDCAVsPIn23deg","p #times DCA versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times DCA (GeV #times cm)",2*fNPBins,fPRange[0],fPRange[1], deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h2, kPDCAVsPIn23deg);
  h2 = new TH2F("hPDCAVsPIn310deg","p #times DCA versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);p #times DCA (GeV #times cm)",2*fNPBins,fPRange[0],fPRange[1], deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h2, kPDCAVsPIn310deg);
  h2 = new TH2F("hPDCAVsPosAbsEndIn02degMC","p #times DCA versus track position at absorber end for tracks with MC angle < 2 degrees;position (cm);p #times DCA (GeV #times cm)",1000,0.,100.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h2, kPDCAVsPosAbsEndIn02degMC);
  h2 = new TH2F("hPDCAVsPosAbsEndIn23degMC","p #times DCA}versus track position at absorber end for tracks with MC angle in [2,3[ degrees;position (cm);p #times DCA (GeV #times cm)",1000,0.,100.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h2, kPDCAVsPosAbsEndIn23degMC);
  h2 = new TH2F("hPDCAVsPosAbsEndIn310degMC","p #times DCA versus track position at absorber end for tracks with MC angle in [3,10[ degrees;position (cm);p #times DCA (GeV #times cm)",1000,0.,100.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h2, kPDCAVsPosAbsEndIn310degMC);
  h2 = new TH2F("hPDCAVsAngleAtAbsEnd","p #times DCA versus track position at absorber end converted to degrees;angle (Deg);p #times DCA (GeV #times cm)",10,0.,10.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h2, kPDCAVsAngleAtAbsEnd);
  h2 = new TH2F("hPDCAVsMCAngle","p #times DCA versus MC angle;MC angle (Deg);p #times DCA (GeV #times cm)",10,0.,10.,deltaPDCANBins, deltaPDCAEdges[0], deltaPDCAEdges[1]);
  fTrackerList->AddAt(h2, kPDCAVsMCAngle);
  
  // MCS angular dispersion
  h2 = new TH2F("hPMCSAngVsPIn23deg","p #times #Delta#theta_{MCS} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #Delta#theta_{MCS} (GeV)",2*fNPBins,fPRange[0],fPRange[1], deltaPDCANBins, deltaPMCSAngEdges[0], deltaPMCSAngEdges[1]);
  fTrackerList->AddAt(h2, kPMCSAngVsPIn23deg);
  h2 = new TH2F("hPMCSAngVsPIn310deg","p #times #Delta#theta_{MCS} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);p #times #Delta#theta_{MCS} (GeV)",2*fNPBins,fPRange[0],fPRange[1], deltaPDCANBins, deltaPMCSAngEdges[0], deltaPMCSAngEdges[1]);
  fTrackerList->AddAt(h2, kPMCSAngVsPIn310deg);
  
  // cluster resolution
  Int_t nCh = AliMUONConstants::NTrackingCh();
  const Int_t clusterResNBins = 5000;
  Double_t clusterResMaxX = 10.*fClusterMaxRes[0];
  Double_t clusterResMaxY = 10.*fClusterMaxRes[1];
  
  h2 = new TH2F("hResClXVsCh", "cluster-track residual-X distribution per chamber;chamber ID;#Delta_{X} (cm)", nCh, 0.5, nCh+0.5, clusterResNBins, -clusterResMaxX, clusterResMaxX);
  fTrackerList->AddAt(h2, kResClXVsCh);
  h2 = new TH2F("hResClYVsCh", "cluster-track residual-Y distribution per chamber;chamber ID;#Delta_{Y} (cm)", nCh, 0.5, nCh+0.5, clusterResNBins, -clusterResMaxY, clusterResMaxY);
  fTrackerList->AddAt(h2, kResClYVsCh);
  h2 = new TH2F("hResClXVsDE", "cluster-track residual-X distribution per DE;DE ID;#Delta_{X} (cm)", fNDE, 0.5, fNDE+0.5, clusterResNBins, -clusterResMaxX, clusterResMaxX);
  for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
  fTrackerList->AddAt(h2, kResClXVsDE);
  h2 = new TH2F("hResClYVsDE", "cluster-track residual-Y distribution per DE;DE ID;#Delta_{Y} (cm)", fNDE, 0.5, fNDE+0.5, clusterResNBins, -clusterResMaxY, clusterResMaxY);
  for (Int_t i = 1; i <= fNDE; i++) h2->GetXaxis()->SetBinLabel(i, Form("%d",fDEIds[i]));
  fTrackerList->AddAt(h2, kResClYVsDE);
  
  // Disable printout of AliMCEvent
  AliLog::SetClassDebugLevel("AliMCEvent",-1);
  
  PostData(1, fCFContainer);
  PostData(3, fTriggerList);
  PostData(4, fTrackerList);
}

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::UserExec(Option_t * /*option*/) 
{
  //
  /// Main loop
  /// Called for each event
  //

  // check that OCDB objects have been properly set
  if (fSigmaCut < 0) return;
  
  // Load ESD event
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if ( ! esd ) {
    AliError ("ESD event not found. Nothing done!");
    return;
  }

  // Load MC event 
  AliMCEventHandler* mcH = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if ( ! mcH ) {
    AliError ("MCH event handler not found. Nothing done!");
    return;
  }
  
  // Get reference tracks
  AliMUONRecoCheck rc(esd,mcH);
  AliMUONVTriggerTrackStore* triggerTrackRefStore = rc.TriggerableTracks(-1);
  AliMUONVTrackStore* reconstructibleStore = rc.ReconstructibleTracks(-1, fRequestedStationMask, fRequest2ChInSameSt45);
  
  Double_t containerInput[kNvars];
  AliMUONTrackParam *trackParam;
  Double_t x1,y1,z1,slopex1,slopey1,pX1,pY1,pZ1,p1,pT1,eta1,phi1;
  Double_t x2,y2,z2,slopex2,slopey2,pX2,pY2,pZ2,p2,pT2,eta2,phi2;
  Double_t dPhi;
  Double_t xAbs,yAbs,dAbs,aAbs,aMCS,aMC;
  Double_t xDCA,yDCA,dca,pU;
  
  // ------ Loop over reconstructed tracks ------
  AliESDMuonTrack *esdTrack = 0x0;
  Int_t nMuTracks = esd->GetNumberOfMuonTracks();
  for (Int_t iMuTrack = 0; iMuTrack < nMuTracks; ++iMuTrack) {
    
    esdTrack = esd->GetMuonTrack(iMuTrack);

    AliMUONTrack* matchedTrackRef = 0x0;
    AliMUONTriggerTrack* matchedTrigTrackRef = 0x0;
    containerInput[kVarMatchMC] = static_cast<Double_t>(kNoMatch);
    
    // Tracker
    if (esdTrack->ContainTrackerData()) {
    
      // convert ESD track to MUON track (without recomputing track parameters at each clusters)
      AliMUONTrack muonTrack;
      AliMUONESDInterface::ESDToMUON(*esdTrack, muonTrack, kFALSE);
      
      // try to match the reconstructed track with a simulated one
      Int_t nMatchClusters = 0;
      matchedTrackRef = rc.FindCompatibleTrack(muonTrack, *reconstructibleStore, nMatchClusters, kFALSE, fSigmaCut);
      if (matchedTrackRef) {
	
	containerInput[kVarMatchMC] = static_cast<Double_t>(kTrackerOnly);
	
	// simulated track parameters at vertex
        trackParam = matchedTrackRef->GetTrackParamAtVertex();
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
	
	// reconstructed track parameters at the end of the absorber
        AliMUONTrackParam trackParamAtAbsEnd(*((AliMUONTrackParam*)muonTrack.GetTrackParamAtCluster()->First()));
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
	
	// reconstructed track parameters at vertex
	trackParam = muonTrack.GetTrackParamAtVertex();
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
        
	// reconstructed track parameters at DCA
        AliMUONTrackParam trackParamAtDCA(*((AliMUONTrackParam*) muonTrack.GetTrackParamAtCluster()->First()));
	pU = trackParamAtDCA.P();
	AliMUONTrackExtrap::ExtrapToVertexWithoutBranson(&trackParamAtDCA, z2);
        xDCA = trackParamAtDCA.GetNonBendingCoor();
        yDCA = trackParamAtDCA.GetBendingCoor();
	dca = TMath::Sqrt(xDCA*xDCA + yDCA*yDCA);
	
	dPhi = phi2-phi1;
	if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();
	else if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
	
	// fill histograms
	static_cast<TH1*>(fTrackerList->At(kResPAtVtx))->Fill(p2-p1);
	static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsP))->Fill(p1,p2-p1);
	static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsAngleAtAbsEnd))->Fill(aAbs,p2-p1);
	static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsMCAngle))->Fill(aMC,p2-p1);
	static_cast<TH3*>(fTrackerList->At(kResPAtVtxVsAngleAtAbsEndVsP))->Fill(p1,aAbs,p2-p1);
	static_cast<TH2*>(fTrackerList->At(kResPtAtVtxVsPt))->Fill(pT1,pT2-pT1);
	
	static_cast<TH1*>(fTrackerList->At(kResSlopeXAtVtx))->Fill(slopex2-slopex1);
	static_cast<TH1*>(fTrackerList->At(kResSlopeYAtVtx))->Fill(slopey2-slopey1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeXAtVtxVsP))->Fill(p1,slopex2-slopex1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeYAtVtxVsP))->Fill(p1,slopey2-slopey1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeXAtVtxVsAngleAtAbsEnd))->Fill(aAbs,slopex2-slopex1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeYAtVtxVsAngleAtAbsEnd))->Fill(aAbs,slopey2-slopey1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeXAtVtxVsMCAngle))->Fill(aMC,slopex2-slopex1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeYAtVtxVsMCAngle))->Fill(aMC,slopey2-slopey1);
	
	static_cast<TH1*>(fTrackerList->At(kResEtaAtVtx))->Fill(eta2-eta1);
	static_cast<TH2*>(fTrackerList->At(kResEtaAtVtxVsP))->Fill(p1,eta2-eta1);
	static_cast<TH2*>(fTrackerList->At(kResEtaAtVtxVsAngleAtAbsEnd))->Fill(aAbs,eta2-eta1);
	static_cast<TH2*>(fTrackerList->At(kResEtaAtVtxVsMCAngle))->Fill(aMC,eta2-eta1);
	
	static_cast<TH1*>(fTrackerList->At(kResPhiAtVtx))->Fill(dPhi);
	static_cast<TH2*>(fTrackerList->At(kResPhiAtVtxVsP))->Fill(p1,dPhi);
	static_cast<TH2*>(fTrackerList->At(kResPhiAtVtxVsAngleAtAbsEnd))->Fill(aAbs,dPhi);
	static_cast<TH2*>(fTrackerList->At(kResPhiAtVtxVsMCAngle))->Fill(aMC,dPhi);
	
	static_cast<TH1*>(fTrackerList->At(kPDCA))->Fill(0.5*(p2+pU)*dca);
	static_cast<TH2*>(fTrackerList->At(kPDCAVsAngleAtAbsEnd))->Fill(aAbs,0.5*(p2+pU)*dca);
	static_cast<TH2*>(fTrackerList->At(kPDCAVsMCAngle))->Fill(aMC,0.5*(p2+pU)*dca);
	
	if (aAbs > 2. && aAbs < 3.) {
	  
	  static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsPIn23deg))->Fill(p1,p2-p1);
	  static_cast<TH2*>(fTrackerList->At(kPDCAVsPIn23deg))->Fill(p1,0.5*(p2+pU)*dca);
	  static_cast<TH2*>(fTrackerList->At(kPMCSAngVsPIn23deg))->Fill(p1,0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad());
	  
	} else if (aAbs >= 3. && aAbs < 10.) {
	  
	  static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsPIn310deg))->Fill(p1,p2-p1);
	  static_cast<TH2*>(fTrackerList->At(kPDCAVsPIn310deg))->Fill(p1,0.5*(p2+pU)*dca);
	  static_cast<TH2*>(fTrackerList->At(kPMCSAngVsPIn310deg))->Fill(p1,0.5*(p2+pU)*(aMCS-aMC)*TMath::DegToRad());
	  
	}
	
	if (aMC < 2.) {
	  
	  static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsPIn02degMC))->Fill(p1,p2-p1);
	  static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsPosAbsEndIn02degMC))->Fill(dAbs,p2-p1);
	  static_cast<TH2*>(fTrackerList->At(kResSlopeXAtVtxVsPosAbsEndIn02degMC))->Fill(dAbs,slopex2-slopex1);
	  static_cast<TH2*>(fTrackerList->At(kResSlopeYAtVtxVsPosAbsEndIn02degMC))->Fill(dAbs,slopey2-slopey1);
	  static_cast<TH2*>(fTrackerList->At(kResEtaAtVtxVsPosAbsEndIn02degMC))->Fill(dAbs,eta2-eta1);
	  static_cast<TH2*>(fTrackerList->At(kResPhiAtVtxVsPosAbsEndIn02degMC))->Fill(dAbs,dPhi);
	  static_cast<TH2*>(fTrackerList->At(kPDCAVsPosAbsEndIn02degMC))->Fill(dAbs,0.5*(p2+pU)*dca);
	  
	} else if (aMC >= 2. && aMC < 3) {
	  
	  static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsPosAbsEndIn23degMC))->Fill(dAbs,p2-p1);
	  static_cast<TH2*>(fTrackerList->At(kResSlopeXAtVtxVsPosAbsEndIn23degMC))->Fill(dAbs,slopex2-slopex1);
	  static_cast<TH2*>(fTrackerList->At(kResSlopeYAtVtxVsPosAbsEndIn23degMC))->Fill(dAbs,slopey2-slopey1);
	  static_cast<TH2*>(fTrackerList->At(kResEtaAtVtxVsPosAbsEndIn23degMC))->Fill(dAbs,eta2-eta1);
	  static_cast<TH2*>(fTrackerList->At(kResPhiAtVtxVsPosAbsEndIn23degMC))->Fill(dAbs,dPhi);
	  static_cast<TH2*>(fTrackerList->At(kPDCAVsPosAbsEndIn23degMC))->Fill(dAbs,0.5*(p2+pU)*dca);
	  
	} else if (aMC >= 3. && aMC < 10.) {
	  
	  static_cast<TH2*>(fTrackerList->At(kResPAtVtxVsPosAbsEndIn310degMC))->Fill(dAbs,p2-p1);
	  static_cast<TH2*>(fTrackerList->At(kResSlopeXAtVtxVsPosAbsEndIn310degMC))->Fill(dAbs,slopex2-slopex1);
	  static_cast<TH2*>(fTrackerList->At(kResSlopeYAtVtxVsPosAbsEndIn310degMC))->Fill(dAbs,slopey2-slopey1);
	  static_cast<TH2*>(fTrackerList->At(kResEtaAtVtxVsPosAbsEndIn310degMC))->Fill(dAbs,eta2-eta1);
	  static_cast<TH2*>(fTrackerList->At(kResPhiAtVtxVsPosAbsEndIn310degMC))->Fill(dAbs,dPhi);
	  static_cast<TH2*>(fTrackerList->At(kPDCAVsPosAbsEndIn310degMC))->Fill(dAbs,0.5*(p2+pU)*dca);
	  
	}
	
	// simulated track parameters at vertex
        trackParam = (AliMUONTrackParam*) matchedTrackRef->GetTrackParamAtCluster()->First();
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
	
	// reconstructed track parameters at vertex
        trackParam = (AliMUONTrackParam*) muonTrack.GetTrackParamAtCluster()->First();
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
        
	// fill histograms
	static_cast<TH1*>(fTrackerList->At(kResPAt1stCl))->Fill(p2-p1);
	static_cast<TH2*>(fTrackerList->At(kResPAt1stClVsP))->Fill(p1,p2-p1);
	static_cast<TH2*>(fTrackerList->At(kResPtAt1stClVsPt))->Fill(pT1,pT2-pT1);
	static_cast<TH1*>(fTrackerList->At(kResSlopeXAt1stCl))->Fill(slopex2-slopex1);
	static_cast<TH1*>(fTrackerList->At(kResSlopeYAt1stCl))->Fill(slopey2-slopey1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeXAt1stClVsP))->Fill(p1,slopex2-slopex1);
	static_cast<TH2*>(fTrackerList->At(kResSlopeYAt1stClVsP))->Fill(p1,slopey2-slopey1);
	
	// clusters residuals
	for (Int_t iCl1 = 0; iCl1 < muonTrack.GetNClusters(); iCl1++) {
	  
	  AliMUONVCluster* cluster1 = static_cast<AliMUONTrackParam*>(muonTrack.GetTrackParamAtCluster()->UncheckedAt(iCl1))->GetClusterPtr();
	  Int_t chId = cluster1->GetChamberId();
	  Int_t deId = cluster1->GetDetElemId();
	  
	  for (Int_t iCl2 = 0; iCl2 < matchedTrackRef->GetNClusters(); iCl2++) {
	    
	    AliMUONVCluster* cluster2 = static_cast<AliMUONTrackParam*>(matchedTrackRef->GetTrackParamAtCluster()->UncheckedAt(iCl2))->GetClusterPtr();
	    
	    if (cluster2->GetDetElemId() == deId) {
	      
	      // fill histograms
	      static_cast<TH2*>(fTrackerList->At(kResClXVsCh))->Fill(chId+1, cluster1->GetX()-cluster2->GetX());
	      static_cast<TH2*>(fTrackerList->At(kResClYVsCh))->Fill(chId+1, cluster1->GetY()-cluster2->GetY());
	      static_cast<TH2*>(fTrackerList->At(kResClXVsDE))->Fill(fDEIndices[deId], cluster1->GetX()-cluster2->GetX());
	      static_cast<TH2*>(fTrackerList->At(kResClYVsDE))->Fill(fDEIndices[deId], cluster1->GetY()-cluster2->GetY());
	      
	      break;
	      
	    }
	    
	  }
	  
	}
	
      } // end if (matchedTrackRef)
      
    } // end if (esdTrack->ContainTrackerData())
    
    // Trigger
    if (esdTrack->ContainTriggerData()) {

      // Convert ESD track to trigger track
      AliMUONLocalTrigger locTrg;
      AliMUONESDInterface::ESDToMUON(*esdTrack, locTrg);
      AliMUONTriggerTrack trigTrack;
      rc.TriggerToTrack(locTrg, trigTrack);
      
      // try to match the trigger track with the same MC track if any
      if (matchedTrackRef) {
	AliMUONTriggerTrack* triggerTrackRef = static_cast<AliMUONTriggerTrack*>(triggerTrackRefStore->FindObject(matchedTrackRef->GetUniqueID()));
        if (triggerTrackRef && trigTrack.Match(*triggerTrackRef, fSigmaCutTrig)) matchedTrigTrackRef = triggerTrackRef;
	if (matchedTrigTrackRef) containerInput[kVarMatchMC] = static_cast<Double_t>(kMatchedSame);
      }
      
      // or try to match with any triggerable track
      if (!matchedTrigTrackRef) {
	matchedTrigTrackRef = rc.FindCompatibleTrack(trigTrack, *triggerTrackRefStore, fSigmaCutTrig);
	if (matchedTrigTrackRef) containerInput[kVarMatchMC] = (matchedTrackRef) ? static_cast<Double_t>(kMatchedDiff) : static_cast<Double_t>(kTriggerOnly);
      }
      
      // fill histograms
      if (matchedTrigTrackRef) {
	static_cast<TH1*>(fTriggerList->At(kResTrigX11))->Fill(trigTrack.GetX11() - matchedTrigTrackRef->GetX11());
	static_cast<TH1*>(fTriggerList->At(kResTrigY11))->Fill(trigTrack.GetY11() - matchedTrigTrackRef->GetY11());
	static_cast<TH1*>(fTriggerList->At(kResTrigSlopeY))->Fill(trigTrack.GetSlopeY() - matchedTrigTrackRef->GetSlopeY());
      }
      
    }
    
    // get MC particle ID
    Int_t mcID = -1;
    if (matchedTrackRef) mcID = static_cast<Int_t>(matchedTrackRef->GetUniqueID());
    else if (matchedTrigTrackRef) mcID = static_cast<Int_t>(matchedTrigTrackRef->GetUniqueID());
    
    // fill particle container
    FillContainerInfo(containerInput, esdTrack, mcID);
    fCFContainer->Fill(containerInput, kStepReconstructed);

  }

  // ------ Loop over reconstructible tracks ------
  AliMUONTrack* trackRef = 0x0;
  TIter next(reconstructibleStore->CreateIterator());
  while ((trackRef = static_cast<AliMUONTrack*>(next()))) {
    
    // find the corresponding triggerable track if any
    AliMUONTriggerTrack* trigTrackRef = static_cast<AliMUONTriggerTrack*>(triggerTrackRefStore->FindObject(trackRef->GetUniqueID()));
    
    // fill particle container
    FillContainerInfo(containerInput, 0x0, static_cast<Int_t>(trackRef->GetUniqueID()));
    containerInput[kVarHasTracker] = 1.;
    if (trigTrackRef) containerInput[kVarTrigger] = static_cast<Double_t>(kAllPtTrig);
    containerInput[kVarMatchMC] = static_cast<Double_t>(kNoMatch);
    fCFContainer->Fill(containerInput, kStepGeneratedMC);
    
  }

  // ------ Loop over triggerable ghosts ------
  AliMUONTriggerTrack* trigTrackRef = 0x0;
  TIter nextTrig(triggerTrackRefStore->CreateIterator());
  while ((trigTrackRef = static_cast<AliMUONTriggerTrack*>(nextTrig()))) {
    
    // discard tracks also reconstructible
    if (reconstructibleStore->FindObject(trigTrackRef->GetUniqueID())) continue;

    // fill particle container
    FillContainerInfo(containerInput, 0x0, static_cast<Int_t>(trigTrackRef->GetUniqueID()));
    containerInput[kVarTrigger] = static_cast<Double_t>(kAllPtTrig);
    containerInput[kVarMatchMC] = static_cast<Double_t>(kNoMatch);
    fCFContainer->Fill(containerInput, kStepGeneratedMC);
    
  }
  
  // Post final data
  PostData(1, fCFContainer);
  PostData(3, fTriggerList);
  PostData(4, fTrackerList);
}

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::Terminate(Option_t *)
{
  //
  /// Draw some histogram at the end.
  //

  // do it only locally
  //if (gROOT->IsBatch()) return;
  
  // get output containers
  fCFContainer = dynamic_cast<AliCFContainer*>(GetOutputData(1));
  fTriggerList = dynamic_cast<TObjArray*>(GetOutputData(3));
  fTrackerList = dynamic_cast<TObjArray*>(GetOutputData(4));
  if (!fCFContainer || !fTriggerList || !fTrackerList) {
    AliWarning("Output containers not found: summary histograms are not created"); 
    return;
  }
  
  // --- compute efficiencies ---
  fEfficiencyList = new TObjArray(100);
  fEfficiencyList->SetOwner();
  
  AliCFEffGrid* efficiency = new AliCFEffGrid("eff","",*fCFContainer);
  efficiency->CalculateEfficiency(kStepReconstructed, kStepGeneratedMC);
  Double_t totalEff = 0., totalEffErr = 0.;
  Int_t sumEffBin = 0;
  TH1* auxHisto = 0x0;
  
  // add histogram summarizing global efficiencies
  TH1D* effSummary = new TH1D("effSummary", "Efficiency summary", 1, 0., 0.);
  effSummary->GetYaxis()->SetTitle("Efficiency");
  fEfficiencyList->AddLast(effSummary);
  
  // Tracker efficiency using all reconstructed tracks
  efficiency->GetNum()->SetRangeUser(kVarHasTracker, 1., 1.);
  efficiency->GetDen()->SetRangeUser(kVarHasTracker, 1., 1.);
  auxHisto = efficiency->Project(kVarPt);
  auxHisto->SetName("effVsPt_trackerTracks");
  fEfficiencyList->AddLast(auxHisto);
  auxHisto = efficiency->Project(kVarEta);
  auxHisto->SetName("effVsEta_trackerTracks");
  fEfficiencyList->AddLast(auxHisto);
  GetEfficiency(efficiency, totalEff, totalEffErr);
  sumEffBin++;
  effSummary->Fill("Tracker_all", totalEff);
  effSummary->SetBinError(sumEffBin, totalEffErr);
  printf("Tracker efficiency using all reconstructed tracks = %f +- %f\n", totalEff, totalEffErr);
  
  // Tracker efficiency using tracks matched with reconstructible ones
  efficiency->GetNum()->SetRangeUser(kVarMatchMC, kTrackerOnly, kMatchedDiff);
  auxHisto = efficiency->Project(kVarPt);
  auxHisto->SetName("effVsPt_trackerTracksMatchMC");
  fEfficiencyList->AddLast(auxHisto);
  auxHisto = efficiency->Project(kVarEta);
  auxHisto->SetName("effVsEta_trackerTracksMatchMC");
  fEfficiencyList->AddLast(auxHisto);
  GetEfficiency(efficiency, totalEff, totalEffErr);
  sumEffBin++;
  effSummary->Fill("Tracker_MCId", totalEff);
  effSummary->SetBinError(sumEffBin, totalEffErr);
  printf("Tracker efficiency using reconstructed tracks matching MC = %f +- %f\n", totalEff, totalEffErr);
  
  // Matched efficiency using all reconstructed tracks
  efficiency->GetNum()->SetRangeUser(kVarTrigger, kAllPtTrig, kHighPtTrig);
  efficiency->GetDen()->SetRangeUser(kVarTrigger, kAllPtTrig, kAllPtTrig);
  efficiency->GetNum()->GetAxis(kVarMatchMC)->SetRange();
  auxHisto = efficiency->Project(kVarPt);
  auxHisto->SetName("effVsPt_matchedTracks");
  fEfficiencyList->AddLast(auxHisto);
  auxHisto = efficiency->Project(kVarEta);
  auxHisto->SetName("effVsEta_matchedTracks");
  fEfficiencyList->AddLast(auxHisto);
  GetEfficiency(efficiency, totalEff, totalEffErr);
  sumEffBin++;
  effSummary->Fill("Matched_all", totalEff);
  effSummary->SetBinError(sumEffBin, totalEffErr);
  printf("Matched efficiency using all reconstructed tracks = %f +- %f\n", totalEff, totalEffErr);
  
  // Matched efficiency using tracks matched with reconstructible & triggerable ones
  efficiency->GetNum()->SetRangeUser(kVarMatchMC, kMatchedSame, kMatchedSame);
  auxHisto = efficiency->Project(kVarPt);
  auxHisto->SetName("effVsPt_matchedTracksMatchMC");
  fEfficiencyList->AddLast(auxHisto);
  auxHisto = efficiency->Project(kVarEta);
  auxHisto->SetName("effVsEta_matchedTracksMatchMC");
  fEfficiencyList->AddLast(auxHisto);
  GetEfficiency(efficiency, totalEff, totalEffErr);
  sumEffBin++;
  effSummary->Fill("Matched_MCId", totalEff);
  effSummary->SetBinError(sumEffBin, totalEffErr);
  printf("Matched efficiency using reconstructed tracks matching MC = %f +- %f\n", totalEff, totalEffErr);
  
  // Trigger efficiency using all reconstructed tracks
  efficiency->GetNum()->GetAxis(kVarHasTracker)->SetRange();
  efficiency->GetDen()->GetAxis(kVarHasTracker)->SetRange();
  efficiency->GetNum()->GetAxis(kVarMatchMC)->SetRange();
  auxHisto = efficiency->Project(kVarPt);
  auxHisto->SetName("effVsPt_triggerTracks");
  fEfficiencyList->AddLast(auxHisto);
  auxHisto = efficiency->Project(kVarEta);
  auxHisto->SetName("effVsEta_triggerTracks");
  fEfficiencyList->AddLast(auxHisto);
  GetEfficiency(efficiency, totalEff, totalEffErr);
  sumEffBin++;
  effSummary->Fill("Trigger_all", totalEff);
  effSummary->SetBinError(sumEffBin, totalEffErr);
  printf("Trigger efficiency using all reconstructed tracks = %f +- %f\n", totalEff, totalEffErr);
  
  // Trigger efficiency using tracks matched with triggerable ones
  efficiency->GetNum()->SetRangeUser(kVarMatchMC, kMatchedSame, kTriggerOnly);
  auxHisto = efficiency->Project(kVarPt);
  auxHisto->SetName("effVsPt_triggerTracksMatchMC");
  fEfficiencyList->AddLast(auxHisto);
  auxHisto = efficiency->Project(kVarEta);
  auxHisto->SetName("effVsEta_triggerTracksMatchMC");
  fEfficiencyList->AddLast(auxHisto);
  GetEfficiency(efficiency, totalEff, totalEffErr);
  sumEffBin++;
  effSummary->Fill("Trigger_MCId", totalEff);
  effSummary->SetBinError(sumEffBin, totalEffErr);
  printf("Trigger efficiency using reconstructed tracks matching MC = %f +- %f\n", totalEff, totalEffErr);
  
  // plot histogram summarizing global efficiencies
  TCanvas* cEffSummary = new TCanvas("cEffSummary","Efficiency summary",20,20,310,310);
  cEffSummary->SetFillColor(10); cEffSummary->SetHighLightColor(10);
  cEffSummary->SetLeftMargin(0.15); cEffSummary->SetBottomMargin(0.15);
  effSummary->DrawCopy("etext");
  
  // --- plot trigger resolution ---
  TCanvas* cTriggerResolution = new TCanvas("cTriggerResolution","Trigger resolution",10,10,310,310);
  cTriggerResolution->SetFillColor(10); cTriggerResolution->SetHighLightColor(10);
  cTriggerResolution->SetLeftMargin(0.15); cTriggerResolution->SetBottomMargin(0.15);
  cTriggerResolution->Divide(2,2);
  cTriggerResolution->cd(1);
  static_cast<TH1*>(fTriggerList->UncheckedAt(kResTrigX11))->DrawCopy();
  cTriggerResolution->cd(2);
  static_cast<TH1*>(fTriggerList->UncheckedAt(kResTrigY11))->DrawCopy();
  cTriggerResolution->cd(3);
  static_cast<TH1*>(fTriggerList->UncheckedAt(kResTrigSlopeY))->DrawCopy();
  
  // --- compute momentum resolution at vertex ---
  fPAtVtxList = new TObjArray(100);
  fPAtVtxList->SetOwner();
  
  // define graphs
  TGraphAsymmErrors* gMeanResPAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResPAtVtxVsP->SetName("gMeanResPAtVtxVsP");
  gMeanResPAtVtxVsP->SetTitle("<#Delta_{p}> at vertex versus p;p (GeV/c);<#Delta_{p}> (GeV/c)");
  fPAtVtxList->AddAt(gMeanResPAtVtxVsP, kMeanResPAtVtxVsP);
  TGraphAsymmErrors* gMostProbResPAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gMostProbResPAtVtxVsP->SetName("gMostProbResPAtVtxVsP");
  gMostProbResPAtVtxVsP->SetTitle("Most probable #Delta_{p} at vertex versus p;p (GeV/c);Most prob. #Delta_{p} (GeV/c)");
  fPAtVtxList->AddAt(gMostProbResPAtVtxVsP, kMostProbResPAtVtxVsP);
  TGraphAsymmErrors* gSigmaResPAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResPAtVtxVsP->SetName("gSigmaResPAtVtxVsP");
  gSigmaResPAtVtxVsP->SetTitle("#sigma_{p}/p at vertex versus p;p (GeV/c);#sigma_{p}/p (%)");
  fPAtVtxList->AddAt(gSigmaResPAtVtxVsP, kSigmaResPAtVtxVsP);
  
  // fit histo and fill graphs
  TH2* h = static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsP));
  FitLandauGausResVsP(h, "momentum residuals at vertex", gMeanResPAtVtxVsP, gMostProbResPAtVtxVsP, gSigmaResPAtVtxVsP);
  
  // convert resolution into relative resolution
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/fNPBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    Double_t x,y;
    gSigmaResPAtVtxVsP->GetPoint(i/rebinFactorX-1, x, y);
    gSigmaResPAtVtxVsP->SetPoint(i/rebinFactorX-1, x, 100.*y/x);
    gSigmaResPAtVtxVsP->SetPointEYlow(i/rebinFactorX-1, 100.*gSigmaResPAtVtxVsP->GetErrorYlow(i/rebinFactorX-1)/x);
    gSigmaResPAtVtxVsP->SetPointEYhigh(i/rebinFactorX-1, 100.*gSigmaResPAtVtxVsP->GetErrorYhigh(i/rebinFactorX-1)/x);
  }
  
  // --- compute momentum resolution at first cluster ---
  fPAt1stClList = new TObjArray(100);
  fPAt1stClList->SetOwner();
  
  // define graphs
  TGraphAsymmErrors* gMeanResPAt1stClVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResPAt1stClVsP->SetName("gMeanResPAt1stClVsP");
  gMeanResPAt1stClVsP->SetTitle("<#Delta_{p}> at first cluster versus p;p (GeV/c);<#Delta_{p}> (GeV/c)");
  fPAt1stClList->AddAt(gMeanResPAt1stClVsP, kMeanResPAt1stClVsP);
  TGraphAsymmErrors* gSigmaResPAt1stClVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResPAt1stClVsP->SetName("gSigmaResPAt1stClVsP");
  gSigmaResPAt1stClVsP->SetTitle("#sigma_{p}/p at first cluster versus p;p (GeV/c);#sigma_{p}/p (%)");
  fPAt1stClList->AddAt(gSigmaResPAt1stClVsP, kSigmaResPAt1stClVsP);
  
  // fit histo and fill graphs
  h = static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAt1stClVsP));
  FitGausResVsMom(h, 0., 1., "momentum residuals at first cluster", gMeanResPAt1stClVsP, gSigmaResPAt1stClVsP);
  
  // convert resolution into relative resolution
  rebinFactorX = TMath::Max(h->GetNbinsX()/fNPBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    Double_t x,y;
    gSigmaResPAt1stClVsP->GetPoint(i/rebinFactorX-1, x, y);
    gSigmaResPAt1stClVsP->SetPoint(i/rebinFactorX-1, x, 100.*y/x);
    gSigmaResPAt1stClVsP->SetPointEYlow(i/rebinFactorX-1, 100.*gSigmaResPAt1stClVsP->GetErrorYlow(i/rebinFactorX-1)/x);
    gSigmaResPAt1stClVsP->SetPointEYhigh(i/rebinFactorX-1, 100.*gSigmaResPAt1stClVsP->GetErrorYhigh(i/rebinFactorX-1)/x);
  }
  
  // --- compute slope resolution at vertex ---
  fSlopeAtVtxList = new TObjArray(100);
  fSlopeAtVtxList->SetOwner();
  
  // define graphs
  TGraphAsymmErrors* gMeanResSlopeXAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResSlopeXAtVtxVsP->SetName("gMeanResSlopeXAtVtxVsP");
  gMeanResSlopeXAtVtxVsP->SetTitle("<#Delta_{slope_{X}}> at vertex versus p;p (GeV/c);<#Delta_{slope_{X}}>");
  fSlopeAtVtxList->AddAt(gMeanResSlopeXAtVtxVsP, kMeanResSlopeXAtVtxVsP);
  TGraphAsymmErrors* gMeanResSlopeYAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResSlopeYAtVtxVsP->SetName("gMeanResSlopeYAtVtxVsP");
  gMeanResSlopeYAtVtxVsP->SetTitle("<#Delta_{slope_{Y}}> at vertex versus p;p (GeV/c);<#Delta_{slope_{Y}}>");
  fSlopeAtVtxList->AddAt(gMeanResSlopeYAtVtxVsP, kMeanResSlopeYAtVtxVsP);
  TGraphAsymmErrors* gSigmaResSlopeXAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResSlopeXAtVtxVsP->SetName("gSigmaResSlopeXAtVtxVsP");
  gSigmaResSlopeXAtVtxVsP->SetTitle("#sigma_{slope_{X}} at vertex versus p;p (GeV/c);#sigma_{slope_{X}}");
  fSlopeAtVtxList->AddAt(gSigmaResSlopeXAtVtxVsP, kSigmaResSlopeXAtVtxVsP);
  TGraphAsymmErrors* gSigmaResSlopeYAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResSlopeYAtVtxVsP->SetName("gSigmaResSlopeYAtVtxVsP");
  gSigmaResSlopeYAtVtxVsP->SetTitle("#sigma_{slope_{Y}} at vertex versus p;p (GeV/c);#sigma_{slope_{Y}}");
  fSlopeAtVtxList->AddAt(gSigmaResSlopeYAtVtxVsP, kSigmaResSlopeYAtVtxVsP);
  
  // fit histo and fill graphs
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeXAtVtxVsP)), 0., 2.e-3,
		  "slopeX residuals at vertex", gMeanResSlopeXAtVtxVsP, gSigmaResSlopeXAtVtxVsP);
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeYAtVtxVsP)), 0., 2.e-3,
		  "slopeY residuals at vertex", gMeanResSlopeYAtVtxVsP, gSigmaResSlopeYAtVtxVsP);
  
  // --- compute slope resolution at first cluster ---
  fSlopeAt1stClList = new TObjArray(100);
  fSlopeAt1stClList->SetOwner();
  
  // define graphs
  TGraphAsymmErrors* gMeanResSlopeXAt1stClVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResSlopeXAt1stClVsP->SetName("gMeanResSlopeXAt1stClVsP");
  gMeanResSlopeXAt1stClVsP->SetTitle("<#Delta_{slope_{X}}> at first cluster versus p;p (GeV/c);<#Delta_{slope_{X}}>");
  fSlopeAt1stClList->AddAt(gMeanResSlopeXAt1stClVsP, kMeanResSlopeXAt1stClVsP);
  TGraphAsymmErrors* gMeanResSlopeYAt1stClVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResSlopeYAt1stClVsP->SetName("gMeanResSlopeYAt1stClVsP");
  gMeanResSlopeYAt1stClVsP->SetTitle("<#Delta_{slope_{Y}}> at first cluster versus p;p (GeV/c);<#Delta_{slope_{Y}}>");
  fSlopeAt1stClList->AddAt(gMeanResSlopeYAt1stClVsP, kMeanResSlopeYAt1stClVsP);
  TGraphAsymmErrors* gSigmaResSlopeXAt1stClVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResSlopeXAt1stClVsP->SetName("gSigmaResSlopeXAt1stClVsP");
  gSigmaResSlopeXAt1stClVsP->SetTitle("#sigma_{slope_{X}} at first cluster versus p;p (GeV/c);#sigma_{slope_{X}}");
  fSlopeAt1stClList->AddAt(gSigmaResSlopeXAt1stClVsP, kSigmaResSlopeXAt1stClVsP);
  TGraphAsymmErrors* gSigmaResSlopeYAt1stClVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResSlopeYAt1stClVsP->SetName("gSigmaResSlopeYAt1stClVsP");
  gSigmaResSlopeYAt1stClVsP->SetTitle("#sigma_{slope_{Y}} at first cluster versus p;p (GeV/c);#sigma_{slope_{Y}}");
  fSlopeAt1stClList->AddAt(gSigmaResSlopeYAt1stClVsP, kSigmaResSlopeYAt1stClVsP);
  
  // fit histo and fill graphs
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeXAt1stClVsP)), 0., 3.e-4,
		  "slopeX residuals at first cluster", gMeanResSlopeXAt1stClVsP, gSigmaResSlopeXAt1stClVsP);
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeYAt1stClVsP)), 0., 3.e-4,
		  "slopeY residuals at first cluster", gMeanResSlopeYAt1stClVsP, gSigmaResSlopeYAt1stClVsP);
  
  // --- compute eta resolution at vertex ---
  fEtaAtVtxList = new TObjArray(100);
  fEtaAtVtxList->SetOwner();
  
  // define graphs
  TGraphAsymmErrors* gMeanResEtaAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResEtaAtVtxVsP->SetName("gMeanResEtaAtVtxVsP");
  gMeanResEtaAtVtxVsP->SetTitle("<#Delta_{eta}> at vertex versus p;p (GeV/c);<#Delta_{eta}>");
  fEtaAtVtxList->AddAt(gMeanResEtaAtVtxVsP, kMeanResEtaAtVtxVsP);
  TGraphAsymmErrors* gSigmaResEtaAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResEtaAtVtxVsP->SetName("gSigmaResEtaAtVtxVsP");
  gSigmaResEtaAtVtxVsP->SetTitle("#sigma_{eta} at vertex versus p;p (GeV/c);#sigma_{eta}");
  fEtaAtVtxList->AddAt(gSigmaResEtaAtVtxVsP, kSigmaResEtaAtVtxVsP);
  
  // fit histo and fill graphs
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kResEtaAtVtxVsP)), 0., 0.1,
		  "eta residuals at vertex", gMeanResEtaAtVtxVsP, gSigmaResEtaAtVtxVsP);
  
  // --- compute phi resolution at vertex ---
  fPhiAtVtxList = new TObjArray(100);
  fPhiAtVtxList->SetOwner();
  
  // define graphs
  TGraphAsymmErrors* gMeanResPhiAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gMeanResPhiAtVtxVsP->SetName("gMeanResPhiAtVtxVsP");
  gMeanResPhiAtVtxVsP->SetTitle("<#Delta_{phi}> at vertex versus p;p (GeV/c);<#Delta_{phi}>");
  fPhiAtVtxList->AddAt(gMeanResPhiAtVtxVsP, kMeanResPhiAtVtxVsP);
  TGraphAsymmErrors* gSigmaResPhiAtVtxVsP = new TGraphAsymmErrors(fNPBins);
  gSigmaResPhiAtVtxVsP->SetName("gSigmaResPhiAtVtxVsP");
  gSigmaResPhiAtVtxVsP->SetTitle("#sigma_{phi} at vertex versus p;p (GeV/c);#sigma_{phi}");
  fPhiAtVtxList->AddAt(gSigmaResPhiAtVtxVsP, kSigmaResPhiAtVtxVsP);
  
  // fit histo and fill graphs
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kResPhiAtVtxVsP)), 0., 0.01,
		  "phi residuals at vertex", gMeanResPhiAtVtxVsP, gSigmaResPhiAtVtxVsP);
  
  // --- compute DCA resolution and MCS dispersion ---
  fDCAList = new TObjArray(100);
  fDCAList->SetOwner();
  
  // define graphs
  TGraphAsymmErrors* gMeanPDCAVsPIn23deg = new TGraphAsymmErrors(fNPBins);
  gMeanPDCAVsPIn23deg->SetName("gMeanPDCAVsPIn23deg");
  gMeanPDCAVsPIn23deg->SetTitle("<p #times DCA> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times DCA> (GeV #times cm)");
  fDCAList->AddAt(gMeanPDCAVsPIn23deg, kMeanPDCAVsPIn23deg);
  TGraphAsymmErrors* gSigmaPDCAVsPIn23deg = new TGraphAsymmErrors(fNPBins);
  gSigmaPDCAVsPIn23deg->SetName("gSigmaPDCAVsPIn23deg");
  gSigmaPDCAVsPIn23deg->SetTitle("#sigma_{p #times DCA} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCA} (GeV #times cm)");
  fDCAList->AddAt(gSigmaPDCAVsPIn23deg, kSigmaPDCAVsPIn23deg);
  TGraphAsymmErrors* gMeanPDCAVsPIn310deg = new TGraphAsymmErrors(fNPBins);
  gMeanPDCAVsPIn310deg->SetName("gMeanPDCAVsPIn310deg");
  gMeanPDCAVsPIn310deg->SetTitle("<p #times DCA> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times DCA> (GeV #times cm)");
  fDCAList->AddAt(gMeanPDCAVsPIn310deg, kMeanPDCAVsPIn310deg);
  TGraphAsymmErrors* gSigmaPDCAVsPIn310deg = new TGraphAsymmErrors(fNPBins);
  gSigmaPDCAVsPIn310deg->SetName("gSigmaPDCAVsPIn310deg");
  gSigmaPDCAVsPIn310deg->SetTitle("#sigma_{p #times DCA} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times DCA} (GeV #times cm)");
  fDCAList->AddAt(gSigmaPDCAVsPIn310deg, kSigmaPDCAVsPIn310deg);
  
  TGraphAsymmErrors* gMeanPMCSAngVsPIn23deg = new TGraphAsymmErrors(fNPBins);
  gMeanPMCSAngVsPIn23deg->SetName("gMeanPMCSAngVsPIn23deg");
  gMeanPMCSAngVsPIn23deg->SetTitle("<p #times #Delta#theta_{MCS}> versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);<p #times #Delta#theta_{MCS}> (GeV)");
  fDCAList->AddAt(gMeanPMCSAngVsPIn23deg, kMeanPMCSAngVsPIn23deg);
  TGraphAsymmErrors* gSigmaPMCSAngVsPIn23deg = new TGraphAsymmErrors(fNPBins);
  gSigmaPMCSAngVsPIn23deg->SetName("gSigmaPMCSAngVsPIn23deg");
  gSigmaPMCSAngVsPIn23deg->SetTitle("#sigma_{p #times #Delta#theta_{MCS}} versus p for tracks within [2,3[ degrees at absorber end;p (GeV/c);#sigma_{p #times #Delta#theta_{MCS}} (GeV)");
  fDCAList->AddAt(gSigmaPMCSAngVsPIn23deg, kSigmaPMCSAngVsPIn23deg);
  TGraphAsymmErrors* gMeanPMCSAngVsPIn310deg = new TGraphAsymmErrors(fNPBins);
  gMeanPMCSAngVsPIn310deg->SetName("gMeanPMCSAngVsPIn310deg");
  gMeanPMCSAngVsPIn310deg->SetTitle("<p #times #Delta#theta_{MCS}> versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);<p #times #Delta#theta_{MCS}> (GeV)");
  fDCAList->AddAt(gMeanPMCSAngVsPIn310deg, kMeanPMCSAngVsPIn310deg);
  TGraphAsymmErrors* gSigmaPMCSAngVsPIn310deg = new TGraphAsymmErrors(fNPBins);
  gSigmaPMCSAngVsPIn310deg->SetName("gSigmaPMCSAngVsPIn310deg");
  gSigmaPMCSAngVsPIn310deg->SetTitle("#sigma_{p #times #Delta#theta_{MCS}} versus p for tracks within [3,10[ degrees at absorber end;p (GeV/c);#sigma_{p #times #Delta#theta_{MCS}} (GeV)");
  fDCAList->AddAt(gSigmaPMCSAngVsPIn310deg, kSigmaPMCSAngVsPIn310deg);
  
  // fit histo and fill graphs
  FitPDCAVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kPDCAVsPIn23deg)),
	       "p*DCA (tracks in [2,3] deg.)", gMeanPDCAVsPIn23deg, gSigmaPDCAVsPIn23deg);
  FitPDCAVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kPDCAVsPIn310deg)),
	       "p*DCA (tracks in [3,10] deg.)", gMeanPDCAVsPIn310deg, gSigmaPDCAVsPIn310deg);
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kPMCSAngVsPIn23deg)), 0., 2.e-3,
		  "p*MCSAngle (tracks in [2,3] deg.)", gMeanPMCSAngVsPIn23deg, gSigmaPMCSAngVsPIn23deg);
  FitGausResVsMom(static_cast<TH2*>(fTrackerList->UncheckedAt(kPMCSAngVsPIn310deg)), 0., 2.e-3,
		  "p*MCSAngle (tracks in [3,10] deg.)", gMeanPMCSAngVsPIn310deg, gSigmaPMCSAngVsPIn310deg);
  
  // --- compute cluster resolution ---
  fClusterList = new TObjArray(100);
  fClusterList->SetOwner();
  
  // define graphs per chamber
  TGraphErrors* gMeanResClXVsCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gMeanResClXVsCh->SetName("gMeanResClXVsCh");
  gMeanResClXVsCh->SetTitle("cluster-trackRef residual-X per Ch: mean;chamber ID;<#Delta_{X}> (cm)");
  gMeanResClXVsCh->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gMeanResClXVsCh, kMeanResClXVsCh);
  TGraphErrors* gMeanResClYVsCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gMeanResClYVsCh->SetName("gMeanResClYVsCh");
  gMeanResClYVsCh->SetTitle("cluster-trackRef residual-Y per Ch: mean;chamber ID;<#Delta_{Y}> (cm)");
  gMeanResClYVsCh->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gMeanResClYVsCh, kMeanResClYVsCh);
  TGraphErrors* gSigmaResClXVsCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gSigmaResClXVsCh->SetName("gSigmaResClXVsCh");
  gSigmaResClXVsCh->SetTitle("cluster-trackRef residual-X per Ch: sigma;chamber ID;#sigma_{X} (cm)");
  gSigmaResClXVsCh->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gSigmaResClXVsCh, kSigmaResClXVsCh);
  TGraphErrors* gSigmaResClYVsCh = new TGraphErrors(AliMUONConstants::NTrackingCh());
  gSigmaResClYVsCh->SetName("gSigmaResClYVsCh");
  gSigmaResClYVsCh->SetTitle("cluster-trackRef residual-Y per Ch: sigma;chamber ID;#sigma_{Y} (cm)");
  gSigmaResClYVsCh->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gSigmaResClYVsCh, kSigmaResClYVsCh);
  
  // define graphs per DE
  TGraphErrors* gMeanResClXVsDE = new TGraphErrors(fNDE);
  gMeanResClXVsDE->SetName("gMeanResClXVsDE");
  gMeanResClXVsDE->SetTitle("cluster-trackRef residual-X per DE: mean;DE ID;<#Delta_{X}> (cm)");
  gMeanResClXVsDE->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gMeanResClXVsDE, kMeanResClXVsDE);
  TGraphErrors* gMeanResClYVsDE = new TGraphErrors(fNDE);
  gMeanResClYVsDE->SetName("gMeanResClYVsDE");
  gMeanResClYVsDE->SetTitle("cluster-trackRef residual-Y per dE: mean;DE ID;<#Delta_{Y}> (cm)");
  gMeanResClYVsDE->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gMeanResClYVsDE, kMeanResClYVsDE);
  TGraphErrors* gSigmaResClXVsDE = new TGraphErrors(fNDE);
  gSigmaResClXVsDE->SetName("gSigmaResClXVsDE");
  gSigmaResClXVsDE->SetTitle("cluster-trackRef residual-X per DE: sigma;DE ID;#sigma_{X} (cm)");
  gSigmaResClXVsDE->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gSigmaResClXVsDE, kSigmaResClXVsDE);
  TGraphErrors* gSigmaResClYVsDE = new TGraphErrors(fNDE);
  gSigmaResClYVsDE->SetName("gSigmaResClYVsDE");
  gSigmaResClYVsDE->SetTitle("cluster-trackRef residual-Y per DE: sigma;DE ID;#sigma_{Y} (cm)");
  gSigmaResClYVsDE->SetMarkerStyle(kFullDotLarge);
  fClusterList->AddAt(gSigmaResClYVsDE, kSigmaResClYVsDE);
  
  // fit histo and fill graphs per chamber
  Double_t clusterResPerCh[10][2];
  for (Int_t i = 0; i < AliMUONConstants::NTrackingCh(); i++) {
    TH1D *tmp = static_cast<TH2*>(fTrackerList->UncheckedAt(kResClXVsCh))->ProjectionY("tmp",i+1,i+1,"e");
    FitClusterResidual(tmp, i, clusterResPerCh[i][0], gMeanResClXVsCh, gSigmaResClXVsCh);
    delete tmp;
    tmp = static_cast<TH2*>(fTrackerList->UncheckedAt(kResClYVsCh))->ProjectionY("tmp",i+1,i+1,"e");
    FitClusterResidual(tmp, i, clusterResPerCh[i][1], gMeanResClYVsCh, gSigmaResClYVsCh);
    delete tmp;
  }
  
  // fit histo and fill graphs per DE
  Double_t clusterResPerDE[200][2];
  for (Int_t i = 0; i < fNDE; i++) {
    TH1D *tmp = static_cast<TH2*>(fTrackerList->UncheckedAt(kResClXVsDE))->ProjectionY("tmp",i+1,i+1,"e");
    FitClusterResidual(tmp, i, clusterResPerDE[i][0], gMeanResClXVsDE, gSigmaResClXVsDE);
    delete tmp;
    tmp = static_cast<TH2*>(fTrackerList->UncheckedAt(kResClYVsDE))->ProjectionY("tmp",i+1,i+1,"e");
    FitClusterResidual(tmp, i, clusterResPerDE[i][1], gMeanResClYVsDE, gSigmaResClYVsDE);
    delete tmp;
  }
  
  // set DE graph labels
  TAxis* xAxis = static_cast<TH2*>(fTrackerList->UncheckedAt(kResClXVsDE))->GetXaxis();
  gMeanResClXVsDE->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
  gMeanResClYVsDE->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
  gSigmaResClXVsDE->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
  gSigmaResClYVsDE->GetXaxis()->Set(fNDE, 0.5, fNDE+0.5);
  for (Int_t i = 1; i <= fNDE; i++) {
    const char* label = xAxis->GetBinLabel(i);
    gMeanResClXVsDE->GetXaxis()->SetBinLabel(i, label);
    gMeanResClYVsDE->GetXaxis()->SetBinLabel(i, label);
    gSigmaResClXVsDE->GetXaxis()->SetBinLabel(i, label);
    gSigmaResClYVsDE->GetXaxis()->SetBinLabel(i, label);
  }
  
  // --- diplay momentum residuals at vertex ---
  TCanvas* cResPAtVtx = DrawVsAng("cResPAtVtx", "momentum residual at vertex in 3 angular regions",
				  static_cast<TH1*>(fTrackerList->UncheckedAt(kResPAtVtx)),
				  static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsAngleAtAbsEnd)));
  fPAtVtxList->AddAt(cResPAtVtx, kcResPAtVtx);
  TCanvas* cResPAtVtxMC = DrawVsAng("cResPAtVtxMC", "momentum residual at vertex in 3 MC angular regions",
				    static_cast<TH1*>(fTrackerList->UncheckedAt(kResPAtVtx)),
				    static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsMCAngle)));
  fPAtVtxList->AddAt(cResPAtVtxMC, kcResPAtVtxMC);
  TCanvas* cResPAtVtxVsPosAbsEndMC = DrawVsPos("cResPAtVtxVsPosAbsEndMC", "momentum residual at vertex versus position at absorber end in 3 MC angular regions",
					       static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsPosAbsEndIn02degMC)),
					       static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsPosAbsEndIn23degMC)),
					       static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsPosAbsEndIn310degMC)));
  fPAtVtxList->AddAt(cResPAtVtxVsPosAbsEndMC, kcResPAtVtxVsPosAbsEndMC);
  TCanvas* cResPAtVtxVsPIn23deg = DrawFitLandauGausResPVsP("cResPAtVtxVsPIn23deg", "momentum residual for tracks between 2 and 3 degrees",
						      static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsPIn23deg)),
						      10, "momentum residuals at vertex (tracks in [2,3] deg.)");
  fPAtVtxList->AddAt(cResPAtVtxVsPIn23deg, kcResPAtVtxVsPIn23deg);
  TCanvas* cResPAtVtxVsPIn310deg = DrawFitLandauGausResPVsP("cResPAtVtxVsPIn310deg", "momentum residual for tracks between 3 and 10 degrees",
						       static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsPIn310deg)),
						       10, "momentum residuals at vertex (tracks in [3,10] deg.)");
  fPAtVtxList->AddAt(cResPAtVtxVsPIn310deg, kcResPAtVtxVsPIn310deg);
  TCanvas* cResPAtVtxVsPIn02degMC = DrawResPVsP("cResPAtVtxVsPIn02degMC", "momentum residuals for tracks with MC angle < 2 degrees",
					   static_cast<TH2*>(fTrackerList->UncheckedAt(kResPAtVtxVsPIn02degMC)), 5);
  fPAtVtxList->AddAt(cResPAtVtxVsPIn02degMC, kcResPAtVtxVsPIn02degMC);
  
  // --- diplay slopeX residuals at vertex ---
  TCanvas* cResSlopeXAtVtx = DrawVsAng("cResSlopeXAtVtx", "slope_{X} residual at vertex in 3 angular regions",
				       static_cast<TH1*>(fTrackerList->UncheckedAt(kResSlopeXAtVtx)),
				       static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeXAtVtxVsAngleAtAbsEnd)));
  fSlopeAtVtxList->AddAt(cResSlopeXAtVtx, kcResSlopeXAtVtx);
  TCanvas* cResSlopeYAtVtx = DrawVsAng("cResSlopeYAtVtx", "slope_{Y} residual at vertex in 3 angular regions",
				       static_cast<TH1*>(fTrackerList->UncheckedAt(kResSlopeYAtVtx)),
				       static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeYAtVtxVsAngleAtAbsEnd)));
  fSlopeAtVtxList->AddAt(cResSlopeYAtVtx, kcResSlopeYAtVtx);
  TCanvas* cResSlopeXAtVtxMC = DrawVsAng("cResSlopeXAtVtxMC", "slope_{X} residual at vertex in 3 MC angular regions",
					 static_cast<TH1*>(fTrackerList->UncheckedAt(kResSlopeXAtVtx)),
					 static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeXAtVtxVsMCAngle)));
  fSlopeAtVtxList->AddAt(cResSlopeXAtVtxMC, kcResSlopeXAtVtxMC);
  TCanvas* cResSlopeYAtVtxMC = DrawVsAng("cResSlopeYAtVtxMC", "slope_{Y} residual at vertex in 3 MC angular regions",
					 static_cast<TH1*>(fTrackerList->UncheckedAt(kResSlopeYAtVtx)),
					 static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeYAtVtxVsMCAngle)));
  fSlopeAtVtxList->AddAt(cResSlopeYAtVtxMC, kcResSlopeYAtVtxMC);
  TCanvas* cResSlopeXAtVtxVsPosAbsEndMC = DrawVsPos("cResSlopeXAtVtxVsPosAbsEndMC", "slope_{X} residual at vertex versus position at absorber end in 3 MC angular regions",
						    static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeXAtVtxVsPosAbsEndIn02degMC)),
						    static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeXAtVtxVsPosAbsEndIn23degMC)),
						    static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeXAtVtxVsPosAbsEndIn310degMC)));
  fSlopeAtVtxList->AddAt(cResSlopeXAtVtxVsPosAbsEndMC, kcResSlopeXAtVtxVsPosAbsEndMC);
  TCanvas* cResSlopeYAtVtxVsPosAbsEndMC = DrawVsPos("cResSlopeYAtVtxVsPosAbsEndMC", "slope_{Y} residual at vertex versus position at absorber end in 3 MC angular regions",
						    static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeYAtVtxVsPosAbsEndIn02degMC)),
						    static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeYAtVtxVsPosAbsEndIn23degMC)),
						    static_cast<TH2*>(fTrackerList->UncheckedAt(kResSlopeYAtVtxVsPosAbsEndIn310degMC)));
  fSlopeAtVtxList->AddAt(cResSlopeYAtVtxVsPosAbsEndMC, kcResSlopeYAtVtxVsPosAbsEndMC);
  
  // --- diplay eta residuals at vertex ---
  TCanvas* cResEtaAtVtx = DrawVsAng("cResEtaAtVtx", "eta residual at vertex in 3 angular regions",
				    static_cast<TH1*>(fTrackerList->UncheckedAt(kResEtaAtVtx)),
				    static_cast<TH2*>(fTrackerList->UncheckedAt(kResEtaAtVtxVsAngleAtAbsEnd)));
  fEtaAtVtxList->AddAt(cResEtaAtVtx, kcResEtaAtVtx);
  TCanvas* cResEtaAtVtxMC = DrawVsAng("cResEtaAtVtxMC", "eta residual at vertex in 3 MC angular regions",
				      static_cast<TH1*>(fTrackerList->UncheckedAt(kResEtaAtVtx)),
				      static_cast<TH2*>(fTrackerList->UncheckedAt(kResEtaAtVtxVsMCAngle)));
  fEtaAtVtxList->AddAt(cResEtaAtVtxMC, kcResEtaAtVtxMC);
  TCanvas* cResEtaAtVtxVsPosAbsEndMC = DrawVsPos("cResEtaAtVtxVsPosAbsEndMC", "eta residual at vertex versus position at absorber end in 3 MC angular regions",
						 static_cast<TH2*>(fTrackerList->UncheckedAt(kResEtaAtVtxVsPosAbsEndIn02degMC)),
						 static_cast<TH2*>(fTrackerList->UncheckedAt(kResEtaAtVtxVsPosAbsEndIn23degMC)),
						 static_cast<TH2*>(fTrackerList->UncheckedAt(kResEtaAtVtxVsPosAbsEndIn310degMC)));
  fEtaAtVtxList->AddAt(cResEtaAtVtxVsPosAbsEndMC, kcResEtaAtVtxVsPosAbsEndMC);
  
  // --- diplay phi residuals at vertex ---
  TCanvas* cResPhiAtVtx = DrawVsAng("cResPhiAtVtx", "phi residual at vertex in 3 angular regions",
				    static_cast<TH1*>(fTrackerList->UncheckedAt(kResPhiAtVtx)),
				    static_cast<TH2*>(fTrackerList->UncheckedAt(kResPhiAtVtxVsAngleAtAbsEnd)));
  fPhiAtVtxList->AddAt(cResPhiAtVtx, kcResPhiAtVtx);
  TCanvas* cResPhiAtVtxMC = DrawVsAng("cResPhiAtVtxMC", "phi residual at vertex in 3 MC angular regions",
				      static_cast<TH1*>(fTrackerList->UncheckedAt(kResPhiAtVtx)),
				      static_cast<TH2*>(fTrackerList->UncheckedAt(kResPhiAtVtxVsMCAngle)));
  fPhiAtVtxList->AddAt(cResPhiAtVtxMC, kcResPhiAtVtxMC);
  TCanvas* cResPhiAtVtxVsPosAbsEndMC = DrawVsPos("cResPhiAtVtxVsPosAbsEndMC", "phi residual at vertex versus position at absorber end in 3 MC angular regions",
						 static_cast<TH2*>(fTrackerList->UncheckedAt(kResPhiAtVtxVsPosAbsEndIn02degMC)),
						 static_cast<TH2*>(fTrackerList->UncheckedAt(kResPhiAtVtxVsPosAbsEndIn23degMC)),
						 static_cast<TH2*>(fTrackerList->UncheckedAt(kResPhiAtVtxVsPosAbsEndIn310degMC)));
  fPhiAtVtxList->AddAt(cResPhiAtVtxVsPosAbsEndMC, kcResPhiAtVtxVsPosAbsEndMC);
  
  // --- diplay P*DCA ---
  TCanvas* cPDCA = DrawVsAng("cPDCA", "p #times DCA in 3 angular regions",
			     static_cast<TH1*>(fTrackerList->UncheckedAt(kPDCA)),
			     static_cast<TH2*>(fTrackerList->UncheckedAt(kPDCAVsAngleAtAbsEnd)));
  fDCAList->AddAt(cPDCA, kcPDCA);
  TCanvas* cPDCAMC = DrawVsAng("cPDCAMC", "p #times DCA in 3 MC angular regions",
			       static_cast<TH1*>(fTrackerList->UncheckedAt(kPDCA)),
			       static_cast<TH2*>(fTrackerList->UncheckedAt(kPDCAVsMCAngle)));
  fDCAList->AddAt(cPDCAMC, kcPDCAMC);
  TCanvas* cPDCAVsPosAbsEndMC = DrawVsPos("cPDCAVsPosAbsEndMC", "p #times DCA versus position at absorber end in 3 MC angular regions",
				  static_cast<TH2*>(fTrackerList->UncheckedAt(kPDCAVsPosAbsEndIn02degMC)),
				  static_cast<TH2*>(fTrackerList->UncheckedAt(kPDCAVsPosAbsEndIn23degMC)),
				  static_cast<TH2*>(fTrackerList->UncheckedAt(kPDCAVsPosAbsEndIn310degMC)));
  fDCAList->AddAt(cPDCAVsPosAbsEndMC, kcPDCAVsPosAbsEndMC);
  
  // post param containers
  PostData(2, fEfficiencyList);
  PostData(5, fPAtVtxList);
  PostData(6, fSlopeAtVtxList);
  PostData(7, fEtaAtVtxList);
  PostData(8, fPhiAtVtxList);
  PostData(9, fPAt1stClList);
  PostData(10, fSlopeAt1stClList);
  PostData(11, fDCAList);
  PostData(12, fClusterList);
}


//________________________________________________________________________
Bool_t AliAnalysisTaskMuonPerformance::GetEfficiency(AliCFEffGrid* efficiency, Double_t& calcEff, Double_t& calcEffErr)
{
  //
  /// Calculate the efficiency when cuts on the THnSparse are applied
  //

  Bool_t isGood = kTRUE;

  TH1* histo = 0x0;
  Double_t sum[2] = {0., 0.};
  for ( Int_t ihisto=0; ihisto<2; ihisto++ ) {
    histo = ( ihisto==0 ) ? efficiency->GetNum()->Project(kVarCharge) : efficiency->GetDen()->Project(kVarCharge);
    sum[ihisto] = histo->Integral();
    delete histo;
  }

  if ( sum[1] == 0. ) isGood = kFALSE;
  
  calcEff = ( isGood ) ? sum[0]/sum[1] : 0.;
  if ( calcEff > 1. ) isGood = kFALSE;

  calcEffErr = ( isGood ) ? TMath::Sqrt(calcEff*(1-calcEff)/sum[1]) : 0.;

  return isGood;
}

//________________________________________________________________________
Int_t AliAnalysisTaskMuonPerformance::RecoTrackMother(AliMCParticle* mcParticle)
{
  //
  /// Find track mother from kinematics
  //

  if ( ! mcParticle )
    return kUnknownPart;

  Int_t recoPdg = mcParticle->PdgCode();

  // Track is not a muon
  if ( TMath::Abs(recoPdg) != 13 ) return kRecoHadron;

  Int_t imother = mcParticle->GetMother();

  Bool_t isFirstMotherHF = kFALSE;
  Int_t step = 0;

  //printf("\n"); // REMEMBER TO CUT
  while ( imother >= 0 ) {
    TParticle* part = static_cast<AliMCParticle*>(fMCEvent->GetTrack(imother))->Particle();
    //printf("Particle %s  -> %s\n", part->GetName(), TMCProcessName[part->GetUniqueID()]); // REMEMBER TO CUT

    if ( part->GetUniqueID() == kPHadronic ) 
      return kSecondaryMu;

    Int_t absPdg = TMath::Abs(part->GetPdgCode());

    step++;
    if ( step == 1 ) // only for first mother
      isFirstMotherHF = ( ( absPdg >= 400  && absPdg < 600 ) || 
			  ( absPdg >= 4000 && absPdg < 6000 ) );
      
    if ( absPdg < 4 )
      return kPrimaryMu;

    if ( isFirstMotherHF) {
      if ( absPdg == 4 )
	return kCharmMu;
      else if ( absPdg == 5 )
	return kBeautyMu;
    }

    imother = part->GetFirstMother();
  }

  return kPrimaryMu;

}

//________________________________________________________________________
Float_t AliAnalysisTaskMuonPerformance::GetBinThetaAbsEnd(Float_t RAtAbsEnd, Bool_t isTheta)
{
  //
  /// Get bin of theta at absorber end region
  //
  Float_t thetaDeg = ( isTheta ) ? RAtAbsEnd : TMath::ATan( RAtAbsEnd / 505. );
  thetaDeg *= TMath::RadToDeg();
  if ( thetaDeg < 2. )
    return 0.;
  else if ( thetaDeg < 3. )
    return 1.;
  else if ( thetaDeg < 9. )
    return 2.;

  return 3.;
}

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::FillContainerInfo(Double_t* containerInput, AliESDMuonTrack* esdTrack, Int_t mcID)
{
  //
  /// Fill container info (except kVarMatchMC)
  //

  AliMCParticle* mcPart = 0x0;
  if (mcID >= 0) mcPart = static_cast<AliMCParticle*>(fMCEvent->GetTrack(mcID));

  AliVParticle* track = esdTrack;
  if (!track) track = mcPart;

  containerInput[kVarPt] = track->Pt();
  containerInput[kVarEta] = track->Eta();
  containerInput[kVarPhi] = track->Phi();
  containerInput[kVarThetaZones] = (esdTrack) ? GetBinThetaAbsEnd(esdTrack->GetRAtAbsorberEnd()) : GetBinThetaAbsEnd(TMath::Pi()-track->Theta(),kTRUE);
  containerInput[kVarCharge] = (esdTrack) ? static_cast<Double_t>(track->Charge()) : static_cast<Double_t>(track->Charge())/3.;
  containerInput[kVarHasTracker] = (esdTrack) ? static_cast<Double_t>(esdTrack->ContainTrackerData()) : 0.;
  containerInput[kVarTrigger] = (esdTrack) ? static_cast<Double_t>(esdTrack->GetMatchTrigger()) : 0.;
  containerInput[kVarMotherType] = static_cast<Double_t>(RecoTrackMother(mcPart));

  if (esdTrack) esdTrack->SetLabel(mcID);
}

//________________________________________________________________________
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

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::FitLandauGausResVsP(TH2* h, const char* fitting, TGraphAsymmErrors* gMean,
							 TGraphAsymmErrors* gMostProb, TGraphAsymmErrors* gSigma)
{
  /// generic function to fit residuals versus momentum with a landau convoluted with a gaussian
  
  static TF1* fLandauGaus = 0x0;
  if (!fLandauGaus) fLandauGaus = new TF1("fLandauGaus",langaufun,h->GetYaxis()->GetBinLowEdge(1),h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1),4);
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/fNPBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,fNPBins)<<flush;
    
    TH1D *tmp = h->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    
    // first fit
    fLandauGaus->SetParameters(0.2,0.,(Double_t)tmp->GetEntries(),1.);
    tmp->Fit("fLandauGaus","WWNQ");
    
    // rebin histo
    Double_t fwhm = fLandauGaus->GetParameter(0);
    Double_t sigma = fLandauGaus->GetParameter(3);
    Double_t sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Int_t rebin = TMath::Max(static_cast<Int_t>(0.5*sigmaP/tmp->GetBinWidth(1)),1);
    while (tmp->GetNbinsX()%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    
    // second fit
    tmp->Fit("fLandauGaus","NQ");
    
    // get fit results and fill histograms
    fwhm = fLandauGaus->GetParameter(0);
    sigma = fLandauGaus->GetParameter(3);
    sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Double_t fwhmErr = fLandauGaus->GetParError(0);
    Double_t sigmaErr = fLandauGaus->GetParError(3);
    Double_t sigmaPErr = TMath::Sqrt(sigma*sigma*sigmaErr*sigmaErr + fwhm*fwhm*fwhmErr*fwhmErr/(64.*log(2.)*log(2.))) / sigmaP;
    h->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = (tmp->GetEntries() > 0) ? h->GetMean() : 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetBinLowEdge(i-rebinFactorX+1), h->GetBinLowEdge(i+1)-p};
    gMean->SetPoint(i/rebinFactorX-1, p, tmp->GetMean());
    gMean->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], tmp->GetMeanError(), tmp->GetMeanError());
    gMostProb->SetPoint(i/rebinFactorX-1, p, -fLandauGaus->GetParameter(1));
    gMostProb->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fLandauGaus->GetParError(1), fLandauGaus->GetParError(1));
    gSigma->SetPoint(i/rebinFactorX-1, p, sigmaP);
    gSigma->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], sigmaPErr, sigmaPErr);
    
    // clean memory
    delete tmp;
  }
  
  cout<<Form("\rFitting %s... %d/%d",fitting,fNPBins,fNPBins)<<endl;
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::FitGausResVsMom(TH2* h, const Double_t mean0,
						     const Double_t sigma0, const char* fitting,
						     TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma)
{
  /// generic function to fit residuals versus momentum with a gaussian
  
  static TF1* fGaus = 0x0;
  if (!fGaus) fGaus = new TF1("fGaus","gaus");
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/fNPBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,fNPBins)<<flush;
    
    TH1D *tmp = h->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    
    // first fit
    fGaus->SetParameters(tmp->GetEntries(), mean0, sigma0);
    tmp->Fit("fGaus","WWNQ");
    
    // rebin histo
    Int_t rebin = TMath::Max(static_cast<Int_t>(0.5*fGaus->GetParameter(2)/tmp->GetBinWidth(1)),1);
    while (tmp->GetNbinsX()%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    
    // second fit
    tmp->Fit("fGaus","NQ");
    
    // get fit results and fill histograms
    h->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = (tmp->GetEntries() > 0) ? h->GetMean() : 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetBinLowEdge(i-rebinFactorX+1), h->GetBinLowEdge(i+1)-p};
    gMean->SetPoint(i/rebinFactorX-1, p, fGaus->GetParameter(1));
    gMean->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fGaus->GetParError(1), fGaus->GetParError(1));
    gSigma->SetPoint(i/rebinFactorX-1, p, fGaus->GetParameter(2));
    gSigma->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fGaus->GetParError(2), fGaus->GetParError(2));
    
    // clean memory
    delete tmp;
  }
  
  cout<<Form("\rFitting %s... %d/%d",fitting,fNPBins,fNPBins)<<endl;
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::FitPDCAVsMom(TH2* h, const char* fitting,
						  TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma)
{
  /// generic function to fit p*DCA distributions
  
  static TF1* fPGaus = 0x0;
  if (!fPGaus) fPGaus = new TF1("fPGaus","x*gaus");
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/fNPBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,fNPBins)<<flush;
    
    TH1D *tmp = h->ProjectionY("tmp",i-rebinFactorX+1,i,"e");
    
    // rebin histo
    Int_t rebin = static_cast<Int_t>(25*(tmp->GetNbinsX()/(tmp->GetBinLowEdge(tmp->GetNbinsX()+1)-tmp->GetBinLowEdge(1))));
    while (tmp->GetNbinsX()%rebin!=0) rebin--;
    tmp->Rebin(rebin);
    
    // fit
    fPGaus->SetParameters(1.,0.,80.);
    tmp->Fit("fPGaus","NQ");
    
    // get fit results and fill histograms
    h->GetXaxis()->SetRange(i-rebinFactorX+1,i);
    Double_t p = (tmp->GetEntries() > 0) ? h->GetMean() : 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    h->GetXaxis()->SetRange();
    Double_t pErr[2] = {p-h->GetBinLowEdge(i-rebinFactorX+1), h->GetBinLowEdge(i+1)-p};
    gMean->SetPoint(i/rebinFactorX-1, p, fPGaus->GetParameter(1));
    gMean->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fPGaus->GetParError(1), fPGaus->GetParError(1));
    gSigma->SetPoint(i/rebinFactorX-1, p, fPGaus->GetParameter(2));
    gSigma->SetPointError(i/rebinFactorX-1, pErr[0], pErr[1], fPGaus->GetParError(2), fPGaus->GetParError(2));
    
    // clean memory
    delete tmp;
  }
  
  cout<<Form("\rFitting %s... %d/%d",fitting,fNPBins,fNPBins)<<endl;
  
}

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::FitClusterResidual(TH1* h, Int_t i, Double_t& sigma,
							TGraphErrors* gMean, TGraphErrors* gSigma)
{
  /// fill graphs with residual mean and sigma
  
  static TF1* fRGaus = 0x0;
  Double_t mean, meanErr, sigmaErr;
  
  if (fFitResiduals) {
    
    if (!fRGaus) fRGaus = new TF1("fRGaus","gaus");
    
    // first fit
    fRGaus->SetParameters(h->GetEntries(), 0., 0.1);
    h->Fit("fRGaus", "WWNQ");
    
    // rebin histo
    Int_t rebin = TMath::Max(static_cast<Int_t>(0.2*fRGaus->GetParameter(2)/h->GetBinWidth(1)),1);
    while (h->GetNbinsX()%rebin!=0) rebin--;
    h->Rebin(rebin);
    
    // second fit
    h->Fit("fRGaus","NQ");
    
    mean = fRGaus->GetParameter(1);
    meanErr = fRGaus->GetParError(1);
    sigma = fRGaus->GetParameter(2);
    sigmaErr = fRGaus->GetParError(2);
    
  } else {
    
    Zoom(h);
    mean = h->GetMean();
    meanErr = h->GetMeanError();
    sigma = h->GetRMS();
    sigmaErr = h->GetRMSError();
    h->GetXaxis()->SetRange(0,0);
    
  }
  
  gMean->SetPoint(i, i+1, mean);
  gMean->SetPointError(i, 0., meanErr);
  
  if (fCorrectForSystematics) {
    Double_t s = TMath::Sqrt(mean*mean + sigma*sigma);
    sigmaErr = (s>0.) ? TMath::Sqrt(sigma*sigma*sigmaErr*sigmaErr + mean*mean*meanErr*meanErr) / s : 0.;
    sigma = s;
  }
  
  gSigma->SetPoint(i, i+1, sigma);
  gSigma->SetPointError(i, 0., sigmaErr);
  
}

//________________________________________________________________________
TCanvas* AliAnalysisTaskMuonPerformance::DrawVsAng(const char* name, const char* title, TH1* h1, TH2* h2)
{
  /// generic function to draw histograms versus absorber angular region
  TCanvas* c = new TCanvas(name, title);
  c->cd();
  h1->DrawCopy();
  TH1D *proj1 = h2->ProjectionY(Form("%s_proj_0_2",h2->GetName()),1,2);
  proj1->SetLineColor(2);
  proj1->Draw("sames");
  TH1D *proj2 = h2->ProjectionY(Form("%s_proj_2_3",h2->GetName()),3,3);
  proj2->SetLineColor(4);
  proj2->Draw("sames");
  TH1D *proj3 = h2->ProjectionY(Form("%s__proj_3_10",h2->GetName()),4,10);
  proj3->SetLineColor(3);
  proj3->Draw("sames");
  return c;
}

//________________________________________________________________________
TCanvas* AliAnalysisTaskMuonPerformance::DrawVsPos(const char* name, const char* title, TH2* h1, TH2* h2, TH2* h3)
{
  /// generic function to draw histograms versus position at absorber end
  TCanvas* c = new TCanvas(name, title);
  c->cd();
  h1->SetMarkerColor(2);
  h1->DrawCopy();
  h2->SetMarkerColor(4);
  h2->DrawCopy("sames");
  h3->SetMarkerColor(3);
  h3->DrawCopy("sames");
  return c;
}

//________________________________________________________________________
TCanvas* AliAnalysisTaskMuonPerformance::DrawFitLandauGausResPVsP(const char* name, const char* title,
								  TH2* h, const Int_t nBins, const char* fitting)
{
  /// generic function to draw and fit momentum residuals versus momentum
  
  static TF1* fLandauGaus = 0x0;
  if (!fLandauGaus) fLandauGaus = new TF1("fLandauGaus",langaufun,h->GetYaxis()->GetBinLowEdge(1),h->GetYaxis()->GetBinLowEdge(h->GetNbinsY()+1),4);
  
  TLegend* l = new TLegend(0.15,0.25,0.3,0.85);
  TCanvas* c = new TCanvas(name, title);
  c->cd();
  
  h->Sumw2();
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/nBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    
    cout<<Form("\rFitting %s... %d/%d",fitting,i/rebinFactorX,nBins)<<flush;
    
    // draw projection
    TH1D* proj = h->ProjectionY(Form("%s_%d",h->GetName(),i/rebinFactorX),i-rebinFactorX+1,i);
    if (proj->GetEntries() > 0) proj->Scale(1./proj->GetEntries());
    proj->Draw((i==rebinFactorX)?"hist":"histsames");
    proj->SetLineColor(i/rebinFactorX);
    
    // first fit
    fLandauGaus->SetParameters(0.2,0.,1.,1.);
    fLandauGaus->SetLineColor(i/rebinFactorX);
    proj->Fit("fLandauGaus","WWNQ","sames");
    
    // rebin histo
    Double_t fwhm = fLandauGaus->GetParameter(0);
    Double_t sigma = fLandauGaus->GetParameter(3);
    Double_t sigmaP = TMath::Sqrt(sigma*sigma + fwhm*fwhm/(8.*log(2.)));
    Int_t rebin = TMath::Max(static_cast<Int_t>(0.5*sigmaP/proj->GetBinWidth(1)),1);
    while (proj->GetNbinsX()%rebin!=0) rebin--;
    proj->Rebin(rebin);
    proj->Scale(1./rebin);
    
    // second fit
    proj->Fit("fLandauGaus","Q","sames");
    
    // set label
    Double_t p = 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    l->AddEntry(proj,Form("%5.1f GeV",p));
    
  }
  
  cout<<Form("\rFitting %s... %d/%d",fitting,nBins,nBins)<<endl;
  
  l->Draw("same");
  
  return c;
}

//________________________________________________________________________
TCanvas* AliAnalysisTaskMuonPerformance::DrawResPVsP(const char* name, const char* title, TH2* h, const Int_t nBins)
{
  /// generic function to draw momentum residuals versus momentum
  
  TLegend* l = new TLegend(0.15,0.25,0.3,0.85);
  TCanvas* c = new TCanvas(name, title);
  c->cd();
  
  h->Sumw2();
  
  Int_t rebinFactorX = TMath::Max(h->GetNbinsX()/nBins, 1);
  for (Int_t i = rebinFactorX; i <= h->GetNbinsX(); i+=rebinFactorX) {
    
    // draw projection
    TH1D* proj = h->ProjectionY(Form("%s_%d",h->GetName(),i/rebinFactorX),i-rebinFactorX+1,i);
    if (proj->GetEntries() > 0) proj->Scale(1./proj->GetEntries());
    proj->Draw((i==rebinFactorX)?"hist":"histsames");
    proj->SetLineColor(i/rebinFactorX);
    proj->SetLineWidth(2);
    
    // set label
    Double_t p = 0.5 * (h->GetBinLowEdge(i-rebinFactorX+1) + h->GetBinLowEdge(i+1));
    l->AddEntry(proj,Form("%5.1f GeV",p));
    
  }
  
  l->Draw("same");
  
  return c;
}

//________________________________________________________________________
void AliAnalysisTaskMuonPerformance::Zoom(TH1* h, Double_t fractionCut)
{
  /// Reduce the range of the histogram by removing a given fration of the statistic at each edge
  
  Double_t maxEventsCut = fractionCut * h->GetEntries();
  Int_t nBins = h->GetNbinsX();
  
  // set low edge  
  Int_t minBin;
  Double_t eventsCut = 0.;
  for (minBin = 1; minBin <= nBins; minBin++) {
    eventsCut += h->GetBinContent(minBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set high edge
  Int_t maxBin;
  eventsCut = 0.;
  for (maxBin = nBins; maxBin >= 1; maxBin--) {
    eventsCut += h->GetBinContent(maxBin);
    if (eventsCut > maxEventsCut) break;
  }
  
  // set new axis range
  h->GetXaxis()->SetRange(--minBin, ++maxBin);
}
