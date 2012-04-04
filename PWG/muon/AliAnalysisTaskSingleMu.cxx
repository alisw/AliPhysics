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

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliAnalysisTaskSingleMu
/// Analysis task for single muons in the spectrometer.
/// The output is a list of histograms and CF containers.
/// The macro class can run on AODs or ESDs.
/// In the latter case a flag can be activated to produce a tree as output.
/// If Monte Carlo information is present, some basics checks are performed.
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

#define AliAnalysisTaskSingleMu_cxx

#include "AliAnalysisTaskSingleMu.h"

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TStyle.h"
//#include "TMCProcess.h"
#include "TArrayI.h"
#include "TPaveStats.h"
#include "TFitResultPtr.h"

// STEER includes
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliVHeader.h"
#include "AliAODMCHeader.h"
#include "AliStack.h"

// ANALYSIS includes
#include "AliAnalysisManager.h"

// CORRFW includes
#include "AliCFContainer.h"
#include "AliCFGridSparse.h"
#include "AliCFEffGrid.h"

// PWG includes
#include "AliVAnalysisMuon.h"
#include "AliMergeableCollection.h"
#include "AliCounterCollection.h"
#include "AliMuonTrackCuts.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSingleMu) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu() :
  AliVAnalysisMuon(),
  fThetaAbsKeys(0x0)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu(const char *name, const AliMuonTrackCuts& cuts) :
  AliVAnalysisMuon(name, cuts),
  fThetaAbsKeys(0x0)
{
  //
  /// Constructor.
  //
  TString thetaAbsKeys = "ThetaAbs23 ThetaAbs310";
  fThetaAbsKeys = thetaAbsKeys.Tokenize(" ");
}


//________________________________________________________________________
AliAnalysisTaskSingleMu::~AliAnalysisTaskSingleMu()
{
  //
  /// Destructor
  //

  delete fThetaAbsKeys;
}

//___________________________________________________________________________
void AliAnalysisTaskSingleMu::MyUserCreateOutputObjects()
{

  TH1* histo = 0x0;
  TString histoName = "", histoTitle = "";
  
  Int_t nVzBins = 40;
  Double_t vzMin = -20., vzMax = 20.;
  TString vzName("Vz"), vzTitle("Vz"), vzUnits("cm");  
  
  histoName = "hIpVtx";
  histo = new TH1F(histoName.Data(), histoName.Data(), nVzBins, vzMin, vzMax);
  histo->SetXTitle("v_{z} (cm)");
  AddObjectToCollection(histo, kIPVz);

  
  Int_t nPtBins = 80;
  Double_t ptMin = 0., ptMax = 80.;
  TString ptName("Pt"), ptTitle("p_{t}"), ptUnits("GeV/c");
  
  Int_t nEtaBins = 25;
  Double_t etaMin = -4.5, etaMax = -2.;
  TString etaName("Eta"), etaTitle("#eta"), etaUnits("");
  
  Int_t nPhiBins = 36;
  Double_t phiMin = 0.; Double_t phiMax = 2.*TMath::Pi();
  TString phiName("Phi"), phiTitle("#phi"), phiUnits("rad");
    
  Int_t nChargeBins = 2;
  Double_t chargeMin = -2, chargeMax = 2.;
  TString chargeName("Charge"), chargeTitle("charge"), chargeUnits("e");
  
  Int_t nThetaAbsEndBins = 2;
  Double_t thetaAbsEndMin = -0.5, thetaAbsEndMax = 1.5;
  TString thetaAbsEndName("ThetaAbsEnd"), thetaAbsEndTitle("#theta_{abs}"), thetaAbsEndUnits("a.u.");    
  
  Int_t nMotherTypeBins = kNtrackSources;
  Double_t motherTypeMin = -0.5, motherTypeMax = (Double_t)kNtrackSources - 0.5;
  TString motherType("MotherType"), motherTypeTitle("motherType"), motherTypeUnits("");
    
  Int_t nbins[kNvars] = {nPtBins, nEtaBins, nPhiBins, nVzBins, nChargeBins, nThetaAbsEndBins, nMotherTypeBins};
  Double_t xmin[kNvars] = {ptMin, etaMin, phiMin, vzMin, chargeMin, thetaAbsEndMin, motherTypeMin};
  Double_t xmax[kNvars] = {ptMax, etaMax, phiMax, vzMax, chargeMax, thetaAbsEndMax, motherTypeMax};
  TString axisTitle[kNvars] = {ptTitle, etaTitle, phiTitle, vzTitle, chargeTitle, thetaAbsEndTitle, motherTypeTitle};
  TString axisUnits[kNvars] = {ptUnits, etaUnits, phiUnits, vzUnits, chargeUnits, thetaAbsEndUnits, motherTypeUnits};

  AliCFContainer* cfContainer = new AliCFContainer("SingleMuContainer","Container for tracks",kNsteps,kNvars,nbins);
  
  for ( Int_t idim = 0; idim<kNvars; idim++){
    histoTitle = Form("%s (%s)", axisTitle[idim].Data(), axisUnits[idim].Data());
    histoTitle.ReplaceAll("()","");
    
    cfContainer->SetVarTitle(idim, histoTitle.Data());
    cfContainer->SetBinLimits(idim, xmin[idim], xmax[idim]);
  }
  
  TString stepTitle[kNsteps] = {"reconstructed", "generated"};

  TAxis* currAxis = 0x0;
  for (Int_t istep=0; istep<kNsteps; istep++){
    cfContainer->SetStepTitle(istep, stepTitle[istep].Data());
    AliCFGridSparse* gridSparse = cfContainer->GetGrid(istep);
        
    currAxis = gridSparse->GetAxis(kHvarMotherType);
    for ( Int_t ibin=0; ibin<fSrcKeys->GetEntries(); ibin++ ) {
      currAxis->SetBinLabel(ibin+1, fSrcKeys->At(ibin)->GetName());
    }
  }
  
  AddObjectToCollection(cfContainer, kTrackContainer);
  
  fMuonTrackCuts->Print("mask");
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality)
{
  //
  /// Fill output objects
  //

  if ( GetVertexSPD()->GetNContributors() < fMinNvtxContirbutors ) return;

  Double_t ipVz = GetVertexSPD()->GetZ();
  Double_t ipVzMC = 0;
  if ( IsMC() ) {
    if ( fMCEvent ) ipVzMC = fMCEvent->GetPrimaryVertex()->GetZ();
    else if ( fAODEvent ) {
      AliAODMCHeader* aodMCHeader = (AliAODMCHeader *)fAODEvent->FindListObject(AliAODMCHeader::StdBranchName());
      if ( aodMCHeader ) ipVzMC = aodMCHeader->GetVtxZ();
    }
  }
  
  for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
    TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
    ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, "hIpVtx"))->Fill(ipVz);
  }

  Bool_t isPileupFromSPD = ( fAODEvent && ! fAODEvent->GetTracklets() ) ? InputEvent()->IsPileupFromSPD(3, 0.8, 3., 2., 5.) : InputEvent()->IsPileupFromSPDInMultBins(); // Avoid break when reading Muon AODs (tracklet info is not present and IsPileupFromSPDInMultBins crashes
  if ( isPileupFromSPD ) return;
  
  Double_t containerInput[kNvars];
  AliVParticle* track = 0x0;

  for ( Int_t istep = 0; istep<2; ++istep ) {
    Int_t nTracks = ( istep == kStepReconstructed ) ? GetNTracks() : GetNMCTracks();
    for (Int_t itrack = 0; itrack < nTracks; itrack++) {
      track = ( istep == kStepReconstructed ) ? GetTrack(itrack) : GetMCTrack(itrack);
      
      Bool_t isSelected = ( istep == kStepReconstructed ) ? fMuonTrackCuts->IsSelected(track) : ( TMath::Abs(track->PdgCode()) == 13 );
      if ( ! isSelected ) continue;
      
      // In W simulations with Pythia, sometimes muon is stored twice.
      // Remove muon in case it has another muon as daugther
      if ( istep == kStepGeneratedMC ) {
        Int_t firstDaughter = GetDaughterIndex(track, 0);
        if ( firstDaughter >= 0 ) {
          Bool_t hasMuonDaughter = kFALSE;
          Int_t lastDaughter = GetDaughterIndex(track, 1);
          for ( Int_t idaugh=firstDaughter; idaugh<=lastDaughter; idaugh++ ) {
            AliVParticle* currTrack = GetMCTrack(idaugh);
            if ( currTrack->PdgCode() == track->PdgCode() ) {
              hasMuonDaughter = kTRUE;
              break;
            }
          }
          if ( hasMuonDaughter ) {
            AliDebug(1, Form("Current muon (%i) has muon daughter: rejecting it", itrack));
            continue;
          }
        }
      }      
      
      Int_t trackSrc = ( istep == kStepReconstructed ) ? GetParticleType(track) : RecoTrackMother(track);
      
      Double_t thetaAbsEndDeg = 0;
      if ( istep == kStepReconstructed ) {
        Double_t rAbsEnd =  ( fAODEvent ) ? ((AliAODTrack*)track)->GetRAtAbsorberEnd(): ((AliESDMuonTrack*)track)->GetRAtAbsorberEnd();
        thetaAbsEndDeg = TMath::ATan( rAbsEnd / 505. ) * TMath::RadToDeg();
      }
      else {
        thetaAbsEndDeg = ( TMath::Pi()-track->Theta() ) * TMath::RadToDeg();
      }
      Int_t thetaAbsBin = ( thetaAbsEndDeg < 3. ) ? kThetaAbs23 : kThetaAbs310;

      containerInput[kHvarPt]         = track->Pt();
      containerInput[kHvarEta]        = track->Eta();
      containerInput[kHvarPhi]        = track->Phi();
      containerInput[kHvarVz]         = ( istep == kStepReconstructed ) ? ipVz : ipVzMC;
      containerInput[kHvarCharge]     = track->Charge()/3.;
      containerInput[kHvarThetaAbs]   = (Double_t)thetaAbsBin;
      containerInput[kHvarMotherType] = (Double_t)trackSrc;
      
      for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
        TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
        if ( istep == kStepReconstructed && ! TrackPtCutMatchTrigClass(track, trigClassName) ) continue;
        ((AliCFContainer*)GetMergeableObject(physSel, trigClassName, centrality, "SingleMuContainer"))->Fill(containerInput,istep);
      } // loop on selected trigger classes
    } // loop on tracks
  } // loop on container steps
}


//________________________________________________________________________
void AliAnalysisTaskSingleMu::Terminate(Option_t *) {
  //
  /// Draw some histograms at the end.
  //

  AliVAnalysisMuon::Terminate("");

  if ( ! fMergeableCollection ) return;
  
  TString physSel = fTerminateOptions->At(0)->GetName();
  TString trigClassName = fTerminateOptions->At(1)->GetName();
  TString centralityRange = fTerminateOptions->At(2)->GetName();
  TString furtherOpt = fTerminateOptions->At(3)->GetName();
  
  TString minBiasTrig = "";
  TObjArray* optArr = furtherOpt.Tokenize(" ");
  TString currName = "";
  for ( Int_t iopt=0; iopt<optArr->GetEntries(); iopt++ ) {
    currName = optArr->At(iopt)->GetName();
    if ( currName.Contains("-B-") ) minBiasTrig = currName;
  }
  delete optArr;

  furtherOpt.ToUpper();
  
  AliCFContainer* cfContainer = static_cast<AliCFContainer*> ( GetSum(physSel,trigClassName,centralityRange,"SingleMuContainer") );
  if ( ! cfContainer ) return;
  
  AliCFEffGrid* effSparse = new AliCFEffGrid(Form("eff%s", cfContainer->GetName()),Form("Efficiency %s", cfContainer->GetTitle()),*cfContainer);
  effSparse->CalculateEfficiency(kStepReconstructed, kStepGeneratedMC);
  
  AliCFGridSparse* gridSparseArray[3] = {effSparse->GetNum(), effSparse->GetDen(), effSparse};
  TString gridSparseName[3] = {cfContainer->GetStepTitle(kStepReconstructed), cfContainer->GetStepTitle(kStepGeneratedMC), "Efficiency"};

  Int_t srcColors[kNtrackSources] = {kBlack, kRed, kSpring, kTeal, kBlue, kViolet, kMagenta, kOrange};
//  TString allSrcNames = "";
//  for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
//    if ( ! allSrcNames.IsNull() ) allSrcNames.Append(" ");
//    allSrcNames += fSrcKeys->At(isrc)->GetName();
//  }

  TCanvas* can = 0x0;
  Int_t xshift = 100;
  Int_t yshift = 100;
  Int_t igroup1 = -1;
  Int_t igroup2 = 0;
  
  Bool_t isMC = furtherOpt.Contains("MC");
  Int_t firstSrc = ( isMC ) ? 0 : kUnidentified;
  Int_t lastSrc  = ( isMC ) ? kNtrackSources - 1 : kUnidentified;
  if ( ! isMC ) srcColors[kUnidentified] = 1;

  TString histoName = "", histoPattern = "", drawOpt = "";
  ////////////////
  // Kinematics //
  ////////////////
  TCanvas* canKine[3] = {0x0, 0x0, 0x0};
  TLegend* legKine[3] = {0x0, 0x0, 0x0};
  for ( Int_t isrc = firstSrc; isrc <= lastSrc; ++isrc ) {
    for ( Int_t icharge=0; icharge<2; ++icharge ) {        
      for ( Int_t igrid=0; igrid<3; ++igrid ) {
        if ( gridSparseArray[igrid]->GetEntries() == 0. ) break;
        if ( gridSparseArray[igrid]->IsA() != AliCFEffGrid::Class() ) {
          SetSparseRange(gridSparseArray[igrid], kHvarEta, "", -3.999, -2.501);
          SetSparseRange(gridSparseArray[igrid], kHvarMotherType, "", isrc+1, isrc+1, "USEBIN");
          SetSparseRange(gridSparseArray[igrid], kHvarCharge, "", icharge+1, icharge+1, "USEBIN");
        }
        if ( ! canKine[igrid] ) {
          igroup1++;
          igroup2 = 0;
          currName = Form("%s_proj_%s", GetName(), gridSparseName[igrid].Data());
          canKine[igrid] = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
          canKine[igrid]->Divide(2,2);
          legKine[igrid] = new TLegend(0.6, 0.6, 0.8, 0.8);
          igroup2++;
        }
        for ( Int_t iproj=0; iproj<4; ++iproj ) {
          canKine[igrid]->cd(iproj+1);
          if ( ( iproj == kHvarPt || iproj == kHvarVz ) && gridSparseArray[igrid]->IsA() != AliCFEffGrid::Class() ) gPad->SetLogy();
          TH1* projHisto = gridSparseArray[igrid]->Project(iproj);
          projHisto->SetName(Form("proj%i_%s_src%i_charge%i", iproj, gridSparseName[igrid].Data(), isrc, icharge));
          if ( projHisto->GetEntries() == 0 ) continue;
          Bool_t isFirst = ( gPad->GetListOfPrimitives()->GetEntries() == 0 );
          drawOpt = isFirst ? "e" : "esames";
          //if ( isrc == kUnidentified && ! drawOpt.Contains("same") ) isMC = kFALSE;
          //if ( ! isMC ) srcColors[kUnidentified] = 1;
          projHisto->SetLineColor(srcColors[isrc]);
          projHisto->SetMarkerColor(srcColors[isrc]);
          projHisto->SetMarkerStyle(20+4*icharge);
          projHisto->Draw(drawOpt.Data());
          gPad->Update();
          TPaveStats* paveStats = (TPaveStats*)projHisto->FindObject("stats");
          if ( paveStats ) paveStats->SetTextColor(srcColors[isrc]);
          if ( iproj == 0 ) {
            TString legEntry = fChargeKeys->At(icharge)->GetName();
            if ( isMC ) legEntry += Form(" %s", fSrcKeys->At(isrc)->GetName());
            legKine[igrid]->AddEntry(projHisto,legEntry.Data(), "lp");
          }
        } // loop on grid sparse
      } // loop on projections
    } // loop on mu charge
  } // loop on track sources
  
  
  for ( Int_t igrid=0; igrid<3; igrid++ ) {
    if ( ! canKine[igrid] ) continue;
    canKine[igrid]->cd(1);
    legKine[igrid]->Draw("same");
    if ( gridSparseArray[igrid]->IsA() == AliCFEffGrid::Class() ) continue;
    SetSparseRange(gridSparseArray[igrid], kHvarCharge, "", 1, gridSparseArray[igrid]->GetAxis(kHvarCharge)->GetNbins(), "USEBIN"); // Reset range
  } // loop on container steps
  
  
  //////////////////////
  // Event statistics //
  //////////////////////
  printf("\nTotal analyzed events:\n");
  TString evtSel = Form("trigger:%s", trigClassName.Data());
  fEventCounters->PrintSum(evtSel.Data());
  printf("Physics selected analyzed events:\n");
  evtSel = Form("trigger:%s/selected:yes", trigClassName.Data());
  fEventCounters->PrintSum(evtSel.Data());
  
  TString countPhysSel = "any";
  if ( physSel.Contains(fPhysSelKeys->At(kPhysSelPass)->GetName()) ) countPhysSel = "yes";
  else if ( physSel.Contains(fPhysSelKeys->At(kPhysSelReject)->GetName()) ) countPhysSel="no";
  countPhysSel.Prepend("selected:");
  printf("Analyzed events vs. centrality:\n");
  evtSel = Form("trigger:%s/%s", trigClassName.Data(), countPhysSel.Data());
  fEventCounters->Print("centrality",evtSel.Data(),kTRUE);
  
  
  ///////////////////
  // Vertex method //
  ///////////////////
  if ( ! furtherOpt.Contains("VERTEX") ) return;
  Int_t firstMother = kUnidentified, lastMother = kUnidentified;
  igroup1++;
  TH1* eventVertex = (TH1*)GetSum(physSel, minBiasTrig, centralityRange, "hIpVtx");
  if ( ! eventVertex ) return;
  Double_t minZ = -9.99, maxZ = 9.99;
  Double_t meanZ = 0., sigmaZ = 4.;
  Double_t nSigma = 2.;
  TString fitOpt = "R0S";
  Bool_t fixFitRange = kFALSE;
  TString fitFormula = Form("[0]+[1]*(x+[2])");
    
  // Get vertex shape    
  if ( eventVertex->GetSumw2N() == 0 ) eventVertex->Sumw2();
  Double_t eventVtxIntegral = eventVertex->Integral(0,eventVertex->GetNbinsX()+1); // Include under/overflow
  printf("Event vertex integral %.0f\n\n", eventVtxIntegral);
  if ( eventVtxIntegral <= 0. ) return;
  eventVertex->Scale(1./eventVtxIntegral);
  printf("\nFit MB vertex\n");
  eventVertex->Fit("gaus",fitOpt.Data(),"",minZ,maxZ);
  TF1* vtxFit = (TF1*)eventVertex->GetListOfFunctions()->FindObject("gaus");
  currName = "vtxIntegrated";
  can = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
  can->SetLogy();
  eventVertex->Draw();
  vtxFit->Draw("same");

  
  enum {kRecoHF, kRecoBkg, kInputHF, kInputDecay, kRecoAll, kNrecoHistos};
  TString baseRecoName[kNrecoHistos] = {"RecoHF", "RecoBkg", "InputHF", "InputDecay", "RecoAll"};
  TArrayI sumMothers[kNrecoHistos];
  sumMothers[kRecoHF].Set(0);
  sumMothers[kRecoBkg].Set(0);
  sumMothers[kInputHF].Set(3);
  sumMothers[kInputHF][0] = kCharmMu;
  sumMothers[kInputHF][1] = kBeautyMu;
  sumMothers[kInputHF][2] = kQuarkoniumMu;
  sumMothers[kInputDecay].Set(1);
  sumMothers[kInputDecay][0] = kDecayMu;
  sumMothers[kRecoAll].Set(kNtrackSources);
  for ( Int_t isrc=0; isrc<kNtrackSources; ++isrc ) {
    sumMothers[kRecoAll][isrc] = isrc;
  }
  
  meanZ = vtxFit->GetParameter(1);
  sigmaZ = vtxFit->GetParameter(2);
  
  Double_t minZfit = ( fixFitRange ) ? minZ : meanZ - nSigma*sigmaZ;
  Double_t maxZfit = ( fixFitRange ) ? maxZ : meanZ + nSigma*sigmaZ;
  
  TF1* fitFunc = new TF1("fitFunc", fitFormula.Data(), minZ, maxZ);
  fitFunc->SetLineColor(2);
  fitFunc->SetParNames("Line norm", "Line slope", "Free path");
  const Double_t kFreePath = 153.; // 150.; // 130.; // cm
  //fitFunc->SetParameters(0.,1.);
  fitFunc->FixParameter(2, kFreePath);

  AliCFGridSparse* gridSparse = cfContainer->GetGrid(kStepReconstructed);
  TAxis* ptAxis = gridSparse->GetAxis(kHvarPt);
  
  Double_t slope = 0.;
  Double_t limitNorm = 0., limitSlope = 0.;
  Int_t firstPtBin = 0, lastPtBin = 0;
  
  gStyle->SetOptFit(1111);

  for ( Int_t itheta=0; itheta<kNthetaAbs; ++itheta ) {
    igroup2++;
    SetSparseRange(gridSparse, kHvarThetaAbs, "", itheta+1, itheta+1, "USEBIN");
    SetSparseRange(gridSparse, kHvarPt, "", 1, ptAxis->GetNbins(), "USEBIN");
    TH1* recoHisto[kNrecoHistos];
    for ( Int_t ireco=0; ireco<kNrecoHistos; ++ireco ) {
      recoHisto[ireco] = gridSparse->Project(kHvarPt);
      histoName = Form("%sMuon_%s", baseRecoName[ireco].Data(), fThetaAbsKeys->At(itheta)->GetName());
      recoHisto[ireco]->SetName(histoName.Data());
      recoHisto[ireco]->SetTitle(histoName.Data());
      recoHisto[ireco]->Reset();
      recoHisto[ireco]->Sumw2();
      for ( Int_t isrc=0; isrc<sumMothers[ireco].GetSize(); ++isrc ) {
        SetSparseRange(gridSparse, kHvarMotherType, "", sumMothers[ireco][isrc]+1, sumMothers[ireco][isrc]+1, "USEBIN");
        TH1* auxHisto = gridSparse->Project(kHvarPt);
        recoHisto[ireco]->Add(auxHisto);
        delete auxHisto;
      }
    }
    SetSparseRange(gridSparse, kHvarMotherType, "", firstMother+1, lastMother+1, "USEBIN");
    Int_t currDraw = 0;

    for ( Int_t ibinpt=0; ibinpt<=ptAxis->GetNbins(); ++ibinpt ) {
      firstPtBin = ibinpt;
      lastPtBin = ( ibinpt == 0 ) ? ptAxis->GetNbins() : ibinpt;
      SetSparseRange(gridSparse, kHvarPt, "", firstPtBin, lastPtBin, "USEBIN");
      TH1* histo = gridSparse->Project(kHvarVz);
      histo->SetName(Form("hVtx_%s_%s_ptBin%i", cfContainer->GetStepTitle(kStepReconstructed), fThetaAbsKeys->At(itheta)->GetName(), ibinpt));
      if ( histo->GetEntries() < 100. ) break;
      printf("\nFit %.2f < pt < %.2f (entries %g)\n", ptAxis->GetBinLowEdge(firstPtBin), ptAxis->GetBinUpEdge(lastPtBin), histo->GetEntries());
      histo->Divide(eventVertex);
      Double_t norm = histo->GetBinContent(histo->FindBin(0.));
      histo->GetYaxis()->SetTitle("#frac{dN_{#mu}}{dv_{z}} / #left(#frac{1}{N_{MB}}#frac{dN_{MB}}{dv_{z}}#right)");
      slope = ( histo->GetBinContent(histo->FindBin(meanZ+sigmaZ)) - 
               histo->GetBinContent(histo->FindBin(meanZ-sigmaZ)) ) / ( 2. * sigmaZ );
      
      if ( slope < 0. ) slope = norm/kFreePath;
      
      // Try to fit twice: it fit fails the first time
      // set some limits on parameters
      for ( Int_t itry=0; itry<2; itry++ ) {
        fitFunc->SetParameter(0, norm);
        fitFunc->SetParameter(1, slope);
        if ( itry > 0 ) {
          limitNorm = 2.*histo->Integral();
          limitSlope = 2.*histo->Integral()/kFreePath;
          //fitFunc->SetParLimits(0, 0., limitNorm); // REMEMBER TO CHECK
          fitFunc->SetParLimits(1, 0., limitSlope); // REMEMBER TO CHECK
          printf("Norm 0. < %f < %f  slope  0. < %f < %f\n", norm, limitNorm, slope, limitSlope);
        }
        TFitResultPtr fitRes = histo->Fit(fitFunc, fitOpt.Data(), "", minZfit, maxZfit);
        
//      if ( gMinuit->fCstatu.Contains("CONVERGED") && 
        if ( ((Int_t)fitRes) == 0 &&
            fitFunc->GetParameter(0) > 0. && 
            fitFunc->GetParameter(1) > 0. )
          break;
        else if ( furtherOpt.Contains("REFIT") ) printf("Re-fit with limits\n");
        else {
          printf("Warning: fit problems !!!\n");
          break;
        }
      }
      
      Double_t p0 = fitFunc->GetParameter(0);
      Double_t p0err = fitFunc->GetParError(0);
      Double_t p1 = fitFunc->GetParameter(1);
      Double_t p1err = fitFunc->GetParError(1);
      
      Double_t nullVz = ( p1 != 0. ) ? -p0/p1 : 0.;
      Double_t nullVzErr = ( p0 != 0. && p1 != 0. ) ? TMath::Abs(nullVz) * TMath::Sqrt(p0err*p0err/(p0*p0) + p1err*p1err/(p1*p1) ) : 0.;
      
      printf("Null value at %f +- %f\n", nullVz - kFreePath, nullVzErr);
      
      recoHisto[kRecoHF]->SetBinContent(ibinpt, p0);
      recoHisto[kRecoHF]->SetBinError(ibinpt, p0err);
      recoHisto[kRecoBkg]->SetBinContent(ibinpt, ( kFreePath + meanZ ) * p1);
      recoHisto[kRecoBkg]->SetBinError(ibinpt, ( kFreePath + meanZ ) * p1err);
      if ( currDraw%4 == 0 ){
        currName = Form("vtx_%s_PtBin%i",fThetaAbsKeys->At(itheta)->GetName(), ibinpt);
        can = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
        can->Divide(2,2);
      }
      can->cd( currDraw%4 + 1 );
      can->SetLogy();
      histo->Draw();
      fitFunc->DrawCopy("same");
      currDraw++;
    } // loop on pt bins
    SetSparseRange(gridSparse, kHvarMotherType, "", firstMother+1, lastMother+1, "USEBIN");
    currName = Form("recoPt_%s",fThetaAbsKeys->At(itheta)->GetName());
    can = new TCanvas(currName.Data(),currName.Data(),(igroup1+1)*xshift,igroup2*yshift,600,600);
    TLegend* leg = new TLegend(0.6, 0.6, 0.8, 0.8);
    drawOpt = "e";
    for ( Int_t ireco=0; ireco<kNrecoHistos-1; ++ireco ) {
      if ( recoHisto[ireco]->GetEntries() == 0. ) continue;
      TH1* ratio = (TH1*)recoHisto[ireco]->Clone(Form("%s_ratio", recoHisto[ireco]->GetName()));
      ratio->Divide(recoHisto[kRecoAll]);
      ratio->SetLineColor(srcColors[ireco]);
      ratio->SetMarkerColor(srcColors[ireco]);
      ratio->SetMarkerStyle(20+ireco);
      ratio->GetYaxis()->SetTitle("fraction of total");
      ratio->Draw(drawOpt.Data());
      leg->AddEntry(ratio,baseRecoName[ireco].Data(), "lp");
      drawOpt = "esame";
    }
    leg->Draw("same");
  } // loop on theta abs
}
