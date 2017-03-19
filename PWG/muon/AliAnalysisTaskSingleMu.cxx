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
#include "TArrayD.h"
#include "TPaveStats.h"
#include "TFitResultPtr.h"
#include "TFile.h"
#include "THashList.h"

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
#include "AliMuonEventCuts.h"
#include "AliMuonTrackCuts.h"
#include "AliAnalysisMuonUtility.h"
#include "AliMuonAnalysisOutput.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSingleMu) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu() :
  AliVAnalysisMuon(),
  fThetaAbsKeys(0x0),
  fCutOnDimu(kFALSE),
  fUseMCKineForRecoTracks(kFALSE)
{
  /// Default ctor.
}

//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu(const char *name, const AliMuonTrackCuts& cuts) :
  AliVAnalysisMuon(name, cuts),
  fThetaAbsKeys(0x0),
  fCutOnDimu(kFALSE),
  fUseMCKineForRecoTracks(kFALSE)
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
  TString ptName("Pt"), ptTitle("p_{T}"), ptUnits("GeV/c");
  
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
  
  histoName = "hRecoDimu";
  TH2* histoDimu = new TH2F(histoName.Data(), histoName.Data(), 24, 0., 120., 30, 0., 150.);
  histoDimu->SetXTitle("min. p_{T}^{#mu} (GeV/c)");
  histoDimu->SetYTitle("M_{#mu#mu} (GeV/c^{2})");
  AddObjectToCollection(histoDimu, kNobjectTypes);
  
  histoName = "hGeneratedZ";
  TH2* histoGenZ = static_cast<TH2*>(histoDimu->Clone(histoName.Data()));
  histoGenZ->SetTitle(histoName.Data());
  AddObjectToCollection(histoGenZ, kNobjectTypes+1);
  
  
  histoName = "hZmuEtaCorr";
  histoDimu = new TH2F(histoName.Data(), histoName.Data(), 160, -8., 8., 160, -8., 8.);
  histoDimu->SetXTitle("#eta_{#mu-}");
  histoDimu->SetYTitle("#eta_{#mu+}");
  AddObjectToCollection(histoDimu, kNobjectTypes+2);
  
  fMuonTrackCuts->Print("mask");
  
  AliInfo(Form("Apply cut on dimuon (60<M_mumu<120 GeV/c^2) to reject Z contribution: %i", fCutOnDimu));
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality)
{
  //
  /// Fill output objects
  //

  Double_t ipVz = InputEvent()->GetPrimaryVertexSPD()->GetZ();
  Double_t ipVzMC = MCEvent() ? AliAnalysisMuonUtility::GetMCVertexZ(InputEvent(),MCEvent()) : 0.;
  
  for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
    TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
    ((TH1*)GetMergeableObject(physSel, trigClassName, centrality, "hIpVtx"))->Fill(ipVz);
  }
  
  // Bool_t isPileupFromSPD = ( fAODEvent && ! fAODEvent->GetTracklets() ) ? InputEvent()->IsPileupFromSPD(3, 0.8, 3., 2., 5.) : InputEvent()->IsPileupFromSPDInMultBins(); // Avoid break when reading Muon AODs (tracklet info is not present and IsPileupFromSPDInMultBins crashes // UNCOMMENT TO REJECT PILEUP
  // if ( isPileupFromSPD ) return; // UNCOMMENT TO REJECT PILEUP
  
  Double_t containerInput[kNvars];
  AliVParticle* track = 0x0;

  Int_t nSteps = MCEvent() ? 2 : 1;
  for ( Int_t istep = 0; istep<nSteps; ++istep ) {
    Int_t nTracks = ( istep == kStepReconstructed ) ? AliAnalysisMuonUtility::GetNTracks(InputEvent()) : MCEvent()->GetNumberOfTracks();
    
    TObjArray selectedTracks(nTracks);
    TArrayI trackSources(nTracks);
    TArrayD trackWgt(nTracks);
    trackWgt.Reset(1.);
    for (Int_t itrack = 0; itrack < nTracks; itrack++) {
      track = ( istep == kStepReconstructed ) ? AliAnalysisMuonUtility::GetTrack(itrack,InputEvent()) : MCEvent()->GetTrack(itrack);
      
      // In case of MC we usually ask that the particle is a muon
      // However, in W or Z simulations, Pythia stores both the initial muon
      // (before ISR, FSR and kt kick) and the final state one.
      // The first muon is of course there only for information and should be rejected.
      // The Pythia code for initial state particles is 21
      // When running with POWHEG, Pythia puts the hard process input of POWHEG in the stack
      // with state 21, and then re-add it to stack before applying ISR, FSR and kt kick.
      // This muon produces the final state muon, and its status code is 11
      // To avoid all problems, keep only final state muon (status code <10)
      // FIXME: is the convention valid for other generators as well?
      Bool_t isSelected = ( istep == kStepReconstructed ) ? fMuonTrackCuts->IsSelected(track) : ( TMath::Abs(track->PdgCode()) == 13 && AliAnalysisMuonUtility::GetStatusCode(track) < 10 );
      if ( ! isSelected ) continue;
      
      selectedTracks.AddAt(track,itrack);
      trackSources[itrack] = GetParticleType(track);
      
      TObject* wgtObj = GetWeight(fSrcKeys->At(trackSources[itrack])->GetName());
      
      if ( wgtObj  ) {
        AliVParticle* mcTrack = ( istep == kStepReconstructed ) ? MCEvent()->GetTrack(track->GetLabel()) : track;
        if ( wgtObj->IsA() == TF1::Class() ) trackWgt[itrack] = static_cast<TF1*>(wgtObj)->Eval(mcTrack->Pt());
        else if ( wgtObj->IsA()->InheritsFrom(TH1::Class()) ) {
          TH1* wgtHisto = static_cast<TH1*>(wgtObj);
          trackWgt[itrack] = wgtHisto->GetBinContent(wgtHisto->GetXaxis()->FindBin(mcTrack->Pt()));
        }
        AliDebug(3,Form("Apply weights %s:  pt %g  gen pt %g  weight %g",wgtObj->GetName(),track->Pt(),mcTrack->Pt(),trackWgt[itrack]));
//        Int_t iancestor = fUtilityMuonAncestor->GetAncestor(track,MCEvent());
//        AliVParticle* motherPart = MCEvent()->GetTrack(iancestor);
//        trackWgt[itrack] = beautyMuWgt->GetBinContent(beautyMuWgt->GetXaxis()->FindBin(motherPart->Pt()));
      }
    } // loop on tracks
    
    // Loop on selected tracks
    TArrayI rejectTrack(nTracks);
    for ( Int_t itrack=0; itrack<nTracks; itrack++) {
      track = static_cast<AliVParticle*>(selectedTracks.At(itrack));
      if ( ! track ) continue;
      
      // Check dimuons
      for ( Int_t jtrack=itrack+1; jtrack<nTracks; jtrack++ ) {
        AliVParticle* auxTrack = static_cast<AliVParticle*>(selectedTracks.At(jtrack));
        if ( ! auxTrack ) continue;
        if ( track->Charge() * auxTrack->Charge() >= 0 ) continue;
          
        TLorentzVector dimuPair = AliAnalysisMuonUtility::GetTrackPair(track,auxTrack);
        Double_t ptMin = TMath::Min(track->Pt(),auxTrack->Pt());
        Double_t invMass = dimuPair.M();
        if ( fCutOnDimu && invMass > 60. && invMass < 120. ) {
          rejectTrack[itrack] = 1;
          rejectTrack[jtrack] = 1;
        }
        Double_t dimuWgt = trackWgt[itrack] * trackWgt[jtrack];
        if ( istep == kStepReconstructed ) {
          for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
            TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
            if ( ! fMuonTrackCuts->TrackPtCutMatchTrigClass(track, fMuonEventCuts->GetTrigClassPtCutLevel(trigClassName)) || ! fMuonTrackCuts->TrackPtCutMatchTrigClass(auxTrack, fMuonEventCuts->GetTrigClassPtCutLevel(trigClassName)) ) continue;
            ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, "hRecoDimu"))->Fill(ptMin,invMass,dimuWgt);
          } // loop on triggers
        }
        else {
          if ( trackSources[itrack] == kZbosonMu && trackSources[jtrack] == kZbosonMu ) {
            Bool_t isAccepted = kTRUE;
            AliVParticle* muDaughter[2] = {0x0, 0x0};
            for ( Int_t imu=0; imu<2; imu++ ) {
              AliVParticle* currPart = ( imu == 0 ) ? track : auxTrack;
              if ( currPart->Charge() < 0. ) muDaughter[0] = currPart;
              else muDaughter[1] = currPart;
              if ( currPart->Eta() < -4.5 || currPart->Eta() > -2. ) {
                isAccepted = kFALSE;
              }
            } // loop on muons in the pair
              
            Double_t pairRapidity = dimuPair.Rapidity();
            if ( pairRapidity < -4. || pairRapidity > -2.5 ) isAccepted = kFALSE;
            //printf("Rapidity Z %g  pair %g\n",track->Y(), pairRapidity);
              
            for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
              TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
              if ( isAccepted ) ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, "hGeneratedZ"))->Fill(ptMin,invMass,dimuWgt);
              ((TH2*)GetMergeableObject(physSel, trigClassName, centrality, "hZmuEtaCorr"))->Fill(muDaughter[0]->Eta(),muDaughter[1]->Eta(),dimuWgt);
            } // loop on selected trig
          } // both muons from Z
        } // kStepGeneratedMC
      } // loop on auxiliary tracks
      if ( rejectTrack[itrack] > 0 ) continue;
      
      Double_t thetaAbsEndDeg = 0;
      if ( istep == kStepReconstructed ) {
        thetaAbsEndDeg = AliAnalysisMuonUtility::GetThetaAbsDeg(track);
      }
      else {
        thetaAbsEndDeg = ( TMath::Pi()-track->Theta() ) * TMath::RadToDeg();
      }
      Int_t thetaAbsBin = ( thetaAbsEndDeg < 3. ) ? kThetaAbs23 : kThetaAbs310;

      AliVParticle* recoTrack = 0x0;
      if ( fUseMCKineForRecoTracks && istep == kStepReconstructed && track->GetLabel() >= 0 ) {
        // Switch to the matching MC track in order to use its kinematic parameters
        recoTrack = track;
        track = MCEvent()->GetTrack(track->GetLabel());
      }

      containerInput[kHvarPt]         = track->Pt();
      containerInput[kHvarEta]        = track->Eta();
      containerInput[kHvarPhi]        = track->Phi();
      containerInput[kHvarVz]         = ( istep == kStepReconstructed && !fUseMCKineForRecoTracks ) ? ipVz : ipVzMC;
      containerInput[kHvarCharge]     = track->Charge()/3.;
      containerInput[kHvarThetaAbs]   = (Double_t)thetaAbsBin;
      containerInput[kHvarMotherType] = (Double_t)trackSources[itrack];

      // Switch back to the reconstructed track (to correctly fill the information for different trigger classes)
      if ( recoTrack ) track = recoTrack;
      
      for ( Int_t itrig=0; itrig<selectTrigClasses.GetEntries(); ++itrig ) {
        TString trigClassName = ((TObjString*)selectTrigClasses.At(itrig))->GetString();
        if ( istep == kStepReconstructed && ! fMuonTrackCuts->TrackPtCutMatchTrigClass(track, fMuonEventCuts->GetTrigClassPtCutLevel(trigClassName)) ) continue;
        ((AliCFContainer*)GetMergeableObject(physSel, trigClassName, centrality, "SingleMuContainer"))->Fill(containerInput,istep,trackWgt[itrack]);
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

  AliMuonAnalysisOutput muonOut(fOutputList);
  
  TString physSel = fTerminateOptions->At(0)->GetName();
  TString trigClassName = fTerminateOptions->At(1)->GetName();
  TString centralityRange = fTerminateOptions->At(2)->GetName();
  TString furtherOpt = fTerminateOptions->At(3)->GetName();
  
  TString minBiasTrig = "";
  TObjArray* optArr = furtherOpt.Tokenize(" ");
  TString currName = "";
  for ( Int_t iopt=0; iopt<optArr->GetEntries(); iopt++ ) {
    currName = optArr->At(iopt)->GetName();
    if ( currName.Contains("-B-") || currName.Contains("ANY") ) minBiasTrig = currName;
  }
  delete optArr;
  
  furtherOpt.ToUpper();
  Bool_t plotChargeAsymmetry = furtherOpt.Contains("ASYM");
  
  AliCFContainer* inContainer = static_cast<AliCFContainer*> ( muonOut.GetSum(physSel,trigClassName,centralityRange,"SingleMuContainer") );
  if ( ! inContainer ) return;
  
  AliCFContainer* cfContainer = inContainer;
  
  if ( ! furtherOpt.Contains("GENINTRIGCLASS") && trigClassName != "ANY" ) {
    // The output container contains both the reconstructed and (in MC)
    // the generated muons in a specific trigger class.
    // Since the trigger pt level of the track is matched to the trigger class,
    // analyzing the MUHigh trigger (for example) is equivalent of analyzing
    // Hpt tracks.
    // However, in this way, the generated muons distribution depend
    // on a condition on the reconstructed track.
    // If we calculate the Acc.xEff. as the ratio of reconstructed and
    // generated muons IN A TRIGGER CLASS, we are biasing the final result.
    // The correct value of Acc.xEff. is hence obtained as the distribution
    // of reconstructed tracks in a muon trigger class, divided by the
    // total number of generated class (which is in the "class" ANY).
    // The following code sets the generated muons as the one in the class ANY.
    // The feature is the default. If you want the generated muons in the same
    // trigger class as the generated tracks, use option "GENINTRIGCLASS"
    AliCFContainer* fullMCcontainer = static_cast<AliCFContainer*> ( muonOut.GetSum(physSel,"ANY",centralityRange,"SingleMuContainer") );
    if ( fullMCcontainer ) {
      cfContainer = static_cast<AliCFContainer*>(cfContainer->Clone("SingleMuContainerCombo"));
      cfContainer->SetGrid(kStepGeneratedMC,fullMCcontainer->GetGrid(kStepGeneratedMC));
    }
  }
  
  AliCFEffGrid* effSparse = new AliCFEffGrid(Form("eff%s", cfContainer->GetName()),Form("Efficiency %s", cfContainer->GetTitle()),*cfContainer);
  effSparse->CalculateEfficiency(kStepReconstructed, kStepGeneratedMC);
  
  AliCFGridSparse* gridSparseArray[3] = {effSparse->GetNum(), effSparse->GetDen(), effSparse};
  TString gridSparseName[3] = {cfContainer->GetStepTitle(kStepReconstructed), cfContainer->GetStepTitle(kStepGeneratedMC), "Efficiency"};
  
  Int_t srcColors[kNtrackSources] = {kBlack, kRed, kSpring, kTeal, kBlue, kViolet, kMagenta, kOrange, kGray};
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
  
//  Bool_t isMC = furtherOpt.Contains("MC");
  Bool_t isMC = ( gridSparseArray[1]->GetEntries() > 0. );
  Int_t nGrids = isMC ? 3 : 1;
  
  TAxis* srcAxis = gridSparseArray[0]->GetAxis(kHvarMotherType);
  Int_t unIdBin = srcAxis->FindBin(fSrcKeys->At(kUnidentified)->GetName());
  if ( unIdBin < 1 ) unIdBin = srcAxis->GetNbins();
  
  Int_t firstSrcBin = ( isMC ) ? 1 : unIdBin;
  Int_t lastSrcBin  = ( isMC ) ? srcAxis->GetNbins() - 1 : unIdBin;
  if ( ! isMC ) srcColors[unIdBin-1] = 1;
  
  ////////////////
  // Kinematics //
  ////////////////
  TString chargeNames[3] = {fChargeKeys->At(0)->GetName(), fChargeKeys->At(1)->GetName(), "Total"};
  THashList histoList[3];
  for ( Int_t icharge=0; icharge<3; icharge++ ) {
    histoList[icharge].SetName(chargeNames[icharge].Data());
  }
  for ( Int_t isrc = firstSrcBin; isrc <= lastSrcBin; ++isrc ) {
    for ( Int_t icharge=0; icharge<3; ++icharge ) {
      Int_t icharge1 = ( icharge < 2 ) ? icharge : 0;
      Int_t icharge2 = ( icharge < 2 ) ? icharge : 1;
      for ( Int_t igrid=0; igrid<nGrids; ++igrid ) {
        if ( gridSparseArray[igrid]->GetEntries() == 0. ) continue;
        if ( gridSparseArray[igrid]->IsA() != AliCFEffGrid::Class() ) {
          SetSparseRange(gridSparseArray[igrid], kHvarEta, "", -3.999, -2.501);
          SetSparseRange(gridSparseArray[igrid], kHvarMotherType, "", isrc, isrc, "USEBIN");
          SetSparseRange(gridSparseArray[igrid], kHvarCharge, "", icharge1+1, icharge2+1, "USEBIN");
        }
        TH1* histo = gridSparseArray[igrid]->Project(kHvarPt, kHvarEta);
        histo->SetName(Form("hPtEta_%s_%s_%s", gridSparseName[igrid].Data(), srcAxis->GetBinLabel(isrc), chargeNames[icharge].Data()));
        histo->SetDirectory(0);
        if ( histo->Integral() > 0 ) histoList[icharge].Add(histo);
        for ( Int_t iproj=0; iproj<4; ++iproj ) {
          histo = gridSparseArray[igrid]->Project(iproj);
          histo->SetName(Form("proj%i_%s_%s_%s", iproj, gridSparseName[igrid].Data(), srcAxis->GetBinLabel(isrc), chargeNames[icharge].Data()));
          histo->SetDirectory(0);
          if ( histo->Integral() > 0 ) histoList[icharge].Add(histo);
        } // loop on projections
      } // loop on grid sparse
    } // loop on charge
  } // loop on track sources
  
  // Get charge asymmetry or mu+/mu-
  THashList histoListRatio;
  TString basePlotName = plotChargeAsymmetry ? "ChargeAsymmetry" : "ChargeRatio";
  histoListRatio.SetName(basePlotName.Data());
  Int_t baseCharge = 1;
  Int_t auxCharge = 1-baseCharge;
  for ( Int_t ihisto=0; ihisto<histoList[baseCharge].GetEntries(); ihisto++ ) {
    TObject* obj = histoList[baseCharge].At(ihisto);
    TString histoName = obj->GetName();
    if ( histoName.Contains(gridSparseName[2].Data()) ) continue;
    TString searchName = histoName;
    searchName.ReplaceAll(fChargeKeys->At(baseCharge)->GetName(), fChargeKeys->At(auxCharge)->GetName());
    TH1* auxHisto = static_cast<TH1*> (histoList[auxCharge].FindObject(searchName.Data()));
    if ( ! auxHisto ) continue;
    histoName.ReplaceAll(fChargeKeys->At(baseCharge)->GetName(),basePlotName.Data());
    TH1* histo = static_cast<TH1*> (obj->Clone(histoName.Data()));
    if ( plotChargeAsymmetry ) {
      histo->Add(auxHisto, -1.);
      // h2 + h1 = 2xh2 + (h1-h2)
      auxHisto->Add(auxHisto, histo, 2.);
    }
    histo->Divide(auxHisto);
    TString axisTitle = plotChargeAsymmetry ? Form("(%s - %s) / (%s + %s)", fChargeKeys->At(1)->GetName(), fChargeKeys->At(0)->GetName(), fChargeKeys->At(1)->GetName(), fChargeKeys->At(0)->GetName()) : Form("%s / %s", fChargeKeys->At(1)->GetName(), fChargeKeys->At(0)->GetName());
    axisTitle.ReplaceAll("MuPlus","#mu^{+}");
    axisTitle.ReplaceAll("MuMinus","#mu^{-}");
    histo->GetYaxis()->SetTitle(axisTitle.Data());
    histo->SetStats(kFALSE);
    histo->SetDirectory(0);
    histoListRatio.Add(histo);
  }
  
  // Plot kinematics
  TString histoName = "", drawOpt = "";
  for ( Int_t itype=0; itype<3; itype++ ) {
    THashList* currList = 0x0;
    Int_t nCharges = 1;
    if ( itype == 1 ) currList = &histoListRatio;
    else if ( itype == 2 ) currList = &histoList[2];
    else nCharges = 2;
    for ( Int_t igrid=0; igrid<nGrids; ++igrid ) {
      igroup1 = igrid;
      TCanvas* canKine = 0x0;
      TLegend* legKine = 0x0;
      for ( Int_t iproj=0; iproj<4; ++iproj ) {
        for ( Int_t isrc = firstSrcBin; isrc <= lastSrcBin; ++isrc ) {
          for ( Int_t icharge=0; icharge<nCharges; ++icharge ) {
            if ( itype == 0 ) currList = &histoList[icharge];
            histoName = Form("proj%i_%s_%s_%s", iproj, gridSparseName[igrid].Data(), srcAxis->GetBinLabel(isrc), currList->GetName());
            TH1* histo = static_cast<TH1*>(currList->FindObject(histoName.Data()));
            if ( ! histo ) continue;
            if ( ! canKine ) {
              igroup2 = itype;
              igroup1 = igrid;
              currName = Form("%s_%s_%s", GetName(), currList->GetName(), gridSparseName[igrid].Data());
              canKine = new TCanvas(currName.Data(),currName.Data(),igroup1*xshift,igroup2*yshift,600,600);
              canKine->Divide(2,2);
              legKine = new TLegend(0.6, 0.6, 0.8, 0.8);
            }
            canKine->cd(iproj+1);
            if ( itype != 1 ) {
              if ( ( iproj == kHvarPt || iproj == kHvarVz ) && gridSparseArray[igrid]->IsA() != AliCFEffGrid::Class() ) gPad->SetLogy();
            }
            Bool_t isFirst = ( gPad->GetListOfPrimitives()->GetEntries() == 0 );
            drawOpt = isFirst ? "e" : "esames";
            histo->SetLineColor(srcColors[isrc-1]);
            histo->SetMarkerColor(srcColors[isrc-1]);
            histo->SetMarkerStyle(20+4*icharge);
            histo->Draw(drawOpt.Data());
            TPaveStats* paveStats = (TPaveStats*)histo->FindObject("stats");
            if ( paveStats ) paveStats->SetTextColor(srcColors[isrc-1]);
            if ( iproj == 0 ) {
              TString legEntry = ( itype == 0 ) ? fChargeKeys->At(icharge)->GetName() : "";
              if ( isMC ) legEntry += Form(" %s", srcAxis->GetBinLabel(isrc));
              if ( ! legEntry.IsNull() ) legKine->AddEntry(histo,legEntry.Data(), "lp");
            }
          } // loop on mu charge
        } // loop on track sources
      } // loop on projections
      if ( legKine && legKine->GetNRows() > 0 ) {
        canKine->cd(1);
        legKine->Draw("same");
      }
    } // loop on grid sparse
  } // loop on types
  
  
  for ( Int_t igrid=0; igrid<nGrids; igrid++ ) {
    if ( gridSparseArray[igrid]->IsA() == AliCFEffGrid::Class() ) continue;
    SetSparseRange(gridSparseArray[igrid], kHvarCharge, "", 1, gridSparseArray[igrid]->GetAxis(kHvarCharge)->GetNbins(), "USEBIN"); // Reset range
  } // loop on container steps
  
  //////////////////////
  // Event statistics //
  //////////////////////
  printf("\nTotal analyzed events:\n");
  TString evtSel = Form("trigger:%s", trigClassName.Data());
  muonOut.GetCounterCollection()->PrintSum(evtSel.Data());
  printf("Physics selected analyzed events:\n");
  evtSel = Form("trigger:%s/selected:yes", trigClassName.Data());
  muonOut.GetCounterCollection()->PrintSum(evtSel.Data());
  
  TString countPhysSel = "any";
  if ( physSel.Contains(fPhysSelKeys->At(kPhysSelPass)->GetName()) ) countPhysSel = "yes";
  else if ( physSel.Contains(fPhysSelKeys->At(kPhysSelReject)->GetName()) ) countPhysSel="no";
  countPhysSel.Prepend("selected:");
  printf("Analyzed events vs. centrality:\n");
  evtSel = Form("trigger:%s/%s", trigClassName.Data(), countPhysSel.Data());
  muonOut.GetCounterCollection()->Print("centrality",evtSel.Data(),kTRUE);
  
  TString outFilename = Form("/tmp/out%s.root", GetName());
  printf("\nCreating file %s\n", outFilename.Data());
  TFile* outFile = new TFile(outFilename.Data(),"RECREATE");
  for ( Int_t icharge=0; icharge<3; icharge++ ) {
    histoList[icharge].Write();
  }
  histoListRatio.Write();
  outFile->Close();
  
  ///////////////////
  // Vertex method //
  ///////////////////
  if ( ! furtherOpt.Contains("VERTEX") ) return;
  igroup1++;
  TH1* eventVertex = static_cast<TH1*>(muonOut.GetSum(physSel, minBiasTrig, centralityRange, "hIpVtx"));
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
  
  sumMothers[kInputHF][0] = srcAxis->FindBin(fSrcKeys->At(kCharmMu)->GetName());
  sumMothers[kInputHF][1] = srcAxis->FindBin(fSrcKeys->At(kBeautyMu)->GetName());
  sumMothers[kInputHF][2] = srcAxis->FindBin(fSrcKeys->At(kQuarkoniumMu)->GetName());
  sumMothers[kInputDecay].Set(1);
  sumMothers[kInputDecay][0] = srcAxis->FindBin(fSrcKeys->At(kDecayMu)->GetName());
  sumMothers[kRecoAll].Set(srcAxis->GetNbins());
  for ( Int_t isrc=0; isrc<srcAxis->GetNbins(); ++isrc ) {
    sumMothers[kRecoAll][isrc] = isrc+1;
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
        SetSparseRange(gridSparse, kHvarMotherType, "", sumMothers[ireco][isrc], sumMothers[ireco][isrc], "USEBIN");
        TH1* auxHisto = gridSparse->Project(kHvarPt);
        recoHisto[ireco]->Add(auxHisto);
        delete auxHisto;
      }
    }
    SetSparseRange(gridSparse, kHvarMotherType, "", firstPtBin, lastSrcBin, "USEBIN");
    Int_t currDraw = 0;
    
    for ( Int_t ibinpt=0; ibinpt<=ptAxis->GetNbins(); ++ibinpt ) {
      firstPtBin = ibinpt;
      lastPtBin = ( ibinpt == 0 ) ? ptAxis->GetNbins() : ibinpt;
      SetSparseRange(gridSparse, kHvarPt, "", firstPtBin, lastPtBin, "USEBIN");
      TH1* histo = gridSparse->Project(kHvarVz);
      histo->SetName(Form("hVtx_%s_%s_ptBin%i", cfContainer->GetStepTitle(kStepReconstructed), fThetaAbsKeys->At(itheta)->GetName(), ibinpt));
      if ( histo->Integral() < 100. ) break;
      printf("\nFit %.2f < pt < %.2f (entries %g)\n", ptAxis->GetBinLowEdge(firstPtBin), ptAxis->GetBinUpEdge(lastPtBin), histo->GetEntries());
      histo->Divide(eventVertex);
      Double_t norm = histo->GetBinContent(histo->FindBin(0.));
      histo->GetYaxis()->SetTitle("#frac{dN_{#mu}}{dv_{z}} / #left(#frac{1}{N_{MB}}#frac{dN_{MB}}{dv_{z}}#right)");
      histo->SetTitle(Form("%g < p_{T} (GeV/c) < %g",ptAxis->GetBinLowEdge(firstPtBin), ptAxis->GetBinUpEdge(lastPtBin)));
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
    SetSparseRange(gridSparse, kHvarMotherType, "", firstSrcBin, lastSrcBin, "USEBIN");
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
