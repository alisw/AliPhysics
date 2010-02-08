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
/// \class AliAnalysisTaskSingleMu
/// Analysis task for single muons in the spectrometer.
/// The output is a list of histograms.
/// The macro class can run on AODs or ESDs.
/// In the latter case a flag can be activated to produce a tree as output.
/// If Monte Carlo information is present, some basics checks are performed.
///
/// \author Diego Stocco
//-----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//    Implementation of the single muon analysis class
// An example of usage can be found in the macro RunSingleMuAnalysisFromAOD.C.
//----------------------------------------------------------------------------

#define AliAnalysisTaskSingleMu_cxx

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
#include "TList.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
#include "TTimeStamp.h"

// STEER includes
#include "AliLog.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"

#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliAnalysisManager.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSingleMu.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSingleMu) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu(const char *name, Bool_t fillTree, Bool_t keepAll) :
  AliAnalysisTaskSE(name),
  fUseMC(0),
  fFillTree(fillTree),
  fKeepAll(keepAll),
  fHistoList(0),
  fHistoListMC(0),
  fTreeSingleMu(0),
  fTreeSingleMuMC(0),
  fVarFloat(0),
  fVarInt(0),
  fVarChar(0),
  fVarUInt(0),
  fVarFloatMC(0),
  fVarIntMC(0)
{
  //
  /// Constructor.
  //
  if ( ! fFillTree )
    fKeepAll = kFALSE;

  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  if ( fFillTree )
    DefineOutput(3, TTree::Class());
}


//________________________________________________________________________
AliAnalysisTaskSingleMu::~AliAnalysisTaskSingleMu()
{
  //
  /// Destructor
  //
  delete fHistoList;
  delete fHistoListMC;
  delete fTreeSingleMu;
  delete fTreeSingleMuMC;
  delete fVarFloat;
  delete fVarInt;
  delete [] fVarChar;
  delete fVarUInt;
  delete fVarFloatMC;
  delete fVarIntMC;
}


//___________________________________________________________________________
void AliAnalysisTaskSingleMu::NotifyRun()
{
  //
  /// Use the event handler information to correctly fill the analysis flags:
  /// - check if Monte Carlo information is present
  //  
  
  if ( MCEvent() ) {
    fUseMC = kTRUE;
    AliInfo("Monte Carlo information is present");
  }
  else {
    AliInfo("No Monte Carlo information in run");
  }
}


//___________________________________________________________________________
void AliAnalysisTaskSingleMu::UserCreateOutputObjects() 
{
  //
  /// Create output histograms
  //
  AliInfo(Form("   CreateOutputObjects of task %s\n", GetName()));

  if ( fFillTree ) OpenFile(1);

  // initialize histogram lists
  if(!fHistoList) fHistoList = new TList();
  if(!fHistoListMC) fHistoListMC = new TList();

  // Init variables
  fVarFloat = new Float_t [kNvarFloat];
  fVarInt = new Int_t [kNvarInt];
  fVarChar = new Char_t *[kNvarChar];
  fVarUInt = new UInt_t [kNvarUInt];
  fVarFloatMC = new Float_t [kNvarFloatMC];
  fVarIntMC = new Int_t [kNvarIntMC];
  
  Int_t nPtBins = 60;
  Float_t ptMin = 0., ptMax = 30.;
  TString ptName("Pt"), ptTitle("p_{t}"), ptUnits("GeV/c");
  
  Int_t nDcaBins = 100;
  Float_t dcaMin = 0., dcaMax = 200.;
  TString dcaName("DCA"), dcaTitle("DCA"), dcaUnits("cm");

  Int_t nVzBins = 60;
  Float_t vzMin = -30, vzMax = 30;
  TString vzName("Vz"), vzTitle("Vz"), vzUnits("cm");
  
  Int_t nEtaBins = 25;
  Float_t etaMin = -4.5, etaMax = -2.;
  TString etaName("Eta"), etaTitle("#eta"), etaUnits("a.u.");
  
  Int_t nRapidityBins = 25;
  Float_t rapidityMin = -4.5, rapidityMax = -2.;
  TString rapidityName("Rapidity"), rapidityTitle("rapidity"), rapidityUnits("a.u.");

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
  srcName[kUnknownPart] = "Unknown";

  TH1F* histo1D = 0x0;
  TH2F* histo2D = 0x0;
  TString histoName, histoTitle;
  Int_t histoIndex = 0;

  // 1D histos
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%sDistrib%sTrig", ptName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s distribution. Trigger: %s", ptTitle.Data(), trigName[itrig].Data());
    histo1D = new TH1F(histoName.Data(), histoTitle.Data(), nPtBins, ptMin, ptMax);
    histo1D->GetXaxis()->SetTitle(Form("%s (%s)", ptTitle.Data(), ptUnits.Data()));
    histoIndex = GetHistoIndex(kHistoPt, itrig);
    fHistoList->AddAt(histo1D, histoIndex);
  }
  
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%sDistrib%sTrig", dcaName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s distribution. Trigger: %s", dcaTitle.Data(), trigName[itrig].Data());
    histo1D = new TH1F(histoName.Data(), histoTitle.Data(), nDcaBins, dcaMin, dcaMax);
    histo1D->GetXaxis()->SetTitle(Form("%s (%s)", dcaTitle.Data(), dcaUnits.Data()));
    histoIndex = GetHistoIndex(kHistoDCA, itrig);
    fHistoList->AddAt(histo1D, histoIndex);
  }
  
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%sDistrib%sTrig", vzName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s distribution. Trigger: %s", vzTitle.Data(), trigName[itrig].Data());
    histo1D = new TH1F(histoName.Data(), histoTitle.Data(), nVzBins, vzMin, vzMax);
    histo1D->GetXaxis()->SetTitle(Form("%s (%s)", vzTitle.Data(), vzUnits.Data()));
    histoIndex = GetHistoIndex(kHistoVz, itrig);
    fHistoList->AddAt(histo1D, histoIndex);
  }
  
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%sDistrib%sTrig", etaName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s distribution. Trigger: %s", etaTitle.Data(), trigName[itrig].Data());
    histo1D = new TH1F(histoName.Data(), histoTitle.Data(), nEtaBins, etaMin, etaMax);
    histo1D->GetXaxis()->SetTitle(Form("%s (%s)", etaTitle.Data(), etaUnits.Data()));
    histoIndex = GetHistoIndex(kHistoEta, itrig);
    fHistoList->AddAt(histo1D, histoIndex);
  }
  
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%sDistrib%sTrig", rapidityName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s distribution. Trigger: %s", rapidityTitle.Data(), trigName[itrig].Data());
    histo1D = new TH1F(histoName.Data(), histoTitle.Data(), nRapidityBins, rapidityMin, rapidityMax);
    histo1D->GetXaxis()->SetTitle(Form("%s (%s)", rapidityTitle.Data(), rapidityUnits.Data()));
    histoIndex = GetHistoIndex(kHistoRapidity, itrig);
    fHistoList->AddAt(histo1D, histoIndex);
  }

  // 2D histos
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%s%sDistrib%sTrig", dcaName.Data(), ptName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s vs. %s distribution. Trigger: %s", dcaTitle.Data(), ptTitle.Data(), trigName[itrig].Data());
    histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
		       nPtBins, ptMin, ptMax,
		       nDcaBins, dcaMin, dcaMax);
    histo2D->GetXaxis()->SetTitle(Form("%s (%s)", ptTitle.Data(), ptUnits.Data()));
    histo2D->GetYaxis()->SetTitle(Form("%s (%s)", dcaTitle.Data(), dcaUnits.Data()));
    histoIndex = GetHistoIndex(kHistoPtDCA, itrig);
    fHistoList->AddAt(histo2D, histoIndex);
  }
  
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%s%sDistrib%sTrig", vzName.Data(), ptName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s vs. %s distribution. Trigger: %s", vzTitle.Data(), ptTitle.Data(), trigName[itrig].Data());
    histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
		       nPtBins, ptMin, ptMax,
		       nVzBins, vzMin, vzMax);
    histo2D->GetXaxis()->SetTitle(Form("%s (%s)", ptTitle.Data(), ptUnits.Data()));
    histo2D->GetYaxis()->SetTitle(Form("%s (%s)", vzTitle.Data(), vzUnits.Data()));
    histoIndex = GetHistoIndex(kHistoPtVz, itrig);
    fHistoList->AddAt(histo2D, histoIndex);
  }
  
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    histoName = Form("%s%sDistrib%sTrig", rapidityName.Data(), ptName.Data(), trigName[itrig].Data());
    histoTitle = Form("%s vs. %s distribution. Trigger: %s", rapidityTitle.Data(), ptTitle.Data(), trigName[itrig].Data());
    histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
		     nPtBins, ptMin, ptMax,
		     nRapidityBins, rapidityMin, rapidityMax);
    histo2D->GetXaxis()->SetTitle(Form("%s (%s)", ptTitle.Data(), ptUnits.Data()));
    histo2D->GetYaxis()->SetTitle(Form("%s (%s)", rapidityTitle.Data(), rapidityUnits.Data()));
    histoIndex = GetHistoIndex(kHistoPtRapidity, itrig);
    fHistoList->AddAt(histo2D, histoIndex);
  }

  // Summary histos
  histoName = "histoMuonMultiplicity";
  histoTitle = "Muon track multiplicity";
  histo1D = new TH1F(histoName.Data(), histoTitle.Data(), 10, -0.5, 10-0.5);
  //histo1D->GetXaxis()->SetBinLabel(1, "All events");
  //histo1D->GetXaxis()->SetBinLabel(2, "Muon events");
  histo1D->GetXaxis()->SetTitle("# of muons");
  histoIndex = GetHistoIndex(kNhistoTypes, 0) + kHistoMuonMultiplicity;
  fHistoList->AddAt(histo1D, histoIndex);

  // MC histos
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    for (Int_t isrc = 0; isrc < kNtrackSources; isrc++) {
      histoName = Form("%sResolution%sTrig%s", ptName.Data(), trigName[itrig].Data(), srcName[isrc].Data());
      histoTitle = Form("%s resolution. Trigger: %s (%s)", ptTitle.Data(), trigName[itrig].Data(), srcName[isrc].Data());
      histo1D = new TH1F(histoName.Data(), histoTitle.Data(), 100, -5., 5.);
      histo1D->GetXaxis()->SetTitle(Form("%s^{reco} - %s^{MC} (%s)", ptTitle.Data(), ptTitle.Data(), ptUnits.Data()));
      histoIndex = GetHistoIndex(kHistoPtResolutionMC, itrig, isrc);
      fHistoListMC->AddAt(histo1D, histoIndex);
    }
  }

  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    for (Int_t isrc = 0; isrc < kNtrackSources; isrc++) {
      histoName = Form("%s%sDistrib%sTrig%s", dcaName.Data(), ptName.Data(),
		       trigName[itrig].Data(), srcName[isrc].Data());
      histoTitle = Form("%s vs. %s distribution. Trigger: %s (%s)", dcaTitle.Data(), ptTitle.Data(),
			trigName[itrig].Data(), srcName[isrc].Data());
      histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
			 nPtBins, ptMin, ptMax,
			 nDcaBins, dcaMin, dcaMax);
      histo2D->GetXaxis()->SetTitle(Form("%s (%s)", ptTitle.Data(), ptUnits.Data()));
      histo2D->GetYaxis()->SetTitle(Form("%s (%s)", dcaTitle.Data(), dcaUnits.Data()));
      histoIndex = GetHistoIndex(kHistoPtDCAMC, itrig, isrc);
      fHistoListMC->AddAt(histo2D, histoIndex);
    }
  }
  
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    for (Int_t isrc = 0; isrc < kNtrackSources; isrc++) {
      histoName = Form("%s%sDistrib%sTrig%s", vzName.Data(), ptName.Data(), 
		       trigName[itrig].Data(), srcName[isrc].Data());
      histoTitle = Form("%s vs. %s distribution. Trigger: %s (%s)", vzTitle.Data(), ptTitle.Data(),
			trigName[itrig].Data(), srcName[isrc].Data());
      histo2D = new TH2F(histoName.Data(), histoTitle.Data(),
			 nPtBins, ptMin, ptMax,
			 nVzBins, vzMin, vzMax);
      histo2D->GetXaxis()->SetTitle(Form("%s (%s)", ptTitle.Data(), ptUnits.Data()));
      histo2D->GetYaxis()->SetTitle(Form("%s (%s)", vzTitle.Data(), vzUnits.Data()));
      histoIndex = GetHistoIndex(kHistoPtVzMC, itrig, isrc);
      fHistoListMC->AddAt(histo2D, histoIndex);
    }
  }

  // Trees
  if ( fFillTree ){
    TString leavesFloat[kNvarFloat] = {"Px", "Py", "Pz", "Pt",
				       "PxAtDCA", "PyAtDCA", "PzAtDCA", "PtAtDCA",
				       "PxUncorrected", "PyUncorrected", "PzUncorrected", "PtUncorrected",
				       "XUncorrected", "YUncorrected", "ZUncorrected",
				       "XatDCA", "YatDCA", "DCA",
				       "Eta", "Rapidity", "Charge",
				       "IPVx", "IPVy", "IPVz"};
    TString leavesInt[kNvarInt] = {"MatchTrig", "IsMuon", "IsGhost", "PassPhysicsSelection"};
    TString leavesChar[kNvarChar] = {"FiredTrigClass"};
    const Int_t charWidth[kNvarChar] = {255};
    TString leavesUInt[kNvarUInt] = {"BunchCrossNum", "OrbitNum", "PeriodNum", "RunNum"};
    TString leavesFloatMC[kNvarFloatMC] = {"PxMC", "PyMC", "PzMC", "PtMC", "EtaMC", "RapidityMC", "Vx", "Vy", "Vz"};
    TString leavesIntMC[kNvarIntMC] = {"Pdg", "MotherType"};

    if ( ! fTreeSingleMu ) fTreeSingleMu = new TTree("fTreeSingleMu", "Single Mu");
    if ( ! fTreeSingleMuMC ) fTreeSingleMuMC = new TTree("fTreeSingleMuMC", "Single Mu MC");

    TTree* currTree[2] = {fTreeSingleMu, fTreeSingleMuMC};
    //TList* currList[2] = {fHistoList, fHistoListMC};

    for(Int_t itree=0; itree<2; itree++){
      for(Int_t ivar=0; ivar<kNvarFloat; ivar++){
	currTree[itree]->Branch(leavesFloat[ivar].Data(), &fVarFloat[ivar], Form("%s/F", leavesFloat[ivar].Data()));
      }
      for(Int_t ivar=0; ivar<kNvarInt; ivar++){
	currTree[itree]->Branch(leavesInt[ivar].Data(), &fVarInt[ivar], Form("%s/I", leavesInt[ivar].Data()));
      }
      for(Int_t ivar=0; ivar<kNvarChar; ivar++){
	if ( itree == 0 ){
	  fVarChar[ivar] = new Char_t [charWidth[ivar]];
	}
	TString addString = leavesChar[ivar] + "/C";
	currTree[itree]->Branch(leavesChar[ivar].Data(), fVarChar[ivar], addString.Data());
      }
      for(Int_t ivar=0; ivar<kNvarUInt; ivar++){
	currTree[itree]->Branch(leavesUInt[ivar].Data(), &fVarUInt[ivar], Form("%s/i", leavesUInt[ivar].Data()));
      }
      if ( itree==1 ){
	for(Int_t ivar=0; ivar<kNvarFloatMC; ivar++){
	  currTree[itree]->Branch(leavesFloatMC[ivar].Data(), &fVarFloatMC[ivar], Form("%s/F", leavesFloatMC[ivar].Data()));
	}
	for(Int_t ivar=0; ivar<kNvarIntMC; ivar++){
	  currTree[itree]->Branch(leavesIntMC[ivar].Data(), &fVarIntMC[ivar], Form("%s/I", leavesIntMC[ivar].Data()));
	}
      }
    } // loop on trees
  } // fillNTuple
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::UserExec(Option_t * /*option*/) 
{
  //
  /// Main loop
  /// Called for each event
  //

  AliESDEvent* esdEvent = 0x0;
  AliAODEvent* aodEvent = 0x0;
  AliMCEvent*  mcEvent  = 0x0;

  esdEvent = dynamic_cast<AliESDEvent*> (InputEvent());
  if ( ! esdEvent ){
    fFillTree = kFALSE;
    aodEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  }
  //else
  //aodEvent = AODEvent();

  if ( ! aodEvent && ! esdEvent ) {
    AliError ("AOD or ESD event not found. Nothing done!");
    return;
  }

  if ( ! fUseMC && InputEvent()->GetEventType() !=7 ) return; // Run only on physics events!

  if ( fFillTree ){
    Reset(kFALSE);
    fVarInt[kVarPassPhysicsSelection] = ((AliESDInputHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())->IsEventSelected();
    strcpy(fVarChar[kVarTrigMask], esdEvent->GetFiredTriggerClasses().Data());

    // Small workaround: in MC the bunch ID are not properly set and the timestamp is in seconds
    // So fill bunchCrossing with the read timestamp
    //    fill the orbit and period number with a timestamp created while reading the run
    TTimeStamp ts;
    fVarUInt[kVarBunchCrossNumber] = ( fUseMC ) ? esdEvent->GetTimeStamp() : esdEvent->GetBunchCrossNumber();
    fVarUInt[kVarOrbitNumber] = ( fUseMC ) ? ts.GetNanoSec() : esdEvent->GetOrbitNumber();
    fVarUInt[kVarPeriodNumber] = ( fUseMC ) ? ts.GetTime() : esdEvent->GetPeriodNumber();
    fVarUInt[kVarRunNumber] = esdEvent->GetRunNumber();

    fVarFloat[kVarIPVx] = esdEvent->GetPrimaryVertex()->GetX();
    fVarFloat[kVarIPVy] = esdEvent->GetPrimaryVertex()->GetY();
    fVarFloat[kVarIPVz] = esdEvent->GetPrimaryVertex()->GetZ();
  }

  if ( fUseMC ) mcEvent = MCEvent();

  // Object declaration
  AliMCParticle* mcTrack = 0x0;

  Int_t trackLabel = -1;

  Int_t nTracks = ( esdEvent ) ? esdEvent->GetNumberOfMuonTracks() : aodEvent->GetNTracks();

  Bool_t isGhost = kFALSE;
  Int_t nGhosts = 0, nMuons = 0;

  AliVParticle* track = 0x0;

  for (Int_t itrack = 0; itrack < nTracks; itrack++) {

    // Get variables
    if ( esdEvent ){
      if (itrack>0) Reset(kTRUE);

      track = esdEvent->GetMuonTrack(itrack);
      isGhost = ( ((AliESDMuonTrack*)track)->ContainTriggerData() && ! ((AliESDMuonTrack*)track)->ContainTrackerData() );

      // ESD only info
      fVarFloat[kVarPxUncorrected] = ( isGhost ) ? -TMath::Tan(((AliESDMuonTrack*)track)->GetThetaXUncorrected()) : ((AliESDMuonTrack*)track)->PxUncorrected();
      fVarFloat[kVarPyUncorrected] = ( isGhost ) ? -TMath::Tan(((AliESDMuonTrack*)track)->GetThetaYUncorrected()) : ((AliESDMuonTrack*)track)->PyUncorrected();
      fVarFloat[kVarPzUncorrected] = ( isGhost ) ? -1 : ((AliESDMuonTrack*)track)->PzUncorrected();
      fVarFloat[kVarPtUncorrected] = TMath::Sqrt(
	fVarFloat[kVarPxUncorrected] * fVarFloat[kVarPxUncorrected] + 
	fVarFloat[kVarPyUncorrected] * fVarFloat[kVarPyUncorrected]);

      fVarFloat[kVarXUncorrected] = ((AliESDMuonTrack*)track)->GetNonBendingCoorUncorrected();
      fVarFloat[kVarYUncorrected] = ((AliESDMuonTrack*)track)->GetBendingCoorUncorrected();
      fVarFloat[kVarZUncorrected] = ((AliESDMuonTrack*)track)->GetZUncorrected();

      fVarInt[kVarMatchTrig] = ((AliESDMuonTrack*)track)->GetMatchTrigger();

      // If is ghost fill only a partial information
      if ( isGhost ){
	fVarInt[kVarIsMuon] = 0;
	nGhosts++;
	fVarInt[kVarIsGhost] = nGhosts;

	if ( ! fUseMC ) fTreeSingleMu->Fill();
	else fTreeSingleMuMC->Fill();

	continue;
      }
    }
    else {
      track = aodEvent->GetTrack(itrack);
      if ( ! ((AliAODTrack*)track)->IsMuonTrack() )
	continue;

      //Reset(kTRUE);
    }

    // Information for tracks in tracker
    nMuons++;
    fVarInt[kVarIsMuon] = nMuons;
    fVarInt[kVarIsGhost] = 0;

    fVarFloat[kVarPt] = track->Pt();
    fVarFloat[kVarXatDCA] = (esdEvent ) ? ((AliESDMuonTrack*)track)->GetNonBendingCoorAtDCA() : ((AliAODTrack*)track)->XAtDCA();
    fVarFloat[kVarYatDCA] = (esdEvent ) ? ((AliESDMuonTrack*)track)->GetBendingCoorAtDCA() : ((AliAODTrack*)track)->YAtDCA();
    fVarFloat[kVarDCA] = TMath::Sqrt( fVarFloat[kVarXatDCA] * fVarFloat[kVarXatDCA] + fVarFloat[kVarYatDCA] * fVarFloat[kVarYatDCA] );
    fVarFloat[kVarEta] = track->Eta();
    fVarFloat[kVarRapidity] = track->Y();
    trackLabel = track->GetLabel();
    fVarInt[kVarMatchTrig] = (esdEvent ) ? ((AliESDMuonTrack*)track)->GetMatchTrigger() : ((AliAODTrack*)track)->GetMatchTrigger();

    // Fill histograms
    FillTriggerHistos(kHistoPt,       fVarInt[kVarMatchTrig], -1, fVarFloat[kVarPt]);
    FillTriggerHistos(kHistoDCA,      fVarInt[kVarMatchTrig], -1, fVarFloat[kVarDCA]);
    FillTriggerHistos(kHistoVz,       fVarInt[kVarMatchTrig], -1, fVarFloat[kVarIPVz]);
    FillTriggerHistos(kHistoEta,      fVarInt[kVarMatchTrig], -1, fVarFloat[kVarEta]);
    FillTriggerHistos(kHistoRapidity, fVarInt[kVarMatchTrig], -1, fVarFloat[kVarRapidity]);

    FillTriggerHistos(kHistoPtDCA,      fVarInt[kVarMatchTrig], -1, fVarFloat[kVarPt], fVarFloat[kVarDCA]);
    FillTriggerHistos(kHistoPtVz,       fVarInt[kVarMatchTrig], -1, fVarFloat[kVarPt], fVarFloat[kVarIPVz]);
    FillTriggerHistos(kHistoPtRapidity, fVarInt[kVarMatchTrig], -1, fVarFloat[kVarPt], fVarFloat[kVarRapidity]);

    if ( fFillTree ){
      fVarFloat[kVarPx] = track->Px();
      fVarFloat[kVarPy] = track->Py();
      fVarFloat[kVarPz] = track->Pz();
      fVarFloat[kVarPxAtDCA] = (esdEvent ) ? ((AliESDMuonTrack*)track)->PxAtDCA() : ((AliAODTrack*)track)->PxAtDCA();
      fVarFloat[kVarPyAtDCA] = (esdEvent ) ? ((AliESDMuonTrack*)track)->PyAtDCA() : ((AliAODTrack*)track)->PyAtDCA();
      fVarFloat[kVarPzAtDCA] = (esdEvent ) ? ((AliESDMuonTrack*)track)->PzAtDCA() : ((AliAODTrack*)track)->PzAtDCA();
      fVarFloat[kVarPtAtDCA] = TMath::Sqrt(fVarFloat[kVarPxAtDCA]*fVarFloat[kVarPxAtDCA] + fVarFloat[kVarPyAtDCA]*fVarFloat[kVarPyAtDCA]);

      fVarFloat[kVarCharge] = track->Charge();

      if ( ! fUseMC ) fTreeSingleMu->Fill();
    }

    // Monte Carlo part
    if ( ! fUseMC ) continue;

    fVarIntMC[kVarMotherType] = kUnknownPart;

    AliMCParticle* matchedMCTrack = 0x0;

    Int_t nMCtracks = mcEvent->GetNumberOfTracks();
    for(Int_t imctrack=0; imctrack<nMCtracks; imctrack++){
      mcTrack = (AliMCParticle*)mcEvent->GetTrack(imctrack);
      if ( trackLabel == mcTrack->GetLabel() ) {
	matchedMCTrack = mcTrack;
	break;
      }
    } // loop on MC tracks

    if ( matchedMCTrack )
      fVarIntMC[kVarMotherType] = RecoTrackMother(matchedMCTrack, mcEvent);

    AliDebug(1, Form("Found mother %i", fVarIntMC[kVarMotherType]));


    FillTriggerHistos(kHistoPtDCAMC, fVarInt[kVarMatchTrig], fVarIntMC[kVarMotherType], fVarFloat[kVarPt], fVarFloat[kVarDCA]);
    FillTriggerHistos(kHistoPtVzMC,  fVarInt[kVarMatchTrig], fVarIntMC[kVarMotherType], fVarFloat[kVarPt], fVarFloat[kVarIPVz]);

    if ( fVarIntMC[kVarMotherType] != kUnknownPart) {
      fVarFloat[kVarPtMC] = mcTrack->Pt();
      FillTriggerHistos(kHistoPtResolutionMC, fVarInt[kVarMatchTrig], fVarIntMC[kVarMotherType], fVarFloat[kVarPt] - fVarFloatMC[kVarPtMC]);
      if ( fFillTree ){
	fVarFloatMC[kVarPxMC] = mcTrack->Px();
	fVarFloatMC[kVarPyMC] = mcTrack->Py();
	fVarFloatMC[kVarPzMC] = mcTrack->Pz();
	fVarFloatMC[kVarEtaMC] = mcTrack->Eta();
	fVarFloatMC[kVarRapidityMC] = mcTrack->Y();
	fVarFloatMC[kVarVxMC] = mcTrack->Xv();
	fVarFloatMC[kVarVyMC] = mcTrack->Yv();
	fVarFloatMC[kVarVzMC] = mcTrack->Zv();
	fVarIntMC[kVarPdg] = mcTrack->PdgCode();

	fTreeSingleMuMC->Fill();
      }
    }
  } // loop on tracks

  if ( fKeepAll &&  ( ( fVarInt[kVarIsMuon] + fVarInt[kVarIsGhost] ) == 0 ) ) {
    // Fill event also if there is not muon (when explicitely required)
    if ( ! fUseMC ) fTreeSingleMu->Fill();
    else fTreeSingleMuMC->Fill();
  }

  if( strstr(fVarChar[kVarTrigMask],"MUON") && fVarInt[kVarIsMuon]==0 ) 
    printf("WARNING: Muon trigger does not match tracker!\n");

  Int_t histoIndex = GetHistoIndex(kNhistoTypes, 0) + kHistoMuonMultiplicity;
  ((TH1F*)fHistoList->At(histoIndex))->Fill(fVarInt[kVarIsMuon]);
  //if ( fVarInt[kVarIsMuon]>0 ) ((TH1F*)fHistoList->At(histoIndex))->Fill(2);
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fHistoList);
  if ( fUseMC ) PostData(2, fHistoListMC);
  if ( fFillTree ){
    if ( fUseMC ) 
      PostData(3, fTreeSingleMuMC);
    else 
      PostData(3, fTreeSingleMu);
  }
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::Terminate(Option_t *) {
  //
  /// Draw some histogram at the end.
  //
  if (!gROOT->IsBatch()) {
    fHistoList = dynamic_cast<TList*> (GetOutputData(1));
    TCanvas *c1_SingleMu = new TCanvas("c1_SingleMu","Vz vs Pt",10,10,310,310);
    c1_SingleMu->SetFillColor(10); c1_SingleMu->SetHighLightColor(10);
    c1_SingleMu->SetLeftMargin(0.15); c1_SingleMu->SetBottomMargin(0.15);
    Int_t histoIndex = GetHistoIndex(kHistoPtVz, kAllPtTrig);
    ((TH2F*)fHistoList->At(histoIndex))->DrawCopy("COLZ");
  }
}

//________________________________________________________________________
 void AliAnalysisTaskSingleMu::FillTriggerHistos(Int_t histoIndex, Int_t matchTrig, Int_t motherType,
						Float_t var1, Float_t var2, Float_t var3)
{
  //
  /// Fill all histograms passing the trigger cuts
  //

  Int_t nTrigs = TMath::Max(1, matchTrig);
  TArrayI trigToFill(nTrigs);
  trigToFill[0] = matchTrig;
  for(Int_t itrig = 1; itrig < matchTrig; itrig++){
    trigToFill[itrig] = itrig;
  }

  TList* histoList = (motherType < 0 ) ? fHistoList : fHistoListMC;

  TString className;
  for(Int_t itrig = 0; itrig < nTrigs; itrig++){
    Int_t currIndex = GetHistoIndex(histoIndex, trigToFill[itrig], motherType);
    className = histoList->At(currIndex)->ClassName();
    AliDebug(3, Form("matchTrig %i  Fill %s", matchTrig, histoList->At(currIndex)->GetName()));
    if (className.Contains("1"))
      ((TH1F*)histoList->At(currIndex))->Fill(var1);
    else if (className.Contains("2"))
      ((TH2F*)histoList->At(currIndex))->Fill(var1, var2);
    else if (className.Contains("3"))
      ((TH3F*)histoList->At(currIndex))->Fill(var1, var2, var3);
    else
      AliWarning("Histogram not filled");
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskSingleMu::RecoTrackMother(AliMCParticle* mcTrack, AliMCEvent* mcEvent)
{
  //
  /// Find track mother from kinematics
  //

  Int_t recoPdg = mcTrack->PdgCode();

  Bool_t isMuon = (TMath::Abs(recoPdg) == 13) ? kTRUE : kFALSE;

  // Track is not a muon
  if ( ! isMuon ) return kRecoHadron;

  Int_t imother = mcTrack->GetMother();

  if ( imother<0 ) 
    return kPrimaryMu; // Drell-Yan Muon

  Int_t igrandma = imother;

  AliMCParticle* motherPart = (AliMCParticle*)mcEvent->GetTrack(imother);
  Int_t motherPdg = motherPart->PdgCode();

  // Track is an heavy flavor muon
  Int_t absPdg = TMath::Abs(motherPdg);
  if(absPdg/100==5 || absPdg/1000==5) {
    return kBeautyMu;
  }
  if(absPdg/100==4 || absPdg/1000==4){
    Int_t newMother = -1;
    igrandma = imother;
    AliInfo("\nFound candidate B->C->mu. History:\n");
    mcTrack->Print();
    printf("- %2i   ", igrandma);
    motherPart->Print();
    Int_t absGrandMotherPdg = TMath::Abs(motherPart->PdgCode());
    while ( absGrandMotherPdg > 10 ) {
      igrandma = ((AliMCParticle*)mcEvent->GetTrack(igrandma))->GetMother();
      if( igrandma < 0 ) break;
      printf("- %2i   ", igrandma);
      mcEvent->GetTrack(igrandma)->Print();
      absGrandMotherPdg = TMath::Abs(((AliMCParticle*)mcEvent->GetTrack(igrandma))->PdgCode());
    }

    if (absGrandMotherPdg==5) newMother = kBeautyMu; // Charm from beauty
    else if (absGrandMotherPdg==4) newMother = kCharmMu;

    if(newMother<0) {
      AliWarning("Mother not correctly found! Set to charm!\n");
      newMother = kCharmMu;
    }

    return newMother;
  }

  Int_t nPrimaries = mcEvent->Stack()->GetNprimary();

  // Track is a bkg. muon
  if (imother<nPrimaries) {
    return kPrimaryMu;
  }
  else {
    return kSecondaryMu;
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskSingleMu::GetHistoIndex(Int_t histoTypeIndex, Int_t trigIndex, Int_t srcIndex)
 {
   //
   /// Get histogram index in the list
   //
   if ( srcIndex < 0 ) return kNtrigCuts * histoTypeIndex + trigIndex;

   return 
     kNtrackSources * kNtrigCuts * histoTypeIndex  + 
     kNtrackSources * trigIndex  + 
     srcIndex;
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::Reset(Bool_t keepGlobal)
{
  //
  /// Reset variables
  //
  Int_t lastVarFloat = ( keepGlobal ) ? kVarIPVx : kNvarFloat;
  for(Int_t ivar=0; ivar<lastVarFloat; ivar++){
    fVarFloat[ivar] = 0.;
  }
  Int_t lastVarInt = ( keepGlobal ) ? kVarPassPhysicsSelection : kNvarInt;
  for(Int_t ivar=0; ivar<lastVarInt; ivar++){
    fVarInt[ivar] = -1;
  }
  fVarInt[kVarIsMuon] = 0;
  fVarInt[kVarIsGhost] = 0;

  if ( ! keepGlobal ){
    for(Int_t ivar=0; ivar<kNvarChar; ivar++){
      //sprintf(fVarChar[ivar],"%253s","");
      sprintf(fVarChar[ivar]," ");
    }
    for(Int_t ivar=0; ivar<kNvarUInt; ivar++){
      fVarUInt[ivar] = 0;
    }
  }
  if ( fUseMC ){
    for(Int_t ivar=0; ivar<kNvarFloatMC; ivar++){
      fVarFloatMC[ivar] = 0.;
    }
    for(Int_t ivar=0; ivar<kNvarIntMC; ivar++){
      fVarIntMC[ivar] = -1;
    }
  }  
}
