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
/// The macro class can run on AOD or in the train with the ESD filter.
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

// STEER includes
#include "AliLog.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"


#include "AliAnalysisTaskSE.h"
//#include "AliAnalysisDataSlot.h"
//#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSingleMu.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSingleMu) // Class implementation in ROOT context
/// \endcond


//________________________________________________________________________
AliAnalysisTaskSingleMu::AliAnalysisTaskSingleMu(const char *name) :
  AliAnalysisTaskSE(name),
  fUseMC(0),
  fHistoList(0),
  fHistoListMC(0)
{
  //
  /// Constructor.
  //
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskSingleMu::~AliAnalysisTaskSingleMu()
{
  //
  /// Destructor
  //
  if ( fHistoList ) delete fHistoList;
  if ( fHistoListMC ) delete fHistoListMC;
}


//___________________________________________________________________________
void AliAnalysisTaskSingleMu::SetFlagsFromHandler()
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

  // initialize histogram lists
  if(!fHistoList) fHistoList = new TList();
  if(!fHistoListMC) fHistoListMC = new TList();
  
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

  // MC histos
  for (Int_t itrig=0; itrig < kNtrigCuts; itrig++) {
    for (Int_t isrc = 0; isrc < kNtrackSources; isrc++) {
      histoName = Form("%sResolution%sTrig%s", ptName.Data(), trigName[itrig].Data(), srcName[isrc].Data());
      histoTitle = Form("%s resolution. Trigger: %s (%s)", ptTitle.Data(), trigName[itrig].Data(), srcName[isrc].Data());
      histo1D = new TH1F(histoName.Data(), histoTitle.Data(), 100, -5., 5.);
      histo1D->GetXaxis()->SetTitle(Form("%s^{reco} - %s^{MC} (%s)", ptTitle.Data(), ptTitle.Data(), ptUnits.Data()));
      histoIndex = GetHistoIndex(kHistoPtResolution, itrig, isrc);
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
      histoIndex = GetHistoIndex(kHistoPtDCAType, itrig, isrc);
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
      histoIndex = GetHistoIndex(kHistoPtVzType, itrig, isrc);
      fHistoListMC->AddAt(histo2D, histoIndex);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::UserExec(Option_t * /*option*/) 
{
  //
  /// Main loop
  /// Called for each event
  //

  AliAODEvent* aodEvent = 0x0;
  AliMCEvent*  mcEvent  = 0x0;

  if ( Entry()==0 ) 
    SetFlagsFromHandler();

  aodEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! aodEvent ) aodEvent = AODEvent();
  if ( ! aodEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }

  Int_t nTracks = aodEvent->GetNTracks();

  if ( fUseMC ) mcEvent = MCEvent();

  // Object declaration
  AliAODTrack *aodTrack = 0x0;
  AliMCParticle* mcTrack = 0x0;

  const Float_t kDefaultFloat = -999.;

  Float_t pt       = kDefaultFloat;
  Float_t dca      = kDefaultFloat;
  Float_t xAtDca   = kDefaultFloat;
  Float_t yAtDca   = kDefaultFloat;
  Float_t eta      = kDefaultFloat;
  Float_t rapidity = kDefaultFloat;
  Float_t vz       = kDefaultFloat;
  Int_t trackLabel = -1;
  Int_t matchTrig = -1;

  for (Int_t itrack = 0; itrack < nTracks; itrack++) {

    // Get variables
    aodTrack = aodEvent->GetTrack(itrack);
    if ( aodTrack->GetMostProbablePID() != AliAODTrack::kMuon ) continue;
    pt = aodTrack->Pt();
    xAtDca = aodTrack->XAtDCA();
    yAtDca = aodTrack->YAtDCA();
    dca = TMath::Sqrt( xAtDca * xAtDca + yAtDca * yAtDca );
    eta = aodTrack->Eta();
    rapidity = aodTrack->Y();
    vz = aodTrack->Zv();
    trackLabel = aodTrack->GetLabel();
    matchTrig = aodTrack->GetMatchTrigger();

    // Fill histograms
    FillTriggerHistos(kHistoPt,       matchTrig, -1, pt);
    FillTriggerHistos(kHistoDCA,      matchTrig, -1, dca);
    FillTriggerHistos(kHistoVz,       matchTrig, -1, vz);
    FillTriggerHistos(kHistoEta,      matchTrig, -1, eta);
    FillTriggerHistos(kHistoRapidity, matchTrig, -1, rapidity);

    FillTriggerHistos(kHistoPtDCA,      matchTrig, -1, pt, dca);
    FillTriggerHistos(kHistoPtVz,       matchTrig, -1, pt, vz);
    FillTriggerHistos(kHistoPtRapidity, matchTrig, -1, pt, rapidity);

    // Monte Carlo part
    if (! fUseMC ) continue;

    Int_t motherType = kUnknownPart;

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
      motherType = RecoTrackMother(matchedMCTrack, mcEvent);

    AliDebug(1, Form("Found mother %i", motherType));

    if ( motherType != kUnknownPart) {
      FillTriggerHistos(kHistoPtResolution, matchTrig, motherType, pt - mcTrack->Pt());
    }
    FillTriggerHistos(kHistoPtDCAType, matchTrig, motherType, pt, dca);
    FillTriggerHistos(kHistoPtVzType,  matchTrig, motherType, pt, vz);

  } // loop on tracks
  
  // Post final data. It will be written to a file with option "RECREATE"
  PostData(1, fHistoList);
  if ( fUseMC ) PostData(2, fHistoListMC);
}

//________________________________________________________________________
void AliAnalysisTaskSingleMu::Terminate(Option_t *) {
  //
  /// Draw some histogram at the end.
  //
  if (!gROOT->IsBatch()) {
    TCanvas *c1 = new TCanvas("c1","Vz vs Pt",10,10,310,310);
    c1->SetFillColor(10); c1->SetHighLightColor(10);
    c1->SetLeftMargin(0.15); c1->SetBottomMargin(0.15);
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
   if ( srcIndex < 0 ) return kNtrigCuts * histoTypeIndex + trigIndex;

   return 
     kNtrackSources * kNtrigCuts * histoTypeIndex  + 
     kNtrackSources * trigIndex  + 
     srcIndex;
}
