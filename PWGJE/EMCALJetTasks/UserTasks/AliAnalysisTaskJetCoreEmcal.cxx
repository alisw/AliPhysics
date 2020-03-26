/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <TClonesArray.h>
#include <TList.h>
#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "AliLog.h"

#include <AliAnalysisManager.h>
#include <AliVEventHandler.h>
#include <AliVEvent.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalPythiaInfo.h"

#include "AliAnalysisTaskJetCoreEmcal.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetCoreEmcal);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskJetCoreEmcal::AliAnalysisTaskJetCoreEmcal() : 
  AliAnalysisTaskEmcalJet(),
  fHistManager(),
	fJetShapeType(AliAnalysisTaskJetCoreEmcal::kData),
	fCentMin(0.),
	fCentMax(100.),
	fTTLowRef(8.),
	fTTUpRef(9.),
	fTTLowSig(20.),
	fTTUpSig(50.),
	fNRPBins(50),
	fFrac(0.8),
	fJetHadronDeltaPhi(0.6),
	fMinFractionSharedPt(0.5),
	fMinEmbJetPt(15.),
	fJetContName(""),
	fJetContTrueName(""),
	fJetContPartName(""),
	fFillTrackHistograms(kTRUE),
	fFillJetHistograms(kTRUE),
	fFillRecoilTHnSparse(kTRUE),
	fFillInclusiveTree(kFALSE),
	fFillRecoilTree(kFALSE),
	fMoreTreeVars(kFALSE),
	fPtHardBin(0.),
	fRejectionFactorInclusiveJets(1),
	fRandom(0),
	fHistEvtSelection(0x0), 
	fHJetSpec(0x0),
	fh1TrigRef(0x0),
	fh1TrigSig(0x0),
	fh2Ntriggers(0x0),
	fhRhoCentSig(0x0),
	fhRhoCentRef(0x0),
	fhDphiPtSigPi(0x0),
	fhDphiPtSig(0x0),
	fhDphiPtRefPi(0x0),
	fhDphiPtRef(0x0),
	fhPtDetPart(0x0),
	fhPtHybrDet(0x0),
	fhPtHybrPart(0x0),
	fhPtHybrPartCor(0x0),
	fhPhiHybrPartCor(0x0),
	fhPtDet(0x0),
	fhPtDetMatchedToPart(0x0),
	fhPtPartMatched(0x0),
	fhPtPartMatchedCent(0x0),
	fhPtPartMatchedWrongCent(0x0),
	fhResidual(0x0),
	fhPtResidual(0x0),
	fhPhiResidual(0x0),
	fhPhiPhiResidual(0x0),
	fhPtDetPartRecoil(0x0),
	fhPtHybrDetRecoil(0x0),
	fhPtHybrPartRecoil(0x0),
	fhPtHybrPartCorRecoil(0x0),
	fhPtDetRecoil(0x0),
	fhPtDetMatchedToPartRecoil(0x0),
	fhResidualRecoil(0x0),
	fhPtResidualRecoil(0x0),
	fhDphiResidualRecoil(0x0),
	fhDphiphiResidualRecoil(0x0),
	fhTTPtDetMatchedToPart(0x0),
	fhTTPhiDetMatchedToPart(0x0),
	fhDPhiHybrPartCorRecoil(0x0),
	fhSelectedTrigger(0x0),
	fhFractionSharedPtInclusive(0x0),
	fhFractionSharedPtRecoil(0x0),
	fTreeEmbInclusive(0x0),
	fTreeEmbRecoil(0x0)
{
  SetMakeGeneralHistograms(kTRUE);

	DefineOutput(1, TList::Class());
	//if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart && fFillInclusiveTree) 
	DefineOutput(2, TTree::Class());
//	if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart && fFillRecoilTree)
	DefineOutput(3, TTree::Class());
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskJetCoreEmcal::AliAnalysisTaskJetCoreEmcal(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fHistManager(name),
	fJetShapeType(AliAnalysisTaskJetCoreEmcal::kData),
	fCentMin(0.),
	fCentMax(100.),
	fTTLowRef(8),
	fTTUpRef(9.),
	fTTLowSig(20.),
	fTTUpSig(50.),
	fNRPBins(50),
	fFrac(0.8),
	fJetHadronDeltaPhi(0.6),
	fMinFractionSharedPt(0.5),
	fMinEmbJetPt(15.),
	fJetContName(""),
	fJetContTrueName(""),
	fJetContPartName(""),
	fFillTrackHistograms(kTRUE),
	fFillJetHistograms(kTRUE),
	fFillRecoilTHnSparse(kTRUE),
	fFillInclusiveTree(kFALSE),
	fFillRecoilTree(kFALSE),
	fMoreTreeVars(kFALSE),
	fPtHardBin(0.),
	fRejectionFactorInclusiveJets(1),
	fRandom(0),
	fHistEvtSelection(0x0), 
	fHJetSpec(0x0),
	fh1TrigRef(0x0),
	fh1TrigSig(0x0),
	fh2Ntriggers(0x0),
	fhRhoCentSig(0x0),
	fhRhoCentRef(0x0),
	fhDphiPtSigPi(0x0),
	fhDphiPtSig(0x0),
	fhDphiPtRefPi(0x0),
	fhDphiPtRef(0x0),
	fhPtDetPart(0x0),
	fhPtHybrDet(0x0),
	fhPtHybrPart(0x0),
	fhPtHybrPartCor(0x0),
	fhPhiHybrPartCor(0x0),
	fhPtDet(0x0),
	fhPtDetMatchedToPart(0x0),
	fhPtPartMatched(0x0),
	fhPtPartMatchedCent(0x0),
	fhPtPartMatchedWrongCent(0x0),
	fhResidual(0x0),
	fhPtResidual(0x0),
	fhPhiResidual(0x0),
	fhPhiPhiResidual(0x0),
	fhPtDetPartRecoil(0x0),
	fhPtHybrDetRecoil(0x0),
	fhPtHybrPartRecoil(0x0),
	fhPtHybrPartCorRecoil(0x0),
	fhPtDetRecoil(0x0),
	fhPtDetMatchedToPartRecoil(0x0),
	fhResidualRecoil(0x0),
	fhPtResidualRecoil(0x0),
	fhDphiResidualRecoil(0x0),
	fhDphiphiResidualRecoil(0x0),
	fhTTPtDetMatchedToPart(0x0),
	fhTTPhiDetMatchedToPart(0x0),
	fhDPhiHybrPartCorRecoil(0x0),
	fhSelectedTrigger(0x0),
	fhFractionSharedPtInclusive(0x0),
	fhFractionSharedPtRecoil(0x0),
	fTreeEmbInclusive(0x0),
	fTreeEmbRecoil(0x0)
{
  SetMakeGeneralHistograms(kTRUE);

	DefineOutput(1, TList::Class());
	DefineOutput(2, TTree::Class());
	DefineOutput(3, TTree::Class());
}

/**
 * Destructor
 */
AliAnalysisTaskJetCoreEmcal::~AliAnalysisTaskJetCoreEmcal()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskJetCoreEmcal::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  if(fFillTrackHistograms) AllocateTrackHistograms();
  if(fFillJetHistograms) AllocateJetHistograms();
  //AllocateClusterHistograms();
  //AllocateCellHistograms();
  AllocateJetCoreHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
	if((fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) && fFillInclusiveTree) PostData(2, fTreeEmbInclusive); // Post data for ALL output slots > 0 here.
  if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart && fFillRecoilTree)    PostData(3, fTreeEmbRecoil); // Post data for ALL output slots > 0 here.
}

/*
 * This function allocates the histograms for basic EMCal cluster QA.
 * A set of histograms (energy, eta, phi, number of cluster) is allocated
 * per each cluster container and per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateClusterHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliClusterContainer* clusCont = 0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The cluster containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{exotic} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{non-lin.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{E}_{cluster}^{had.corr.} (GeV);counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{custer};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNClusters_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of clusters;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNClusters";
  histtitle = TString::Format("%s;Sum of n clusters;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 3000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/*
 * This function allocates the histograms for basic EMCal QA.
 * One 2D histogram with the cell energy spectra and the number of cells
 * per event is allocated per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateCellHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname(fCaloCellsName);

  fHistManager.CreateHistoGroup(groupname);
  for (Int_t cent = 0; cent < fNcentBins; cent++) {
    histname = TString::Format("%s/histCellEnergy_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;#it{E}_{cell} (GeV);counts", histname.Data());
    fHistManager.CreateTH1(histname, histtitle, 300, 0, 150);

    histname = TString::Format("%s/histNCells_%d", groupname.Data(), cent);
    histtitle = TString::Format("%s;number of cells;events", histname.Data());
    if (fForceBeamType != kpp) {
      fHistManager.CreateTH1(histname, histtitle, 500, 0, 6000);
    }
    else {
      fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
    }
  }
}

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateTrackHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();

    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The track containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt / 2);

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      if (TClass(partCont->GetClassName()).InheritsFrom("AliVTrack")) {
        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#eta}_{track}^{vertex} - #it{#eta}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 50, -0.5, 0.5);

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{#phi}_{track}^{vertex} - #it{#phi}_{track}^{EMCal};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, 200, -2, 2);

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,track}^{vertex} (GeV/#it{c});#it{p}_{T,track}^{vertex} - #it{p}_{T,track}^{EMCal} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, -fMaxBinPt/2, fMaxBinPt/2);

        histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{P}_{track} (GeV/#it{c});#it{E}_{cluster} / #it{P}_{track} #it{c};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins / 2, fMinBinPt, fMaxBinPt, fNbins / 2, 0, 4);
      }

      histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of tracks;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }
  }

  histname = "fHistSumNTracks";
  histtitle = TString::Format("%s;Sum of n tracks;events", histname.Data());
  if (fForceBeamType != kpp) {
    fHistManager.CreateTH1(histname, histtitle, 500, 0, 5000);
  }
  else {
    fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
  }
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskJetCoreEmcal::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    // Protect against creating the histograms twice
    if (fHistManager.FindObject(groupname)) {
      AliWarning(TString::Format("%s: Found groupname %s in hist manager. The jet containers will be filled into the same histograms.", GetName(), groupname.Data()));
      continue;
    }
    fHistManager.CreateHistoGroup(groupname);
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, fMinBinPt, fMaxBinPt);

      histname = TString::Format("%s/histJetPtLow_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, 50, 0., 0.5);

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{A}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, 3);

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 2, 0, TMath::TwoPi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{jet};counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins / 6, -1, 1);

      histname = TString::Format("%s/histNJets_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of jets;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 500, 0, 500);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 100, 0, 100);
      }

      if (!jetCont->GetRhoName().IsNull()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);

        histname = TString::Format("%s/histJetCorrPtLeadingTrackPt_%d", groupname.Data(), cent);
        histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});#it{p}_{T,leading} (GeV/#it{c});counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2, 1000,0.,100.);

      }
    }
  }
}

void AliAnalysisTaskJetCoreEmcal::AllocateJetCoreHistograms() 
{

	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);

	// set seed
	fRandom = new TRandom3(0);

	fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
	fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
	fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
	fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
	fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
	fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
	fHistEvtSelection->GetXaxis()->SetBinLabel(6,"multiplicity (rejected)");
	
	fh1TrigRef=new TH1D("Trig Ref","",10,0.,10);
	fh1TrigSig=new TH1D("Trig Sig","",10,0.,10);  
	fh2Ntriggers=new TH2F("# of triggers","",100,0.,100.,50,0.,50.);

  fhRhoCentSig=new TH2F("hRhoCentSig","rho vs centrality signal",1500,0,300,100,0,100);
  fhRhoCentRef=new TH2F("hRhoCentRef","rho vs centrality reference",1500,0,300,100,0,100);

	fOutput->Add(fHistEvtSelection);

	fOutput->Add(fh1TrigRef);
	fOutput->Add(fh1TrigSig); 
	fOutput->Add(fh2Ntriggers);
  fOutput->Add(fhRhoCentSig);
  fOutput->Add(fhRhoCentRef);

	if(fFillRecoilTHnSparse) {
		const Int_t dimSpec = 6;
		const Int_t nBinsSpec[dimSpec]     = {10,10, 280, 50,200, fNRPBins};
		const Double_t lowBinSpec[dimSpec] = {0,0,-80, 0,0, 0};
		const Double_t hiBinSpec[dimSpec]  = {100,1, 200, 50,2*TMath::Pi(),  static_cast<Double_t>(fNRPBins)};
		fHJetSpec = new THnSparseF("fHJetSpec","Recoil jet spectrum",dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);

	//change binning in jet area
    Double_t *xArea = new Double_t[11];
    xArea[0]=0.;
    xArea[1]=0.06;
    xArea[2]=0.07;
    xArea[3]=0.08;
    xArea[4]=0.36;
    xArea[5]=0.4;
    xArea[6]=0.44;
    xArea[7]=0.55;
    xArea[8]=0.6;
    xArea[9]=0.65;
    xArea[10]=1.;
 //   Double_t xArea[11] = {0., 0.06, 0.07, 0.08, 0.36, 0.4, 0.44, 0.55, 0.6, 0.65, 1.};
    fHJetSpec->SetBinEdges(1,xArea);
    delete [] xArea;

	}

	fOutput->Add(fHJetSpec);  

	// azimuthal correlation

	fhDphiPtSigPi = new TH2F("hDphiPtSPi","recoil #Delta #phi vs jet pT signal",200,0,2*TMath::Pi(),25000,-50,200);  
	fhDphiPtSigPi->GetXaxis()->SetTitle("#Delta #phi"); 
	fhDphiPtSigPi->GetYaxis()->SetTitle("p^{reco,ch}_{T,jet} (GeV/c)"); 
	fhDphiPtRefPi = new TH2F("hDphiPtRPi","recoil #Delta #phi vs jet pT reference",200,0,2*TMath::Pi(),25000,-50,200);  
	fhDphiPtRefPi->GetXaxis()->SetTitle("#Delta #phi"); 
	fhDphiPtRefPi->GetYaxis()->SetTitle("p^{reco,ch}_{T,jet} (GeV/c)"); 

	fhDphiPtSig = new TH2F("hDphiPtS","recoil #Delta #phi vs jet pT signal",100,-2,5,250,-50,200);  
	fhDphiPtSig->GetXaxis()->SetTitle("#Delta #phi"); 
	fhDphiPtSig->GetYaxis()->SetTitle("p^{reco,ch}_{T,jet} (GeV/c)"); 
	fhDphiPtRef = new TH2F("hDphiPtR","recoil #Delta #phi vs jet pT reference",100,-2,5,250,-50,200);  
	fhDphiPtRef->GetXaxis()->SetTitle("#Delta #phi"); 
	fhDphiPtRef->GetYaxis()->SetTitle("p^{reco,ch}_{T,jet} (GeV/c)"); 


	fOutput->Add(fhDphiPtRef);  
	fOutput->Add(fhDphiPtSig);  
	fOutput->Add(fhDphiPtRefPi);  
	fOutput->Add(fhDphiPtSigPi);  


	fhPtHybrDet= new TH2F("hPtHybrDet","pT response Pb-Pb+PYTHIA vs PYTHIA",200,0,200,200,0,200);
	fhPtHybrDet->GetXaxis()->SetTitle("p^{Pb-Pb+PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtHybrDet->GetYaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtDetPart= new TH2F("hPtDetPart","pT response PYTHIA vs part",200,0,200,200,0,200);
	fhPtDetPart->GetXaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtDetPart->GetYaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 
	fhPtHybrPart = new TH2F("hPtHybrPart",Form("pT response Pb-Pb+PYTHIA vs part, min shared pT > %.0f",fMinFractionSharedPt*100),200,0,200,200,0,200);
	fhPtHybrPart->GetXaxis()->SetTitle("p^{Pb-Pb+PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtHybrPart->GetYaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 
	fhPtHybrPartCor = new TH2F("hPtHybrPartCor",Form("pT response Pb-Pb+PYTHIA corrected vs part, min shared pT > %.0f",fMinFractionSharedPt*100),200,0,200,200,0,200);
	fhPtHybrPartCor->GetXaxis()->SetTitle("p^{Pb-Pb+PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtHybrPartCor->GetYaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 

	fhPhiHybrPartCor = new TH2F("hPhiHybrPartCor",Form("phi response Pb-Pb+PYTHIA corrected vs part, min shared pT > %.0f",fMinFractionSharedPt*100),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-0.5*TMath::Pi(),1.5*TMath::Pi());
	fhPhiHybrPartCor->GetXaxis()->SetTitle("#phi^{Pb-Pb+PYTHIA,ch}"); 
	fhPhiHybrPartCor->GetYaxis()->SetTitle("#phi^{true}"); 
	fhResidual = new TH1F("hResidual","residual",50,-1,1);
	fhResidual->GetXaxis()->SetTitle("p^{reco}_{T} - p^{true}_{T} / p^{true}_{T}"); 
	fhPtResidual= new TH2F("hPtResidual","pT vs residual",200,0,200,50,-1,1);
	fhPtResidual->GetXaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 
	fhPtResidual->GetYaxis()->SetTitle("p^{reco}_{T} - p^{true}_{T} / p^{true}_{T}"); 
	fhPhiResidual = new TH1F("hPhiResidual","residual phi",50,-1,1);
	fhPhiResidual->GetXaxis()->SetTitle("#phi^{reco} - #phi^{true} / #phi^{true}"); 
	fhPhiPhiResidual= new TH2F("hPhiPhiResidual","pT vs residual",600,-3,3,200,-0.5*TMath::Pi(),1.5*TMath::Pi());
	fhPhiPhiResidual->GetXaxis()->SetTitle("#phi^{true}"); 
	fhPhiPhiResidual->GetYaxis()->SetTitle("#phi^{reco} - #phi^{true} / #phi^{true}"); 

	fhPtDet= new TH1F("hPtDet","pT detector level",200,0,200);
	fhPtDet->GetXaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtDetMatchedToPart = new TH1F("hPtDetMatchedToPart","pT detector level matched to particle level jet",200,0,200);
	fhPtDetMatchedToPart->GetXaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtPartMatched= new TH1F("hPtPartMatched","pT particle level matched",200,0,200);
	fhPtPartMatched->GetXaxis()->SetTitle("p^{part}_{T} (GeV/c)"); 
	fhPtPartMatchedCent= new TH2F("hPtPartMatchedCent","pT particle level matched",200,0,200,100,0,100);
	fhPtPartMatchedCent->GetXaxis()->SetTitle("p^{part}_{T} (GeV/c)"); 
	fhPtPartMatchedCent->GetYaxis()->SetTitle("centrality (%)"); 
	fhPtPartMatchedWrongCent= new TH2F("hPtPartMatchedWrongCent","pT particle level matched incorrectly",200,0,200,100,0,100);
	fhPtPartMatchedWrongCent->GetXaxis()->SetTitle("p^{part}_{T} (GeV/c)"); 
	fhPtPartMatchedWrongCent->GetYaxis()->SetTitle("centrality (%)"); 

	fhPtHybrDetRecoil= new TH2F("hPtHybrDetRecoil","pT response Pb-Pb+PYTHIA vs PYTHIA",200,0,200,200,0,200);
	fhPtHybrDetRecoil->GetXaxis()->SetTitle("p^{Pb-Pb+PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtHybrDetRecoil->GetYaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtDetPartRecoil= new TH2F("hPtDetPartRecoil","pT response PYTHIA vs part",200,0,200,200,0,200);
	fhPtDetPartRecoil->GetXaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtDetPartRecoil->GetYaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 
	fhPtHybrPartRecoil = new TH2F("hPtHybrPartRecoil",Form("pT response Pb-Pb+PYTHIA vs part, min shared pT > %.0f",fMinFractionSharedPt*100),200,0,200,200,0,200);
	fhPtHybrPartRecoil->GetXaxis()->SetTitle("p^{Pb-Pb+PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtHybrPartRecoil->GetYaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 
	fhPtHybrPartCorRecoil = new TH2F("hPtHybrPartCorRecoil",Form("pT response Pb-Pb+PYTHIA corrected vs part, min shared pT > %.0f",fMinFractionSharedPt*100),200,0,200,200,0,200);
	fhPtHybrPartCorRecoil->GetXaxis()->SetTitle("p^{Pb-Pb+PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtHybrPartCorRecoil->GetYaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 

	fhDPhiHybrPartCorRecoil = new TH2F("hDPhiHybrPartCorRecoil",Form("#Delta#phi response Pb-Pb+PYTHIA corrected vs part, min shared pT > %.0f",fMinFractionSharedPt*100),200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-0.5*TMath::Pi(),1.5*TMath::Pi());
	fhDPhiHybrPartCorRecoil->GetXaxis()->SetTitle("#Delta#phi^{Pb-Pb+PYTHIA,ch}"); 
	fhDPhiHybrPartCorRecoil->GetYaxis()->SetTitle("#Delta#phi^{true}"); 
	fhResidualRecoil = new TH1F("hResidualRecoil","residual",50,-1,1);
	fhResidualRecoil->GetXaxis()->SetTitle("p^{reco}_{T} - p^{true}_{T} / p^{true}_{T}"); 
	fhPtResidualRecoil= new TH2F("hPtResidualRecoil","pT vs residual",200,0,200,50,-1,1);
	fhPtResidualRecoil->GetXaxis()->SetTitle("p^{true}_{T} (GeV/c)"); 
	fhPtResidualRecoil->GetYaxis()->SetTitle("p^{reco}_{T} - p^{true}_{T} / p^{true}_{T}"); 

	fhDphiResidualRecoil = new TH1F("hDphiResidualRecoil","residual",50,-1,1);
	fhDphiResidualRecoil->GetXaxis()->SetTitle("#Delta#phi^{reco} - #Delta#phi^{true} / #Delta#phi^{true}"); 
	fhDphiphiResidualRecoil= new TH2F("hDphiphiResidualRecoil","#Delta#phi vs residual in phi",400,-3,3,200,-0.5*TMath::Pi(),1.5*TMath::Pi());
	fhDphiphiResidualRecoil->GetXaxis()->SetTitle("#Delta#phi^{true}"); 
	fhDphiphiResidualRecoil->GetYaxis()->SetTitle("#Delta#phi^{reco} - #Delta#phi^{true} / #Delta#phi^{true}"); 
	fhPtDetRecoil= new TH1F("hPtDetRecoil","pT detector level",200,0,200);
	fhPtDetRecoil->GetXaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 
	fhPtDetMatchedToPartRecoil = new TH1F("hPtDetMatchedToPartRecoil","pT detector level matched to particle level jet",200,0,200);
	fhPtDetMatchedToPartRecoil->GetXaxis()->SetTitle("p^{reco,PYTHIA,ch}_{T} (GeV/c)"); 

	fhTTPtDetMatchedToPart = new TH2F("hTTPtDetMatchedToPart","trigger track pT response reco vs partice",140,0,70,140,0,70);
	fhTTPtDetMatchedToPart->GetXaxis()->SetTitle("p^{TT,reco}_{T} (GeV/c)"); 
	fhTTPtDetMatchedToPart->GetYaxis()->SetTitle("p^{TT,part}_{T} (GeV/c)"); 
	fhTTPhiDetMatchedToPart = new TH2F("hTTPhiDetMatchedToPart","trigger track #phi response reco vs partice",200,-0.5*TMath::Pi(),1.5*TMath::Pi(),200,-0.5*TMath::Pi(),1.5*TMath::Pi());
	fhTTPhiDetMatchedToPart->GetXaxis()->SetTitle("p^{TT,reco}_{T} (GeV/c)"); 
	fhTTPhiDetMatchedToPart->GetYaxis()->SetTitle("p^{TT,part}_{T} (GeV/c)"); 

	fOutput->Add(fhPtHybrDet);
	fOutput->Add(fhPtDetPart);
	fOutput->Add(fhPtHybrPart);
	fOutput->Add(fhPtHybrPartCor);
	fOutput->Add(fhPhiHybrPartCor);
	fOutput->Add(fhPtDet);
	fOutput->Add(fhPtDetMatchedToPart);
	fOutput->Add(fhPtPartMatched);
	fOutput->Add(fhPtPartMatchedCent);
	fOutput->Add(fhPtPartMatchedWrongCent);
	fOutput->Add(fhResidual);
	fOutput->Add(fhPtResidual);
	fOutput->Add(fhPhiResidual);
	fOutput->Add(fhPhiPhiResidual);

	fOutput->Add(fhPtHybrDetRecoil);
	fOutput->Add(fhPtDetPartRecoil);
	fOutput->Add(fhPtHybrPartRecoil);
	fOutput->Add(fhPtHybrPartCorRecoil);
	fOutput->Add(fhDPhiHybrPartCorRecoil);
	fOutput->Add(fhPtDetRecoil);
	fOutput->Add(fhPtDetMatchedToPartRecoil);
	fOutput->Add(fhResidualRecoil);
	fOutput->Add(fhPtResidualRecoil);
	fOutput->Add(fhDphiResidualRecoil);
	fOutput->Add(fhDphiphiResidualRecoil);

	fOutput->Add(fhTTPtDetMatchedToPart);
	fOutput->Add(fhTTPhiDetMatchedToPart);

  fhSelectedTrigger= new TH2F("hSelectedTrigger","ID of selected trigger",2,0,2,200,0,100);
  fhSelectedTrigger->GetXaxis()->SetBinLabel(1,"Pb-Pb trigger");
  fhSelectedTrigger->GetXaxis()->SetBinLabel(2,"pp trigger");
  fhSelectedTrigger->GetYaxis()->SetTitle("p^{TT}_{T} (GeV/c)"); 
	fOutput->Add(fhSelectedTrigger);


  fhFractionSharedPtInclusive = new TH2F("hFractionSharedPtInclusive","fraction of shared pT",200,-50,150,50,0,1); 
  fhFractionSharedPtInclusive ->GetXaxis()->SetTitle("p_{T}^{Pb-Pb}"); 
  fhFractionSharedPtInclusive ->GetYaxis()->SetTitle("f"); 
  fhFractionSharedPtRecoil = new TH2F("hFractionSharedPtRecoil","fraction of shared pT",200,-50,150,50,0,1); 
  fhFractionSharedPtRecoil ->GetXaxis()->SetTitle("p_{T}^{Pb-Pb}"); 
  fhFractionSharedPtRecoil ->GetYaxis()->SetTitle("f"); 
	fOutput->Add(fhFractionSharedPtInclusive);
	fOutput->Add(fhFractionSharedPtRecoil);

  TString varNamesInclusive[9]={"centrality","ptRawRec","areaRec","ptCorrRec","phiRec","ptPart","phiPart","ptLeadingTrackRec","ptLeadingTrackPart"};
  TString varNamesInclusiveMoreVars[13]={"centrality","ptRawRec","areaRec","ptCorrRec","phiRec","ptPart","phiPart","ptLeadingTrackRec","ptLeadingTrackPart","ptDet","phiDet","matchedJetDistanceRec","matchedJetDistancePart"};
  TString varNamesRecoil[8]={"centrality","ptTT","ptRawRec","areaRec","ptCorrRec","DPhiRec","ptPart","DPhiPart"};
  TString varNamesRecoilMoreVars[12]={"centrality","ptTT","ptRawRec","areaRec","ptCorrRec","DPhiRec","ptPart","DPhiPart","ptDet","DPhiDet","matchedJetDistanceRec","matchedJetDistancePart"};
	if((fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) && fFillInclusiveTree) {
		const char* nameEmbInclusive = GetOutputSlot(2)->GetContainer()->GetName();
		fTreeEmbInclusive = new TTree(nameEmbInclusive, nameEmbInclusive);
    if(fMoreTreeVars) {
      for(Int_t ivar=0; ivar < 13; ivar++){
        fTreeEmbInclusive->Branch(varNamesInclusiveMoreVars[ivar].Data(), &fTreeVarsInclusiveMoreVars[ivar], Form("%s/F", varNamesInclusiveMoreVars[ivar].Data()));
      }
    }
    else {
      for(Int_t ivar=0; ivar < 9; ivar++){
        fTreeEmbInclusive->Branch(varNamesInclusive[ivar].Data(), &fTreeVarsInclusive[ivar], Form("%s/F", varNamesInclusive[ivar].Data()));
      }
    }
	}

	if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart && fFillRecoilTree) {
		const char* nameEmbRecoil= GetOutputSlot(3)->GetContainer()->GetName();
		fTreeEmbRecoil = new TTree(nameEmbRecoil, nameEmbRecoil);
    if(fMoreTreeVars) {
      for(Int_t ivar=0; ivar < 12; ivar++){
        fTreeEmbRecoil->Branch(varNamesRecoilMoreVars[ivar].Data(), &fTreeVarsRecoilMoreVars[ivar], Form("%s/F", varNamesRecoilMoreVars[ivar].Data()));
      }
    }
    else {
      for(Int_t ivar=0; ivar < 8; ivar++){
        fTreeEmbRecoil->Branch(varNamesRecoil[ivar].Data(), &fTreeVarsRecoil[ivar], Form("%s/F", varNamesRecoil[ivar].Data()));
      }
    }
	}

	// =========== Switch on Sumw2 for all histos ===========
	for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
		TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
		if (h1){
			h1->Sumw2();
			continue;
		}
		THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
		if (hn){
			hn->Sumw2();
		}	  
	}

	// add QA plots from fEventCuts
	fEventCuts.AddQAplotsToList(fOutput);

	TH1::AddDirectory(oldStatus);

	PostData(1, fOutput);
	if((fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) && fFillInclusiveTree) PostData(2, fTreeEmbInclusive);

  if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart && fFillRecoilTree)    PostData(3, fTreeEmbRecoil);
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskJetCoreEmcal::FillHistograms()
{

	fHistEvtSelection->Fill(1); // number of events before event selection
	AliVEvent *ev = InputEvent();
	if (!fEventCuts.AcceptEvent(ev)) {
		fHistEvtSelection->Fill(2);
		return kTRUE;
	}

	// centrality selection 
	if(fDebug) Printf("centrality: %f\n", fCent);
	if (fCent>fCentMax || fCent<fCentMin) {
		fHistEvtSelection->Fill(4);
		return kTRUE;
	}

  if(fFillJetHistograms) DoJetLoop();
  if(fFillTrackHistograms) DoTrackLoop();
  //DoClusterLoop();
  //DoCellLoop();
	DoJetCoreLoop();
	if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart || 
      fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart || 
      fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) DoMatchingLoop();


  return kTRUE;
}

void AliAnalysisTaskJetCoreEmcal::DoJetCoreLoop()
{
	// Do jet core analysis and fill histograms.

	AliJetContainer *jetCont = GetJetContainer(fJetContName);
	AliJetContainer *jetContPart = 0x0; 
	AliJetContainer *jetContTrue = 0x0; 
	if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart){
 	 jetContTrue = GetJetContainer(fJetContTrueName);
 	 jetContPart = GetJetContainer(fJetContPartName);
	}

	if(!jetCont ||
			(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart && !jetContPart)) {
		AliError(Form("jet container not found - check name %s",fJetContName.Data()));
		TIter next(&fJetCollArray);
		while ((jetCont = static_cast<AliJetContainer*>(next())))
			AliError(Form("%s",jetCont->GetName()));
		AliFatal("Exit...");
		return;
	}

	fHistEvtSelection->Fill(0); 

	// Background
  // If rho exists get it, otherwise it is set to 0 and ptJet_corr = ptJet_raw
	Double_t rho = 0;
	if (jetCont->GetRhoParameter()) rho = jetCont->GetRhoVal(); 
	if(fDebug) Printf("rho = %f",rho);

	// get MC particle container in case running embedding, to match
	// reconstructed and MC-level trigger tracks
	AliParticleContainer *partCont = 0x0; 
	if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart) partCont = GetParticleContainer(2);
//	AliParticleContainer *partCont = 0x0;
//	TIter next(&fParticleCollArray);
//	while ((partCont = static_cast<AliParticleContainer*>(next()))) {
//		TString groupname = partCont->GetName();
//		Printf("particle name = %s",groupname.Data());
//	}

	// Choose trigger track
	Int_t nT=0;
	TList ParticleList;
	Double_t minT=0;
	Double_t maxT=0;
	Int_t number=0;
	Double_t dice=fRandom->Uniform(0,1);
	Bool_t isSignal = kFALSE;
	if(dice>fFrac){ 
		minT=fTTLowRef;
		maxT=fTTUpRef;
	}
	if(dice<=fFrac){
		isSignal = kTRUE;
		minT=fTTLowSig;
		maxT=fTTUpSig;
	} 
	nT=SelectTrigger(&ParticleList,minT,maxT,number);
	if(fDebug) Printf("%s class ---> n triggers between %f and %f = %i, index of trigger chosen = %i",dice>fFrac?"ref.":"sig.",minT,maxT,number,nT);
	if(nT<0) return;

	if(isSignal) {
    fh1TrigSig->Fill(number);
    fhRhoCentSig->Fill(rho,fCent);
  }
  else         {
    fh1TrigRef->Fill(number);
    fhRhoCentRef->Fill(rho,fCent);
  }


	// particle loop - 
	for(Int_t tt=0;tt<ParticleList.GetEntries();tt++){
		// histogram 
		//if(fHardest==0||fHardest==1){if(tt!=nT) continue;}
		if(tt!=nT) continue;
		AliVParticle *partback = (AliVParticle*)ParticleList.At(tt);     
		if(!partback) continue;
		if(fDebug) Printf("trigger particle pt = %f \teta = %f \t phi = %f",partback->Pt(),partback->Eta(),partback->Phi());
		//     if(partback->Pt()<8) continue;

    fh2Ntriggers->Fill(fCent,partback->Pt());
    Double_t phiBinT = RelativePhi(partback->Phi(),fEPV0);

		Double_t etabig=0;
		Double_t ptbig=0;
		Double_t areabig=0;
		Double_t phibig=0.;
		//   Double_t areasmall=0;

		TString histname;
		TString groupname;
		groupname = jetCont->GetName();
		UInt_t count = 0;
		for(auto jetbig : jetCont->accepted()) {
			if (!jetbig) continue;
			count++;
			ptbig   = jetbig->Pt();
			etabig  = jetbig->Eta();
			phibig  = jetbig->Phi();
			if(ptbig==0) continue; 
			Double_t phiBin = RelativePhi(phibig,fEPV0); //relative phi between jet and ev. plane
			areabig = jetbig->Area();
			Double_t ptcorr=ptbig-rho*areabig;
			Double_t dphi=RelativePhi(partback->Phi(),phibig); 
			if(fDebug) Printf("jet properties...\n\teta = %f \t phi = %f \t pt = %f \t relativephi = %f\t area = %f\t rho = %f",etabig,phibig,ptbig,dphi,areabig,rho);

			// do azimuthal correlation analysis
			// dPhi between -0.5 < dPhi < 1.5
			Double_t dPhiShift=phibig-partback->Phi();
			if(dPhiShift>2*TMath::Pi()) dPhiShift -= 2*TMath::Pi();
			if(dPhiShift<-2*TMath::Pi()) dPhiShift += 2*TMath::Pi();
			if(dPhiShift<-0.5*TMath::Pi()) dPhiShift += 2*TMath::Pi();
			if(dPhiShift>1.5*TMath::Pi()) dPhiShift -= 2*TMath::Pi();

			// dPhi between 0 < dPhi < 2pi
			Double_t dPhiShiftPi=phibig-partback->Phi();
			if(dPhiShiftPi>2*TMath::Pi()) dPhiShiftPi -= 2*TMath::Pi();
			if(dPhiShiftPi<-2*TMath::Pi()) dPhiShiftPi += 2*TMath::Pi();
			if(dPhiShiftPi<0)              dPhiShiftPi += 2*TMath::Pi();


			if(isSignal) {
        fhDphiPtSigPi->Fill(dPhiShiftPi,ptcorr);
        fhDphiPtSig->Fill(dPhiShift,ptcorr);
      }
			else         {
        fhDphiPtRefPi->Fill(dPhiShiftPi,ptcorr);
        fhDphiPtRef->Fill(dPhiShift,ptcorr);
      }

			// selection on relative phi
			if(fJetHadronDeltaPhi>0. &&
					TMath::Abs(dphi)<TMath::Pi()-fJetHadronDeltaPhi) continue;

			if(fFillRecoilTHnSparse) {
				Float_t phitt=partback->Phi();
				if(phitt<0)phitt+=TMath::Pi()*2.; 
				Int_t phiBintt = GetPhiBin(phitt-fEPV0);

				Double_t fillspec[] = {fCent,areabig,ptcorr,partback->Pt(),dPhiShiftPi, static_cast<Double_t>(phiBintt)};
				fHJetSpec->Fill(fillspec);
			}

			if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart) {
				//
				// embedding for recoil jets
				// get MC info
				//
				Double_t ptTTMC = 0;
				Double_t phiTTMC = 0;
				Int_t TTmatched = 0;
        Double_t distanceClosestJet1=-1, distanceClosestJet2=-1;
				for(auto partMC : partCont->accepted()) {
					Int_t labtr = partback->GetLabel();
					Int_t labpa = partMC->GetLabel();
					if(labtr==labpa) {
						ptTTMC = partMC->Pt();
						phiTTMC = partMC->Phi();
						TTmatched++;
						break;
					}
				}
				if(TTmatched!=1) continue;
				Double_t ptTTreco = partback->Pt();
				Double_t phiTTreco = partback->Phi();
				if(fDebug) Printf("found corresponding truth-level particle, pt reco = %f pt MC = %f, phi reco = %f phi MC = %f",ptTTreco,ptTTMC,phiTTreco,phiTTMC);
				fhTTPtDetMatchedToPart->Fill(ptTTreco,ptTTMC);
				fhTTPhiDetMatchedToPart->Fill(phiTTreco,phiTTMC);

				auto jet2 = jetbig->ClosestJet();
				if(!jet2) {
					//Printf("jet 2 cant be found");
					continue;}
        distanceClosestJet1 = jetbig->ClosestJetDistance();
				Double_t ptJet2 = jet2->Pt();
				Double_t phiJet2 = jet2->Phi();
				fhPtDetRecoil->Fill(ptJet2);
				auto jet3 = jet2->ClosestJet();
				if(!jet3) {
					//Printf("jet3 can't be found");
					continue;
				}
        distanceClosestJet2 = jet2->ClosestJetDistance();
				fhPtDetMatchedToPartRecoil->Fill(ptJet2);
				Double_t ptJet3 = jet3->Pt();
				Double_t phiJet3 = jet3->Phi();


				Double_t dPhiPart=phiJet3-phiTTMC;
        Double_t dPhiPartShiftPi=phiJet3-phiTTMC;
				if(dPhiPart>2*TMath::Pi()) dPhiPart -= 2*TMath::Pi();
				if(dPhiPart<-2*TMath::Pi()) dPhiPart += 2*TMath::Pi();
				if(dPhiPart<-0.5*TMath::Pi()) dPhiPart += 2*TMath::Pi();
				if(dPhiPart>1.5*TMath::Pi()) dPhiPart -= 2*TMath::Pi();

        // dPhi between 0 < dPhi < 2pi
        if(dPhiPartShiftPi>2*TMath::Pi()) dPhiPartShiftPi -= 2*TMath::Pi();
        if(dPhiPartShiftPi<-2*TMath::Pi()) dPhiPartShiftPi += 2*TMath::Pi();
        if(dPhiPartShiftPi<0)              dPhiPartShiftPi += 2*TMath::Pi();

        Double_t dPhiDetShiftPi=phiJet2-partback->Phi();
        if(dPhiDetShiftPi>2*TMath::Pi())  dPhiDetShiftPi -= 2*TMath::Pi();
        if(dPhiDetShiftPi<-2*TMath::Pi()) dPhiDetShiftPi += 2*TMath::Pi();
        if(dPhiDetShiftPi<0)              dPhiDetShiftPi += 2*TMath::Pi();

				if(fDebug) Printf("--- recoil - jet pt = jet hybrid pt = %f\t jet matched det pt = %f\t jet matched particle level pt = %f\t\n\tjet reco phi = %f\t jet particle phi = %f",ptbig,ptJet2,ptJet3,phibig,phiJet3);

				fhPtDetPartRecoil->Fill(ptJet2,ptJet3);
				Double_t fraction = jetCont->GetFractionSharedPt(jetbig);
        fhFractionSharedPtRecoil->Fill(ptcorr,fraction);
				if(fraction < fMinFractionSharedPt) continue;

				Double_t residual = (ptcorr - ptJet3) / ptJet3;
				Double_t residualDphi = (dPhiShift - dPhiPart) / dPhiPart;

				fhPtHybrDetRecoil->Fill(ptbig,ptJet2);
				fhPtHybrPartRecoil->Fill(ptbig,ptJet3);
				fhPtHybrPartCorRecoil->Fill(ptcorr,ptJet3);
				fhDPhiHybrPartCorRecoil->Fill(dPhiShift,dPhiPart);
				fhResidualRecoil->Fill(residual);
				fhPtResidualRecoil->Fill(ptJet3,residual);
				fhDphiResidualRecoil->Fill(residualDphi);
				fhDphiphiResidualRecoil->Fill(dPhiPart,residualDphi);

				if(ptcorr<fMinEmbJetPt) continue;

				if(fFillRecoilTree) {
          //TString varNamesRecoilMoreVars[12]={"centrality","ptTT","ptRawRec","areaRec","ptCorrRec","DPhiRec","ptPart","DPhiPart","ptDet","DPhiDet","matchedJetDistanceRec","matchedJetDistancePart"};
          if(fMoreTreeVars) {
            fTreeVarsRecoilMoreVars[0] = fCent;
            fTreeVarsRecoilMoreVars[1] = partback->Pt();
            fTreeVarsRecoilMoreVars[2] = ptbig;
            fTreeVarsRecoilMoreVars[3] = areabig;
            fTreeVarsRecoilMoreVars[4] = ptcorr;
            fTreeVarsRecoilMoreVars[5] = dPhiShiftPi;
            fTreeVarsRecoilMoreVars[6] = ptJet3;
            fTreeVarsRecoilMoreVars[7] = dPhiPartShiftPi;
            fTreeVarsRecoilMoreVars[8] = ptJet2;
            fTreeVarsRecoilMoreVars[9] = dPhiDetShiftPi;
            fTreeVarsRecoilMoreVars[10] = distanceClosestJet1;
            fTreeVarsRecoilMoreVars[11] = distanceClosestJet2;
          }
          else {
            fTreeVarsRecoil[0] = fCent;
            fTreeVarsRecoil[1] = partback->Pt();
            fTreeVarsRecoil[2] = ptbig;
            fTreeVarsRecoil[3] = areabig;
            fTreeVarsRecoil[4] = ptcorr;
            fTreeVarsRecoil[5] = dPhiShiftPi;
            fTreeVarsRecoil[6] = ptJet3;
            fTreeVarsRecoil[7] = dPhiPartShiftPi;
          }
					fTreeEmbRecoil->Fill();
				}
			}
		}
	}
}

void AliAnalysisTaskJetCoreEmcal::DoMatchingLoop() {

	AliParticleContainer *partCont0 = GetParticleContainer(0);
	AliParticleContainer *partCont1 = GetParticleContainer(1);
	AliJetContainer *jetCont = GetJetContainer(fJetContName);
	AliJetContainer *jetContPart = GetJetContainer(fJetContPartName);
  AliJetContainer *jetContTrue = 0x0;
  if(fJetShapeType==AliAnalysisTaskJetCoreEmcal::kDetEmbPart) {
    AliParticleContainer *partCont2 = GetParticleContainer(2);
    jetContTrue = GetJetContainer(fJetContTrueName);
  }
	//if(fDebug) Printf("particle container 0 entries = %i \t1 entries = %i\t 2 entries = %i",partCont0->GetNParticles(),partCont1->GetNParticles(),partCont2->GetNParticles());

	if((fJetShapeType==AliAnalysisTaskJetCoreEmcal::kDetEmbPart && (!jetCont || !jetContPart || !jetContTrue)) || 
      ((fJetShapeType==AliAnalysisTaskJetCoreEmcal::kDetPart || fJetShapeType==AliAnalysisTaskJetCoreEmcal::kDetEmbDet) && (!jetCont || !jetContPart ))
    )
	{ // if jet containers not found
		AliError(Form("jet container not found - check name %s(base), %s (part) or %s (true)",fJetContName.Data(), fJetContPartName.Data(), fJetContTrueName.Data()));
		TIter next(&fJetCollArray);
		while ((jetCont = static_cast<AliJetContainer*>(next())))
			AliError(Form("%s",jetCont->GetName()));
		AliFatal("Exit...");
		return;
	}

	if(fDebug) {
		Printf("n particle jets = %i",jetContPart->GetNJets());
    Printf("n reco jets = %i",jetCont->GetNJets());
    if(fJetShapeType==AliAnalysisTaskJetCoreEmcal::kDetEmbPart) {
      Printf("n PYTHIA jets = %i",jetContTrue->GetNJets());
    }
	}

	// Background
	Double_t rho = 0;
	if (jetCont->GetRhoParameter()) rho = jetCont->GetRhoVal(); 

	// PYTHIA event weight
	// note - not used
	//  AliGenPythiaEventHeader *pyHeader = 0x0; //!<! Pythia header of the current external event
	//	AliAODEvent *ev = dynamic_cast<AliAODEvent*>(InputEvent());
	//  AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(ev->FindListObject(AliAODMCHeader::StdBranchName()));
	//  if (aodMCH) {
	//    for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
	//      pyHeader= dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
	//      if (pyHeader) break;
	//    }
	//  }
	//
	//	Double_t pythiaCrossSection = 0;
	//	Double_t pythiaTrials = 0;
	//	Double_t pythiaWeight = 0;
	//  if (pyHeader)
	//  {
	//		if(fDebug) Printf("have pythia header - get weight");
	//    pythiaCrossSection = pyHeader->GetXsection();
	//    pythiaTrials = pyHeader->Trials();
	//    //fPythiaPtHard = fPythiaHeader->GetPtHard();
	//		pythiaWeight = pythiaCrossSection / pythiaTrials;
	//	}
	//	if(fDebug) Printf("pythia weight = %f",pythiaWeight);

	for(auto jet1 : jetCont->accepted()) { // loop over hybrid jets

		Double_t ptJet1 = jet1->Pt();
		Double_t phiJet1 = jet1->Phi();
		Double_t area = jet1->Area();
		Double_t ptCorr = ptJet1-rho*area;
		Double_t ptLeadingTrackJet1 = jet1->GetLeadingTrack()->Pt();
    // closest jet
    Double_t distanceClosestJet1=-1, distanceClosestJet2=-1;
    Double_t ptJet2=0, phiJet2=0;
    Double_t ptJet3=0, phiJet3=0;
    Double_t ptLeadingTrackJet3=0;
		if(fDebug) Printf("--- jet pt hybrid = %f\t ",ptJet1);
    if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart) {
      auto jet2 = jet1->ClosestJet();
      if(!jet2) {
        //Printf("jet 2 cant be found");
        continue;}
      distanceClosestJet1 = jet1->ClosestJetDistance();
      ptJet2 = jet2->Pt();
      phiJet2 = jet2->Phi();
      fhPtDet->Fill(ptJet2);

      auto jet3 = jet2->ClosestJet();
      if(!jet3) {
        //Printf("jet3 can't be found");
        continue;
      }
      distanceClosestJet2 = jet2->ClosestJetDistance();
      fhPtDetMatchedToPart->Fill(ptJet2);
      ptJet3 = jet3->Pt();
      phiJet3 = jet3->Phi();
      ptLeadingTrackJet3 = jet3->GetLeadingTrack()->Pt();
      if(ptJet3==ptJet1) {
        fhPtPartMatchedWrongCent->Fill(ptJet3,fCent);
        continue;
      }
      fhPtPartMatched->Fill(ptJet3);
      fhPtPartMatchedCent->Fill(ptJet3,fCent);
    }

    if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) { // loop over detector jets
      auto jet3 = jet1->ClosestJet();
      if(!jet3) {
        if(fDebug) Printf("jet3 can't be found");
        continue;
      }
      distanceClosestJet1 = jet1->ClosestJetDistance();
      ptJet3 = jet3->Pt();
      phiJet3 = jet3->Phi();
      fhPtDetMatchedToPart->Fill(ptJet3);
    }

		if(fDebug) Printf("--- jet pt = jet hybrid pt = %f\t jet matched det pt = %f\t jet matched particle level pt = %f\t",ptJet1,ptJet2,ptJet3);

		fhPtDetPart->Fill(ptJet2,ptJet3);
		Double_t fraction = 1.;
    if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) fraction = jetCont->GetFractionSharedPt(jet1); 
    if(fDebug) Printf("FRACTION shared pT = %f",fraction);
    fhFractionSharedPtInclusive->Fill(ptCorr,fraction);
		if(fraction < fMinFractionSharedPt) continue;

		Double_t residual = (ptCorr - ptJet3) / ptJet3;
		Double_t residualPhi = (phiJet1 - phiJet3) / phiJet3;

		fhPtHybrDet->Fill(ptJet1,ptJet2);
		fhPtHybrPart->Fill(ptJet1,ptJet3);
		fhPtHybrPartCor->Fill(ptCorr,ptJet3);
		fhPhiHybrPartCor->Fill(phiJet1,phiJet3);

		fhResidual->Fill(residual);
		fhPtResidual->Fill(ptJet3,residual);
		fhPhiResidual->Fill(residualPhi);
		fhPhiPhiResidual->Fill(phiJet3,residualPhi);

		if(ptCorr<fMinEmbJetPt) continue;

		if(fFillInclusiveTree && fRandom->Integer(fRejectionFactorInclusiveJets)==0 ) {
      if(fMoreTreeVars) {
        fTreeVarsInclusiveMoreVars[0] = fCent;
        fTreeVarsInclusiveMoreVars[1] = ptJet1;
        fTreeVarsInclusiveMoreVars[2] = area;
        fTreeVarsInclusiveMoreVars[3] = ptCorr;
        fTreeVarsInclusiveMoreVars[4] = phiJet1;
        fTreeVarsInclusiveMoreVars[5] = ptJet3;
        fTreeVarsInclusiveMoreVars[6] = phiJet3;
        fTreeVarsInclusiveMoreVars[7] = ptLeadingTrackJet1;
        fTreeVarsInclusiveMoreVars[8] = ptLeadingTrackJet3;
        fTreeVarsInclusiveMoreVars[9] = ptJet2;
        fTreeVarsInclusiveMoreVars[10] = phiJet2;
        fTreeVarsInclusiveMoreVars[11] = distanceClosestJet1;
        fTreeVarsInclusiveMoreVars[12] = distanceClosestJet2;
        fTreeEmbInclusive->Fill();
      }
      else {
        fTreeVarsInclusive[0] = fCent;
        fTreeVarsInclusive[1] = ptJet1;
        fTreeVarsInclusive[2] = area;
        fTreeVarsInclusive[3] = ptCorr;
        fTreeVarsInclusive[4] = phiJet1;
        fTreeVarsInclusive[5] = ptJet3;
        fTreeVarsInclusive[6] = phiJet3;
        fTreeVarsInclusive[7] = ptLeadingTrackJet1;
        fTreeVarsInclusive[8] = ptLeadingTrackJet3;
        fTreeEmbInclusive->Fill();
      }
		}
	}
//	for(auto jettrue : jetContPart->accepted()) {
//		//				jettrue
//		Double_t ptTrue = jettrue->Pt();
//		Double_t phiTrue = jettrue->Phi();
//		Double_t etaTrue= jettrue->Eta();
//		auto jetmatched = jettrue->ClosestJet();
//		if(!jetmatched) continue;
//		Double_t ptMatched = jetmatched->Pt();
//		Double_t phiMatched = jetmatched->Phi();
//		Double_t etaMatched = jetmatched->Eta();
//		Double_t residual = (ptMatched - ptTrue) / ptTrue;
//		fhPtHybrTrue->Fill(ptTrue,ptMatched);
//		fhResidual->Fill(residual);
//		fhPtResidual->Fill(ptTrue,residual);
//	}
}



/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoJetLoop()
{
  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    UInt_t count = 0;
    for(auto jet : jetCont->accepted()) {
      if (!jet) continue;
      count++;

      histname = TString::Format("%s/histJetPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Pt());

      histname = TString::Format("%s/histJetPtLow_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Pt());

      histname = TString::Format("%s/histJetArea_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Area());

      histname = TString::Format("%s/histJetPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Phi());

      histname = TString::Format("%s/histJetEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, jet->Eta());

      if (jetCont->GetRhoParameter()) {
        histname = TString::Format("%s/histJetCorrPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area());

        histname = TString::Format("%s/histJetCorrPtLeadingTrackPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH2(histname, jet->Pt() - jetCont->GetRhoVal() * jet->Area(), jet->GetLeadingTrack()->Pt());

      }
    }
    histname = TString::Format("%s/histNJets_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoTrackLoop()
{
  AliClusterContainer* clusCont = GetClusterContainer(0);

  TString histname;
  TString groupname;
  UInt_t sumAcceptedTracks = 0;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    UInt_t count = 0;
    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;

      histname = TString::Format("%s/histTrackPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Pt());

      histname = TString::Format("%s/histTrackPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Phi());

      histname = TString::Format("%s/histTrackEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, part->Eta());

      if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
        const AliVTrack* track = static_cast<const AliVTrack*>(part);

        histname = TString::Format("%s/fHistDeltaEtaPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Eta() - track->GetTrackEtaOnEMCal());

        histname = TString::Format("%s/fHistDeltaPhiPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Phi() - track->GetTrackPhiOnEMCal());

        histname = TString::Format("%s/fHistDeltaPtvsPt_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, track->Pt(), track->Pt() - track->GetTrackPtOnEMCal());

        if (clusCont) {
          Int_t iCluster = track->GetEMCALcluster();
          if (iCluster >= 0) {
            AliVCluster* cluster = clusCont->GetAcceptCluster(iCluster);
            if (cluster) {
              histname = TString::Format("%s/fHistEoverPvsP_%d", groupname.Data(), fCentBin);
              fHistManager.FillTH2(histname, track->P(), cluster->GetNonLinCorrEnergy() / track->P());
            }
          }
        }
      }
    }
    sumAcceptedTracks += count;

    histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNTracks";
  fHistManager.FillTH1(histname, sumAcceptedTracks);
}

/**
 * This function performs a loop over the reconstructed EMCal clusters
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoClusterLoop()
{
  TString histname;
  TString groupname;
  UInt_t sumAcceptedClusters = 0;
  AliClusterContainer* clusCont = 0;
  TIter next(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(next()))) {
    groupname = clusCont->GetName();

    for(auto cluster : clusCont->all()) {
      if (!cluster) continue;

      if (cluster->GetIsExotic()) {
        histname = TString::Format("%s/histClusterEnergyExotic_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH1(histname, cluster->E());
      }
    }

    UInt_t count = 0;
    for(auto cluster : clusCont->accepted()) {
      if (!cluster) continue;
      count++;

      AliTLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);

      histname = TString::Format("%s/histClusterEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->E());

      histname = TString::Format("%s/histClusterNonLinCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetNonLinCorrEnergy());

      histname = TString::Format("%s/histClusterHadCorrEnergy_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, cluster->GetHadCorrEnergy());

      histname = TString::Format("%s/histClusterPhi_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Phi_0_2pi());

      histname = TString::Format("%s/histClusterEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH1(histname, nPart.Eta());
    }
    sumAcceptedClusters += count;

    histname = TString::Format("%s/histNClusters_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }

  histname = "fHistSumNClusters";
  fHistManager.FillTH1(histname, sumAcceptedClusters);
}

/**
 * This function performs a loop over the reconstructed EMCal cells
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskJetCoreEmcal::DoCellLoop()
{
  if (!fCaloCells) return;

  TString histname;

  const Short_t ncells = fCaloCells->GetNumberOfCells();

  histname = TString::Format("%s/histNCells_%d", fCaloCellsName.Data(), fCentBin);
  fHistManager.FillTH1(histname, ncells);

  histname = TString::Format("%s/histCellEnergy_%d", fCaloCellsName.Data(), fCentBin);
  for (Short_t pos = 0; pos < ncells; pos++) {
    Double_t amp   = fCaloCells->GetAmplitude(pos);
    
    fHistManager.FillTH1(histname, amp);
  }
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskJetCoreEmcal::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskJetCoreEmcal::Run()
{
  return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskJetCoreEmcal::Terminate(Option_t *) 
{
}

Int_t  AliAnalysisTaskJetCoreEmcal::SelectTrigger(TList *list,Double_t minT,Double_t maxT,Int_t &number){

	Int_t index=-1;
	Int_t triggers[100];

	for(Int_t cr=0;cr<100;cr++) triggers[cr]=-1;

	Int_t im=0;

	AliParticleContainer* partCont = 0x0;
	AliParticleContainer* partContDet = 0x0;
	if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPartCorr || fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbPart) partCont = GetParticleContainer(1);
  else if(fJetShapeType == AliAnalysisTaskJetCoreEmcal::kDetEmbDet) {
    partCont = GetParticleContainer(0);
    partContDet = GetParticleContainer(1);
  }
	else partCont = GetParticleContainer(0);
	UInt_t iCount = 0;
  // loop over first container
	for(auto part : partCont->accepted()) {
		if (!part) continue;
		list->Add(part);
		iCount++;
		if(part->Pt()>=minT && part->Pt()<maxT){
			triggers[im]=iCount-1;
			im=im+1;
      fhSelectedTrigger->Fill(0.5,part->Pt());
      //Printf("Pb-Pb data trigger added - pt = %f, number = %i",part->Pt(),im);
		}
	}
  // loop over second container (embedded) if requested
  if(partContDet) {
    for(auto part : partContDet->accepted()) {
      if (!part) continue;
      list->Add(part);
      iCount++;
      if(part->Pt()>=minT && part->Pt()<maxT){
        triggers[im]=iCount-1;
        im=im+1;
        fhSelectedTrigger->Fill(1.5,part->Pt());
        //Printf("embedded trigger added - pt = %f, number = %i",part->Pt(),im);
      }
    }
  }
	number=im;
	Int_t rd=0;
	if(im>0) rd=fRandom->Integer(im);
	index=triggers[rd];
	return index;
}

Double_t AliAnalysisTaskJetCoreEmcal::RelativePhi(Double_t mphi,Double_t vphi){

	// get relative DeltaPhi from -pi < DeltaPhi < pi

  if (vphi < -1*TMath::Pi()) vphi += (2*TMath::Pi());
  else if (vphi > TMath::Pi()) vphi -= (2*TMath::Pi());
  if (mphi < -1*TMath::Pi()) mphi += (2*TMath::Pi());
  else if (mphi > TMath::Pi()) mphi -= (2*TMath::Pi());
  double dphi = mphi-vphi;
  if (dphi < -1*TMath::Pi()) dphi += (2*TMath::Pi());
  else if (dphi > TMath::Pi()) dphi -= (2*TMath::Pi());
  return dphi;//dphi in [-Pi, Pi]
}

Int_t AliAnalysisTaskJetCoreEmcal::GetPhiBin(Double_t phi)
{
    Int_t phibin=-1;
    if(!(TMath::Abs(phi)<=2*TMath::Pi())) return -1;
    Double_t phiwrtrp=TMath::ACos(TMath::Abs(TMath::Cos(phi)));
    phibin=Int_t(fNRPBins*phiwrtrp/(0.5*TMath::Pi()));
    //if(phibin<0||phibin>=fNRPBins){AliError("Phi Bin not defined");}
    return phibin;
}
