#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TVector2.h>
#include "TClonesArray.h"
#include <TList.h>
#include <TLorentzVector.h>
#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"
#include "AliHFJetTaggingIP.h"
#include "AliHFJetsTagging.h"
#include "TClonesArray.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDv0.h"
#include "AliAODv0.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliRDHFJetsCuts.h"
#include <AliPID.h>
#include <AliPIDCombined.h>
#include <AliPIDResponse.h>
// Includes to add invariant mass cross check fits

#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexerTracks.h"
#include "AliKFParticle.h"
#include "AliAODVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDVertex.h"
#include "AliKFVertex.h"
#include "AliStack.h"
#include "TRandom3.h"

#include "AliAnalysisTaskEmcalJetBJetTaggingIP.h"

// Changelog  MC Jet Container retrieval fixed
ClassImp(AliAnalysisTaskEmcalJetBJetTaggingIP)
    ////////////////////////////////////////////////////////////////////////////////
    /// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::AliAnalysisTaskEmcalJetBJetTaggingIP
    ///
    AliAnalysisTaskEmcalJetBJetTaggingIP::AliAnalysisTaskEmcalJetBJetTaggingIP()
    : AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetBJetTaggingIP", kTRUE)
    , fJetsCont(NULL)
    , fJetsContMC(NULL)
    , fTracksCont(NULL)
    , fTaggingHFClass(new AliHFJetsTagging("fTaggingHFClass"))
    , fTrackCountingTagger(new AliHFJetTaggingIP())
    , fRandom(new TRandom3(0))
    , fUtils(new AliAnalysisUtils())
    , fJetCutsHF(new AliRDHFJetsCuts())
    , fPIDResponse(NULL)
    , fPIDCombined(NULL)
    , fMCHeader(NULL)
    , fSelectPtHard(kFALSE)
    , fPtHardBin("")
    , fPtHardFlavour("")
    , fMCparticles(NULL)
    , fCurrentNProngVertex(NULL)
    , fCurrentNProngs(0)
    , fIsMC(kFALSE)
    , fIsTrackQA(kFALSE)
    , fIsTrackQAConstituent(kFALSE)
    , fUseCorrectedJetPt(kFALSE)
    , fDoRandomCones(kFALSE)
    , fVetexingMassFitTest(kFALSE)
    , fVetexingPtRelTest(kFALSE)
    , fUseAtlasSignTCCalculation(kFALSE)
    , fUseEventSelection(0)
    , fUseJetSelection(0)
    , fUseMCTagger(0)
    , fSigmaTPCElectronLow(-0.5)
    , fSigmaTPCElectronHigh(3.)
    , fSigmaTOFElectronLow(-3)
    , fSigmaTOFElectronHigh(3.)
    , fSigmaITSElectronLow(-5)
    , fSigmaITSElectronHigh(1.)
    , fpTRelSoftElectronPt(1.)
    , fpTRelSoftElectronJetPt(5.)
    , fCurrentTrackContainerSecVtx(NULL)
    , fhist_Events(NULL)
    , fhist_Jets(NULL)
    , fhists_SPD_cluster_vs_tracklet_correlation(NULL)
    , fhists_SPD_cluster_vs_tracklet_correlation_PostSelection(NULL)
    , fhist_MonteCarloFlavour(NULL)
    , fhist_Tracks_Eta_Phi(NULL)
    , fhist_Tracks_Eta_Phi_Bit4(NULL)
    , fhist_Tracks_Eta_Phi_Bit9(NULL)
    , fhist_QualityClasses(NULL)
    , fhist_Jet_Eta_Phi(NULL)
    , fhist_Jet_Nconst_Pt(NULL)
    , fhist_Jet_Pt(NULL)
    , fhist_Jet_Background_Fluctuation(NULL)
    , fhist_Rho(NULL)
    , fhist_parton_genjet_dR(NULL)
    , fhist_parton_genjet_pT(NULL)
    , fhist_parton_genjet_Eta(NULL)
    , fhist_parton_genjet_Phi(NULL)
    , fHistTPCnSigmaElectron(NULL)
    , fHistTPCnSigmaPion(NULL)
{
    memset(fhist_TC_sIP_Pt, 0, sizeof fhist_TC_sIP_Pt);
    memset(fhist_TC_Eta_Phi, 0, sizeof fhist_TC_sIP_Pt);
    memset(fhist_QualityClasses_sIP, 0, sizeof fhist_QualityClasses_sIP);
    memset(fhist_QualityClasses_Eta_Phi, 0, sizeof fhist_QualityClasses_Eta_Phi);
    memset(fhist_momentum_response, 0, sizeof fhist_momentum_response);
    SetMakeGeneralHistograms(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::AliAnalysisTaskEmcalJetBJetTaggingIP
/// \param name
///
AliAnalysisTaskEmcalJetBJetTaggingIP::AliAnalysisTaskEmcalJetBJetTaggingIP(const char* name)
    : AliAnalysisTaskEmcalJet(name, kTRUE)
    , fJetsCont(NULL)
    , fJetsContMC(NULL)
    , fTracksCont(NULL)
    , fTaggingHFClass(new AliHFJetsTagging("fTaggingHFClass"))
    , fTrackCountingTagger(new AliHFJetTaggingIP())
    , fRandom(new TRandom3(0))
    , fUtils(new AliAnalysisUtils())
    , fJetCutsHF(new AliRDHFJetsCuts())
    , fPIDResponse(NULL)
    , fPIDCombined(NULL)
    , fMCHeader(NULL)
    , fSelectPtHard(kFALSE)
    , fPtHardBin("")
    , fPtHardFlavour("")
    , fMCparticles(NULL)
    , fCurrentNProngVertex(NULL)
    , fCurrentNProngs(0)
    , fIsMC(kFALSE)
    , fIsTrackQA(kFALSE)
    , fIsTrackQAConstituent(kFALSE)
    , fUseCorrectedJetPt(kFALSE)
    , fDoRandomCones(kFALSE)
    , fVetexingMassFitTest(kFALSE)
    , fVetexingPtRelTest(kFALSE)
    , fUseAtlasSignTCCalculation(kFALSE)
    , fUseEventSelection(0)
    , fUseJetSelection(0)
    , fUseMCTagger(0)
    , fSigmaTPCElectronLow(-0.5)
    , fSigmaTPCElectronHigh(3.)
    , fSigmaTOFElectronLow(-3)
    , fSigmaTOFElectronHigh(3.)
    , fSigmaITSElectronLow(-5)
    , fSigmaITSElectronHigh(1.)
    , fpTRelSoftElectronPt(1.)
    , fpTRelSoftElectronJetPt(5.)
    , fCurrentTrackContainerSecVtx(NULL)
    , fhist_Events(NULL)
    , fhist_Jets(NULL)
    , fhists_SPD_cluster_vs_tracklet_correlation(NULL)
    , fhists_SPD_cluster_vs_tracklet_correlation_PostSelection(NULL)
    , fhist_MonteCarloFlavour(NULL)
    , fhist_Tracks_Eta_Phi(NULL)
    , fhist_Tracks_Eta_Phi_Bit4(NULL)
    , fhist_Tracks_Eta_Phi_Bit9(NULL)
    , fhist_QualityClasses(NULL)
    , fhist_Jet_Eta_Phi(NULL)
    , fhist_Jet_Nconst_Pt(NULL)
    , fhist_Jet_Pt(NULL)
    , fhist_Jet_Background_Fluctuation(NULL)
    , fhist_Rho(NULL)
    , fhist_parton_genjet_dR(NULL)
    , fhist_parton_genjet_pT(NULL)
    , fhist_parton_genjet_Eta(NULL)
    , fhist_parton_genjet_Phi(NULL)
    , fHistTPCnSigmaElectron(NULL)
    , fHistTPCnSigmaPion(NULL)
{
    // Standard constructor.
    memset(fhist_TC_sIP_Pt, 0, sizeof fhist_TC_sIP_Pt);
    memset(fhist_TC_Eta_Phi, 0, sizeof fhist_TC_sIP_Pt);
    memset(fhist_QualityClasses_sIP, 0, sizeof fhist_QualityClasses_sIP);
    memset(fhist_QualityClasses_Eta_Phi, 0, sizeof fhist_QualityClasses_Eta_Phi);
    memset(fhist_momentum_response, 0, sizeof fhist_momentum_response);
    SetMakeGeneralHistograms(kTRUE);
}

////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::Run
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::Run()
{
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse)
	AliFatal("This Task needs the PID response attached to the inputHandler");

    // Event selection called automatically

    if(fIsMC)
	fMCparticles = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));

    AliEmcalJet* curjet = NULL;
    AliEmcalJet* curjetMatched = NULL;
    Double_t tagvalue[3] = { 0., 0., 0. };
    Bool_t istagged[3] = { kFALSE, kFALSE, kFALSE };
    Int_t flavourtag = 0;
    fTrackCountingTagger->SetEvent(InputEvent());
    if(!fJetsCont) {
	AliError("Missing jet container!");
	return kFALSE;
    }

    if(fDoRandomCones) {
	Double_t randomConePt = GetDeltaPtRandomCone();
	fhist_Jet_Background_Fluctuation->Fill(randomConePt);
    }

    fhist_Rho->Fill(fJetsCont->GetRhoVal());

    // Loop over all available jets
    for(Int_t ijet = 0; ijet < fJetsCont->GetNJets(); ++ijet) {
	memset(tagvalue, 0, sizeof tagvalue);
	memset(istagged, 0, sizeof istagged);
	curjet = NULL;
	curjetMatched = NULL;
	flavourtag = 0;

	fCurrentNProngs = 0;
	curjet = fJetsCont->GetJet(ijet);
	if(!IsJetSelected(curjet))
	    continue;
	Double_t jetPt = 0.;
	fUseCorrectedJetPt ? jetPt = GetPtCorrected(curjet) : jetPt = curjet->Pt();

	if(fIsMC) {
	    curjetMatched = curjet->MatchedJet();
	    if(curjetMatched) {
		fhist_Jets->Fill("Matched", jetPt, 1.);
		AddTagJet(curjetMatched);
		if(curjetMatched->TestFlavourTag(kLFgJet))
		    flavourtag = 1;
		else if(curjetMatched->TestFlavourTag(kBeautyJet))
		    flavourtag = 2;
		else if(curjetMatched->TestFlavourTag(kCharmJet))
		    flavourtag = 3;
		Double_t ptm = curjetMatched->Pt();
		Double_t ptmcorr = ptm - fJetsContMC->GetRhoVal() * curjetMatched->Area();
		Double_t ptr = curjet->Pt();
		Double_t ptrcorr = GetPtCorrected(curjet);
		fhist_momentum_response[0][0]->Fill(ptm, ptr);
		fhist_momentum_response[0][1]->Fill(ptmcorr, ptrcorr);
		if(flavourtag > 0) {
		    fhist_momentum_response[flavourtag][0]->Fill(ptm, ptr);
		    fhist_momentum_response[flavourtag][1]->Fill(ptmcorr, ptrcorr);
		}
	    }
	}

	// Standard track counting algorithm

	fhist_Jet_Pt->Fill(jetPt);
	fhist_Jet_Eta_Phi->Fill(curjet->Eta(), curjet->Phi());
	fhist_Jet_Nconst_Pt->Fill(jetPt, curjet->GetNumberOfTracks());

	Double_t n3value = -9999.;

	if(fTrackCountingTagger->GetJetDiscriminator(curjet, tagvalue, istagged)) {
	    if(istagged[0]) {
		fhist_Jets->Fill("TaggingN1", jetPt, 1.);
		fhist_TC_sIP_Pt[0][0][0]->Fill(tagvalue[0], jetPt);
		fhist_TC_Eta_Phi[0][0][0]->Fill(curjet->Eta(), curjet->Phi());
		if(fIsMC)
		    if(flavourtag > 0) {
			fhist_TC_sIP_Pt[0][flavourtag][0]->Fill(tagvalue[0], jetPt);
			fhist_TC_Eta_Phi[0][flavourtag][0]->Fill(curjet->Eta(), curjet->Phi());
		    }
	    }
	    if(istagged[1]) {
		fhist_Jets->Fill("TaggingN2", jetPt, 1.);
		fhist_TC_sIP_Pt[1][0][0]->Fill(tagvalue[1], jetPt);
		fhist_TC_Eta_Phi[1][0][0]->Fill(curjet->Eta(), curjet->Phi());
		if(fIsMC)
		    if(flavourtag > 0) {
			fhist_TC_Eta_Phi[1][flavourtag][0]->Fill(curjet->Eta(), curjet->Phi());
			fhist_TC_sIP_Pt[1][flavourtag][0]->Fill(tagvalue[1], jetPt);
		    }
	    }
	    if(istagged[2]) {
		n3value = tagvalue[2];
		fhist_Jets->Fill("TaggingN3", jetPt, 1.);
		fhist_TC_sIP_Pt[2][0][0]->Fill(tagvalue[2], jetPt);
		fhist_TC_Eta_Phi[2][0][0]->Fill(curjet->Eta(), curjet->Phi());

		if(fIsMC)
		    if(flavourtag > 0) {
			fhist_TC_sIP_Pt[2][flavourtag][0]->Fill(tagvalue[2], jetPt);
			fhist_TC_Eta_Phi[2][flavourtag][0]->Fill(curjet->Eta(), curjet->Phi());
		    }
	    }
	}

	if(fVetexingPtRelTest)
	    ProcessPtRelTemplateAnalysis(curjet, jetPt, flavourtag, n3value);
	if(fVetexingMassFitTest)
	    ProcessVtxMassTemplateAnalysis(curjet, jetPt, flavourtag, n3value);

	// Track quality track counting
	for(int i = 1; i < 5; ++i) {
	    memset(tagvalue, 0, sizeof tagvalue);
	    memset(istagged, 0, sizeof istagged);

	    if(fTrackCountingTagger->GetJetDiscriminatorQualityClass(i, curjet, tagvalue, istagged)) {
		if(istagged[0]) {
		    fhist_TC_sIP_Pt[0][0][i]->Fill(tagvalue[0], jetPt);
		    fhist_TC_Eta_Phi[0][0][i]->Fill(curjet->Eta(), curjet->Phi());
		    if(fIsMC)
			if(flavourtag > 0) {
			    fhist_TC_sIP_Pt[0][flavourtag][i]->Fill(tagvalue[0], jetPt);
			    fhist_TC_Eta_Phi[0][flavourtag][i]->Fill(curjet->Eta(), curjet->Phi());
			}
		}
		if(istagged[1]) {
		    fhist_TC_sIP_Pt[1][0][i]->Fill(tagvalue[1], jetPt);
		    fhist_TC_Eta_Phi[1][0][i]->Fill(curjet->Eta(), curjet->Phi());

		    if(fIsMC)
			if(flavourtag > 0) {
			    fhist_TC_sIP_Pt[1][flavourtag][i]->Fill(tagvalue[1], jetPt);
			    fhist_TC_Eta_Phi[0][flavourtag][i]->Fill(curjet->Eta(), curjet->Phi());
			}
		}
		if(istagged[2]) {
		    fhist_TC_sIP_Pt[2][0][i]->Fill(tagvalue[2], jetPt);
		    fhist_TC_Eta_Phi[2][0][i]->Fill(curjet->Eta(), curjet->Phi());
		    if(fIsMC)
			if(flavourtag > 0) {
			    fhist_TC_Eta_Phi[2][flavourtag][i]->Fill(curjet->Eta(), curjet->Phi());
			    fhist_TC_sIP_Pt[2][flavourtag][i]->Fill(tagvalue[2], jetPt);
			}
		}
	    }
	}

	// Constituent QA
	if(fIsTrackQAConstituent) {
	    RunQATracksJet(curjet);
	}
    }
    // Run Event Track QA
    if(fIsTrackQA) {
	RunQATracksEvent();
    }

    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::RunQATracksJet
/// \param jet
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::RunQATracksJet(const AliEmcalJet* jet)
{
    Double_t dv[2] = { 0., 0. };
    Double_t covv[3] = { 0., 0. };
    Double_t sip = 0.;
    AliAODTrack* track = 0x0;
    AliAODMCParticle* part = 0x0;

    if(!jet)
	return kFALSE;
    Int_t ntracks = (Int_t)jet->GetNumberOfTracks();

    for(Int_t itrack = 0; itrack < ntracks; ++itrack) {

	memset(dv, 0, sizeof dv);
	memset(covv, 0, sizeof covv);
	sip = 0.;

	track = (AliAODTrack*)(((AliPicoTrack*)(fTracksCont->GetParticle(jet->TrackAt(itrack))))->GetTrack());
	if(!track)
	    continue;
	if(fabs(track->Eta()) > 0.9)
	    continue;
	if(track->Pt() < 1)
	    continue;

	Bool_t isPrimary = kFALSE;
	if(fIsMC && fMCparticles) {
	    part = dynamic_cast<AliAODMCParticle*>(fMCparticles->At(TMath::Abs(track->GetLabel())));
	    if(part) {
		isPrimary = part->IsPhysicalPrimary();
	    }
	    fTaggingHFClass->GetSignedRPhiImpactParameter(InputEvent(), track, jet, sip, dv, covv);
	}

	fhist_QualityClasses->Fill("all", 1.);
	if(IsQuality(track, kQtyVeryGood)) {
	    fhist_QualityClasses->Fill("kQtyVeryGood", 1.);
	    if(fIsMC) {
		if(isPrimary) {
		    fhist_QualityClasses_sIP[0][0]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[0][0]->Fill(track->Eta(), track->Phi());
		} else {
		    fhist_QualityClasses_sIP[0][1]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[0][1]->Fill(track->Eta(), track->Phi());
		}
	    }
	}
	if(IsQuality(track, kQtyGood)) {
	    fhist_QualityClasses->Fill("kQtyGood", 1.);
	    if(fIsMC) {
		if(isPrimary) {
		    fhist_QualityClasses_sIP[1][0]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[1][0]->Fill(track->Eta(), track->Phi());
		} else {
		    fhist_QualityClasses_sIP[1][1]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[1][1]->Fill(track->Eta(), track->Phi());
		}
	    }
	}
	if(IsQuality(track, kQtyMedium)) {
	    fhist_QualityClasses->Fill("kQtyMedium", 1.);
	    if(fIsMC) {
		if(isPrimary) {
		    fhist_QualityClasses_sIP[2][0]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[2][0]->Fill(track->Eta(), track->Phi());
		} else {
		    fhist_QualityClasses_sIP[2][1]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[2][1]->Fill(track->Eta(), track->Phi());
		}
	    }
	}
	if(IsQuality(track, kQtyBad)) {
	    fhist_QualityClasses->Fill("kQtyBad", 1.);
	    if(fIsMC) {
		if(isPrimary) {
		    fhist_QualityClasses_sIP[3][0]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[3][0]->Fill(track->Eta(), track->Phi());
		} else {
		    fhist_QualityClasses_sIP[3][1]->Fill(sip, track->Pt());
		    fhist_QualityClasses_Eta_Phi[3][1]->Fill(track->Eta(), track->Phi());
		}
	    }
	}
    }
    return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::RunQATracksEvent
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::RunQATracksEvent()
{

    Int_t ntracks = (Int_t)InputEvent()->GetNumberOfTracks();
    AliAODTrack* track = 0x0;
    for(Int_t itrack = 0; itrack < ntracks; ++itrack) {
	track = (AliAODTrack*)InputEvent()->GetTrack(itrack);
	if(fabs(track->Eta()) > 0.9)
	    continue;
	if(track->Pt() < 0.150)
	    continue;
	if(!((track)->TestFilterBit(1 << 4) || ((track)->TestFilterBit(1 << 9))))
	    continue;
	fhist_Tracks_Eta_Phi->Fill(track->Eta(), track->Phi());
	if((track)->TestFilterBit(1 << 4) && !((track)->TestFilterBit(1 << 9)))
	    fhist_Tracks_Eta_Phi_Bit4->Fill(track->Eta(), track->Phi());
	if(!((track)->TestFilterBit(1 << 4)) && (track)->TestFilterBit(1 << 9))
	    fhist_Tracks_Eta_Phi_Bit9->Fill(track->Eta(), track->Phi());
    }

    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelected
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelected()
{
    AliAODEvent* aev = dynamic_cast<AliAODEvent*>(InputEvent());

    if(fIsMC && fSelectPtHard) {
	fMCHeader = NULL;
	fMCHeader = dynamic_cast<AliAODMCHeader*>(aev->FindListObject(AliAODMCHeader::StdBranchName()));
	TString title = fMCHeader->GetGeneratorName();
	if(!title.Contains(fPtHardBin) || !title.Contains(fPtHardFlavour))
	    return kFALSE;
    }

    // SPD Cluster vs Tracklet plot to estimate pileup effect
    Int_t nClustersLayer0 = aev->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = aev->GetNumberOfITSClusters(1);
    Int_t nTracklets = aev->GetMultiplicity()->GetNumberOfTracklets();
    fhists_SPD_cluster_vs_tracklet_correlation->Fill(nTracklets, nClustersLayer0 + nClustersLayer1);

    switch(fUseEventSelection) {
    case 0:
	if(!IsEventSelectedLegacy(aev))
	    return kFALSE;
	break;
    case 1:
	if(!IsEventSelectedpp(aev))
	    return kFALSE;
	break;
    case 2:
	if(!IsEventSelectedpA(aev))
	    return kFALSE;
	break;
    case 3:
	if(!IsEventSelectedHF(aev))
	    return kFALSE;
	break;
    }
    fhists_SPD_cluster_vs_tracklet_correlation_PostSelection->Fill(nTracklets, nClustersLayer0 + nClustersLayer1);
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedLegacy
/// \param aev
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedLegacy(AliAODEvent* aev)
{
    UInt_t res = 0;
    if(aev)
	res = ((AliVAODHeader*)aev->GetHeader())->GetOfflineTrigger();
    if((res & AliVEvent::kAny) == 0) {
	return kFALSE;
    }
    fhist_Events->Fill("AliVEvent::kAny", 1.);
    if((res & fOffTrigger) == 0) {
	return kFALSE;
    }
    fhist_Events->Fill("AliVEvent::kMB", 1.);

    if(!aev->GetPrimaryVertex() || aev->GetPrimaryVertex()->GetNContributors() < 1)
	fhist_Events->Fill("NoVertex", 1.);

    if(aev->IsPileupFromSPD(5, 0.8, 3.0, 2.0, 5.0)) {
	return kFALSE;
    }

    fhist_Events->Fill("!Pileup SPD", 1.);
    if(fabs(aev->GetPrimaryVertex()->GetZ()) > 10.) {
	return kFALSE;
    }
    fhist_Events->Fill("Vertex  z < 10 cm", 1.);
    const AliVVertex* trkVtx = dynamic_cast<const AliVVertex*>(aev->GetPrimaryVertex());
    const AliVVertex* spdVtx = dynamic_cast<const AliVVertex*>(aev->GetPrimaryVertexSPD());
    TString vtxTtl = trkVtx->GetTitle();
    if(!vtxTtl.Contains("VertexerTracks")) {
	return kFALSE;
    }
    fhist_Events->Fill("No VertexerTracks", 1.);
    Double_t cov[6] = { 0 };
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if(spdVtx->IsFromVertexerZ() && (zRes > 0.25)) {
	return kFALSE;
    }
    if((TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ()) > 0.5)) {
	return kFALSE;
    }
    fhist_Events->Fill("Vertex Z Resoulution", 1.);
    if(trkVtx->GetNContributors() < 2)
	return kFALSE;
    fhist_Events->Fill(">1 contributors", 1.);
    if(spdVtx->GetNContributors() < 1)
	return kFALSE;
    fhist_Events->Fill(">0 contributors SPD", 1.);
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedpA
/// \param aev
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedpA(AliAODEvent* aev)
{
    if(fUtils->IsFirstEventInChunk(aev))
	return kFALSE;
    if(!fUtils->IsVertexSelected2013pA(aev))
	return kFALSE;
    if(!fUtils->IsPileUpEvent(aev))
	return kFALSE;
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedpp
/// \param aev
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedpp(AliAODEvent* aev)
{
    if(!fUtils->IsVertexSelected2013pA(aev))
	return kFALSE;
    if(!fUtils->IsPileUpEvent(aev))
	return kFALSE;
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedHF
/// \param aev
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedHF(AliAODEvent* aev)
{
    if(fJetCutsHF->IsEventSelected(aev))
	return kTRUE;
    else
	return kFALSE;
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelectedHF
/// \param jet
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelectedHF(const AliEmcalJet* jet)
{
    if(fJetCutsHF->IsJetSelected(jet))
	return kTRUE;
    else
	return kFALSE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelectedLegacy
/// \param jet
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelectedLegacy(const AliEmcalJet* jet)
{
    Double_t jetPt = 0.;
    fUseCorrectedJetPt ? jetPt = GetPtCorrected(jet) : jetPt = jet->Pt();
    if(!(jetPt > 0))
	return kFALSE;
    fhist_Jets->Fill("In", jetPt, 1.);
    if(!jet)
	return kFALSE;
    fhist_Jets->Fill("Pointer", jetPt, 1.);
    Double_t jetradius = fJetsCont->GetJetRadius();
    jetradius *= jetradius;
    if(jetPt < 1.0)
	return kFALSE;
    fhist_Jets->Fill("Pt", jetPt, 1.);
    if(fabs(jet->Eta()) > 0.5)
	return kFALSE;
    fhist_Jets->Fill("Eta", jetPt, 1.);
    if(jet->Area() < 0.6 * jetradius * TMath::Pi())
	return kFALSE;
    fhist_Jets->Fill("Area", jetPt, 1.);
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelected
/// \param jet
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelected(const AliEmcalJet* jet)
{
    switch(fUseJetSelection) {
    case 0:
	return IsJetSelectedLegacy(jet);
	break;
    case 1:
	return IsJetSelectedHF(jet);
	break;
    }
    return kFALSE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::GetJetCutsHF
/// \return
///
AliRDHFJetsCuts* AliAnalysisTaskEmcalJetBJetTaggingIP::GetJetCutsHF()
{
    return fJetCutsHF;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::AddTagJet
/// \param jet
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::AddTagJet(AliEmcalJet* jet)
{
    if(!jet)
	return kFALSE;
    Int_t partonPDG = 0;
    Double_t jetPt = 0.;
    jetPt = jet->Pt();
    if(fUseMCTagger == 0) {
	AliAODMCParticle* parton = NULL;
	parton = fTaggingHFClass->IsMCJetParton(fMCparticles, jet, 0.7);
	if(!parton)
	    return kFALSE;
	partonPDG = abs(parton->PdgCode());
	Double_t deltaR = sqrt((parton->Eta() - jet->Eta()) * (parton->Eta() - jet->Eta()) +
	                       TVector2::Phi_mpi_pi((parton->Phi() - jet->Phi())) *
	                           TVector2::Phi_mpi_pi((parton->Phi() - jet->Phi())));
	fhist_parton_genjet_dR->Fill(deltaR);
	fhist_parton_genjet_pT->Fill(parton->Pt(), jet->Pt());
	fhist_parton_genjet_Eta->Fill(parton->Eta() - jet->Eta());
	fhist_parton_genjet_Phi->Fill(TVector2::Phi_mpi_pi(parton->Phi() - jet->Phi()));
    } else if(fUseMCTagger == 1) {
	AliAODMCParticle* meson = NULL;
	meson = fTaggingHFClass->IsMCJetMeson(fMCparticles, jet, 0.7);
	if(meson) {
	    Int_t mesonPDG = abs(meson->PdgCode());
	    if(mesonPDG >= 500 && mesonPDG < 600)
		partonPDG = 5;
	    else if(mesonPDG >= 400 && mesonPDG < 500)
		partonPDG = 4;
	    else if(mesonPDG >= 5000 && mesonPDG < 6000)
		partonPDG = 5;
	    else if(mesonPDG >= 4000 && mesonPDG < 5000)
		partonPDG = 4;
	} else
	    partonPDG = 1; // Set to 1 for LFg
    }
    fhist_MonteCarloFlavour->Fill("Sum", jetPt, 1.);

    switch(partonPDG) {
    case 21:
	jet->AddFlavourTag(kLFgJet);
	fhist_MonteCarloFlavour->Fill("LFg", jetPt, 1.);
	break;
    case 1:
	jet->AddFlavourTag(kLFgJet);
	fhist_MonteCarloFlavour->Fill("LFg", jetPt, 1.);
	break;
    case 2:
	jet->AddFlavourTag(kLFgJet);
	fhist_MonteCarloFlavour->Fill("LFg", jetPt, 1.);
	break;
    case 3:
	jet->AddFlavourTag(kLFgJet);
	fhist_MonteCarloFlavour->Fill("LFg", jetPt, 1.);
	break;
    case 4:
	jet->AddFlavourTag(kCharmJet);
	fhist_MonteCarloFlavour->Fill("Charm", jetPt, 1.);
	break;
    case 5:
	jet->AddFlavourTag(kBeautyJet);
	fhist_MonteCarloFlavour->Fill("Beauty", jetPt, 1.);
	break;
    default:
	fhist_MonteCarloFlavour->Fill("Untagged", jetPt, 1.);
	break;
    }
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsQuality
/// \param track
/// \param qclass
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsQuality(const AliAODTrack* track, EQualityClass qclass)
{
    if(!track)
	return kFALSE;
    if(track->Pt() < 1.)
	return kFALSE;
    ULong_t status = track->GetStatus();

    if(!(status & AliAODTrack::kTPCrefit))
	return kFALSE;
    if(!(status & AliAODTrack::kITSrefit))
	return kFALSE;

    int nSPDHits = 0;
    if(track->HasPointOnITSLayer(0))
	nSPDHits++;
    if(track->HasPointOnITSLayer(1))
	nSPDHits++;
    int nITSHits = nSPDHits;
    for(int j = 2; j < 6; ++j)
	if(track->HasPointOnITSLayer(j))
	    nITSHits++;
    int nTPCcls = 0;
    nTPCcls = ((AliAODTrack*)track)->GetTPCNcls();
    Float_t cRatioTPC = track->GetTPCNclsF() > 0. ?
                            static_cast<Float_t>(track->GetTPCNcls()) / static_cast<Float_t>(track->GetTPCNclsF()) :
                            1.;
    Bool_t isV0Daughter = kFALSE;
    Double_t v0Radius = 0.;
    isV0Daughter = IsV0DaughterRadius(track, v0Radius);

    switch(qclass) {
    case kQtyVeryGood:
	if(nSPDHits < 2)
	    return kFALSE;
	if(nITSHits < 4)
	    return kFALSE;
	if(nTPCcls < 90)
	    return kFALSE;
	if(cRatioTPC < 0.6)
	    return kFALSE;
	if(isV0Daughter && v0Radius > 2.0)
	    return kFALSE;
	break;
    case kQtyGood:
	if(nSPDHits < 1)
	    return kFALSE;
	if(nITSHits < 4)
	    return kFALSE;
	if(nTPCcls < 90)
	    return kFALSE;
	if(cRatioTPC < 0.6)
	    return kFALSE;
	if(isV0Daughter && v0Radius > 2.0)
	    return kFALSE;
	break;
    case kQtyMedium:
	if(nSPDHits < 1)
	    return kFALSE;
	if(nITSHits < 3)
	    return kFALSE;
	if(nTPCcls < 90)
	    return kFALSE;
	if(cRatioTPC < 0.6)
	    return kFALSE;
	break;
    case kQtyBad:
	if(nITSHits < 3)
	    return kFALSE;
	if(nTPCcls < 80)
	    return kFALSE;
	if(cRatioTPC < 0.6)
	    return kFALSE;
	break;
    default:
	return kFALSE;
	break;
    }
    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsV0DaughterRadius
/// \param track
/// \param Radius
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsV0DaughterRadius(const AliAODTrack* track, Double_t& Radius)
{
    AliAODv0* v0aod = NULL;
    int posid = -1;
    int negid = -1;
    int trackid = -1;
    Double_t P[3];
    for(int i = 0; i < InputEvent()->GetNumberOfV0s(); ++i) {
	memset(P, 0, sizeof P);
	v0aod = ((AliAODEvent*)InputEvent())->GetV0(i);
	posid = v0aod->GetPosID();
	negid = v0aod->GetNegID();
	trackid = track->GetID();
	if(posid == trackid || negid == trackid) {
	    P[0] = v0aod->DecayVertexV0X();
	    P[1] = v0aod->DecayVertexV0Y();
	    P[2] = v0aod->DecayVertexV0Z();
	    Radius = sqrt(P[0] * P[0] + P[1] * P[1]);
	    return kTRUE;
	}
    }
    return kFALSE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::~AliAnalysisTaskEmcalJetBJetTaggingIP
///
AliAnalysisTaskEmcalJetBJetTaggingIP::~AliAnalysisTaskEmcalJetBJetTaggingIP()
{
    // Destructor.
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::UserCreateOutputObjects
///
void AliAnalysisTaskEmcalJetBJetTaggingIP::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

    // Create user output.
    fTrackCountingTagger->InitTrackSelectionParams(NULL);

    // ------- setup PIDCombined
    fPIDCombined = new AliPIDCombined;
    fPIDCombined->SetDefaultTPCPriors();
    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC + AliPIDResponse::kDetTOF);

    if(fUseAtlasSignTCCalculation) {
	Printf("Use ATLAS sign!!!");
	fTrackCountingTagger->SetUseSignDefinitionAtlas(kTRUE);
    }

    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    OpenFile(1);
    if(!fOutput) {
	fOutput = new TList();
	fOutput->SetOwner();
    }
    if(!fOutput) {
	AliError("No container!");
	return;
    }

    const char* tc_track[3] = { "n_1", "n_2", "n_3" };
    const char* tc_flavour[4] = { "inclusive", "lfg", "beauty", "charm" };
    const char* tc_classes[5] = { "default", "verygood", "good", "medium", "bad" };

    TList* glist = new TList();
    glist->SetName("module_EventSelection");
    glist->SetOwner(kTRUE);
    fOutput->Add(glist);

    // Event selection histograms
    AddHistTH1(&fhist_Events, "fhist_Events", "Events", "", "Events", 9, 0., 9., kTRUE, glist);
    fhist_Events->GetXaxis()->SetBinLabel(1, "AliVEvent::kAny");
    fhist_Events->GetXaxis()->SetBinLabel(2, "AliVEvent::kMB");
    fhist_Events->GetXaxis()->SetBinLabel(3, "!Pileup SPD");
    fhist_Events->GetXaxis()->SetBinLabel(4, "Vertex  z < 10 cm");
    fhist_Events->GetXaxis()->SetBinLabel(5, "No VertexerTracks");
    fhist_Events->GetXaxis()->SetBinLabel(6, "Vertex Z Resoulution");
    fhist_Events->GetXaxis()->SetBinLabel(7, ">1 contributors");
    fhist_Events->GetXaxis()->SetBinLabel(8, ">0 contributors SPD");
    fhist_Events->GetXaxis()->SetBinLabel(9, "NoVertex");

    AddHistTH1(&fhist_Rho, "fhist_Rho", "<#rho>", "<#rho> in GeV/c", "a.u.", 500, 0., 10., kTRUE, glist);

    //  Jet selection statistics histograms
    TList* glistjets = new TList();
    glistjets->SetName("module_Jets");
    glistjets->SetOwner(kTRUE);
    fOutput->Add(glistjets);

    AddHistTH2(&fhist_Jets, "fhist_Jets", "Jets", "Cut", "Jets", 9, 0., 9., 500, 0, 250., kTRUE, glistjets);
    fhist_Jets->GetXaxis()->SetBinLabel(1, "In");
    fhist_Jets->GetXaxis()->SetBinLabel(2, "Pointer");
    fhist_Jets->GetXaxis()->SetBinLabel(3, "Pt");
    fhist_Jets->GetXaxis()->SetBinLabel(4, "Eta");
    fhist_Jets->GetXaxis()->SetBinLabel(5, "Area");
    fhist_Jets->GetXaxis()->SetBinLabel(6, "Matched");
    fhist_Jets->GetXaxis()->SetBinLabel(7, "TaggingN1");
    fhist_Jets->GetXaxis()->SetBinLabel(8, "TaggingN2");
    fhist_Jets->GetXaxis()->SetBinLabel(9, "TaggingN3");

    AddHistTH2(&fhist_Jet_Eta_Phi,
               "fhist_Jet_Eta_Phi",
               "jet rec. #eta-#phi",
               "#eta",
               "#phi",
               201,
               -0.9,
               0.9,
               100,
               0.,
               2. * TMath::Pi(),
               kTRUE,
               glistjets);
    AddHistTH2(&fhist_Jet_Nconst_Pt,
               "fhist_Jet_Nconst_Pt",
               "N constituents vs jet pT",
               "p_T",
               "N constituents",
               500,
               0.,
               250.,
               100,
               0.,
               100.,
               kTRUE,
               glistjets);
    AddHistTH1(&fhist_Jet_Pt, "fhist_Jet_Pt", "jet rec. pT", "pT", "counts", 500, 0., 250., kTRUE, glistjets);

    fhists_SPD_cluster_vs_tracklet_correlation = new TH2D(
        "fhists_SPD_cluster_vs_tracklet_correlation", ";SPD Tracklets;SPD Clusters", 200, 0., 200., 1000, 0., 1000.);
    fOutput->Add(fhists_SPD_cluster_vs_tracklet_correlation);
    fhists_SPD_cluster_vs_tracklet_correlation_PostSelection =
        new TH2D("fhists_SPD_cluster_vs_tracklet_correlation_PostSelection",
                 ";SPD Tracklets;SPD Clusters",
                 200,
                 0.,
                 200.,
                 1000,
                 0.,
                 1000.);
    fOutput->Add(fhists_SPD_cluster_vs_tracklet_correlation_PostSelection);

    if(fVetexingPtRelTest)
	fOutput->Add(AddHistsPtRelTemplates());
    if(fVetexingMassFitTest)
	fOutput->Add(AddHistsVtxMassTemplates());

    if(fDoRandomCones) {
	AddHistTH1(&fhist_Jet_Background_Fluctuation,
	           "fhist_Jet_Background_Fluctuation",
	           "random cone #delta #it{p}_{T}",
	           "#delta #it{p}_{T}",
	           "count",
	           500,
	           -125.,
	           125.,
	           kTRUE,
	           glistjets);
    }
    TList* glistunfold = new TList();
    glistunfold->SetName("module_Unfolding");
    glistunfold->SetOwner(kTRUE);
    fOutput->Add(glistunfold);
    const char* uf_flavour[4] = { "inclusive", "lfg", "beauty", "charm" };
    const char* uf_corr[2] = { "raw", "corrected" };
    int uplimit = 0;
    fIsMC ? uplimit = 4 : uplimit = 1;
    for(int itcflavour = 0; itcflavour < 4; itcflavour++)
	for(int itcorr = 0; itcorr < 2; itcorr++) {
	    AddHistTH2(&(fhist_momentum_response[itcflavour][itcorr]),
	               Form("fhist_momentum_response_%s_%s", uf_flavour[itcflavour], uf_corr[itcorr]),
	               Form("Momentum response matrix (%s,%s)", uf_flavour[itcflavour], uf_corr[itcorr]),
	               "#it{p}^{Jet,gen}_{T}",
	               "#it{p}^{Jet,rec}_{T}",
	               500,
	               0.,
	               250.,
	               500,
	               0.,
	               250.,
	               kTRUE,
	               glistunfold);
	}
    // Track counting histograms
    TList* glisttc = new TList();
    glisttc->SetName("module_TrackCounting");
    glisttc->SetOwner(kTRUE);
    fOutput->Add(glisttc);

    for(int itctrack = 0; itctrack < 3; itctrack++) {
	for(int itcflavour = 0; itcflavour < 4; itcflavour++) {
	    if(!fIsMC && itcflavour > 0)
		continue;
	    for(int itcclass = 0; itcclass < 5; itcclass++) {

		AddHistTH2(&(fhist_TC_sIP_Pt[itctrack][itcflavour][itcclass]),
		           Form("fhist_TC_%s_%s_%s", tc_track[itctrack], tc_flavour[itcflavour], tc_classes[itcclass]),
		           Form("Track Counting Output %s_%s_%s",
		                tc_track[itctrack],
		                tc_flavour[itcflavour],
		                tc_classes[itcclass]),
		           "Signed impact parameter (cm)",
		           "#it{p}^{Jet,rec}_{T}",
		           1001,
		           -0.5,
		           0.5,
		           500,
		           0.,
		           250,
		           kTRUE,
		           glisttc);
		AddHistTH2(
		    &(fhist_TC_Eta_Phi[itctrack][itcflavour][itcclass]),
		    Form("fhist_TC_Eta_Phi_%s_%s_%s", tc_track[itctrack], tc_flavour[itcflavour], tc_classes[itcclass]),
		    Form("Track Counting #eta-#phi %s_%s_%s",
		         tc_track[itctrack],
		         tc_flavour[itcflavour],
		         tc_classes[itcclass]),
		    "#eta",
		    "#phi",
		    501,
		    -0.9,
		    0.9,
		    500,
		    0.,
		    2. * TMath::Pi(),
		    kTRUE,
		    glisttc);
	    }
	}
    }
    if(fIsTrackQAConstituent) {
	TList* glisttracksjet = new TList();
	glisttracksjet->SetName("module_TrackQAJet");
	glisttracksjet->SetOwner(kTRUE);
	fOutput->Add(glisttracksjet);
	AddHistTH1(&fhist_QualityClasses,
	           "fhist_QualityClasses",
	           "Tracks in quality class ",
	           "Classname",
	           "count",
	           5,
	           0.,
	           5.,
	           kTRUE,
	           glisttracksjet);
	fhist_QualityClasses->GetXaxis()->SetBinLabel(1, "all");
	fhist_QualityClasses->GetXaxis()->SetBinLabel(2, "kQtyVeryGood");
	fhist_QualityClasses->GetXaxis()->SetBinLabel(3, "kQtyGood");
	fhist_QualityClasses->GetXaxis()->SetBinLabel(4, "kQtyMedium");
	fhist_QualityClasses->GetXaxis()->SetBinLabel(5, "kQtyBad");

	if(fIsMC) {
	    const char* t_classes[4] = { "verygood", "good", "medium", "bad" };
	    const char* t_type[2] = { "Primary", "Secondary" };

	    for(int itcclass = 0; itcclass < 4; itcclass++)
		for(int ittype = 0; ittype < 2; ittype++) {
		    AddHistTH2(
		        &(fhist_QualityClasses_sIP[itcclass][ittype]),
		        Form("fhist_QualityClasses_sIP_%s_%s", t_classes[itcclass], t_type[ittype]),
		        Form("Track Probability PDF sIP distributions (%s,%s)", t_classes[itcclass], t_type[ittype]),
		        "Signed impact parameter (cm)",
		        "#it{p}^{Jet,rec}_{T}",
		        1001,
		        -0.5,
		        0.5,
		        500,
		        0.,
		        250.,
		        kTRUE,
		        glisttracksjet);

		    AddHistTH2(&(fhist_QualityClasses_Eta_Phi[itcclass][ittype]),
		               Form("fhist_QualityClasses_Eta_Phi_%s_%s", t_classes[itcclass], t_type[ittype]),
		               Form("Track Probability classes #eta #phi (%s,%s)", t_classes[itcclass], t_type[ittype]),
		               "#eta",
		               "#phi",
		               501,
		               -0.9,
		               0.9,
		               500,
		               0.,
		               2. * TMath::Pi(),
		               kTRUE,
		               glisttracksjet);
		}
	}
    }
    if(fIsTrackQA) {
	TList* glisttracks = new TList();
	glisttracks->SetName("module_TrackQA");
	glisttracks->SetOwner(kTRUE);
	fOutput->Add(glisttracks);
	AddHistTH2(&fhist_Tracks_Eta_Phi,
	           "fhist_Tracks_Eta_Phi",
	           "Accepted Tracks(9||4)",
	           "Eta",
	           "Phi",
	           201,
	           -1.5,
	           1.5,
	           100,
	           0.,
	           2. * TMath::Pi(),
	           kTRUE,
	           glisttracks);
	AddHistTH2(&fhist_Tracks_Eta_Phi_Bit4,
	           "fhist_Tracks_Eta_Phi_Bit4",
	           "Accepted Tracks(4)",
	           "Eta",
	           "Phi",
	           201,
	           -1.5,
	           1.5,
	           100,
	           0.,
	           2. * TMath::Pi(),
	           kTRUE,
	           glisttracks);
	AddHistTH2(&fhist_Tracks_Eta_Phi_Bit9,
	           "fhist_Tracks_Eta_Phi_Bit9",
	           "Accepted Tracks(9)",
	           "Eta",
	           "Phi",
	           201,
	           -1.5,
	           1.5,
	           100,
	           0.,
	           2. * TMath::Pi(),
	           kTRUE,
	           glisttracks);
    }
    if(fIsMC) {
	// Monte carlo jet flavour statistics histograms
	TList* glistmc = new TList();
	glistmc->SetName("module_MonteCarlo");
	glistmc->SetOwner(kTRUE);
	fOutput->Add(glistmc);
	fhist_MonteCarloFlavour =
	    new TH2D("fhist_MonteCarloFlavour", "fhist_MonteCarloFlavour;pT;Flavour", 5, 0., 5., 500, 0., 250.);
	fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(1, "LFg");
	fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(2, "Charm");
	fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(3, "Beauty");
	fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(4, "Untagged");
	fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(5, "Sum");
	AddHistTH1(&fhist_parton_genjet_dR,
	           "fhist_parton_genjet_dR",
	           "#Delta R parton gen. jet ",
	           "#Delta R",
	           "count",
	           100,
	           0.,
	           1.,
	           kTRUE,
	           glistmc);
	AddHistTH2(&fhist_parton_genjet_pT,
	           "fhist_parton_genjet_pT",
	           "parton gen. jet pT",
	           "pT parton",
	           "pT gen. jet ",
	           500,
	           0.,
	           250.,
	           500,
	           0.,
	           250.,
	           kTRUE,
	           glistmc);
	AddHistTH1(&fhist_parton_genjet_Eta,
	           "fhist_parton_genjet_Eta",
	           "#Delta #eta  parton gen. jet ",
	           "#Delta #eta",
	           "count",
	           500,
	           -1.,
	           1.,
	           kTRUE,
	           glistmc);
	AddHistTH1(&fhist_parton_genjet_Phi,
	           "fhist_parton_genjet_Phi",
	           "#Delta #phi parton gen. jet ",
	           "#Delta #phi",
	           "count",
	           500,
	           -2.,
	           2.,
	           kTRUE,
	           glistmc);
	AddHistTH2(&fhist_MonteCarloFlavour,
	           "fhist_MonteCarloFlavour",
	           "MC flavour",
	           "MC flavour Tag",
	           "#it{p}^{gen}_{T,jet} (GeV/#it{c})",
	           5,
	           0,
	           5.,
	           500,
	           0.,
	           250,
	           kTRUE,
	           glistmc);
    }
    // SetContainer
    fJetsCont = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    fJetsContMC = static_cast<AliJetContainer*>(fJetCollArray.At(1));
    if(fJetsCont)
	fTracksCont = fJetsCont->GetParticleContainer();
    else
	fTracksCont = GetParticleContainer(0);
    if(fTracksCont)
	fTracksCont->SetClassName("AliVTrack");
    fTrackCountingTagger->SetParticleContainer(fTracksCont);
    fTrackCountingTagger->SetAnalysisTypeAOD(kTRUE);
    TH1::AddDirectory(oldStatus);
    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
    return;
}
// Helper functions
void AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistTH1(TH1** hist,
                                                      const char* histname,
                                                      const char* title,
                                                      const char* titlex,
                                                      const char* titley,
                                                      Int_t nBinsX,
                                                      Double_t minX,
                                                      Double_t maxX,
                                                      Bool_t setSumw2,
                                                      TList* container)
{
    if(!container)
	return;
    *hist = new TH1D(histname, Form("%s;%s;%s", title, titlex, titley), nBinsX, minX, maxX);
    if(!(*hist))
	return;
    if(setSumw2)
	(*hist)->Sumw2();
    if(container->FindObject(histname)) {
	AliError(Form("Object with name  %s already exists in %s...returning!", histname, container->GetName()));
	return;
    }
    container->Add(*hist);
    return;
}
void AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistTH2(TH2** hist,
                                                      const char* histname,
                                                      const char* title,
                                                      const char* titlex,
                                                      const char* titley,
                                                      Int_t nBinsX,
                                                      Double_t minX,
                                                      Double_t maxX,
                                                      Int_t nBinsY,
                                                      Double_t minY,
                                                      Double_t maxY,
                                                      Bool_t setSumw2,
                                                      TList* container)
{
    if(!container)
	return;

    *hist = new TH2D(histname, Form("%s;%s;%s", title, titlex, titley), nBinsX, minX, maxX, nBinsY, minY, maxY);
    if(!(*hist))
	return;
    if(setSumw2)
	(*hist)->Sumw2();
    if(container->FindObject(histname)) {
	AliError(Form("Object with name  %s already exists in %s...returning!", histname, container->GetName()));
	return;
    }
    container->Add(*hist);
    return;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::GetDeltaPtRandomCone
/// \return
///
Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetDeltaPtRandomCone()
{
    Double_t deltaPt = -1000.;
    Double_t jetradius = fJetsCont->GetJetRadius();
    Double_t minEta = -0.9 + jetradius;
    Double_t maxEta = 0.9 - jetradius;
    Double_t tmpRandConeEta = minEta + fRandom->Rndm() * (maxEta - minEta);
    Double_t tmpRandConePhi = fRandom->Rndm() * TMath::TwoPi();
    Double_t tmpConePt = -1.;

    for(Int_t i = 0; i < static_cast<AliAODEvent*>(InputEvent())->GetNumberOfTracks(); i++) {
	AliAODTrack* tmpTrack = static_cast<AliAODTrack*>(InputEvent()->GetTrack(i));
	if(fabs(tmpTrack->Eta()) < 0.9) {
	    if(tmpTrack->Pt() > 0.15) {
		if(sqrt((tmpTrack->Eta() - tmpRandConeEta) * (tmpTrack->Eta() - tmpRandConeEta) +
		        TVector2::Phi_mpi_pi((tmpTrack->Phi() - tmpRandConePhi)) *
		            TVector2::Phi_mpi_pi((tmpTrack->Phi() - tmpRandConePhi))) < jetradius) {
		    tmpConePt += tmpTrack->Pt();
		}
	    }
	}
    }
    if(tmpConePt > 0) {
	deltaPt = tmpConePt - 0.4 * 0.4 * TMath::Pi() * fJetsCont->GetRhoVal();
	return deltaPt;
    }
    return deltaPt;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::GetPtCorrected
/// \param jet
/// \return
///
Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetPtCorrected(const AliEmcalJet* jet)
{
    if(jet && fJetsCont)
	return jet->Pt() - fJetsCont->GetRhoVal() * jet->Area();
    return -1.;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::FindVertexNProngSimple
/// \param jet
/// \param vtx
/// \param nProng
/// \return
///
Bool_t
AliAnalysisTaskEmcalJetBJetTaggingIP::FindVertexNProngSimple(const AliEmcalJet* jet, AliAODVertex*& vtx, Int_t& nProng)
{
    // populate array from tracks within jet
    AliVTrack* bTrack = 0x0;
    TList* lTrackInJet = new TList();
    lTrackInJet->SetOwner(0);
    for(Int_t j = 0; j < jet->GetNumberOfTracks(); ++j) {
	bTrack = (AliVTrack*)((AliPicoTrack*)fTracksCont->GetParticle(jet->TrackAt((int)j)))->GetTrack();
	lTrackInJet->Add(bTrack);
    }
    Int_t nprongs = lTrackInJet->GetEntries();
    if(nprongs < 1)
	return kFALSE;
    // Construct Kalmann vertex

    AliKFParticle::SetField(InputEvent()->GetMagneticField());
    AliKFVertex vertexKF;

    for(Int_t i = 0; i < nprongs; i++) {
	AliVTrack* esdTrack = (AliVTrack*)lTrackInJet->At(i);
	AliKFParticle daughterKF(*esdTrack, 211);
	vertexKF.AddDaughter(daughterKF);
    }
    // Printf("nCont %i",vertexKF.GetNContributors());
    AliESDVertex* vertexESD = new AliESDVertex(
        vertexKF.Parameters(), vertexKF.CovarianceMatrix(), vertexKF.GetChi2(), vertexKF.GetNContributors());

    Double_t pos[3], cov[6], chi2perNDF;
    vertexESD->GetXYZ(pos);       // position
    vertexESD->GetCovMatrix(cov); // covariance matrix
    chi2perNDF = vertexESD->GetChi2toNDF();
    Double_t dispersion = vertexESD->GetDispersion();

    AliAODVertex* vertexAOD = new AliAODVertex(pos, cov, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, nprongs);

    for(Int_t i = 0; i < nprongs; i++) {
	vertexAOD->AddDaughter((AliAODTrack*)lTrackInJet->At(i));
    }
    delete vertexESD;
    vertexESD = NULL;
    vtx = vertexAOD;
    nProng = nprongs;

    delete lTrackInJet;

    return kTRUE;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::GetVertexInvariantMass
/// \param vtx
/// \param massParticle
/// \return
///
Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetVertexInvariantMass(AliAODVertex* vtx, Double_t massParticle)
{
    Double_t pxyz[3];
    Double_t pxyzSum[4] = { 0., 0., 0., 0. };
    if(!vtx)
	return -1.;
    for(Int_t jp = 0; jp < vtx->GetNDaughters(); jp++) {
	AliAODTrack* tr = NULL;

	tr = (AliAODTrack*)vtx->GetDaughter(jp);

	if(!tr) {
	    Printf("Error getting Track (GetVertexInvariantMass)");
	    continue;
	}
	tr->GetPxPyPz(pxyz);
	pxyzSum[1] += pxyz[0];
	pxyzSum[2] += pxyz[1];
	pxyzSum[3] += pxyz[2];
	pxyzSum[0] +=
	    TMath::Sqrt(massParticle * massParticle + pxyz[0] * pxyz[0] + pxyz[1] * pxyz[1] + pxyz[2] * pxyz[2]);
    }
    double mass = TMath::Sqrt(pxyzSum[0] * pxyzSum[0] - pxyzSum[1] * pxyzSum[1] - pxyzSum[2] * pxyzSum[2] -
                              pxyzSum[3] * pxyzSum[3]);
    return mass;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::FillVertexingHists
/// \param jet
/// \param vtx
/// \param nProng
/// \param tag
/// \param hist
///
void AliAnalysisTaskEmcalJetBJetTaggingIP::FillVertexingHists(const AliEmcalJet* jet,
                                                              AliAODVertex* vtx,
                                                              Int_t nProng,
                                                              Int_t tag,
                                                              Int_t hist)
{

    //    if(!vtx) return;
    //    Double_t invMass = GetVertexInvariantMass(vtx,0.13957018); //pion mass assumed
    //    if(invMass < 0 ) return;
    //    if(nProng<2) return;
    //    if (nProng>1)fHist_MassVsJetPtHE[tag][hist]->Fill(invMass,jet->Pt());
    //    if (nProng>2)fHist_MassVsJetPtHP[tag][hist]->Fill(invMass,jet->Pt());
    return;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief GetElectronPIDnSigmaTPC
/// \param track
/// \return
///
Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetElectronPIDnSigmaTPC(const AliVTrack* track)
{
    Double_t eleLineDist = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    return eleLineDist;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::GetPionPIDnSigmaTPC
/// \param track
/// \return
///
Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetPionPIDnSigmaTPC(const AliVTrack* track)
{
    Double_t eleLineDist = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
    return eleLineDist;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::JetTrackLoop
/// \param jet
///
void AliAnalysisTaskEmcalJetBJetTaggingIP::JetTrackLoop(const AliEmcalJet* jet, Double_t n3tag, Int_t MCtag)
{

    return;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::GetPtRel
/// \param track
/// \param jet
/// \return
///
Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetPtRel(const AliVTrack* track, const AliEmcalJet* jet)
{
    TVector3 p_vjet(jet->Px(), jet->Py(), jet->Pz());
    TVector3 p_vtrack(track->Px(), track->Py(), track->Pz());
    Double_t theta = p_vjet.Angle(p_vtrack);
    Double_t pTrel = TMath::Sin(theta) * track->Pt();
    return pTrel;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::IsElectronTPC
/// \param nSigmaTPCelectron
/// \param nSigmaTPCpion
/// \return
///
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsElectronTPC(double nSigmaTPCelectron, double nSigmaTOFElectron)
{
    if(nSigmaTPCelectron > fSigmaTPCElectronLow && nSigmaTPCelectron < fSigmaTPCElectronHigh) {
	if(nSigmaTOFElectron > fSigmaTOFElectronLow && nSigmaTOFElectron < fSigmaTOFElectronHigh) {
	    return kTRUE;
	}
    }
    return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistsPtRelTemplates
/// \return
///
TList* AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistsPtRelTemplates()
{
    TList* list_SoftElectron = new TList();
    list_SoftElectron->SetName("module_SoftElectronPtRel");
    list_SoftElectron->SetOwner(kTRUE);
    const char* thesholds[4] = { "ALL", "0", "0.0028", "0.01" };
    const char* tc_flavour[4] = { "inclusive", "lfg", "beauty", "charm" };
    for(int i = 0; i < 4; ++i) {
	fHistpTrelElectron[i] = new TH2D(
	    Form("fHistpTrelElectron_Tagged_%s", thesholds[i]), "jets pTrel;pTrel (GeV)", 120, 0, 3., 400, 0, 200.);
	list_SoftElectron->Add(fHistpTrelElectron[i]);
	for(int j = 0; j < 4; ++j) {
	    fHistpTrelElectronMC[i][j] =
	        new TH2D(Form("fHistpTrelElectronMC_True_Flavour_%s_%s", tc_flavour[j], thesholds[i]),
	                 "jets pTrel;pTrel (GeV)",
	                 120,
	                 0,
	                 3.,
	                 400,
	                 0,
	                 200.);
	    list_SoftElectron->Add(fHistpTrelElectronMC[i][j]);
	}
    }
    fHistTPCnSigmaElectron =
        new TH2D("fHistTPCnSigmaElectron", "fHistTPCnSigmaElectron;pT;nSigma TPC", 500, 0, 250, 500, -10, 10);
    fHistTPCnSigmaPion =
        new TH2D("fHistTPCnSigmaElectron", "fHistTPCnSigmaElectron;pT;nSigma TPC", 500, 0, 250, 500, -10, 10);

    list_SoftElectron->Add(fHistTPCnSigmaElectron);
    list_SoftElectron->Add(fHistTPCnSigmaPion);
    return list_SoftElectron;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistsVtxMassTemplates
/// \return
///
TList* AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistsVtxMassTemplates()
{
    TList* list_VtxMass = new TList();
    list_VtxMass->SetName("module_VertexMass");
    list_VtxMass->SetOwner(kTRUE);
    const char* thesholds[4] = { "ALL", "0", "0.0028", "0.01" };
    const char* tc_flavour[4] = { "inclusive", "lfg", "beauty", "charm" };
    for(int i = 0; i < 4; ++i) {
	fHist_MassVsJetPtHE[i] = new TH2D(
	    Form("fHist_MassVsJetPtHE_%s", thesholds[i]), "inv mass HE;inv. Mass (GeV/c^2)", 100, 0, 20., 400, 0, 200.);
	fHist_MassVsJetPtHP[i] = new TH2D(
	    Form("fHist_MassVsJetPtHP_%s", thesholds[i]), "inv mass HE;inv. Mass (GeV/c^2)", 100, 0, 20., 400, 0, 200.);
	list_VtxMass->Add(fHist_MassVsJetPtHE[i]);
	list_VtxMass->Add(fHist_MassVsJetPtHP[i]);
	for(int j = 0; j < 4; ++j) {
	    fHist_MassVsJetPtHEMC[i][j] =
	        new TH2D(Form("fHist_MassVsJetPtHEMC_True_Flavour_%s_%s", tc_flavour[j], thesholds[i]),
	                 "inv mass HE;inv. Mass (GeV/c^2)",
	                 100,
	                 0,
	                 20.,
	                 400,
	                 0,
	                 200.);
	    fHist_MassVsJetPtHPMC[i][j] =
	        new TH2D(Form("fHist_MassVsJetPtHEMC_True_Flavour_%s_%s", tc_flavour[j], thesholds[i]),
	                 "inv mass HE;inv. Mass (GeV/c^2)",
	                 100,
	                 0,
	                 20.,
	                 400,
	                 0,
	                 200.);
	    list_VtxMass->Add(fHist_MassVsJetPtHEMC[i][j]);
	    list_VtxMass->Add(fHist_MassVsJetPtHPMC[i][j]);
	}
    }
    return list_VtxMass;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::ProcessPtRelTemplateAnalysis
/// \param jet
/// \param mcflavour
/// \param n3value
///
void AliAnalysisTaskEmcalJetBJetTaggingIP::ProcessPtRelTemplateAnalysis(const AliEmcalJet* jet,
                                                                        double corrected_jet_pt,
                                                                        int mcflavour,
                                                                        double n3value)
{
    Bool_t tagged = kFALSE;
    AliVTrack* track = NULL;
    if(n3value > -999.)
	tagged = kTRUE;
    if(!jet)
	return;
    for(Int_t j = 0; j < jet->GetNumberOfTracks(); ++j) {
	track = (AliVTrack*)((AliPicoTrack*)fTracksCont->GetParticle(jet->TrackAt((int)j)))->GetTrack();
	double signal = GetElectronPIDnSigmaTPC(track);
	double pTtrack = track->Pt();
	if(pTtrack > fpTRelSoftElectronPt && jet->Pt() > fpTRelSoftElectronJetPt &&
	   IsElectronTPC(signal, fPIDResponse->NumberOfSigmasTOF(track, AliPID::kElectron))) {
	    fHistTPCnSigmaElectron->Fill(pTtrack, signal);
	    // Calclulate pT rel for identified electrons
	    Double_t pT_rel = GetPtRel(track, jet);
	    // fill inclusive (mens no tagging ) histograms
	    fHistpTrelElectron[0]->Fill(pT_rel, corrected_jet_pt);
	    if(fIsMC) {
		fHistpTrelElectronMC[0][mcflavour]->Fill(pT_rel, corrected_jet_pt);
		fHistpTrelElectronMC[0][0]->Fill(pT_rel, corrected_jet_pt);
	    }
	    if(tagged) {
		if(n3value >= 0.0) {
		    fHistpTrelElectron[1]->Fill(pT_rel, corrected_jet_pt);
		    if(fIsMC) {
			fHistpTrelElectronMC[1][mcflavour]->Fill(pT_rel, corrected_jet_pt);
			fHistpTrelElectronMC[1][0]->Fill(pT_rel, corrected_jet_pt);
		    }
		}
		if(n3value >= 0.0028) {
		    fHistpTrelElectron[2]->Fill(pT_rel, corrected_jet_pt);
		    if(fIsMC) {
			fHistpTrelElectronMC[2][mcflavour]->Fill(pT_rel, corrected_jet_pt);
			fHistpTrelElectronMC[2][0]->Fill(pT_rel, corrected_jet_pt);
		    }
		}
		if(n3value >= 0.01) {
		    fHistpTrelElectron[3]->Fill(pT_rel, corrected_jet_pt);
		    if(fIsMC) {
			fHistpTrelElectronMC[3][mcflavour]->Fill(pT_rel, corrected_jet_pt);
			fHistpTrelElectronMC[3][0]->Fill(pT_rel, corrected_jet_pt);
		    }
		}
	    }
	}
    }
    return;
}
////////////////////////////////////////////////////////////////////////////////
/// \brief AliAnalysisTaskEmcalJetBJetTaggingIP::ProcessVtxMassTemplateAnalysis
/// \param jet
/// \param corrected_jet_pt
/// \param mcflavour
/// \param n3value
///
void AliAnalysisTaskEmcalJetBJetTaggingIP::ProcessVtxMassTemplateAnalysis(const AliEmcalJet* jet,
                                                                          double corrected_jet_pt,
                                                                          int mcflavour,
                                                                          double n3value)
{
    Bool_t tagged = kFALSE;
    AliVTrack* track = NULL;
    if(n3value > -999.)
	tagged = kTRUE;
    if(!jet)
	return;
    AliAODVertex* vertex = NULL;
    Int_t nVertexProng = 0;
    FindVertexNProngSimple(jet, vertex, nVertexProng);
    if(!vertex)
	return;
    if(nVertexProng < 2)
	return;
    Double_t invMass = GetVertexInvariantMass(vertex, 0.13957018); // pion mass assumed
    if(invMass < 0)
	return;
    fHist_MassVsJetPtHE[0]->Fill(invMass, corrected_jet_pt);
    if(nVertexProng > 2)
	fHist_MassVsJetPtHP[0]->Fill(invMass, corrected_jet_pt);

    if(fIsMC) {
	fHist_MassVsJetPtHEMC[0][mcflavour]->Fill(invMass, corrected_jet_pt);
	fHist_MassVsJetPtHEMC[0][0]->Fill(invMass, corrected_jet_pt);
	if(nVertexProng > 2) {
	    fHist_MassVsJetPtHPMC[0][mcflavour]->Fill(invMass, corrected_jet_pt);
	    fHist_MassVsJetPtHPMC[0][0]->Fill(invMass, corrected_jet_pt);
	}
    }
    if(tagged) {
	if(n3value >= 0.0) {
	    fHist_MassVsJetPtHE[1]->Fill(invMass, corrected_jet_pt);
	    if(nVertexProng > 2)
		fHist_MassVsJetPtHP[1]->Fill(invMass, corrected_jet_pt);
	    if(fIsMC) {
		fHist_MassVsJetPtHEMC[1][mcflavour]->Fill(invMass, corrected_jet_pt);
		if(nVertexProng > 2)
		    fHist_MassVsJetPtHPMC[1][mcflavour]->Fill(invMass, corrected_jet_pt);
		fHist_MassVsJetPtHEMC[1][0]->Fill(invMass, corrected_jet_pt);
		if(nVertexProng > 2)
		    fHist_MassVsJetPtHPMC[1][0]->Fill(invMass, corrected_jet_pt);
	    }
	}
	if(n3value >= 0.0028) {
	    fHist_MassVsJetPtHE[2]->Fill(invMass, corrected_jet_pt);
	    if(nVertexProng > 2)
		fHist_MassVsJetPtHP[2]->Fill(invMass, corrected_jet_pt);
	    if(fIsMC) {
		fHist_MassVsJetPtHEMC[2][mcflavour]->Fill(invMass, corrected_jet_pt);
		if(nVertexProng > 2)
		    fHist_MassVsJetPtHPMC[2][mcflavour]->Fill(invMass, corrected_jet_pt);
		fHist_MassVsJetPtHEMC[2][0]->Fill(invMass, corrected_jet_pt);
		if(nVertexProng > 2)
		    fHist_MassVsJetPtHPMC[2][0]->Fill(invMass, corrected_jet_pt);
	    }
	}
	if(n3value >= 0.01) {
	    fHist_MassVsJetPtHE[3]->Fill(invMass, corrected_jet_pt);
	    if(nVertexProng > 2)
		fHist_MassVsJetPtHP[3]->Fill(invMass, corrected_jet_pt);
	    if(fIsMC) {
		fHist_MassVsJetPtHEMC[3][mcflavour]->Fill(invMass, corrected_jet_pt);
		if(nVertexProng > 2)
		    fHist_MassVsJetPtHPMC[3][mcflavour]->Fill(invMass, corrected_jet_pt);
		fHist_MassVsJetPtHEMC[3][0]->Fill(invMass, corrected_jet_pt);
		if(nVertexProng > 2)
		    fHist_MassVsJetPtHPMC[3][0]->Fill(invMass, corrected_jet_pt);
	    }
	}
    }
    return;
}
