#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TVector2.h>

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
#include "AliRDHFJetsCuts.h"
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

//Changelog  MC Jet Container retrieval fixed
ClassImp(AliAnalysisTaskEmcalJetBJetTaggingIP)

//________________________________________________________________________
AliAnalysisTaskEmcalJetBJetTaggingIP::AliAnalysisTaskEmcalJetBJetTaggingIP() : AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetBJetTaggingIP", kTRUE),
    fJetsCont(NULL),
    fJetsContMC(NULL),
    fTracksCont(NULL),
    fTaggingHFClass( new   AliHFJetsTagging("fTaggingHFClass")),
    fTrackCountingTagger(new AliHFJetTaggingIP()),
    fRandom(new TRandom3(0)),
    fUtils(new AliAnalysisUtils()),
    fJetCutsHF(new AliRDHFJetsCuts()),
    fMCparticles(NULL),
    fIsMC(kFALSE),
    fIsTrackQA(kFALSE),
    fIsTrackQAConstituent(kFALSE),
    fUseCorrectedJetPt(kFALSE),
    fDoRandomCones(kFALSE),
    fUseEventSelection(0),
    fUseJetSelection(0),
    fUseMCTagger(0),
    fhist_Events(NULL),
    fhist_Jets(NULL),
    fhists_SPD_cluster_vs_tracklet_correlation(NULL),
    fhists_SPD_cluster_vs_tracklet_correlation_PostSelection(NULL),
    fhist_MonteCarloFlavour(NULL),
    fhist_Tracks_Eta_Phi(NULL),
    fhist_Tracks_Eta_Phi_Bit4(NULL),
    fhist_Tracks_Eta_Phi_Bit9(NULL),
    fhist_QualityClasses(NULL),
    fhist_Jet_Eta_Phi(NULL),
    fhist_Jet_Nconst_Pt(NULL),
    fhist_Jet_Pt(NULL),
    fhist_Jet_Background_Fluctuation(NULL),
    fhist_Rho(NULL),
    fhist_parton_genjet_dR(NULL),
    fhist_parton_genjet_pT(NULL),
    fhist_parton_genjet_Eta(NULL),
    fhist_parton_genjet_Phi(NULL)
{
    memset(fhist_TC_sIP_Pt,0,sizeof fhist_TC_sIP_Pt);
    memset(fhist_TC_Eta_Phi,0,sizeof fhist_TC_sIP_Pt);
    memset(fhist_QualityClasses_sIP,0,sizeof fhist_QualityClasses_sIP);
    memset(fhist_QualityClasses_Eta_Phi,0,sizeof fhist_QualityClasses_Eta_Phi);
    memset(fhist_momentum_response,0,sizeof fhist_momentum_response);
    SetMakeGeneralHistograms(kTRUE);



}

//________________________________________________________________________
AliAnalysisTaskEmcalJetBJetTaggingIP::AliAnalysisTaskEmcalJetBJetTaggingIP(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fJetsCont(NULL),
    fJetsContMC(NULL),
    fTracksCont(NULL),
    fTaggingHFClass(new AliHFJetsTagging("fTaggingHFClass")),
    fTrackCountingTagger(new AliHFJetTaggingIP()),
    fRandom(new TRandom3(0)),
    fUtils(new AliAnalysisUtils()),
    fJetCutsHF(new AliRDHFJetsCuts()),
    fMCparticles(NULL),
    fIsMC(kFALSE),
    fIsTrackQA(kFALSE),
    fIsTrackQAConstituent(kFALSE),
    fUseCorrectedJetPt(kFALSE),
    fDoRandomCones(kFALSE),
    fUseEventSelection(0),
    fUseJetSelection(0),
    fUseMCTagger(0),
    fhist_Events(NULL),
    fhist_Jets(NULL),
    fhists_SPD_cluster_vs_tracklet_correlation(NULL),
    fhists_SPD_cluster_vs_tracklet_correlation_PostSelection(NULL),
    fhist_MonteCarloFlavour(NULL),
    fhist_Tracks_Eta_Phi(NULL),
    fhist_Tracks_Eta_Phi_Bit4(NULL),
    fhist_Tracks_Eta_Phi_Bit9(NULL),
    fhist_QualityClasses(NULL),
    fhist_Jet_Eta_Phi(NULL),
    fhist_Jet_Nconst_Pt(NULL),
    fhist_Jet_Pt(NULL),
    fhist_Jet_Background_Fluctuation(NULL),
    fhist_Rho(NULL),
    fhist_parton_genjet_dR(NULL),
    fhist_parton_genjet_pT(NULL),
    fhist_parton_genjet_Eta(NULL),
    fhist_parton_genjet_Phi(NULL)
{
    // Standard constructor.
    memset(fhist_TC_sIP_Pt,0,sizeof fhist_TC_sIP_Pt);
    memset(fhist_TC_Eta_Phi,0,sizeof fhist_TC_sIP_Pt);
    memset(fhist_QualityClasses_sIP,0,sizeof fhist_QualityClasses_sIP);
    memset(fhist_QualityClasses_Eta_Phi,0,sizeof fhist_QualityClasses_Eta_Phi);
    memset(fhist_momentum_response,0,sizeof fhist_momentum_response);
    SetMakeGeneralHistograms(kTRUE);
}


Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::Run()
{
    //Event selection called automatically

    if(fIsMC)  fMCparticles = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));


    AliEmcalJet *curjet = NULL;
    AliEmcalJet *curjetMatched = NULL;
    Double_t  tagvalue[3] ={0.,0.,0.};
    Bool_t    istagged[3] ={kFALSE,kFALSE,kFALSE};
    Int_t flavourtag = 0;
    fTrackCountingTagger->SetEvent(InputEvent());
    if(!fJetsCont){
        AliError("Missing jet container!");
        return kFALSE;
    }

    if(fDoRandomCones){
        Double_t randomConePt = GetDeltaPtRandomCone();
        fhist_Jet_Background_Fluctuation->Fill(randomConePt);
    }

    fhist_Rho->Fill(fJetsCont->GetRhoVal());

    //Loop over all available jets
    for(Int_t ijet = 0 ; ijet <fJetsCont->GetNJets();++ijet){
        memset(tagvalue,0,sizeof tagvalue);
        memset(istagged,0,sizeof istagged);
        curjet = NULL;
        curjetMatched =NULL;
        flavourtag = 0;


        curjet =  fJetsCont->GetJet(ijet);

        if(!IsJetSelected(curjet)) continue;
        Double_t jetPt = 0.;
        fUseCorrectedJetPt ? jetPt = GetPtCorrected(curjet)  : jetPt = curjet->Pt();
        if(fIsMC) {
            curjetMatched = curjet->MatchedJet();
            if(curjetMatched) {
                fhist_Jets->Fill("Matched",jetPt,1.);
                AddTagJet(curjetMatched);
                if(curjetMatched->TestFlavourTag(kLFgJet)) flavourtag =1;
                else if(curjetMatched->TestFlavourTag(kBeautyJet)) flavourtag=2;
                else if(curjetMatched->TestFlavourTag(kCharmJet)) flavourtag=3;

                Double_t ptm = curjetMatched->Pt();
                Double_t ptmcorr =ptm - fJetsContMC->GetRhoVal()*curjetMatched->Area();
                Double_t ptr = curjet->Pt();
                Double_t ptrcorr = GetPtCorrected(curjet);
                fhist_momentum_response[0][0]->Fill(ptm,ptr);
                fhist_momentum_response[0][1]->Fill(ptmcorr,ptrcorr);
                if(flavourtag>0){
                    fhist_momentum_response[flavourtag][0]->Fill(ptm,ptr);
                    fhist_momentum_response[flavourtag][1]->Fill(ptmcorr,ptrcorr);
                }
            }



        }


        //Standard track counting algorithm

        fhist_Jet_Pt->Fill(jetPt);
        fhist_Jet_Eta_Phi->Fill(curjet->Eta(),curjet->Phi());
        fhist_Jet_Nconst_Pt->Fill(jetPt,curjet->GetNumberOfTracks());

        if(fTrackCountingTagger->GetJetDiscriminator(curjet,tagvalue,istagged)){
            if(istagged[0]){
                fhist_Jets->Fill("TaggingN1",jetPt,1.);
                fhist_TC_sIP_Pt[0][0][0]->Fill(tagvalue[0],jetPt);
                fhist_TC_Eta_Phi[0][0][0]->Fill(curjet->Eta(),curjet->Phi());
                if(fIsMC)
                    if (flavourtag>0) {
                        fhist_TC_sIP_Pt[0][flavourtag][0]->Fill(tagvalue[0],jetPt);
                        fhist_TC_Eta_Phi[0][flavourtag][0]->Fill(curjet->Eta(),curjet->Phi());
                    }
            }
            if(istagged[1]){
                fhist_Jets->Fill("TaggingN2",jetPt,1.);
                fhist_TC_sIP_Pt[1][0][0]->Fill(tagvalue[1],jetPt);
                fhist_TC_Eta_Phi[1][0][0]->Fill(curjet->Eta(),curjet->Phi());
                if(fIsMC)
                    if (flavourtag>0){
                        fhist_TC_Eta_Phi[1][flavourtag][0]->Fill(curjet->Eta(),curjet->Phi());
                        fhist_TC_sIP_Pt[1][flavourtag][0]->Fill(tagvalue[1],jetPt);
                    }
            }
            if(istagged[2]){
                fhist_Jets->Fill("TaggingN3",jetPt,1.);
                fhist_TC_sIP_Pt[2][0][0]->Fill(tagvalue[2],jetPt);
                fhist_TC_Eta_Phi[2][0][0]->Fill(curjet->Eta(),curjet->Phi());
                if(fIsMC)
                    if (flavourtag>0){
                        fhist_TC_sIP_Pt[2][flavourtag][0]->Fill(tagvalue[2],jetPt);
                        fhist_TC_Eta_Phi[2][flavourtag][0]->Fill(curjet->Eta(),curjet->Phi());
                    }
            }
        }

        //Track quality track counting
        for(int i=1;i<5;++i){
            memset(tagvalue,0,sizeof tagvalue);
            memset(istagged,0,sizeof istagged);

            if(fTrackCountingTagger->GetJetDiscriminatorQualityClass(i,curjet,tagvalue,istagged)){
                if(istagged[0]){
                    fhist_TC_sIP_Pt[0][0][i]->Fill(tagvalue[0],jetPt);
                    fhist_TC_Eta_Phi[0][0][i]->Fill(curjet->Eta(),curjet->Phi());
                    if(fIsMC)
                        if (flavourtag>0) {
                            fhist_TC_sIP_Pt[0][flavourtag][i]->Fill(tagvalue[0],jetPt);
                            fhist_TC_Eta_Phi[0][flavourtag][i]->Fill(curjet->Eta(),curjet->Phi());
                        }
                }
                if(istagged[1]){
                    fhist_TC_sIP_Pt[1][0][i]->Fill(tagvalue[1],jetPt);
                    fhist_TC_Eta_Phi[1][0][i]->Fill(curjet->Eta(),curjet->Phi());

                    if(fIsMC)
                        if (flavourtag>0){
                            fhist_TC_sIP_Pt[1][flavourtag][i]->Fill(tagvalue[1],jetPt);
                            fhist_TC_Eta_Phi[0][flavourtag][i]->Fill(curjet->Eta(),curjet->Phi());
                        }
                }
                if(istagged[2]){
                    fhist_TC_sIP_Pt[2][0][i]->Fill(tagvalue[2],jetPt);
                    fhist_TC_Eta_Phi[2][0][i]->Fill(curjet->Eta(),curjet->Phi());
                    if(fIsMC)
                        if (flavourtag>0) {
                            fhist_TC_Eta_Phi[2][flavourtag][i]->Fill(curjet->Eta(),curjet->Phi());
                            fhist_TC_sIP_Pt[2][flavourtag][i]->Fill(tagvalue[2],jetPt);
                        }
                }
            }
        }

        //Constituent QA
        if(fIsTrackQAConstituent){
            RunQATracksJet(curjet);
        }

    }
    //Run Event Track QA
    if(fIsTrackQA){
        RunQATracksEvent();
    }


    return kTRUE;
}
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::RunQATracksJet(const AliEmcalJet * jet)
{
    Double_t dv[2]={0.,0.};
    Double_t covv[3]={0.,0.};
    Double_t sip = 0.;
    AliAODTrack * track = 0x0;
    AliAODMCParticle * part =0x0;

    if(!jet) return kFALSE;
    Int_t  ntracks = (Int_t)jet->GetNumberOfTracks();

    for (Int_t itrack = 0 ; itrack< ntracks; ++itrack){

        memset(dv,0,sizeof dv);
        memset(covv,0,sizeof covv);
        sip = 0.;

        track= (AliAODTrack*) ((( AliPicoTrack*)(fTracksCont->GetParticle(jet->TrackAt(itrack))))->GetTrack());
        if(!track) continue;
        if(fabs(track->Eta()) > 0.9 ) continue;
        if(track->Pt() < 1) continue;

        Bool_t isPrimary=kFALSE;
        if (fIsMC && fMCparticles){
            part = dynamic_cast<AliAODMCParticle*>(fMCparticles->At(TMath::Abs(track->GetLabel())));
            if(part){
                isPrimary = part->IsPhysicalPrimary();
            }
            fTaggingHFClass->GetSignedRPhiImpactParameter(InputEvent(),track,jet,sip,dv,covv);
        }


        fhist_QualityClasses->Fill("all",1.) ;
        if(IsQuality(track, kQtyVeryGood)){
            fhist_QualityClasses->Fill("kQtyVeryGood",1.);
            if(fIsMC){
                if(isPrimary) {
                    fhist_QualityClasses_sIP[0][0]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[0][0]->Fill(track->Eta(),track->Phi());
                }
                else {
                    fhist_QualityClasses_sIP[0][1]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[0][1]->Fill(track->Eta(),track->Phi());
                }
            }

        }
        if(IsQuality(track, kQtyGood)){
            fhist_QualityClasses->Fill("kQtyGood",1.);
            if(fIsMC){
                if(isPrimary) {
                    fhist_QualityClasses_sIP[1][0]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[1][0]->Fill(track->Eta(),track->Phi());
                }
                else {
                    fhist_QualityClasses_sIP[1][1]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[1][1]->Fill(track->Eta(),track->Phi());
                }
            }
        }
        if(IsQuality(track, kQtyMedium)){
            fhist_QualityClasses->Fill("kQtyMedium",1.);
            if(fIsMC){
                if(isPrimary) {
                    fhist_QualityClasses_sIP[2][0]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[2][0]->Fill(track->Eta(),track->Phi());
                }
                else {
                    fhist_QualityClasses_sIP[2][1]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[2][1]->Fill(track->Eta(),track->Phi());
                }
            }
        }
        if(IsQuality(track, kQtyBad)){
            fhist_QualityClasses->Fill("kQtyBad",1.);
            if(fIsMC){
                if(isPrimary) {
                    fhist_QualityClasses_sIP[3][0]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[3][0]->Fill(track->Eta(),track->Phi());
                }
                else {
                    fhist_QualityClasses_sIP[3][1]->Fill(sip,track->Pt());
                    fhist_QualityClasses_Eta_Phi[3][1]->Fill(track->Eta(),track->Phi());
                }
            }
        }
    }
    return kTRUE;
}


Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::RunQATracksEvent(){

    Int_t  ntracks = (Int_t)InputEvent()->GetNumberOfTracks();
    AliAODTrack * track = 0x0;
    for (Int_t itrack = 0 ; itrack< ntracks; ++itrack){
        track = (AliAODTrack*) InputEvent()->GetTrack(itrack);
        if(fabs(track->Eta()) > 0.9 ) continue;
        if(track->Pt() < 0.150) continue;
        if(!((track)->TestFilterBit(1<<4) ||
             ((track)->TestFilterBit(1<<9)))) continue;
        fhist_Tracks_Eta_Phi->Fill(track->Eta(),track->Phi());
        if((track)->TestFilterBit(1<<4) &&
                !((track)->TestFilterBit(1<<9) ))
            fhist_Tracks_Eta_Phi_Bit4->Fill(track->Eta(),track->Phi());
        if(!((track)->TestFilterBit(1<<4)) &&
                (track)->TestFilterBit(1<<9) )
            fhist_Tracks_Eta_Phi_Bit9->Fill(track->Eta(),track->Phi());
    }

    return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelected(){
    AliAODEvent *aev = dynamic_cast< AliAODEvent*>(InputEvent());
    // SPD Cluster vs Tracklet plot to estimate pileup effect
    Int_t nClustersLayer0 = aev->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = aev->GetNumberOfITSClusters(1);
    Int_t nTracklets      = aev->GetMultiplicity()->GetNumberOfTracklets();
    fhists_SPD_cluster_vs_tracklet_correlation->Fill(nTracklets,nClustersLayer0+nClustersLayer1);

    switch(fUseEventSelection){
    case 0:
        if(!IsEventSelectedLegacy(aev))return kFALSE;
        break;
    case 1:
        if(!IsEventSelectedpp(aev))return kFALSE;
        break;
    case 2:
        if(!IsEventSelectedpA(aev))return kFALSE;
        break;
    case 3:
        if(!IsEventSelectedHF(aev))return kFALSE;
        break;
    }
    fhists_SPD_cluster_vs_tracklet_correlation_PostSelection->Fill(nTracklets,nClustersLayer0+nClustersLayer1);
    return kTRUE;
}


Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedLegacy( AliAODEvent * aev){
    UInt_t res = 0;
    if(aev)res = ((AliVAODHeader*)aev->GetHeader())->GetOfflineTrigger();
    if ((res & AliVEvent::kAny) == 0) {
        return kFALSE;
    }
    fhist_Events->Fill("AliVEvent::kAny",1.);
    if ((res & fOffTrigger) == 0) {
        return kFALSE;
    }
    fhist_Events->Fill("AliVEvent::kMB", 1.);
    if(aev->IsPileupFromSPD(5,0.8, 3.0, 2.0, 5.0)){
        return kFALSE;
    }
    fhist_Events->Fill("!Pileup SPD",1.);
    if(fabs(aev->GetPrimaryVertex()->GetZ())>10.){
        return kFALSE;
    }
    fhist_Events->Fill("Vertex  z < 10 cm", 1.);
    const AliVVertex *trkVtx =    dynamic_cast<const AliVVertex*>(aev->GetPrimaryVertex()) ;
    const AliVVertex* spdVtx =     dynamic_cast<const AliVVertex*>(aev->GetPrimaryVertexSPD()) ;
    TString vtxTtl = trkVtx->GetTitle();
    if (!vtxTtl.Contains("VertexerTracks")){
        return kFALSE;
    }
    fhist_Events->Fill("No VertexerTracks",1.);
    Double_t cov[6]={0};
    spdVtx->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (spdVtx->IsFromVertexerZ() && (zRes>0.25)){
        return kFALSE;
    }
    if ((TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5)){
        return kFALSE;
    }
    fhist_Events->Fill("Vertex Z Resoulution",1.);
    if(trkVtx->GetNContributors()<2) return kFALSE;
    fhist_Events->Fill(">1 contributors",1.);
    if(spdVtx->GetNContributors()<1) return kFALSE;
    fhist_Events->Fill(">0 contributors SPD",1.);
    return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedpA( AliAODEvent * aev){
    if(fUtils->IsFirstEventInChunk(aev)) return kFALSE;
    if(!fUtils->IsVertexSelected2013pA(aev)) return kFALSE;
    if(!fUtils->IsPileUpEvent(aev)) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedpp( AliAODEvent * aev){
    if(!fUtils->IsVertexSelected2013pA(aev)) return kFALSE;
    if(!fUtils->IsPileUpEvent(aev)) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsEventSelectedHF( AliAODEvent * aev){
    if(fJetCutsHF->IsEventSelected(aev)) return kTRUE;
    else return kFALSE;
    return kTRUE;
}
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelectedHF(const AliEmcalJet*jet){
    if(fJetCutsHF->IsJetSelected(jet)) return kTRUE;
    else return kFALSE;
}
Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelectedLegacy(const AliEmcalJet*jet){
    Double_t jetPt = 0.;
    fUseCorrectedJetPt ? jetPt = GetPtCorrected(jet)  : jetPt = jet->Pt();
    if(!(jetPt > 0)) return kFALSE;
    fhist_Jets->Fill("In",jetPt,1.);
    if (!jet) return kFALSE;
    fhist_Jets->Fill("Pointer",jetPt,1.);
    Double_t jetradius = fJetsCont->GetJetRadius();
    jetradius *=jetradius;
    if(jetPt < 1.0) return kFALSE;
    fhist_Jets->Fill("Pt",jetPt,1.);
    if(fabs(jet->Eta()) > 0.5) return kFALSE;
    fhist_Jets->Fill("Eta",jetPt,1.);
    if(jet->Area() < 0.6 * jetradius*TMath::Pi())return kFALSE;
    fhist_Jets->Fill("Area",jetPt,1.);
    return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsJetSelected(const AliEmcalJet*jet)
{
    switch (fUseJetSelection){
    case 0:
        return IsJetSelectedLegacy(jet);
        break;
    case 1:
        return IsJetSelectedHF(jet);
        break;

    }
    return kFALSE;
}


AliRDHFJetsCuts * AliAnalysisTaskEmcalJetBJetTaggingIP::GetJetCutsHF(){
    return fJetCutsHF;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::AddTagJet(AliEmcalJet * jet){
    if(!jet) return kFALSE;
    Int_t partonPDG = 0;
    Double_t jetPt = 0.;
    jetPt = jet->Pt();
    if (fUseMCTagger == 0){
        AliAODMCParticle * parton = NULL;
        parton  = fTaggingHFClass->IsMCJetParton(  fMCparticles,jet,0.7);
        if(!parton) return kFALSE;
        partonPDG =  abs(parton->PdgCode());
        Double_t deltaR =sqrt((parton->Eta() -jet->Eta())*(parton->Eta() -jet->Eta()) + TVector2::Phi_mpi_pi((parton->Phi()-jet->Phi()))*TVector2::Phi_mpi_pi((parton->Phi()-jet->Phi())));
        fhist_parton_genjet_dR->Fill(deltaR);
        fhist_parton_genjet_pT->Fill(parton->Pt(),jet->Pt());
        fhist_parton_genjet_Eta->Fill(parton->Eta() - jet->Eta());
        fhist_parton_genjet_Phi->Fill(TVector2::Phi_mpi_pi(parton->Phi() - jet->Phi()));
    }
    else     if (fUseMCTagger == 1){
        AliAODMCParticle * meson = NULL;
        meson  = fTaggingHFClass->IsMCJetMeson(fMCparticles, jet, 0.7);
        if(meson){
            Int_t mesonPDG = abs(meson->PdgCode() );
            if(mesonPDG >=500 && mesonPDG <600) partonPDG=5;
            else if(mesonPDG >=400 && mesonPDG <500) partonPDG=4;
            else if(mesonPDG >=5000 && mesonPDG <6000) partonPDG=5;
            else if(mesonPDG >=4000 && mesonPDG <5000) partonPDG=4;
        }
        else partonPDG=1; // Set to 1 for LFg
    }
    fhist_MonteCarloFlavour->Fill("Sum",jetPt,1.);

    switch(partonPDG){
    case 21:
        jet->AddFlavourTag(kLFgJet);
        fhist_MonteCarloFlavour->Fill("LFg",jetPt,1.);
        break;
    case 1:
        jet->AddFlavourTag(kLFgJet);
        fhist_MonteCarloFlavour->Fill("LFg",jetPt,1.);
        break;
    case 2:
        jet->AddFlavourTag(kLFgJet);
        fhist_MonteCarloFlavour->Fill("LFg",jetPt,1.);
        break;
    case 3:
        jet->AddFlavourTag(kLFgJet);
        fhist_MonteCarloFlavour->Fill("LFg",jetPt,1.);
        break;
    case 4:
        jet->AddFlavourTag(kCharmJet);
        fhist_MonteCarloFlavour->Fill("Charm",jetPt,1.);
        break;
    case 5:
        jet->AddFlavourTag(kBeautyJet);
        fhist_MonteCarloFlavour->Fill("Beauty",jetPt,1.);
        break;
    default:
        fhist_MonteCarloFlavour->Fill("Untagged",jetPt,1.);
        break;
    }
    return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsQuality(const AliAODTrack *track ,EQualityClass qclass )
{
    if(!track) return kFALSE;
    if(track->Pt() <1.) return kFALSE;
    ULong_t status = track->GetStatus();

    if(!(status & AliAODTrack::kTPCrefit)) return kFALSE;
    if(!(status & AliAODTrack::kITSrefit)) return kFALSE;

    int nSPDHits = 0;
    if (track->HasPointOnITSLayer(0)) nSPDHits++;
    if (track->HasPointOnITSLayer(1)) nSPDHits++;
    int nITSHits = nSPDHits;
    for (int j=2;j<6;++j)  if (track->HasPointOnITSLayer(j)) nITSHits++;
    int nTPCcls =0;
    nTPCcls = ((AliAODTrack*)track)->GetTPCNcls();
    Float_t cRatioTPC = track->GetTPCNclsF() > 0. ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t> (track->GetTPCNclsF()) : 1.;
    Bool_t isV0Daughter = kFALSE;
    Double_t v0Radius = 0.;
    isV0Daughter = IsV0DaughterRadius(track, v0Radius);

    switch(qclass)
    {
    case kQtyVeryGood:
        if(nSPDHits < 2 ) return kFALSE;
        if(nITSHits < 4 ) return kFALSE;
        if(nTPCcls < 90 ) return kFALSE;
        if(cRatioTPC < 0.6 ) return kFALSE;
        if(isV0Daughter &&v0Radius >2.0) return kFALSE;
        break;
    case kQtyGood:
        if(nSPDHits < 1 ) return kFALSE;
        if(nITSHits < 4 ) return kFALSE;
        if(nTPCcls < 90 ) return kFALSE;
        if(cRatioTPC < 0.6 ) return kFALSE;
        if(isV0Daughter &&v0Radius >2.0) return kFALSE;
        break;
    case kQtyMedium:
        if(nSPDHits < 1 ) return kFALSE;
        if(nITSHits < 3 ) return kFALSE;
        if(nTPCcls < 90 ) return kFALSE;
        if(cRatioTPC < 0.6 ) return kFALSE;
        break;
    case kQtyBad:
        if(nITSHits < 3 ) return kFALSE;
        if(nTPCcls < 80 ) return kFALSE;
        if(cRatioTPC < 0.6 ) return kFALSE;
        break;
    default:
        return kFALSE;
        break;
    }
    return kTRUE;
}

Bool_t AliAnalysisTaskEmcalJetBJetTaggingIP::IsV0DaughterRadius(const AliAODTrack *track ,Double_t &Radius){
    AliAODv0 * v0aod =NULL;
    int posid = -1;
    int negid = -1;
    int trackid = -1;
    Double_t P[3];
    for (int i =0 ; i< InputEvent()->GetNumberOfV0s() ;++i )
    {
        memset(P,0,sizeof P);
        v0aod = ((AliAODEvent*) InputEvent())->GetV0(i);
        posid = v0aod->GetPosID();
        negid = v0aod->GetNegID();
        trackid = track->GetID();
        if(posid==trackid || negid==trackid){
            P[0] = v0aod->DecayVertexV0X() ;
            P[1] = v0aod->DecayVertexV0Y() ;
            P[2] = v0aod->DecayVertexV0Z() ;
            Radius = sqrt(P[0]*P[0] +P[1]*P[1]);
            return kTRUE;
        }
    }
    return kFALSE;

}
AliAnalysisTaskEmcalJetBJetTaggingIP::~AliAnalysisTaskEmcalJetBJetTaggingIP()
{
    // Destructor.
}


void AliAnalysisTaskEmcalJetBJetTaggingIP::UserCreateOutputObjects()
{
    AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
    // Create user output.
    fTrackCountingTagger->InitTrackSelectionParams(NULL);

    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    OpenFile(1);
    if(!fOutput){
        fOutput  = new TList();
        fOutput->SetOwner();
    }
    if(!fOutput) {
        AliError("No container!");
        return;
    }
    TList * glist = new TList();
    glist->SetName("module_EventSelection");
    glist->SetOwner(kTRUE);
    fOutput->Add(glist);

    //Event selection histograms
    AddHistTH1 (&fhist_Events,"fhist_Events","Events","", "Events", 8, 0.,8., kTRUE,glist);
    fhist_Events->GetXaxis()->SetBinLabel(1,"AliVEvent::kAny");
    fhist_Events->GetXaxis()->SetBinLabel(2,"AliVEvent::kMB");
    fhist_Events->GetXaxis()->SetBinLabel(3,"!Pileup SPD");
    fhist_Events->GetXaxis()->SetBinLabel(4,"Vertex  z < 10 cm");
    fhist_Events->GetXaxis()->SetBinLabel(5,"No VertexerTracks");
    fhist_Events->GetXaxis()->SetBinLabel(6,"Vertex Z Resoulution");
    fhist_Events->GetXaxis()->SetBinLabel(7,">1 contributors");
    fhist_Events->GetXaxis()->SetBinLabel(8,">0 contributors SPD");
    AddHistTH1 (&fhist_Rho,"fhist_Rho","<#rho>","<#rho> in GeV/c", "a.u.", 500, 0.,10., kTRUE,glist);

    //  Jet selection statistics histograms
    TList * glistjets = new TList();
    glistjets->SetName("module_Jets");
    glistjets->SetOwner(kTRUE);
    fOutput->Add(glistjets);

    AddHistTH2 (&fhist_Jets,"fhist_Jets","Jets","Cut", "Jets", 9, 0.,9.,500,0,250.,kTRUE,glistjets);
    fhist_Jets->GetXaxis()->SetBinLabel(1,"In");
    fhist_Jets->GetXaxis()->SetBinLabel(2,"Pointer");
    fhist_Jets->GetXaxis()->SetBinLabel(3,"Pt");
    fhist_Jets->GetXaxis()->SetBinLabel(4,"Eta");
    fhist_Jets->GetXaxis()->SetBinLabel(5,"Area");
    fhist_Jets->GetXaxis()->SetBinLabel(6,"Matched");
    fhist_Jets->GetXaxis()->SetBinLabel(7,"TaggingN1");
    fhist_Jets->GetXaxis()->SetBinLabel(8,"TaggingN2");
    fhist_Jets->GetXaxis()->SetBinLabel(9,"TaggingN3");

    AddHistTH2(&fhist_Jet_Eta_Phi,"fhist_Jet_Eta_Phi","jet rec. #eta-#phi","#eta","#phi",201,-0.9,0.9, 100, 0.,2.*TMath::Pi(), kTRUE,glistjets);
    AddHistTH2(&fhist_Jet_Nconst_Pt,"fhist_Jet_Nconst_Pt","N constituents vs jet pT","p_T","N constituents",500,0.,250.,100,0.,100.,  kTRUE,glistjets);
    AddHistTH1(&fhist_Jet_Pt,"fhist_Jet_Pt","jet rec. pT","pT","counts",500,0.,250., kTRUE,glistjets);

    fhists_SPD_cluster_vs_tracklet_correlation = new TH2D("fhists_SPD_cluster_vs_tracklet_correlation",";SPD Tracklets;SPD Clusters",200,0.,200.,1000,0.,1000.);
    fOutput->Add(fhists_SPD_cluster_vs_tracklet_correlation);
    fhists_SPD_cluster_vs_tracklet_correlation_PostSelection = new TH2D("fhists_SPD_cluster_vs_tracklet_correlation_PostSelection",";SPD Tracklets;SPD Clusters",200,0.,200.,1000,0.,1000.);
    fOutput->Add(fhists_SPD_cluster_vs_tracklet_correlation_PostSelection);

    if(fDoRandomCones){
        AddHistTH1(&fhist_Jet_Background_Fluctuation,"fhist_Jet_Background_Fluctuation","random cone #delta #it{p}_{T}","#delta #it{p}_{T}","count",500,-125.,125.,  kTRUE,glistjets);
    }
    TList * glistunfold = new TList();
    glistunfold->SetName("module_Unfolding");
    glistunfold->SetOwner(kTRUE);
    fOutput->Add(glistunfold);
    const char * uf_flavour[4] = {"inclusive","lfg","beauty","charm"};
    const char * uf_corr[2] = {"raw","corrected"};
    int uplimit=0;
    fIsMC ? uplimit =4 :uplimit=1;
    for(int itcflavour =0; itcflavour<4;itcflavour++)
        for(int itcorr =0; itcorr<2;itcorr++){
            AddHistTH2 (&(fhist_momentum_response[itcflavour][itcorr]),
                        Form("fhist_momentum_response_%s_%s",uf_flavour[itcflavour],uf_corr[itcorr]),
                        Form("Momentum response matrix (%s,%s)",uf_flavour[itcflavour],uf_corr[itcorr]),"#it{p}^{Jet,gen}_{T}", "#it{p}^{Jet,rec}_{T}",500,0.,250.,500,0.,250.,kTRUE,glistunfold);

        }
    //Track counting histograms
    TList * glisttc = new TList();
    glisttc->SetName("module_TrackCounting");
    glisttc->SetOwner(kTRUE);
    fOutput->Add(glisttc);

    const char * tc_track[3] = {"n_1","n_2","n_3"};
    const char * tc_flavour[4] = {"inclusive","lfg","beauty","charm"};
    const char * tc_classes[5] = {"default","verygood","good","medium","bad"};

    for(int itctrack =0; itctrack<3;itctrack++){
        for(int itcflavour =0; itcflavour<4;itcflavour++){
            if(!fIsMC && itcflavour >0) continue;
            for(int itcclass =0; itcclass<5;itcclass++){

                AddHistTH2(&(fhist_TC_sIP_Pt[itctrack][itcflavour][itcclass]),
                           Form("fhist_TC_%s_%s_%s",tc_track[itctrack],tc_flavour[itcflavour],tc_classes[itcclass]),
                           Form("Track Counting Output %s_%s_%s",tc_track[itctrack],tc_flavour[itcflavour],tc_classes[itcclass]),
                           "Signed impact parameter (cm)",
                           "#it{p}^{Jet,rec}_{T}",
                           1001,-0.5,0.5, 500, 0.,250, kTRUE,glisttc);
                AddHistTH2(&(fhist_TC_Eta_Phi[itctrack][itcflavour][itcclass]),
                           Form("fhist_TC_Eta_Phi_%s_%s_%s",tc_track[itctrack],tc_flavour[itcflavour],tc_classes[itcclass]),
                           Form("Track Counting #eta-#phi %s_%s_%s",tc_track[itctrack],tc_flavour[itcflavour],tc_classes[itcclass]),
                           "#eta",
                           "#phi",
                           501,-0.9,0.9, 500, 0.,2.*TMath::Pi(), kTRUE,glisttc);
            }
        }}
    if(fIsTrackQAConstituent){
        TList * glisttracksjet = new TList();
        glisttracksjet->SetName("module_TrackQAJet");
        glisttracksjet->SetOwner(kTRUE);
        fOutput->Add(glisttracksjet);
        AddHistTH1 (&fhist_QualityClasses,"fhist_QualityClasses","Tracks in quality class ", "Classname", "count", 5, 0.,5., kTRUE,glisttracksjet);
        fhist_QualityClasses->GetXaxis()->SetBinLabel(1,"all");
        fhist_QualityClasses->GetXaxis()->SetBinLabel(2,"kQtyVeryGood");
        fhist_QualityClasses->GetXaxis()->SetBinLabel(3,"kQtyGood");
        fhist_QualityClasses->GetXaxis()->SetBinLabel(4,"kQtyMedium");
        fhist_QualityClasses->GetXaxis()->SetBinLabel(5,"kQtyBad");

        if(fIsMC){
            const char * t_classes[4] = {"verygood","good","medium","bad"};
            const char * t_type[2] = {"Primary","Secondary"};

            for(int itcclass =0; itcclass<4;itcclass++)
                for(int ittype =0; ittype<2;ittype++){
                    AddHistTH2 (&(fhist_QualityClasses_sIP[itcclass][ittype]),
                                Form("fhist_QualityClasses_sIP_%s_%s",t_classes[itcclass],t_type[ittype]),
                                Form("Track Probability PDF sIP distributions (%s,%s)",t_classes[itcclass],t_type[ittype]),"Signed impact parameter (cm)", "#it{p}^{Jet,rec}_{T}",1001,-0.5,0.5,500,0.,250.,kTRUE,glisttracksjet);

                    AddHistTH2 (&(fhist_QualityClasses_Eta_Phi[itcclass][ittype]),
                                Form("fhist_QualityClasses_Eta_Phi_%s_%s",t_classes[itcclass],t_type[ittype]),
                                Form("Track Probability classes #eta #phi (%s,%s)",t_classes[itcclass],t_type[ittype]),"#eta", "#phi", 501,-0.9,0.9, 500, 0.,2.*TMath::Pi(),kTRUE,glisttracksjet);
                }

        }

    }
    if(fIsTrackQA){
        TList * glisttracks = new TList();
        glisttracks->SetName("module_TrackQA");
        glisttracks->SetOwner(kTRUE);
        fOutput->Add(glisttracks);
        AddHistTH2 (&fhist_Tracks_Eta_Phi,"fhist_Tracks_Eta_Phi","Accepted Tracks(9||4)","Eta", "Phi",201,-1.5,1.5, 100, 0.,2.*TMath::Pi(), kTRUE,glisttracks);
        AddHistTH2 (&fhist_Tracks_Eta_Phi_Bit4,"fhist_Tracks_Eta_Phi_Bit4","Accepted Tracks(4)","Eta", "Phi",201,-1.5,1.5, 100, 0.,2.*TMath::Pi(), kTRUE,glisttracks);
        AddHistTH2 (&fhist_Tracks_Eta_Phi_Bit9,"fhist_Tracks_Eta_Phi_Bit9","Accepted Tracks(9)","Eta", "Phi",201,-1.5,1.5, 100, 0.,2.*TMath::Pi(), kTRUE,glisttracks);
    }
    if(fIsMC){
        //Monte carlo jet flavour statistics histograms
        TList * glistmc = new TList();
        glistmc->SetName("module_MonteCarlo");
        glistmc->SetOwner(kTRUE);
        fOutput->Add(glistmc);
        fhist_MonteCarloFlavour = new TH2D("fhist_MonteCarloFlavour","fhist_MonteCarloFlavour;pT;Flavour",5,0.,5.,500,0.,250.);
        fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(1,"LFg");
        fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(2,"Charm");
        fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(3,"Beauty");
        fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(4,"Untagged");
        fhist_MonteCarloFlavour->GetXaxis()->SetBinLabel(5,"Sum");
        AddHistTH1 (&fhist_parton_genjet_dR,"fhist_parton_genjet_dR","#Delta R parton gen. jet ", "#Delta R", "count", 100, 0.,1., kTRUE,glistmc);
        AddHistTH2 (&fhist_parton_genjet_pT,"fhist_parton_genjet_pT","parton gen. jet pT", "pT parton", "pT gen. jet ", 500, 0.,250., 500, 0.,250., kTRUE,glistmc);
        AddHistTH1 (&fhist_parton_genjet_Eta,"fhist_parton_genjet_Eta","#Delta #eta  parton gen. jet ", "#Delta #eta", "count", 500, -1.,1., kTRUE,glistmc);
        AddHistTH1 (&fhist_parton_genjet_Phi,"fhist_parton_genjet_Phi","#Delta #phi parton gen. jet ", "#Delta #phi", "count", 500,-2.,2., kTRUE,glistmc);
        AddHistTH2 (&fhist_MonteCarloFlavour,"fhist_MonteCarloFlavour","MC flavour", "MC flavour Tag", "#it{p}^{gen}_{T,jet} (GeV/#it{c})",5, 0,5.,500,0.,250,kTRUE,glistmc);

    }
    //SetContainer
    fJetsCont           = static_cast<AliJetContainer*>(fJetCollArray.At(0));
    fJetsContMC         = static_cast<AliJetContainer*>(fJetCollArray.At(1));
    if(fJetsCont)
        fTracksCont = fJetsCont->GetParticleContainer();
    else  fTracksCont = GetParticleContainer(0);
    if(fTracksCont) fTracksCont->SetClassName("AliVTrack");
    fTrackCountingTagger->SetParticleContainer(fTracksCont);
    fTrackCountingTagger->SetAnalysisTypeAOD(kTRUE);
    TH1::AddDirectory(oldStatus);
    PostData(1, fOutput); // Post data for ALL output slots > 0 here.
    return;
}
// Helper functions
void AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistTH1 (TH1 **hist,const char* histname,const char * title,const char *titlex, const char *titley,Int_t nBinsX, Double_t minX,Double_t maxX,Bool_t setSumw2,TList * container){
    if(!container)return;
    *hist = new TH1D(histname,Form("%s;%s;%s",title,titlex,titley),nBinsX,minX,maxX);
    if (!(*hist)) return;
    if (setSumw2) (*hist)->Sumw2();
    if (container->FindObject(histname)){
        AliError( Form("Object with name  %s already exists in %s...returning!",histname,container->GetName()));
        return;
    }
    container->Add(*hist);
    return;
}
void AliAnalysisTaskEmcalJetBJetTaggingIP::AddHistTH2 (TH2 **hist,const char* histname,const char * title,const char *titlex, const char *titley,Int_t nBinsX, Double_t minX,Double_t maxX,Int_t nBinsY, Double_t minY,Double_t maxY,Bool_t setSumw2,TList * container){
    if(!container)return;

    *hist = new TH2D(histname,Form("%s;%s;%s",title,titlex,titley),nBinsX,minX,maxX,nBinsY,minY,maxY);
    if (!(*hist)) return;
    if (setSumw2) (*hist)->Sumw2();
    if (container->FindObject(histname)){
        AliError( Form("Object with name  %s already exists in %s...returning!",histname,container->GetName()));
        return;
    }
    container->Add(*hist);
    return;
}

Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetDeltaPtRandomCone(){
    Double_t deltaPt =-1000.;
    Double_t jetradius =fJetsCont->GetJetRadius();
    Double_t minEta = -0.9 + jetradius;
    Double_t maxEta =  0.9 - jetradius;
    Double_t tmpRandConeEta = minEta + fRandom->Rndm()*(maxEta-minEta);
    Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();
    Double_t tmpConePt = -1.;

    for (Int_t i = 0; i <  static_cast<AliAODEvent*>(InputEvent())->GetNumberOfTracks(); i++){
        AliAODTrack* tmpTrack = static_cast<AliAODTrack*>(InputEvent()->GetTrack(i));
        if(fabs(tmpTrack->Eta())< 0.9){
            if(tmpTrack->Pt()> 0.15){
                if(sqrt((tmpTrack->Eta() -tmpRandConeEta)*(tmpTrack->Eta() -tmpRandConeEta) + TVector2::Phi_mpi_pi((tmpTrack->Phi()-tmpRandConePhi))*TVector2::Phi_mpi_pi((tmpTrack->Phi()-tmpRandConePhi)))<jetradius){
                    tmpConePt+=tmpTrack->Pt();
                }
            }
        }
    }
    if(tmpConePt>0){
        deltaPt = tmpConePt - 0.4*0.4*TMath::Pi()*fJetsCont->GetRhoVal();
        return deltaPt;
    }
    return deltaPt;
}

Double_t AliAnalysisTaskEmcalJetBJetTaggingIP::GetPtCorrected(const AliEmcalJet * jet){
    if(jet && fJetsCont)
        return  jet->Pt() - fJetsCont->GetRhoVal()*jet->Area();
    return -1.;
}
