#include "AliAnalysisTaskJetAngCorrelations.h"
#include "AliInputEventHandler.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliAODTrackSelection.h"
#include "AliESDtrackCuts.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

using namespace std;
ClassImp(AliAnalysisTaskJetAngCorrelations)
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetAngCorrelations::AliAnalysisTaskJetAngCorrelations():
AliAnalysisTaskSE()
//JET RECONSTRUCTION
,fAODEvent(nullptr)
,fPIDResponse(nullptr)
,fAODTrackCuts(nullptr)
,fESDTrackCuts(nullptr)
,fOutputData(nullptr)
,fQAList(nullptr)
,fTrackBufferProton()
,fTrackBufferAntiproton()
,fTrackBufferDeuteron()
,fTrackBufferAntideuteron()
,fTrackBufferSize()
//
,hNumberOfEvents(nullptr)
,hNumberOfJets(nullptr)
,hEventProtocol(nullptr)
,hTrackProtocol(nullptr)
,hNumberOfParticlesInJet(nullptr)
,hNumberOfJetsInEvent(nullptr)
//
,hITSSignal(nullptr)
,hTPCSignal(nullptr)
,hTOFSignal(nullptr)
,hTPCnsigma(nullptr)
,hTPCnsigmaProton(nullptr)
,hTPCnsigmaDeuteron(nullptr)
,hTOFnsigma(nullptr)
,hTOFnsigmaProton(nullptr)
,hTOFnsigmaDeuteron(nullptr)
,hITSnsigma(nullptr)
,hITSnsigmaProton(nullptr)
,hITSnsigmaDeuteron(nullptr)
//
,hDCAxyFullEvent(nullptr)
,hDCAzFullEvent(nullptr)
,hDCAzJetParticle(nullptr)
//
,hPtFullEvent(nullptr)
,hPtJetParticle(nullptr)
,hJetRapidity(nullptr)
,hPtTotalJet(nullptr)
,hPtSubtractedJet(nullptr)
,hPtJetProtonDeuteron(nullptr)
,hPtJetDeuteron(nullptr)
,hPtDiff(nullptr)
//
,hJetConeRadius(nullptr)
,hEtaFullEvent(nullptr)
//
,hDeltaPhiSE(nullptr)
,hDeltaPhiEtaSE(nullptr)
,hDeltaPhiME(nullptr)
,hDeltaPhiEtaME(nullptr)
//Preliminary Cuts
,fKinkDaugh(kFALSE), fTPCRefit(kFALSE), fITSRefit(kFALSE), fTPCnCr(70), fTPCnCrnCl(0.7), fITSnCl(2), fITSchi2(36), fTPCchi2(4), fEtaCut(0.8), fDCAxy(0.5), fDCAz(2.4)
//Jet Cuts
,fJetRadius(0.4), fMinJetPt(10.0), fMinJetParticlePt(0.0), fMinLeadingPt(5.0)
//Proton Cuts
,fProtonDCAxy(0.5), fProtonDCAz(1.0), fProtonITSTOFpT(0.7), fProtonITSnsig(3.0), fProtonTPCnsigITS(5.0), fProtonTPCnsigTOF(4.0), fProtonTOFnsig(10.0)
//Antiproton Cuts
,fAntiprotonDCAxy(0.5), fAntiprotonDCAz(1.0), fAntiprotonITSTOFpT(0.7), fAntiprotonITSnsig(3.0), fAntiprotonTPCnsigITS(5.0), fAntiprotonTPCnsigTOF(4.0), fAntiprotonTOFnsig(10.0)
//Deuteron Cuts
,fDeuteronDCAxy(0.5), fDeuteronDCAz(1.0), fDeuteronITSTOFpT(0.7), fDeuteronITSnsig(3.0), fDeuteronTPCnsigITS(5.0), fDeuteronTPCnsigTOF(4.0), fDeuteronTOFnsig(10.0)
//Antideuteron Cuts
,fAntideuteronDCAxy(0.5), fAntideuteronDCAz(1.0), fAntideuteronITSTOFpT(0.7), fAntideuteronITSnsig(3.0), fAntideuteronTPCnsigITS(5.0), fAntideuteronTPCnsigTOF(4.0), fAntideuteronTOFnsig(10.0)
{}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetAngCorrelations::AliAnalysisTaskJetAngCorrelations(const char *name):
AliAnalysisTaskSE(name)
//JET RECONSTRUCTION
,fAODEvent(nullptr)
,fPIDResponse(nullptr)
,fAODTrackCuts(nullptr)
,fESDTrackCuts(nullptr)
,fOutputData(nullptr)
,fQAList(nullptr)
,fTrackBufferProton()
,fTrackBufferAntiproton()
,fTrackBufferDeuteron()
,fTrackBufferAntideuteron()
,fTrackBufferSize(2000)
//
,hNumberOfEvents(nullptr)
,hNumberOfJets(nullptr)
,hEventProtocol(nullptr)
,hTrackProtocol(nullptr)
,hNumberOfParticlesInJet(nullptr)
,hNumberOfJetsInEvent(nullptr)
//
,hITSSignal(nullptr)
,hTPCSignal(nullptr)
,hTOFSignal(nullptr)
,hTPCnsigma(nullptr)
,hTPCnsigmaProton(nullptr)
,hTPCnsigmaDeuteron(nullptr)
,hTOFnsigma(nullptr)
,hTOFnsigmaProton(nullptr)
,hTOFnsigmaDeuteron(nullptr)
,hITSnsigma(nullptr)
,hITSnsigmaProton(nullptr)
,hITSnsigmaDeuteron(nullptr)
//
,hDCAxyFullEvent(nullptr)
,hDCAzFullEvent(nullptr)
,hDCAzJetParticle(nullptr)
//
,hPtFullEvent(nullptr)
,hPtJetParticle(nullptr)
,hJetRapidity(nullptr)
,hPtTotalJet(nullptr)
,hPtSubtractedJet(nullptr)
,hPtJetProtonDeuteron(nullptr)
,hPtJetDeuteron(nullptr)
,hPtDiff(nullptr)
//
,hJetConeRadius(nullptr)
,hEtaFullEvent(nullptr)
//
,hDeltaPhiSE(nullptr)
,hDeltaPhiEtaSE(nullptr)
,hDeltaPhiME(nullptr)
,hDeltaPhiEtaME(nullptr)
//Preliminary Cuts
,fKinkDaugh(kFALSE), fTPCRefit(kFALSE), fITSRefit(kFALSE), fTPCnCr(70), fTPCnCrnCl(0.7), fITSnCl(2), fITSchi2(36), fTPCchi2(4), fEtaCut(0.8), fDCAxy(0.5), fDCAz(2.4)
//Jet Cuts
,fJetRadius(0.4), fMinJetPt(10.0), fMinJetParticlePt(0.0), fMinLeadingPt(5.0)
//Proton Cuts
,fProtonDCAxy(0.5), fProtonDCAz(1.0), fProtonITSTOFpT(0.7), fProtonITSnsig(3.0), fProtonTPCnsigITS(5.0), fProtonTPCnsigTOF(4.0), fProtonTOFnsig(10.0)
//Antiproton Cuts
,fAntiprotonDCAxy(0.5), fAntiprotonDCAz(1.0), fAntiprotonITSTOFpT(0.7), fAntiprotonITSnsig(3.0), fAntiprotonTPCnsigITS(5.0), fAntiprotonTPCnsigTOF(4.0), fAntiprotonTOFnsig(10.0)
//Deuteron Cuts
,fDeuteronDCAxy(0.5), fDeuteronDCAz(1.0), fDeuteronITSTOFpT(0.7), fDeuteronITSnsig(3.0), fDeuteronTPCnsigITS(5.0), fDeuteronTPCnsigTOF(4.0), fDeuteronTOFnsig(10.0)
//Antideuteron Cuts
,fAntideuteronDCAxy(0.5), fAntideuteronDCAz(1.0), fAntideuteronITSTOFpT(0.7), fAntideuteronITSnsig(3.0), fAntideuteronTPCnsigITS(5.0), fAntideuteronTPCnsigTOF(4.0), fAntideuteronTOFnsig(10.0)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskJetAngCorrelations::~AliAnalysisTaskJetAngCorrelations()  {
    //JET RECONSTRUCTION
    if (fOutputData) {
        delete fOutputData;
    }
    if (fQAList) {
        delete fQAList;
    }
    delete fAODEvent;
    delete fPIDResponse;
    delete fAODTrackCuts;
    delete fESDTrackCuts;

}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetAngCorrelations::UserCreateOutputObjects()  {
    //Create Output List
    fOutputData  = new TList();
    fQAList      = new TList();
    fOutputData  -> SetOwner();
    fQAList      -> SetOwner();

    //Histograms ------------------------------------------------------------------------
    //Counters
    hNumberOfEvents = new TH1F ("hNumberOfEvents","Total number of events", 1,0,1);
    fOutputData     -> Add(hNumberOfEvents);
    hNumberOfJets   = new TH1F ("hNumberOfJets","Total number of jets", 1,0,1);
    fOutputData     -> Add(hNumberOfJets);
    hEventProtocol  = new TH1F ("hEventProtocol","Event protocol", 20,0,20);
    fOutputData     -> Add(hEventProtocol);
    hTrackProtocol  = new TH1F ("hTrackProtocol","Track protocol", 20,0,20);
    fOutputData     -> Add(hTrackProtocol);
    hNumberOfParticlesInJet = new TH1I ("hNumberOfParticlesInJet","Number of particles in jet", 200,0,200); // number of particles per jet
    fOutputData             -> Add(hNumberOfParticlesInJet);
    hNumberOfJetsInEvent    = new TH1I ("hNumberOfJetsInEvent", "Number of PseudoJets per event", 10, 0, 10); // how many jets FastJet collects per event
    fOutputData             -> Add(hNumberOfJetsInEvent);

    //(Pseudo)Rapidity
    hEtaFullEvent   = new TH1F ("hEtaFullEvent", "Pseudorapidity of full event;#eta", 200,-1,1); // eta of all particles
    hEtaFullEvent   -> Sumw2();
    fOutputData     -> Add(hEtaFullEvent);
    hJetRapidity    = new TH1D ("hJetRapidity", "Rapidity of jet;y", 200,-1,1);
    hJetRapidity    -> Sumw2();
    fOutputData     -> Add(hJetRapidity);

    //pT
    hPtFullEvent    = new TH1D ("hPtFullEvent", "pT of full event;pT [GeV/c]", 1000,0,100); // pT of all particles
    hPtFullEvent    -> Sumw2();
    fOutputData     -> Add(hPtFullEvent);
    hPtJetParticle  = new TH1D ("hPtJetParticle", "pT of jet particles;pT [GeV/c]", 1000,0,100); // pT of jet using FastJet
    hPtJetParticle  -> Sumw2();
    fOutputData     -> Add(hPtJetParticle);
    hPtSubtractedJet = new TH1D ("hPtSubtractedJet","pT of jet particles (background subtracted);pT [GeV/c]", 1000,0,100); // pT of background subtracted jet
    hPtSubtractedJet -> Sumw2();
    fOutputData      -> Add(hPtSubtractedJet);
    hPtJetProtonDeuteron = new TH2D ("hPtJetProtonDeuteron", "pT of (anti)p, (anti)d;[p, ap, d, ad];pT [GeV/c]", 4,1,5, 1000,0,100); // pT of p, ap, d, ad in jet
    hPtJetProtonDeuteron -> Sumw2();
    fOutputData     -> Add(hPtJetProtonDeuteron);
    hPtTotalJet     = new TH1D ("hPtTotalJet", "pT of entire jet;pT [GeV/c]", 2000,0,500); // pT of total jet
    hPtTotalJet     -> Sumw2();
    fOutputData     -> Add(hPtTotalJet);
    hPtDiff         = new TH1D ("hPtDiff","pT difference between PseudoJet and corresponding track;pT [GeV/c]",100,-0.0000005,0.0000005);
    hPtDiff         -> Sumw2();
    fOutputData     -> Add(hPtDiff);

    //nSigma
    hITSSignal         = new TH2F ("hITSSignal", "ITS dE/dx for full event;pT [GeV/c];dE/dx", 1,0,1, 1,0,1);
    hITSSignal         -> Sumw2();
    fOutputData        -> Add(hITSSignal);
    hTPCSignal         = new TH2F ("hTPCSignal", "TPC dE/dx for full event;pT [GeV/c];dE/dx", 1,0,1, 1,0,1);
    hTPCSignal         -> Sumw2();
    fOutputData        -> Add(hTPCSignal);
    hTOFSignal         = new TH2F ("hTOFSignal", "TOF #beta for full event;pT [GeV/c];dE/dx", 1,0,1, 1,0,1);
    hTOFSignal         -> Sumw2();
    fOutputData        -> Add(hTOFSignal);
    hTPCnsigma         = new TH2F ("hTPCnsigma", "TPC n#sigma for full event;pT [GeV/c];n#sigma", 1000,-50,50, 1000,-5,5);
    hTPCnsigma         -> Sumw2();
    fOutputData        -> Add(hTPCnsigma);
    hTPCnsigmaProton   = new TH2F ("hTPCnsigmaProton", "TPC n#sigma for (anti)proton;pT [GeV/c];n#sigma", 1000,-50,50, 1000,-5,5);
    hTPCnsigmaProton   -> Sumw2();
    fOutputData        -> Add(hTPCnsigmaProton);
    hTPCnsigmaDeuteron = new TH2F ("hTPCnsigmaDeuteron", "TPC n#sigma for (anti)deuteron;pT [GeV/c];n#sigma", 1000,-50,50, 1000,-5,5);
    hTPCnsigmaDeuteron -> Sumw2();
    fOutputData        -> Add(hTPCnsigmaDeuteron);
    hTOFnsigma         = new TH2F ("hTOFnsigma", "TOF n#sigma for full event;pT [GeV/c];n#sigma", 200,-50,50, 1000,-10,10);
    hTOFnsigma         -> Sumw2();
    fOutputData        -> Add(hTOFnsigma);
    hTOFnsigmaProton   = new TH2F ("hTOFnsigmaProton", "TOF n#sigma for (anti)proton;pT [GeV/c];n#sigma", 200,-50,50, 1000,-10,10);
    hTOFnsigmaProton   -> Sumw2();
    fOutputData        -> Add(hTOFnsigmaProton);
    hTOFnsigmaDeuteron = new TH2F ("hTOFnsigmaDeuteron", "TOF n#sigma for (anti)deuteron;pT [GeV/c];n#sigma", 200,-50,50, 1000,-10,10);
    hTOFnsigmaDeuteron -> Sumw2();
    fOutputData        -> Add(hTOFnsigmaDeuteron);
    hITSnsigma         = new TH2F ("hITSnsigma", "ITS n#sigma for full event;pT [GeV/c];n#sigma", 200,-5,5, 1000,-5,5);
    hITSnsigma         -> Sumw2();
    fOutputData        -> Add(hITSnsigma);
    hITSnsigmaProton   = new TH2F ("hITSnsigmaProton", "ITS n#sigma for (anti)proton;pT [GeV/c];n#sigma", 200,-50,50, 1000,-10,10);
    hITSnsigmaProton   -> Sumw2();
    fOutputData        -> Add(hITSnsigmaProton);
    hITSnsigmaDeuteron = new TH2F ("hITSnsigmaDeuteron", "ITS n#sigma for (anti)deuteron;pT [GeV/c];n#sigma", 200,-50,50, 1000,-10,10);
    hITSnsigmaDeuteron -> Sumw2();
    fOutputData        -> Add(hITSnsigmaDeuteron);

    //DCA
    hDCAxyFullEvent  = new TH1F ("hDCAxyFullEvent", "DCA xy of full event;DCAxy [cm]", 200,-1,1);
    hDCAxyFullEvent  -> Sumw2();
    fOutputData      -> Add(hDCAxyFullEvent);
    hDCAzFullEvent   = new TH1F ("hDCAzFullEvent", "DCAz of full event;DCAz [cm]", 200,-2.4,2.4);
    hDCAzFullEvent   -> Sumw2();
    fOutputData      -> Add(hDCAzFullEvent);
    hDCAzJetParticle = new TH1F ("hDCAzJetParticle", "hDCAz of charged particles in jets;DCAz [cm]", 200,-2.4,2.4);
    hDCAzJetParticle -> Sumw2();
    fOutputData      -> Add(hDCAzJetParticle);

    hJetConeRadius            = new TH1D ("hJetConeRadius", "Radius of jet cone", 100,0,1);
    hJetConeRadius            -> Sumw2();
    fOutputData               -> Add(hJetConeRadius);

    //Angular Distributions
    hDeltaPhiSE    = new TH2D ("hDeltaPhiSE", "#Delta#varphi of (anti)p, (anti)d in single event;[p, ap, d, ad];#Delta#varphi", 4,1,5, 1000,-2,5);
    hDeltaPhiSE    -> Sumw2();
    fOutputData    -> Add(hDeltaPhiSE);
    hDeltaPhiEtaSE = new TH3D ("hDeltaPhiEtaSE", "#Delta#varphi vs #Delta#eta of (anti)p, (anti)d in single event;[p, ap, d, ad];#Delta#eta", 4,1,5, 1000,-2,5, 1000,-2,2);
    hDeltaPhiEtaSE -> Sumw2();
    fOutputData    -> Add(hDeltaPhiEtaSE);
    hDeltaPhiME    = new TH2D ("hDeltaPhiME", "#Delta#varphi of (anti)p, (anti)d in mixed event;[p, ap, d, ad];#Delta#varphi", 4,1,5, 1000,-2,5);
    hDeltaPhiME    -> Sumw2();
    fOutputData    -> Add(hDeltaPhiME);
    hDeltaPhiEtaME  = new TH3D ("hDeltaPhiEtaME", "#Delta#varphi vs #Delta#eta of (anti)p, (anti)d in mixed event;[p, ap, d, ad];#Delta#eta", 4,1,5, 1000,-2,5, 1000,-2,2);
    hDeltaPhiEtaME -> Sumw2();
    fOutputData    -> Add(hDeltaPhiEtaME);

    //------------------------------------------------------------------------------------
    //Track Selection
    fESDTrackCuts = new AliESDtrackCuts("fESDTrackCuts");
    fESDTrackCuts -> SetAcceptKinkDaughters(fKinkDaugh);
    fESDTrackCuts -> SetRequireTPCRefit(fTPCRefit);
    fESDTrackCuts -> SetRequireITSRefit(fITSRefit);
    fESDTrackCuts -> SetMinNCrossedRowsTPC(fTPCnCr);
    fESDTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(fTPCnCrnCl);
    fESDTrackCuts -> SetMinNClustersITS(fITSnCl);
    fESDTrackCuts -> SetMaxChi2PerClusterITS(fITSchi2);
    fESDTrackCuts -> SetMaxChi2PerClusterTPC(fTPCchi2);

    fAODTrackCuts = new AliAODTrackSelection(fESDTrackCuts, 1);

    //Post Data
    PostData(1, fOutputData);
    PostData(2, fQAList);
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetAngCorrelations::UserExec(Option_t *)  {

    //Get Input Event
    if (!GetEvent()) return;

    RunData();

}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetAngCorrelations::RunData()  {
    //TRACK SELECTION
    AliAODEvent *fAODEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if (!fAODEvent) Fatal("AliAnalysisTaskJetAngCorrelations::RunData", "No AOD event found, check the event handler.");

    //Initialisation
    AliAODTrack* leadingTrack;
    vector<Int_t> particle_ID;
    map<int,AliAODTrack*> particles;
    vector<fastjet::PseudoJet> jetInput;
    vector<fastjet::PseudoJet> jets;
    vector<fastjet::PseudoJet> constituents;
    int index = 0;
    fastjet::PseudoJet hardestJet(0.,0.,0.,0.);
    fastjet::PseudoJet subtractedJet(0.,0.,0.,0.);
    particles.clear();
    jetInput.clear();
    jets.clear();
    constituents.clear();

    Double_t pTbinningTPC[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0};
    Double_t pTbinningTOF[] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0};

    //Preliminary Particle Selection
    for (Int_t i=0; i<fAODEvent->GetNumberOfTracks(); i++) {
        AliAODTrack *track = static_cast<AliAODTrack*>(fAODEvent->GetTrack(i));
        if (!track) continue;
        hTrackProtocol -> Fill(0.5); // # tracks before cuts
        Double_t eta = track->Eta();

        if (TMath::Abs(eta)>fEtaCut) continue;
        hTrackProtocol -> Fill(1.5); // # tracks with |eta| < 0.8

        if (!PassedTrackSelection(track)) continue; // check track cuts
        hTrackProtocol -> Fill(2.5); // # tracks after basic cuts
        hEtaFullEvent  -> Fill(eta);
        hPtFullEvent   -> Fill(track->Pt());

        hDCAzFullEvent  -> Fill(GetDCAtoPrimaryVertex(track,1));
        hDCAxyFullEvent -> Fill(GetDCAtoPrimaryVertex(track,0));

        //Full Event nSigma
        Double_t nSigmaTPC_proton = fPIDResponse -> NumberOfSigmasTPC(track, AliPID::kProton);
        Double_t nSigmaTOF_proton = fPIDResponse -> NumberOfSigmasTOF(track, AliPID::kProton);
        Double_t nSigmaITS_proton = fPIDResponse -> NumberOfSigmasITS(track, AliPID::kProton);
        hTPCnsigma -> Fill((track->Pt()*track->Charge()),nSigmaTPC_proton);
        hTOFnsigma -> Fill((track->Pt()*track->Charge()),nSigmaTOF_proton);
        hITSnsigma -> Fill((track->Pt()*track->Charge()),nSigmaITS_proton);

        if (track->Pt()>fMinLeadingPt) leadingTrack = track;
        particle_ID.emplace_back(index);
        particles[index] = track;
        fastjet::PseudoJet inputPseudoJet(track->Px(),track->Py(),track->Pz(),track->E());
        inputPseudoJet.set_user_index(index);
        jetInput.emplace_back(inputPseudoJet);
        index++;
    } //for (Int_t i=0; i<fAODEvent->GetNumberOfTracks(); i++)

    if (fMinLeadingPt>5.0) hTrackProtocol -> Fill(3.5); // # leading tracks with pT > 5.0

    if ((Int_t)particles.size()<2) return;
    hEventProtocol -> Fill(1.5); // # events with more than 1 particle

    //JET RECONSTRUCTION -------------------------------------------------------------------------------------------------------------------
    double ghost_maxrap = 1.0;
    int ghost_repeat    = 1;
    double ghost_area   = 0.005;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, fJetRadius);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, 0.5); // maybe increase to R=0.6
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap,ghost_repeat,ghost_area));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef);
    jets = sorted_by_pt(clusterSeq.inclusive_jets());
    if (jets.size()==0) return;
    hEventProtocol -> Fill(2.5); // # events with jet
    hNumberOfJetsInEvent -> Fill(jets.size());
    for (int i=0; i<(int)jets.size(); i++) {
        hNumberOfJets -> Fill(0.5);
    }

    hardestJet = jets[0];

    if (hardestJet.pt()<fMinJetPt) return;
    hEventProtocol -> Fill(3.5); // # events with jet with pT > 10 GeV
    hJetRapidity   -> Fill(hardestJet.rap());
    hPtTotalJet    -> Fill(hardestJet.pt());

    if (hardestJet.constituents().size()<2) return;
    hEventProtocol -> Fill(4.5); // # jets with more than 2 particles (should be all of them)
    constituents = hardestJet.constituents();

    for (Int_t j=0; j<(int)constituents.size(); j++) {// jet radius histogram
        hPtJetParticle    -> Fill(constituents[j].pt());
        Double_t DeltaPhi = TVector2::Phi_0_2pi(constituents[j].phi() - hardestJet.phi());
        Double_t DeltaEta = constituents[j].eta() - hardestJet.eta();
        Double_t Delta    = TMath::Sqrt(DeltaPhi*DeltaPhi + DeltaEta*DeltaEta);
        hJetConeRadius    -> Fill(Delta);
    }

    //Background Subtraction
    fastjet::Selector selector = fastjet::SelectorAbsEtaMax(1.0) * (!fastjet::SelectorNHardest(2));
    fastjet::JetMedianBackgroundEstimator bkgEst(selector, jetDefBkg, areaDefBkg);
    fastjet::Subtractor subtractor(&bkgEst);
    bkgEst.set_particles(jetInput);

    subtractedJet = subtractor(hardestJet);
    if (subtractedJet.has_constituents()) {
        for (int i=0; i<(int)subtractedJet.constituents().size(); i++) {
            hPtSubtractedJet -> Fill(subtractedJet.constituents()[i].pt());
        }
    }


    //Jet Proton Selection
    map<TString,vector<AliAODTrack*> > jetNuclei;
    vector<AliAODTrack*> jetProtons;
    vector<AliAODTrack*> jetAntiprotons;
    vector<AliAODTrack*> jetDeuterons;
    vector<AliAODTrack*> jetAntideuterons;

    for (Int_t i=0; i<(Int_t)constituents.size();i++) {
        fastjet::PseudoJet pseudoParticle = constituents[i];
        int id = pseudoParticle.user_index();
        hTrackProtocol -> Fill(4.5); // # jet particle ids

        if (id < 0) continue;
        hTrackProtocol -> Fill(5.5); // # non-negative ids
        AliAODTrack *jetParticle = particles[id];
        if (!jetParticle) continue;

        double jetParticlePt = (double)jetParticle->Pt();
        double diff = (double)pseudoParticle.pt() - jetParticlePt;
        int particleType(0); // 1=p, 2=ap, 3=d, 4=ad
        hPtDiff -> Fill(diff);
        hDCAzJetParticle -> Fill(GetDCAtoPrimaryVertex(jetParticle,1));
        if (jetParticlePt<fMinJetParticlePt) continue;
        if (IsProton(jetParticle) || IsAntiproton(jetParticle)) {// collect (anti)protons in jet
            Double_t nSigmaTPC_proton = fPIDResponse -> NumberOfSigmasTPC(jetParticle, AliPID::kProton);
            Double_t nSigmaTOF_proton = fPIDResponse -> NumberOfSigmasTOF(jetParticle, AliPID::kProton);
            Double_t nSigmaITS_proton = fPIDResponse -> NumberOfSigmasITS(jetParticle, AliPID::kProton);
            hTPCnsigmaProton -> Fill((jetParticle->Pt()*jetParticle->Charge()),nSigmaTPC_proton);
            hTOFnsigmaProton -> Fill((jetParticle->Pt()*jetParticle->Charge()),nSigmaTOF_proton);
            hITSnsigmaProton -> Fill((jetParticle->Pt()*jetParticle->Charge()),nSigmaITS_proton);
            if (IsProton(jetParticle)) {
                particleType = 1;
                hTrackProtocol -> Fill(6.5); // # protons
                hPtJetProtonDeuteron -> Fill(1, jetParticlePt);
                jetProtons.emplace_back(jetParticle);
            } else {
                particleType = 2;
                hTrackProtocol -> Fill(7.5); // # antiprotons
                hPtJetProtonDeuteron -> Fill(2, jetParticlePt);
                jetAntiprotons.emplace_back(jetParticle);
            }
        } else if (IsDeuteron(jetParticle) || IsAntideuteron(jetParticle)) {// collect (anti)deuterons in jet
            Double_t nSigmaTPC_deuteron = fPIDResponse -> NumberOfSigmasTPC(jetParticle, AliPID::kDeuteron);
            Double_t nSigmaTOF_deuteron = fPIDResponse -> NumberOfSigmasTOF(jetParticle, AliPID::kDeuteron);
            Double_t nSigmaITS_deuteron = fPIDResponse -> NumberOfSigmasITS(jetParticle, AliPID::kDeuteron);
            hTPCnsigmaDeuteron -> Fill((jetParticle->Pt()*jetParticle->Charge()),nSigmaTPC_deuteron);
            hTOFnsigmaDeuteron -> Fill((jetParticle->Pt()*jetParticle->Charge()),nSigmaTOF_deuteron);
            hITSnsigmaDeuteron -> Fill((jetParticle->Pt()*jetParticle->Charge()),nSigmaITS_deuteron);
            if (IsDeuteron(jetParticle)) {
                particleType = 3;
                hTrackProtocol -> Fill(8.5); // # deuterons
                hPtJetProtonDeuteron -> Fill(3, jetParticlePt);
                jetDeuterons.emplace_back(jetParticle);
            } else {
                particleType = 4;
                hTrackProtocol -> Fill(9.5); // # antideuterons
                hPtJetProtonDeuteron -> Fill(4, jetParticlePt);
                jetAntideuterons.emplace_back(jetParticle);
            }
        }
    } //for (Int_t i=0; i<(Int_t)constituents.size();i++)

    if ((jetProtons.size()<2) && (jetAntiprotons.size()<2) && (jetDeuterons.size()<2) && (jetAntideuterons.size()<2)) return;
    hEventProtocol -> Fill(5.5); // # events with more than 1 (anti)proton/(anti)deuteron in a jet

    if (jetProtons.size()>1) {
        for (int i=0; i<(int)jetProtons.size(); i++) {
            for (int j=i+1; j<(int)jetProtons.size(); j++) {
                Double_t DeltaPhi = TVector2::Phi_0_2pi(jetProtons[i]->Phi() - jetProtons[j]->Phi());
                Double_t DeltaEta = jetProtons[i] - jetProtons[j];
                if (DeltaPhi>(1.5*TMath::Pi())) {
                    DeltaPhi = DeltaPhi - 2*TMath::Pi();
                }
                hDeltaPhiSE -> Fill(1,DeltaPhi);
                hDeltaPhiEtaSE -> Fill(1,DeltaPhi,DeltaEta);
            }
        }
    }
    if (jetAntiprotons.size()>1) {
        for (int i=0; i<(int)jetAntiprotons.size(); i++) {
            for (int j=i+1; j<(int)jetAntiprotons.size(); j++) {
                Double_t DeltaPhi = TVector2::Phi_0_2pi(jetAntiprotons[i]->Phi() - jetAntiprotons[j]->Phi());
                Double_t DeltaEta = jetAntiprotons[i] - jetAntiprotons[j];
                if (DeltaPhi>(1.5*TMath::Pi())) {
                    DeltaPhi = DeltaPhi - 2*TMath::Pi();
                }
                hDeltaPhiSE -> Fill(2,DeltaPhi);
                hDeltaPhiEtaSE -> Fill(2,DeltaPhi,DeltaEta);
            }
        }
    }
    if (jetDeuterons.size()>1) {
        for (int i=0; i<(int)jetDeuterons.size(); i++) {
            for (int j=i+1; j<(int)jetDeuterons.size(); j++) {
                Double_t DeltaPhi = TVector2::Phi_0_2pi(jetDeuterons[i]->Phi() - jetDeuterons[j]->Phi());
                Double_t DeltaEta = jetDeuterons[i] - jetDeuterons[j];
                if (DeltaPhi>(1.5*TMath::Pi())) {
                    DeltaPhi = DeltaPhi - 2*TMath::Pi();
                }
                hDeltaPhiSE -> Fill(3,DeltaPhi);
                hDeltaPhiEtaSE -> Fill(3,DeltaPhi,DeltaEta);
            }
        }
    }
    if (jetAntideuterons.size()>1) {
        for (int i=0; i<(int)jetAntideuterons.size(); i++) {
            for (int j=i+1; j<(int)jetAntideuterons.size(); j++) {
                Double_t DeltaPhi = TVector2::Phi_0_2pi(jetAntideuterons[i]->Phi() - jetAntideuterons[j]->Phi());
                Double_t DeltaEta = jetAntideuterons[i] - jetAntideuterons[j];
                if (DeltaPhi>(1.5*TMath::Pi())) {
                    DeltaPhi = DeltaPhi - 2*TMath::Pi();
                }
                hDeltaPhiSE -> Fill(4,DeltaPhi);
                hDeltaPhiEtaSE -> Fill(4,DeltaPhi,DeltaEta);
            }
        }
    }


    //Check
    if (fAODEvent->GetNumberOfTracks()==0) hEventProtocol->Fill(18.5); // # events with no tracks

    //Post Output Data
    PostData(1, fOutputData);
    PostData(2, fQAList);
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetAngCorrelations::GetEvent ()  {
    //Get Input Event
    fAODEvent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODEvent) return kFALSE;
    hNumberOfEvents -> Fill(0.5); // # events
    hEventProtocol -> Fill(0.5);

    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    return kTRUE;
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskJetAngCorrelations::GetDCAtoPrimaryVertex (AliAODTrack *track, Int_t index)  {

    Double_t dca[2]; // 0 = DCAxy, 1 = DCAz
    Double_t covMatrix[3];
    if (!track->PropagateToDCA (fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),10000,dca,covMatrix)) return -999;

    return dca[index];
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetAngCorrelations::IsProton(AliAODTrack *track) {
    //Initialisation
    Bool_t isProton=(kFALSE);
    if (track->Charge()<0) return isProton;

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
    Double_t pt = track->Pt();
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));

    if (DCAxy>fProtonDCAxy) return isProton;
    if (DCAz>fProtonDCAz)   return isProton;

    //Selection
    if (pt<fProtonITSTOFpT && TMath::Abs(nsigmaITS)<fProtonITSnsig && TMath::Abs(nsigmaTPC)<fProtonTPCnsigITS) isProton = kTRUE;
    if (pt>fProtonITSTOFpT && TMath::Abs(nsigmaTPC)<fProtonTPCnsigTOF && TMath::Abs(nsigmaTOF)<fProtonTOFnsig) isProton = kTRUE;

    return isProton;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetAngCorrelations::IsAntiproton(AliAODTrack *track) {
    //Initialisation
    Bool_t isAntiproton=(kFALSE);
    if (track->Charge()>0) return isAntiproton;

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
    Double_t pt = track->Pt();
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));

    if (DCAxy>fAntiprotonDCAxy) return isAntiproton;
    if (DCAz>fAntiprotonDCAz)   return isAntiproton;

    //Selection
    if (pt<fAntiprotonITSTOFpT && TMath::Abs(nsigmaITS)<fAntiprotonITSnsig && TMath::Abs(nsigmaTPC)<fAntiprotonTPCnsigITS) isAntiproton = kTRUE;
    if (pt>fAntiprotonITSTOFpT && TMath::Abs(nsigmaTPC)<fAntiprotonTPCnsigTOF && TMath::Abs(nsigmaTOF)<fAntiprotonTOFnsig) isAntiproton = kTRUE;

    return isAntiproton;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetAngCorrelations::IsDeuteron(AliAODTrack *track) {
    Bool_t isDeuteron=(kFALSE);
    if (track->Charge()<0) return isDeuteron;

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kDeuteron);
    Double_t pt = track->Pt();
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));

    if (DCAxy>fDeuteronDCAxy) return isDeuteron;
    if (DCAz>fDeuteronDCAz)   return isDeuteron;

    //Selection
    if (pt<fDeuteronITSTOFpT && TMath::Abs(nsigmaITS)<fDeuteronITSnsig && TMath::Abs(nsigmaTPC)<fDeuteronTPCnsigITS) isDeuteron = kTRUE;
    if (pt>fDeuteronITSTOFpT && TMath::Abs(nsigmaTPC)<fDeuteronTPCnsigTOF && TMath::Abs(nsigmaTOF)<fDeuteronTOFnsig) isDeuteron = kTRUE;

    return isDeuteron;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetAngCorrelations::IsAntideuteron(AliAODTrack *track) {
    Bool_t isAntideuteron=(kFALSE);
    if (track->Charge()>0) return isAntideuteron;

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kDeuteron);
    Double_t pt = track->Pt();
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));

    if (DCAxy>fAntideuteronDCAxy) return isAntideuteron;
    if (DCAz>fAntideuteronDCAz)   return isAntideuteron;

    //Selection
    if (pt<fAntideuteronITSTOFpT && TMath::Abs(nsigmaITS)<fAntideuteronITSnsig && TMath::Abs(nsigmaTPC)<fAntideuteronTPCnsigITS) isAntideuteron = kTRUE;
    if (pt>fAntideuteronITSTOFpT && TMath::Abs(nsigmaTPC)<fAntideuteronTPCnsigTOF && TMath::Abs(nsigmaTOF)<fAntideuteronTOFnsig) isAntideuteron = kTRUE;

    return isAntideuteron;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskJetAngCorrelations::PassedTrackSelection(AliAODTrack *track)  {

    //Initialisation
    Bool_t passedTrkSelection=(kFALSE);

    //Basic Track Selection
    if (!fAODTrackCuts->IsTrackAccepted(track)) return passedTrkSelection;

    //Fixed Cuts
    if (!track->HasPointOnITSLayer(0)) return passedTrkSelection;
    if (!track->HasPointOnITSLayer(1)) return passedTrkSelection;

    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));

    if (DCAxy>fDCAxy) return passedTrkSelection;
    if (DCAz>fDCAz)   return passedTrkSelection;

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetAngCorrelations::SetTrackBuffer(AliAODTrack *track, int particle) {
    switch (particle) {
        case 1:
            if ((int)fTrackBufferProton.size()==fTrackBufferSize) {
                fTrackBufferProton.insert(fTrackBufferProton.begin(), track);
                fTrackBufferProton.resize(fTrackBufferSize);
            } else if ((int)fTrackBufferProton.size()<fTrackBufferSize) {
                fTrackBufferProton.emplace_back(track);
            }
        break;
        case 2:
            if ((int)fTrackBufferAntiproton.size()==fTrackBufferSize) {
                fTrackBufferAntiproton.insert(fTrackBufferAntiproton.begin(), track);
                fTrackBufferAntiproton.resize(fTrackBufferSize);
            } else if ((int)fTrackBufferAntiproton.size()<fTrackBufferSize) {
                fTrackBufferAntiproton.emplace_back(track);
            }
        break;
        case 3:
            if ((int)fTrackBufferDeuteron.size()==fTrackBufferSize) {
                fTrackBufferDeuteron.insert(fTrackBufferDeuteron.begin(), track);
                fTrackBufferDeuteron.resize(fTrackBufferSize);
            } else if ((int)fTrackBufferDeuteron.size()<fTrackBufferSize) {
                fTrackBufferDeuteron.emplace_back(track);
            }
        break;
        case 4:
            if ((int)fTrackBufferAntideuteron.size()==fTrackBufferSize) {
                fTrackBufferAntideuteron.insert(fTrackBufferAntideuteron.begin(), track);
                fTrackBufferAntideuteron.resize(fTrackBufferSize);
            } else if ((int)fTrackBufferAntideuteron.size()<fTrackBufferSize) {
                fTrackBufferAntideuteron.emplace_back(track);
            }
        break;
        default:
        cout << "JetAngCorrelations: invalid particle ID!" << endl;
    }
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskJetAngCorrelations::Terminate(Option_t *)  {

    fOutputData = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputData) return;
}
//__________________________________________________________________________________________________________________________________
