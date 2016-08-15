#include "AliAnalysisTask_Syst_PtDistributionsMC.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "AliESDVertex.h"
#include "TDatabasePDG.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "TParticle.h"
#include "TVector3.h"
#include "AliStack.h"
#include "TRandom.h"
#include "TChain.h"
#include "TList.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliAnalysisTask_Syst_PtDistributionsMC)

//_________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask_Syst_PtDistributionsMC::AliAnalysisTask_Syst_PtDistributionsMC():
AliAnalysisTaskSE(),
fOutputList(0x0),
fESDevent(0x0),
fMCEvent(0x0),
fStack(0x0),
fPIDResponse(0x0),
fESDTrackCuts_Std(0x0),
fESDTrackCuts_Loose(0x0),
fCentralityMin(0),
fCentralityMax(0),
fDCAxy_param0(0),
fDCAxy_param1(0),
fDCAxy_param2(0),
fDCAz_max(0),
fITS_minNcls(0),
fTPC_minNcls(0),
fTPC_nClsdEdx(0),
fTPC_minCr(0),
fMinCrOverFindableCls(0),
fMaxGoldenChi2(0),
fMaxTPCchi2(0),
fMaxITSchi2(0),
fMaxFracSharedCls(0),
fITSreq(0),
fnsigmaTOF_max(0),
fnsigmaITS_max(0),
fnsigmaTPC_min(0),
fnsigmaTPC_max(0),
fMassMin(0),
fMassMax(0),
fPhivLim(0),
fPtMin(0),
fPtMax(0),
fEtaLim(0),
fHisto_Hijing_PizeroWeight(0x0),
fHisto_Hijing_EtaWeight(0x0),
fHisto_Hijing_EtaPrimeWeight(0x0),
fHisto_Hijing_RhoWeight(0x0),
fHisto_Hijing_OmegaWeight(0x0),
fHisto_Hijing_PhiWeight(0x0),
fHistoCentralityBins(0x0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask_Syst_PtDistributionsMC::AliAnalysisTask_Syst_PtDistributionsMC(const char *name):
AliAnalysisTaskSE(name),
fOutputList(0x0),
fESDevent(0x0),
fMCEvent(0x0),
fStack(0x0),
fPIDResponse(0x0),
fESDTrackCuts_Std(0x0),
fESDTrackCuts_Loose(0x0),
fCentralityMin(0),
fCentralityMax(0),
fDCAxy_param0(0),
fDCAxy_param1(0),
fDCAxy_param2(0),
fDCAz_max(0),
fITS_minNcls(0),
fTPC_minNcls(0),
fTPC_nClsdEdx(0),
fTPC_minCr(0),
fMinCrOverFindableCls(0),
fMaxGoldenChi2(0),
fMaxTPCchi2(0),
fMaxITSchi2(0),
fMaxFracSharedCls(0),
fITSreq(0),
fnsigmaTOF_max(0),
fnsigmaITS_max(0),
fnsigmaTPC_min(0),
fnsigmaTPC_max(0),
fMassMin(0),
fMassMax(0),
fPhivLim(0),
fPtMin(0),
fPtMax(0),
fEtaLim(0),
fHisto_Hijing_PizeroWeight(0x0),
fHisto_Hijing_EtaWeight(0x0),
fHisto_Hijing_EtaPrimeWeight(0x0),
fHisto_Hijing_RhoWeight(0x0),
fHisto_Hijing_OmegaWeight(0x0),
fHisto_Hijing_PhiWeight(0x0),
fHistoCentralityBins(0x0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTask_Syst_PtDistributionsMC::~AliAnalysisTask_Syst_PtDistributionsMC()  {
    
    fOutputList->Clear();
    delete fOutputList;
    delete fESDevent;
    delete fMCEvent;
    delete fStack;
    delete fPIDResponse;
    delete fESDTrackCuts_Std;
    delete fESDTrackCuts_Loose;
    delete fHisto_Hijing_PizeroWeight;
    delete fHisto_Hijing_EtaWeight;
    delete fHisto_Hijing_EtaPrimeWeight;
    delete fHisto_Hijing_RhoWeight;
    delete fHisto_Hijing_OmegaWeight;
    delete fHisto_Hijing_PhiWeight;
    delete fHistoCentralityBins;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTask_Syst_PtDistributionsMC::UserCreateOutputObjects()  {
    
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
    //PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    //Track Quality Cuts
    fESDTrackCuts_Std   = new AliESDtrackCuts ("fESDTrackCuts_Std");
    fESDTrackCuts_Loose = new AliESDtrackCuts ("fESDTrackCuts_Loose");

    
    //Statistics
    fHistoEvents = new TH1F ("fHistoEvents","",2,0,2);
    fOutputList -> Add(fHistoEvents);
    
    
    //Pt distributions
    fHistoPtDistribution_Electrons      = new TH1F ("fHistoPtDistribution_Electrons","",5000,0,5);
    fHistoPtDistribution_Positrons      = new TH1F ("fHistoPtDistribution_Positrons","",5000,0,5);
    fHistoPtDistribution_Electrons_Pref = new TH1F ("fHistoPtDistribution_Electrons_Pref","",5000,0,5);
    fHistoPtDistribution_Positrons_Pref = new TH1F ("fHistoPtDistribution_Positrons_Pref","",5000,0,5);

    fOutputList -> Add (fHistoPtDistribution_Electrons);
    fOutputList -> Add (fHistoPtDistribution_Positrons);
    fOutputList -> Add (fHistoPtDistribution_Electrons_Pref);
    fOutputList -> Add (fHistoPtDistribution_Positrons_Pref);

    
    PostData(1, fOutputList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTask_Syst_PtDistributionsMC::UserExec(Option_t *)  {
    
    
    //Get Event
    if ( !GetEvent()) return;
    
    
    //Loop over esdtracks
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {
        
        //Track Selection
        AliESDtrack *esdtrack1 = (AliESDtrack*) fESDevent -> GetTrack(i);
        if ( !esdtrack1 ) continue;
        if ( !IsTrackFromHijing (esdtrack1)) continue;
        if ( !PassedTrackQualityCuts (esdtrack1)) continue;
        if ( !PassedPIDcuts (esdtrack1)) continue;
        
        //MC Particle
        TParticle *particle = (TParticle*) fStack->Particle (TMath::Abs(esdtrack1->GetLabel()));
        if ( !particle ) continue;
        if ( TMath::Abs(particle->GetPdgCode()) != 11 ) continue;//electron selection

        
        //pt of the track
        Double_t pt = esdtrack1->Pt();

        //Pt Distributions (before pre-filter)
        if (esdtrack1->Charge()<0) fHistoPtDistribution_Electrons -> Fill(pt,Weight(particle));
        if (esdtrack1->Charge()>0) fHistoPtDistribution_Positrons -> Fill(pt,Weight(particle));
        
        
        //Pre-filtering
        Bool_t IsConversionCandidate = false;
        for (Int_t j=0 ; j<fESDevent->GetNumberOfTracks() ; j++)  {
            
            if (i==j) continue;//Skip Self-pairing
            
            AliESDtrack *esdtrack2 = (AliESDtrack*) fESDevent -> GetTrack(j);
            if ( !esdtrack2 ) continue;
            if ( !IsTrackFromHijing(esdtrack2)) continue;
            if ( !PassedLooseTrackQualityCuts (esdtrack2)) continue;
            
            Double_t mass = GetMass (esdtrack1,esdtrack2);
            Double_t phiV = GetPhiV (esdtrack1,esdtrack2);
            
            //Pre-filtering Cuts
            if ( mass>=fMassMin && mass<=fMassMax && phiV<=fPhivLim )         { IsConversionCandidate = true; break; }
            if ( mass>=fMassMin && mass<=fMassMax && phiV>=(180.0-fPhivLim) ) { IsConversionCandidate = true; break; }
        }
        
        if (IsConversionCandidate) continue;//Skip Conversion Candidates
        
        //Pt Distributions (after pre-filter)
        if (esdtrack1->Charge()<0) fHistoPtDistribution_Electrons_Pref -> Fill(pt,Weight(particle));
        if (esdtrack1->Charge()>0) fHistoPtDistribution_Positrons_Pref -> Fill(pt,Weight(particle));
    }

    
    PostData(1, fOutputList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTask_Syst_PtDistributionsMC::GetEvent ()  {
    
    fESDevent = dynamic_cast <AliESDEvent*> (InputEvent());
    if (!fESDevent) return false;
    fMCEvent = MCEvent();
    if (!fMCEvent) return false;
    fStack = fMCEvent->Stack();
    if (!fStack) return false;
    
    //Total number of events
    fHistoEvents -> Fill(0.5);
    
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    if ( TMath::Abs(vertex->GetZ() ) > 10.0 ) return false;
    if ( vertex->GetNContributors() < 1 ) return false;
    AliCentrality *centrality = fESDevent->GetCentrality();
    if (!centrality) return false;
    Double_t centr = centrality->GetCentralityPercentile("V0M");
    if (centr<fCentralityMin || centr>=fCentralityMax ) return false;
    
    //Number of selected events
    fHistoEvents -> Fill(1.5);
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTask_Syst_PtDistributionsMC::Weight (TParticle *particle)  {
    
    //Centrality
    AliCentrality *centrality = fESDevent->GetCentrality();
    Double_t centralityPerc = centrality->GetCentralityPercentile("V0M");
    
    //Electron Parent
    TParticle *parent = (TParticle*) fStack->Particle(particle->GetMother(0));
    Int_t pdg = TMath::Abs(parent->GetPdgCode());
    
    //Find Bins
    Int_t nx = fHisto_Hijing_PizeroWeight->GetXaxis()->FindBin(parent->Pt());
    Int_t ny = fHistoCentralityBins->FindBin(centralityPerc);
    
    Double_t weight(1);
    switch (pdg) {
        case 22:  weight = WeightConversions (parent); break;
        case 111: weight = fHisto_Hijing_PizeroWeight   -> GetBinContent(nx,ny); break;
        case 221: weight = fHisto_Hijing_EtaWeight      -> GetBinContent(nx,ny); break;
        case 331: weight = fHisto_Hijing_EtaPrimeWeight -> GetBinContent(nx,ny); break;
        case 113: weight = fHisto_Hijing_RhoWeight      -> GetBinContent(nx,ny); break;
        case 223: weight = fHisto_Hijing_OmegaWeight    -> GetBinContent(nx,ny); break;
        case 333: weight = fHisto_Hijing_PhiWeight      -> GetBinContent(nx,ny); break;
            
        default: weight = 1; break;
    }
    
    return weight;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTask_Syst_PtDistributionsMC::WeightConversions (TParticle *parent)  {
    
    //Centrality
    AliCentrality *centrality = fESDevent->GetCentrality();
    Double_t centralityPerc = centrality->GetCentralityPercentile("V0M");
    
    //Photon Parent
    Int_t lm = parent->GetMother(0);
    if (lm<0) return 1;
    TParticle *gparent = (TParticle*) fStack->Particle(lm);
    Int_t pdg = TMath::Abs(gparent->GetPdgCode());
    
    //Find Bins
    Int_t nx = fHisto_Hijing_PizeroWeight->GetXaxis()->FindBin(gparent->Pt());
    Int_t ny = fHistoCentralityBins->FindBin(centralityPerc);
    
    Double_t weight(1);
    switch (pdg) {
        case 111: weight = fHisto_Hijing_PizeroWeight   -> GetBinContent(nx,ny); break;
        case 221: weight = fHisto_Hijing_EtaWeight      -> GetBinContent(nx,ny); break;
        case 331: weight = fHisto_Hijing_EtaPrimeWeight -> GetBinContent(nx,ny); break;
        case 113: weight = fHisto_Hijing_RhoWeight      -> GetBinContent(nx,ny); break;
        case 223: weight = fHisto_Hijing_OmegaWeight    -> GetBinContent(nx,ny); break;
        case 333: weight = fHisto_Hijing_PhiWeight      -> GetBinContent(nx,ny); break;
            
        default: weight = 1; break;
    }
    
    return weight;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTask_Syst_PtDistributionsMC::IsTrackFromHijing (AliESDtrack *track)  {
    
    Bool_t IsHijing = (false);
    
    Int_t lp = track->GetLabel();
    if ( lp==-1 )  return IsHijing;
    if ( fMCEvent->IsFromBGEvent(fStack->GetPrimary(TMath::Abs(lp))) ) IsHijing = true;
    
    return IsHijing;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTask_Syst_PtDistributionsMC::PassedTrackQualityCuts (AliESDtrack* track)  {
    
    fESDTrackCuts_Std -> SetMaxDCAToVertexZ(fDCAz_max);
    fESDTrackCuts_Std -> SetMaxDCAToVertexXY(fDCAxy_param0 + fDCAxy_param1/TMath::Power(track->Pt(),fDCAxy_param2));
    fESDTrackCuts_Std -> SetAcceptKinkDaughters(false);
    fESDTrackCuts_Std -> SetMinNCrossedRowsTPC(fTPC_minCr);
    fESDTrackCuts_Std -> SetMinNClustersITS(fITS_minNcls);
    fESDTrackCuts_Std -> SetMinNClustersTPC(fTPC_minNcls);
    fESDTrackCuts_Std -> SetRequireTPCRefit(true);
    fESDTrackCuts_Std -> SetRequireITSRefit(true);
    
    //ITS Requirement
    if (strcmp(fITSreq,"kOff")==0)        fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
    if (strcmp(fITSreq,"kNone")==0)       fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
    if (strcmp(fITSreq,"kAny")==0)        fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    if (strcmp(fITSreq,"kFirst")==0)      fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    if (strcmp(fITSreq,"kOnlyFirst")==0)  fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOnlyFirst);
    if (strcmp(fITSreq,"kSecond")==0)     fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kSecond);
    if (strcmp(fITSreq,"kOnlySecond")==0) fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOnlySecond);
    if (strcmp(fITSreq,"kBoth")==0)       fESDTrackCuts_Std -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    
    fESDTrackCuts_Std -> SetMinRatioCrossedRowsOverFindableClustersTPC(fMinCrOverFindableCls);
    fESDTrackCuts_Std -> SetMaxChi2TPCConstrainedGlobal(fMaxGoldenChi2);
    fESDTrackCuts_Std -> SetMaxChi2PerClusterTPC(fMaxTPCchi2);
    fESDTrackCuts_Std -> SetMaxChi2PerClusterITS(fMaxITSchi2);
    fESDTrackCuts_Std -> SetEtaRange(-fEtaLim,fEtaLim);
    fESDTrackCuts_Std -> SetPtRange(fPtMin,fPtMax);
    
    if ( track -> GetTPCsignalN() < fTPC_nClsdEdx ) return false;
    if ( FractionSharedClsITS(track) > fMaxFracSharedCls ) return false;
    if ( fESDTrackCuts_Std->AcceptTrack (track) ) return true;
    return false;
}
//________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTask_Syst_PtDistributionsMC::FractionSharedClsITS (AliESDtrack *track)  {
    
    Double_t nSharedCls(0);
    Double_t fSharedCls(0);
    
    for ( Int_t i=0 ; i<6 ; i++ )
        if ( track->HasPointOnITSLayer(i) && track->HasSharedPointOnITSLayer(i) ) nSharedCls++;
    
    fSharedCls = nSharedCls/(Double_t)track->GetNcls(0);
    
    return fSharedCls;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTask_Syst_PtDistributionsMC::PassedLooseTrackQualityCuts (AliESDtrack *track)  {
        
    fESDTrackCuts_Loose -> SetMinNCrossedRowsTPC(70);
    fESDTrackCuts_Loose -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
    fESDTrackCuts_Loose -> SetMaxChi2PerClusterTPC(4);
    fESDTrackCuts_Loose -> SetAcceptKinkDaughters(false);
    fESDTrackCuts_Loose -> SetRequireTPCRefit(true);
    fESDTrackCuts_Loose -> SetRequireITSRefit(true);
    fESDTrackCuts_Loose -> SetMinNClustersITS(3);
    fESDTrackCuts_Loose -> SetMinNClustersTPC(50);
    fESDTrackCuts_Loose -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    fESDTrackCuts_Loose -> SetMaxDCAToVertexXY(1.0);
    fESDTrackCuts_Loose -> SetMaxDCAToVertexZ(3.0);
    fESDTrackCuts_Loose -> SetDCAToVertex2D(false);
    fESDTrackCuts_Loose -> SetRequireSigmaToVertex(false);
    fESDTrackCuts_Loose -> SetMaxChi2PerClusterITS(36);
    
    if ( fESDTrackCuts_Loose->AcceptTrack (track) ) return true;
    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTask_Syst_PtDistributionsMC::PassedPIDcuts (AliESDtrack *track)  {
    
    //TOF Response
    if ( !fPIDResponse -> CheckPIDStatus(AliPIDResponse::kTOF,track)) return false;
    
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kElectron);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kElectron);
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kElectron);
    
    if ( TMath::Abs(nsigmaTOF ) > fnsigmaTOF_max )            return false;
    if ( nsigmaITS > fnsigmaITS_max )                         return false;
    if ( nsigmaTPC < fnsigmaTPC_min*TMath::Exp(-track->P()) ) return false;
    if ( nsigmaTPC > fnsigmaTPC_max )                         return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTask_Syst_PtDistributionsMC::GetMass (AliESDtrack *track1,AliESDtrack *track2)  {
    
    TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
    TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
    
    Double_t electron_mass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    
    Double_t E1 = TMath::Sqrt( TMath::Power(electron_mass,2) + P1.Mag2());
    Double_t E2 = TMath::Sqrt( TMath::Power(electron_mass,2) + P2.Mag2());
    Double_t mass = TMath::Sqrt( TMath::Power(E1+E2,2) - (P1+P2).Mag2());
    
    return mass;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTask_Syst_PtDistributionsMC::GetPhiV (AliESDtrack *track1,AliESDtrack *track2 ) {
    
    TVector3 P1,P2;
    Double_t r(0);
    
    //Randomization of pair ordering (symmetric peaks)
    do { r = gRandom -> Uniform (0.0,1.0); } while (r==0.5);
    
    if (r < 0.5) { P1.SetXYZ (track1->Px(),track1->Py(),track1->Pz()); P2.SetXYZ (track2->Px(),track2->Py(),track2->Pz()); }
    if (r > 0.5) { P1.SetXYZ (track2->Px(),track2->Py(),track2->Pz()); P2.SetXYZ (track1->Px(),track1->Py(),track1->Pz()); }
    
    //PhiV Calculation
    TVector3 P = P1 + P2;
    TVector3 U ( P.X()/P.Mag(),P.Y()/P.Mag(),P.Z()/P.Mag() );
    TVector3 A ( U.Y()/TMath::Sqrt(U.X()*U.X()+U.Y()*U.Y()),-U.X()/TMath::Sqrt(U.X()*U.X()+ U.Y()*U.Y()),0 );
    TVector3 Vp = P1.Cross(P2);
    TVector3 V (Vp.X()/Vp.Mag(),Vp.Y()/Vp.Mag(),Vp.Z()/Vp.Mag());
    TVector3 W = U.Cross(V);
    
    return (180.0/TMath::Pi())*A.Angle(W);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTask_Syst_PtDistributionsMC::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________

