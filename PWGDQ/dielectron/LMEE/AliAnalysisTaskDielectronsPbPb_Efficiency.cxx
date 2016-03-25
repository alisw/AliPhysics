#include "AliAnalysisTaskDielectronsPbPb_Efficiency.h"
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

ClassImp(AliAnalysisTaskDielectronsPbPb_Efficiency)

//_________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDielectronsPbPb_Efficiency::AliAnalysisTaskDielectronsPbPb_Efficiency():
AliAnalysisTaskSE(),
fOutputList(0x0),
fESDevent(0x0),
fMCEvent(0x0),
fStack(0x0),
fPIDResponse(0x0),
fESDTrackCuts(0x0),
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
fnsigmaTOF_max(0),
fnsigmaITS_max(0),
fnsigmaTPC_min(0),
fnsigmaTPC_max(0),
fMassLim(0),
fPhivLim(0),
fHistoDetResponseMatrix_Momentum(0x0),
fHistoDetResponseMatrix_Theta(0x0),
fHistoDetResponseMatrix_Phi_Electrons(0x0),
fHistoDetResponseMatrix_Phi_Positrons(0x0),
fHisto_Hijing_PizeroWeight(0x0),
fHisto_Hijing_EtaWeight(0x0),
fHisto_Hijing_EtaPrimeWeight(0x0),
fHisto_Hijing_RhoWeight(0x0),
fHisto_Hijing_OmegaWeight(0x0),
fHisto_Hijing_PhiWeight(0x0),
fHistoCentralityBins(0x0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDielectronsPbPb_Efficiency::AliAnalysisTaskDielectronsPbPb_Efficiency(const char *name):
AliAnalysisTaskSE(name),
fOutputList(0x0),
fESDevent(0x0),
fMCEvent(0x0),
fStack(0x0),
fPIDResponse(0x0),
fESDTrackCuts(0x0),
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
fnsigmaTOF_max(0),
fnsigmaITS_max(0),
fnsigmaTPC_min(0),
fnsigmaTPC_max(0),
fMassLim(0),
fPhivLim(0),
fHistoDetResponseMatrix_Momentum(0x0),
fHistoDetResponseMatrix_Theta(0x0),
fHistoDetResponseMatrix_Phi_Electrons(0x0),
fHistoDetResponseMatrix_Phi_Positrons(0x0),
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
AliAnalysisTaskDielectronsPbPb_Efficiency::~AliAnalysisTaskDielectronsPbPb_Efficiency()  {
    
    fOutputList->Clear();
    delete fOutputList;
    delete fESDevent;
    delete fMCEvent;
    delete fStack;
    delete fPIDResponse;
    delete fESDTrackCuts;
    delete fHistoDetResponseMatrix_Momentum;
    delete fHistoDetResponseMatrix_Theta;
    delete fHistoDetResponseMatrix_Phi_Electrons;
    delete fHistoDetResponseMatrix_Phi_Positrons;
    delete fHisto_Hijing_PizeroWeight;
    delete fHisto_Hijing_EtaWeight;
    delete fHisto_Hijing_EtaPrimeWeight;
    delete fHisto_Hijing_RhoWeight;
    delete fHisto_Hijing_OmegaWeight;
    delete fHisto_Hijing_PhiWeight;
    delete fHistoCentralityBins;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDielectronsPbPb_Efficiency::UserCreateOutputObjects()  {
    
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
    //PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    fESDTrackCuts = new AliESDtrackCuts ("AliESDtrackCuts");

    
    //Statistics
    fHistoEvents = new TH1F ("fHistoEvents","",2,0,2);
    fOutputList -> Add(fHistoEvents);
    
   
    //Pair Efficiency
    fHistoInvMass_Gen =                new TH2F ("fHistoInvMass_Gen","",500,0,5,20,0,10);
    fHistoInvMass_Rec =                new TH2F ("fHistoInvMass_Rec","",500,0,5,20,0,10);
    fHistoInvMass_Rec_Pref =           new TH2F ("fHistoInvMass_Rec_Pref","",500,0,5,20,0,10);
    fHistoInvMass_Gen_noweights =      new TH2F ("fHistoInvMass_Gen_noweights","",500,0,5,20,0,10);
    fHistoInvMass_Rec_noweights =      new TH2F ("fHistoInvMass_Rec_noweights","",500,0,5,20,0,10);
    fHistoInvMass_Rec_Pref_noweights = new TH2F ("fHistoInvMass_Rec_Pref_noweights","",500,0,5,20,0,10);

    fOutputList -> Add (fHistoInvMass_Gen);
    fOutputList -> Add (fHistoInvMass_Rec);
    fOutputList -> Add (fHistoInvMass_Rec_Pref);
    fOutputList -> Add (fHistoInvMass_Gen_noweights);
    fOutputList -> Add (fHistoInvMass_Rec_noweights);
    fOutputList -> Add (fHistoInvMass_Rec_Pref_noweights);

    
    //Conversions Residual Contribution
    fHistoInvariantMass_Dielectrons =           new TH2F ("fHistoInvariantMass_Dielectrons","",500,0,5,20,0,10);
    fHistoInvariantMass_Conversions =           new TH2F ("fHistoInvariantMass_Conversions","",500,0,5,20,0,10);
    fHistoInvariantMass_Dielectrons_noweights = new TH2F ("fHistoInvariantMass_Dielectrons_noweights","",500,0,5,20,0,10);
    fHistoInvariantMass_Conversions_noweights = new TH2F ("fHistoInvariantMass_Conversions_noweights","",500,0,5,20,0,10);
    
    fOutputList -> Add (fHistoInvariantMass_Dielectrons);
    fOutputList -> Add (fHistoInvariantMass_Conversions);
    fOutputList -> Add (fHistoInvariantMass_Dielectrons_noweights);
    fOutputList -> Add (fHistoInvariantMass_Conversions_noweights);

        
    PostData(1, fOutputList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDielectronsPbPb_Efficiency::UserExec(Option_t *)  {
    
    
    //Get Event
    if ( !GetEvent()) return;
    
    
    
    
    //Primary Electron Selection
    Int_t nParticles=0;
    TParticle *particle[500];
    
    for ( Int_t i=0 ; i<fStack->GetNtrack() ; i++ )  {
        
        TParticle *MCparticle = fStack->Particle(i);
        if ( !MCparticle ) continue;
        if ( !fMCEvent->IsFromBGEvent(fStack->GetPrimary(i))) continue;
        if ( !IsPrimaryElectron (MCparticle)) continue;
        particle[nParticles] = MCparticle;
        nParticles++;
    }

    
    //Generated Mass Distribution
    for ( Int_t i=0 ; i<nParticles ; i++ )  {
        
        Short_t q1 = -particle[i]->GetPdgCode()/11;
        TVector3 P1 = GetReconstructedMomentum (particle[i]->P(),particle[i]->Theta(),particle[i]->Phi(),q1);
        if ( P1.Pt()<0.4 || P1.Pt()>5.0 ) continue;
        if ( TMath::Abs(P1.Eta())>0.8 ) continue;
        
        for ( Int_t j=i+1 ; j<nParticles ; j++ )  {
            
            if ( !IsCorrelatedPair (particle[i],particle[j])) continue;//Correlated Pairs Selection
            
            Short_t q2 = -particle[j]->GetPdgCode()/11;
            TVector3 P2 = GetReconstructedMomentum (particle[j]->P(),particle[j]->Theta(),particle[j]->Phi(),q2);
            if ( P2.Pt()<0.4 || P2.Pt()>5.0 ) continue;
            if ( TMath::Abs(P2.Eta())>0.8 ) continue;
            
            //Pair Variables
            Double_t mass = GetMass (P1,P2);
            Double_t ptee = (P1+P2).Pt();
            
            //Generated Mass Distribution
            fHistoInvMass_Gen -> Fill(mass,ptee, Weight(particle[i]));
            fHistoInvMass_Gen_noweights -> Fill(mass,ptee);
        }
    }
     
    
    
    
    //Reconstructed Electron Selection
    Int_t nTracks=0;
    AliESDtrack *track[200];
    
    Int_t nPrimTracks=0;
    AliESDtrack *primtrack[200];

    Int_t nTracksLooseCuts=0;
    AliESDtrack *track_looseCuts[5000];
    
    
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {
        
        AliESDtrack *esdtrack = (AliESDtrack*) fESDevent->GetTrack(i);
        if ( !esdtrack ) continue;
        if ( !IsTrackFromHijing(esdtrack)) continue;
        if (  PassedLooseTrackQualityCuts(esdtrack))  { track_looseCuts[nTracksLooseCuts] = esdtrack; nTracksLooseCuts++; }
        if ( !PassedTrackQualityCuts(esdtrack)) continue;
        if ( !PassedPIDCuts(esdtrack)) continue;
       
        TParticle *MCparticle = (TParticle*) fStack->Particle(TMath::Abs(esdtrack->GetLabel()));
        if ( !MCparticle ) continue;
        if ( TMath::Abs(MCparticle->GetPdgCode()) != 11 ) continue;//Electron Selection
        
        track[nTracks] = esdtrack;
        nTracks++;

        //Primary Electron Selection
        if ( !IsPrimaryElectron (MCparticle)) continue;
        primtrack[nPrimTracks] = esdtrack;
        nPrimTracks++;
    }
    
    
    
    
    //Reconstructed Mass Distribution
    for ( Int_t i=0 ; i<nPrimTracks ; i++ )  {
        
        TParticle *particle1 = (TParticle*) fStack->Particle(TMath::Abs(primtrack[i]->GetLabel()));
        
        for ( Int_t j=i+1 ; j<nPrimTracks ; j++ )  {
            
            TParticle *particle2 = (TParticle*) fStack->Particle(TMath::Abs(primtrack[j]->GetLabel()));
            
            if ( !IsCorrelatedPair (particle1,particle2)) continue;//Correlated Pairs Selection
           
            //Pair Variables
            Double_t mass = GetMass (primtrack[i],primtrack[j]);
            Double_t ptee = GetPtee (primtrack[i],primtrack[j]);
            
            //Reconstructed Mass Distribution
            fHistoInvMass_Rec -> Fill(mass,ptee,Weight(particle1));
            fHistoInvMass_Rec_noweights -> Fill(mass,ptee);
        }
    }
    

    
    
    //Pre-Filtering
    Int_t nTracksPref=0;
    AliESDtrack *track_pref[200];
    
    for ( Int_t i=0 ; i<nTracks ; i++ )  {
        
        Bool_t IsConversionCandidate = false;
        for ( Int_t j=0 ; j<nTracksLooseCuts ; j++ )  {
            
            if (track[i]->GetLabel() == track_looseCuts[j]->GetLabel()) continue;//skip self-pairing
            
            Double_t mass = GetMass (track[i],track_looseCuts[j]);
            Double_t phiV = GetPhiV (track[i],track_looseCuts[j]);
            
            if ( mass<=fMassLim && phiV<=fPhivLim )         { IsConversionCandidate=true; break; }
            if ( mass<=fMassLim && phiV>=(180.0-fPhivLim) ) { IsConversionCandidate=true; break; }
        }
        
        if (IsConversionCandidate) continue;
        
        track_pref[nTracksPref] = track[i];
        nTracksPref++;
    }
    
    
    
    
    //Reconstructed Mass Distribution (After Pre-filter)
    for ( Int_t i=0 ; i<nTracksPref ; i++ )  {
        
        TParticle *particle1 = (TParticle*) fStack->Particle(TMath::Abs(track_pref[i]->GetLabel()));
        
        for ( Int_t j=i+1 ; j<nTracksPref ; j++ )  {
            
            TParticle *particle2 = (TParticle*) fStack->Particle(TMath::Abs(track_pref[j]->GetLabel()));

            if ( !IsCorrelatedPair (particle1,particle2)) continue;//Correlated Pairs Selection

            //Pair Variables
            Double_t mass = GetMass (track_pref[i],track_pref[j]);
            Double_t ptee = GetPtee (track_pref[i],track_pref[j]);
            Double_t phiV = GetPhiV (track_pref[i],track_pref[j]);
            
            
            //Conversions: Residual Contribution
            fHistoInvariantMass_Dielectrons -> Fill(mass,ptee,Weight(particle1));
            fHistoInvariantMass_Dielectrons_noweights -> Fill(mass,ptee);
            
            if ( IsFromConversion(particle1) && IsFromConversion(particle2)) {
                
                fHistoInvariantMass_Conversions -> Fill(mass,ptee,Weight(particle1));
                fHistoInvariantMass_Conversions_noweights -> Fill(mass,ptee);
            }
            
            //Primary Electrons
            if ( !IsPrimaryElectron(particle1)) continue;
            if ( !IsPrimaryElectron(particle2)) continue;
            
            
            //Reconstructed Mass Distribution
            fHistoInvMass_Rec_Pref -> Fill(mass,ptee,Weight(particle1));
            fHistoInvMass_Rec_Pref_noweights -> Fill(mass,ptee);
        }
    }

    
    PostData(1, fOutputList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::GetEvent ()  {
    
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
TVector3 AliAnalysisTaskDielectronsPbPb_Efficiency::GetReconstructedMomentum (Double_t p_gen,Double_t theta_gen,Double_t phi_gen,Short_t q)  {
    
    Int_t x1,x2;
    Double_t p_rec(0),theta_rec(0),phi_rec(0);
    
    //Momentum Smearing
    gRandom->SetSeed(0);
    x1 = fHistoDetResponseMatrix_Momentum->GetXaxis()->FindBin(p_gen);
    x2 = fHistoDetResponseMatrix_Momentum->GetXaxis()->FindBin(p_gen);
    TH1F *fHistoDeltaP = (TH1F*) fHistoDetResponseMatrix_Momentum->ProjectionY("fHistoDeltaP",x1,x2,"E");
    p_rec = p_gen + fHistoDeltaP->GetRandom();
    
    //Theta Smearing
    gRandom->SetSeed(0);
    x1 = fHistoDetResponseMatrix_Theta->GetXaxis()->FindBin(p_gen);
    x2 = fHistoDetResponseMatrix_Theta->GetXaxis()->FindBin(p_gen);
    TH1F *fHistoDeltaTheta = (TH1F*) fHistoDetResponseMatrix_Theta->ProjectionY("fHistoDeltaTheta",x1,x2,"E");
    theta_rec = theta_gen + fHistoDeltaTheta->GetRandom();
    
    //Phi Smearing (Electrons)
    if (q<0)  {
        
        gRandom->SetSeed(0);
        x1 = fHistoDetResponseMatrix_Phi_Electrons->GetXaxis()->FindBin(p_gen);
        x2 = fHistoDetResponseMatrix_Phi_Electrons->GetXaxis()->FindBin(p_gen);
        TH1F *fHistoDeltaPhi_Electrons = (TH1F*) fHistoDetResponseMatrix_Phi_Electrons->ProjectionY("fHistoDeltaPhi_Electrons",x1,x2,"E");
        phi_rec = phi_gen + fHistoDeltaPhi_Electrons->GetRandom();
    }
    
    //Phi Smearing (Positrons)
    if (q>0)  {
        
        gRandom->SetSeed(0);
        x1 = fHistoDetResponseMatrix_Phi_Positrons->GetXaxis()->FindBin(p_gen);
        x2 = fHistoDetResponseMatrix_Phi_Positrons->GetXaxis()->FindBin(p_gen);
        TH1F *fHistoDeltaPhi_Positrons = (TH1F*) fHistoDetResponseMatrix_Phi_Positrons->ProjectionY("fHistoDeltaPhi_Positrons",x1,x2,"E");
        phi_rec = phi_gen + fHistoDeltaPhi_Positrons->GetRandom();
    }
    
    //Momentum Vector
    TVector3 Prec;
    Prec.SetX (p_rec*TMath::Sin(theta_rec)*TMath::Cos(phi_rec));
    Prec.SetY (p_rec*TMath::Sin(theta_rec)*TMath::Sin(phi_rec));
    Prec.SetZ (p_rec*TMath::Cos(theta_rec));
    
    return Prec;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::IsPrimaryElectron (TParticle *particle)  {
    
    //Electron Selection
    if ( TMath::Abs(particle->GetPdgCode()) != 11 ) return false;
    
    //Electron Parent
    Int_t lm = particle->GetMother(0);
    if (lm<0) return false;
    TParticle *parent = (TParticle*) fStack->Particle(lm);
    if (!parent) return false;
    
    //Electron Source: LF or HF
    if ( IsLightFlavorParticle(parent) || IsHeavyFlavorParticle(parent) ) return true;
    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::IsLightFlavorParticle (TParticle *parent)  {
    
    Bool_t IsPrimaryPart = (false);
    Bool_t IsFromCocktail = (false);
   
    //Primary Particle
    Double_t xv = fMCEvent->GetPrimaryVertex()->GetX();
    Double_t yv = fMCEvent->GetPrimaryVertex()->GetY();
    Double_t xp = parent->Vx();
    Double_t yp = parent->Vy();
    Double_t Dr = TMath::Sqrt(  (xv-xp)*(xv-xp) + (yv-yp)*(yv-yp) );
    if (Dr < 1.0e-15) IsPrimaryPart=true;
    
    //Cocktail Particle
    Int_t pdgMeson[] = { 111,221,331,113,223,333 };
    const Int_t nPart = sizeof(pdgMeson)/sizeof(Int_t);
    for ( Int_t i=0 ; i<nPart ; i++ )
        if (  TMath::Abs(parent->GetPdgCode()) == pdgMeson[i] ) { IsFromCocktail=true; break; }
    
    if ( IsPrimaryPart && IsFromCocktail ) return true;
    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::IsHeavyFlavorParticle (TParticle *parent)  {
    
    Bool_t IsFromHF = (false);
    Int_t pdg = TMath::Abs(parent->GetPdgCode());
    
    switch (pdg) {
            
        //Charmed Mesons
        case 411:   IsFromHF=true; break;
        case 421:   IsFromHF=true; break;
        case 10411: IsFromHF=true; break;
        case 10421: IsFromHF=true; break;
        case 413:   IsFromHF=true; break;
        case 423:   IsFromHF=true; break;
        case 10413: IsFromHF=true; break;
        case 10423: IsFromHF=true; break;
        case 20413: IsFromHF=true; break;
        case 20423: IsFromHF=true; break;
        case 415:   IsFromHF=true; break;
        case 425:   IsFromHF=true; break;
        case 431:   IsFromHF=true; break;
        case 10431: IsFromHF=true; break;
        case 433:   IsFromHF=true; break;
        case 10433: IsFromHF=true; break;
        case 20433: IsFromHF=true; break;
        case 435:   IsFromHF=true; break;
            
        //Beauty Mesons
        case 511:   IsFromHF=true; break;
        case 521:   IsFromHF=true; break;
        case 10511: IsFromHF=true; break;
        case 10521: IsFromHF=true; break;
        case 513:   IsFromHF=true; break;
        case 523:   IsFromHF=true; break;
        case 10513: IsFromHF=true; break;
        case 10523: IsFromHF=true; break;
        case 20513: IsFromHF=true; break;
        case 20523: IsFromHF=true; break;
        case 515:   IsFromHF=true; break;
        case 525:   IsFromHF=true; break;
        case 531:   IsFromHF=true; break;
        case 10531: IsFromHF=true; break;
        case 533:   IsFromHF=true; break;
        case 10533: IsFromHF=true; break;
        case 20533: IsFromHF=true; break;
        case 535:   IsFromHF=true; break;
        case 541:   IsFromHF=true; break;
        case 10541: IsFromHF=true; break;
        case 543:   IsFromHF=true; break;
        case 10543: IsFromHF=true; break;
        case 20543: IsFromHF=true; break;
        case 545:   IsFromHF=true; break;
            
        default: IsFromHF=false; break;
    }
   
    return IsFromHF;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::IsFromConversion (TParticle *particle)  {
    
    Bool_t IsFromPhotonConv = (false);
    
    TParticle *parent = (TParticle*) fStack->Particle (particle->GetMother(0));
    if (TMath::Abs(parent->GetPdgCode()) == 22) IsFromPhotonConv = true;
    
    return IsFromPhotonConv;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::IsCorrelatedPair (TParticle *particle1,TParticle *particle2)  {
    
    //Skip Like-Sign Pairs
    if ( particle1->GetPdgCode() == particle2->GetPdgCode() ) return false;
    
    //Parents
    Int_t lm1 = particle1->GetMother(0);
    Int_t lm2 = particle2->GetMother(0);
    TParticle *parent1 = (TParticle*) fStack->Particle(lm1);
    TParticle *parent2 = (TParticle*) fStack->Particle(lm2);
    
    //Primary Mothers
    Int_t lfirst_mother1 = fStack->GetPrimary (lm1);
    Int_t lfirst_mother2 = fStack->GetPrimary (lm2);
    TParticle *fm1 = (TParticle*) fStack->Particle(lfirst_mother1);
    TParticle *fm2 = (TParticle*) fStack->Particle(lfirst_mother2);

    //Distance On Stack
    Int_t Dl = TMath::Abs(lfirst_mother1-lfirst_mother2);
    
    //Correlated Pairs Selection
    if ( IsLightFlavorParticle(parent1) && IsLightFlavorParticle(parent2) && lm1==lm2 ) return true;
    if ( IsHeavyFlavorParticle(parent1) && IsHeavyFlavorParticle(parent2) && IsHeavyFlavorParticle(fm1) && IsHeavyFlavorParticle(fm2) && Dl<=2 ) return true;
    if ( lm1==lm2 ) return true;

    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Efficiency::Weight (TParticle *particle)  {
    
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
        case 111: weight = fHisto_Hijing_PizeroWeight   -> GetBinContent(nx,ny); break;
        case 221: weight = fHisto_Hijing_EtaWeight      -> GetBinContent(nx,ny); break;
        case 331: weight = fHisto_Hijing_EtaPrimeWeight -> GetBinContent(nx,ny); break;
        case 113: weight = fHisto_Hijing_RhoWeight      -> GetBinContent(nx,ny); break;
        case 223: weight = fHisto_Hijing_OmegaWeight    -> GetBinContent(nx,ny); break;
        case 333: weight = fHisto_Hijing_PhiWeight      -> GetBinContent(nx,ny); break;
            
        default: weight = 1; break;
    }
    
    
    //Weights For Conversions
    if ( IsFromConversion(particle)) weight = WeightConversions (parent);
    
    return weight;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Efficiency::WeightConversions (TParticle *parent)  {
    
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
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::IsTrackFromHijing (AliESDtrack *track)  {
    
    Bool_t IsHijing = (false);
    
    Int_t lp = track->GetLabel();
    if ( lp==-1 )  return IsHijing;
    if ( fMCEvent->IsFromBGEvent(fStack->GetPrimary(TMath::Abs(lp))) ) IsHijing = true;
    
    return IsHijing;}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::PassedTrackQualityCuts (AliESDtrack* track)  {
    
    fESDTrackCuts -> SetMaxDCAToVertexZ(fDCAz_max);
    fESDTrackCuts -> SetMaxDCAToVertexXY(fDCAxy_param0 + fDCAxy_param1/TMath::Power(track->Pt(),fDCAxy_param2));
    fESDTrackCuts -> SetAcceptKinkDaughters(false);
    fESDTrackCuts -> SetMinNCrossedRowsTPC(fTPC_minCr);
    fESDTrackCuts -> SetMinNClustersITS(fITS_minNcls);
    fESDTrackCuts -> SetMinNClustersTPC(fTPC_minNcls);
    fESDTrackCuts -> SetRequireTPCRefit(true);
    fESDTrackCuts -> SetRequireITSRefit(true);
    fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fESDTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(fMinCrOverFindableCls);
    fESDTrackCuts -> SetMaxChi2TPCConstrainedGlobal(fMaxGoldenChi2);
    fESDTrackCuts -> SetMaxChi2PerClusterTPC(fMaxTPCchi2);
    fESDTrackCuts -> SetMaxChi2PerClusterITS(fMaxITSchi2);
    fESDTrackCuts -> SetEtaRange(-0.8, 0.8);
    fESDTrackCuts -> SetPtRange(0.4,5.0);
    
    if ( track -> GetTPCsignalN() < fTPC_nClsdEdx ) return false;
    if ( FractionSharedClsITS(track) > fMaxFracSharedCls ) return false;
    if ( fESDTrackCuts->AcceptTrack (track) ) return true;
    return false;
}
//________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Efficiency::FractionSharedClsITS (AliESDtrack *track)  {
    
    Double_t nSharedCls(0);
    Double_t fSharedCls(0);
    
    for ( Int_t i=0 ; i<6 ; i++ )
        if ( track->HasPointOnITSLayer(i) && track->HasSharedPointOnITSLayer(i) ) nSharedCls++;
    
    fSharedCls = nSharedCls/(Double_t)track->GetNcls(0);
    
    return fSharedCls;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::PassedLooseTrackQualityCuts (AliESDtrack* track)  {
    
    fESDTrackCuts -> SetMinNCrossedRowsTPC(70);
    fESDTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fESDTrackCuts -> SetMaxChi2PerClusterTPC(4);
    fESDTrackCuts -> SetAcceptKinkDaughters(false);
    fESDTrackCuts -> SetRequireTPCRefit(true);
    fESDTrackCuts -> SetRequireITSRefit(true);
    fESDTrackCuts -> SetMinNClustersITS(3);
    fESDTrackCuts -> SetMinNClustersTPC(50);
    fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    fESDTrackCuts -> SetMaxDCAToVertexXY(1.0);
    fESDTrackCuts -> SetMaxDCAToVertexZ(3.0);
    fESDTrackCuts -> SetDCAToVertex2D(false);
    fESDTrackCuts -> SetRequireSigmaToVertex(false);
    fESDTrackCuts -> SetMaxChi2PerClusterITS(36);
    
    if ( fESDTrackCuts->AcceptTrack (track) ) return true;
    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Efficiency::PassedPIDCuts (AliESDtrack *track)  {
    
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
Double_t AliAnalysisTaskDielectronsPbPb_Efficiency::GetMass (TVector3 P1,TVector3 P2)  {
    
    Double_t electron_mass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    
    Double_t E1 = TMath::Sqrt( TMath::Power(electron_mass,2) + P1.Mag2());
    Double_t E2 = TMath::Sqrt( TMath::Power(electron_mass,2) + P2.Mag2());
    Double_t mass = TMath::Sqrt( TMath::Power(E1+E2,2) - (P1+P2).Mag2());
    
    return mass;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Efficiency::GetMass (AliESDtrack *track1,AliESDtrack *track2)  {
    
    TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
    TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
    
    Double_t electron_mass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    
    Double_t E1 = TMath::Sqrt( TMath::Power(electron_mass,2) + P1.Mag2());
    Double_t E2 = TMath::Sqrt( TMath::Power(electron_mass,2) + P2.Mag2());
    Double_t mass = TMath::Sqrt( TMath::Power(E1+E2,2) - (P1+P2).Mag2());
    
    return mass;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Efficiency::GetPtee (AliESDtrack *track1,AliESDtrack *track2)  {
    
    TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
    TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
    
    Double_t ptee = (P1+P2).Pt();
    
    return ptee;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Efficiency::GetPhiV (AliESDtrack *track1,AliESDtrack *track2 ) {
    
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
void AliAnalysisTaskDielectronsPbPb_Efficiency::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________________

