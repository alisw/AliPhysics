#include "AliAnalysisTaskXiNucleusInteraction.h"
#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "AliESDv0.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1I.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
class AliAnalysisTaskXiNucleusInteraction;
ClassImp(AliAnalysisTaskXiNucleusInteraction)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskXiNucleusInteraction::AliAnalysisTaskXiNucleusInteraction():
AliAnalysisTaskSE(),
fESDeventCuts(),
fESDevent(nullptr),
fPIDResponse(nullptr),
fUtils(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fESDtrackCuts(nullptr),
fTrigger(AliVEvent::kINT7),
fMultLow(0),
fMultHigh(100),
hNumberOfEvents(nullptr),
hNumberOfXi(nullptr),
hNumberOfAntiXi(nullptr)
{}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskXiNucleusInteraction::AliAnalysisTaskXiNucleusInteraction(const char *name):
AliAnalysisTaskSE(name),
fESDeventCuts(),
fESDevent(nullptr),
fPIDResponse(nullptr),
fUtils(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fESDtrackCuts(nullptr),
fTrigger(AliVEvent::kINT7),
fMultLow(0),
fMultHigh(100),
hNumberOfEvents(nullptr),
hNumberOfXi(nullptr),
hNumberOfAntiXi(nullptr)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskXiNucleusInteraction::~AliAnalysisTaskXiNucleusInteraction()  {
    
    fOutputList -> Clear();
    fQAList     -> Clear();
    
    delete fESDevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    delete fQAList;
    delete fESDtrackCuts;
    delete hNumberOfEvents;
    delete hNumberOfXi;
    delete hNumberOfAntiXi;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskXiNucleusInteraction::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();
    
    //QA Plots of Event Selection
    fESDeventCuts.AddQAplotsToList(fQAList);
    fESDeventCuts.OverrideAutomaticTriggerSelection(fTrigger);
    fESDeventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    
    //Event Counter
    hNumberOfEvents = new TH1F ("hNumberOfEvents","",20,0,20);
    hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    
    //Rejection Factor Calculation
    hNumberOfXi     = new TH1F ("hNumberOfXi","",10,0,10);
    hNumberOfAntiXi = new TH1F ("hNumberOfAntiXi","",10,0,10);
    hNumberOfXi     -> Sumw2();
    hNumberOfAntiXi -> Sumw2();
    fOutputList -> Add(hNumberOfXi);
    fOutputList -> Add(hNumberOfAntiXi);

    for (Int_t i=0 ; i<7 ; i++)  {
        
        //Invariant-Mass Plots
        hMassXi_vs_P[i]     = new TH2F (Form("hMassXi_vs_P[%d]",i),"",100,0,10,500,1.30,1.35);
        hMassAntiXi_vs_P[i] = new TH2F (Form("hMassAntiXi_vs_P[%d]",i),"",100,0,10,500,1.30,1.35);
        hMassXi_vs_P[i]     -> Sumw2();
        hMassAntiXi_vs_P[i] -> Sumw2();
        fOutputList -> Add(hMassXi_vs_P[i]);
        fOutputList -> Add(hMassAntiXi_vs_P[i]);
        
        //Pointing Angle Xi
        hXiPointingAngle_vs_P[i]     = new TH2F (Form("hXiPointingAngle_vs_P[%d]",i),"",100,0,10,200,0,20);
        hAntiXiPointingAngle_vs_P[i] = new TH2F (Form("hAntiXiPointingAngle_vs_P[%d]",i),"",100,0,10,200,0,20);
        hXiPointingAngle_vs_P[i]     -> Sumw2();
        hAntiXiPointingAngle_vs_P[i] -> Sumw2();
        fOutputList -> Add(hXiPointingAngle_vs_P[i]);
        fOutputList -> Add(hAntiXiPointingAngle_vs_P[i]);
    }
    
    //Track Cuts Objects
    fESDtrackCuts = new AliESDtrackCuts ("fESDtrackCuts");
    SetTrackSelectionCriteria ();
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//______________________________________________________________________________________________________________________________________
void AliAnalysisTaskXiNucleusInteraction::SetTrackSelectionCriteria()  {
    
    //Track Cuts Cascade and V0 Daughters
    fESDtrackCuts -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts -> SetMinNCrossedRowsTPC(80);
    fESDtrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fESDtrackCuts -> SetMaxChi2PerClusterTPC(4);
    fESDtrackCuts -> SetEtaRange(-0.8,0.8);
    fESDtrackCuts -> SetPtRange(0.1,200);
    fESDtrackCuts -> SetMinDCAToVertexXY(0.05);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskXiNucleusInteraction::UserExec(Option_t *)  {
    
    //Get Input Event
    if ( !GetEvent ()) return;
    
    //Primary Vertex
    AliESDVertex *primaryVertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    Double_t vx = primaryVertex->GetX();
    Double_t vy = primaryVertex->GetY();
    Double_t vz = primaryVertex->GetZ();
    
    //Flags
    Bool_t containsRecXi_R0(kFALSE);
    Bool_t containsRecXi_R1(kFALSE);
    Bool_t containsRecXi_R2(kFALSE);
    Bool_t containsRecXi_R3(kFALSE);
    Bool_t containsRecXi_R4(kFALSE);
    Bool_t containsRecXi_R5(kFALSE);
    Bool_t containsRecXi_R6(kFALSE);
    Bool_t containsRecAntiXi_R0(kFALSE);
    Bool_t containsRecAntiXi_R1(kFALSE);
    Bool_t containsRecAntiXi_R2(kFALSE);
    Bool_t containsRecAntiXi_R3(kFALSE);
    Bool_t containsRecAntiXi_R4(kFALSE);
    Bool_t containsRecAntiXi_R5(kFALSE);
    Bool_t containsRecAntiXi_R6(kFALSE);

    //Loop over Reconstructed Cascades
    for (Int_t icasc=0 ; icasc<fESDevent->GetNumberOfCascades(); icasc++)  {
        
        //Get Reconstructed Cascade
        AliESDcascade *cascade = fESDevent->GetCascade(icasc);
        if (!cascade) continue;
        
        //Get Decay Daughters
        AliESDtrack *positive_track = (AliESDtrack*) fESDevent->GetTrack(TMath::Abs(cascade->GetPindex()));
        AliESDtrack *negative_track = (AliESDtrack*) fESDevent->GetTrack(TMath::Abs(cascade->GetNindex()));
        AliESDtrack *bachelor_track = (AliESDtrack*) fESDevent->GetTrack(TMath::Abs(cascade->GetBindex()));
        if (!positive_track) continue;
        if (!negative_track) continue;
        if (!bachelor_track) continue;

        //Track Quality Cuts
        if (!PassedTrackQualityCuts(positive_track)) continue;
        if (!PassedTrackQualityCuts(negative_track)) continue;
        if (!PassedTrackQualityCuts(bachelor_track)) continue;

        //Cascade Cuts
        if (!PassedCascadeSelectionCuts(cascade)) continue;
        
        //Cascade Radius
        Double_t xCasc(0),yCasc(0),zCasc(0);
        cascade->GetXYZcascade (xCasc,yCasc,zCasc);
        Double_t rC = TMath::Sqrt(xCasc*xCasc + yCasc*yCasc);
        
        //Trigger Requirement: TOF Hits (reduce pile-up)
        Bool_t hasTOFhit_pos = (positive_track->GetStatus() & AliVTrack::kTOFout) && (positive_track->GetStatus() & AliVTrack::kTIME);
        Bool_t hasTOFhit_neg = (negative_track->GetStatus() & AliVTrack::kTOFout) && (negative_track->GetStatus() & AliVTrack::kTIME);
        Bool_t hasTOFhit_bac = (bachelor_track->GetStatus() & AliVTrack::kTOFout) && (bachelor_track->GetStatus() & AliVTrack::kTIME);
        if ((!hasTOFhit_pos) && (!hasTOFhit_neg) && (!hasTOFhit_bac)) continue;
        
        //Candidate Selection
        Bool_t isXi(kFALSE),isAntiXi(kFALSE);
        Double_t m_xi(0),m_antixi(0),p_xi(0),p_antixi(0);
        if (IsXiCandidate(cascade,positive_track,negative_track,bachelor_track,m_xi,p_xi))             isXi     = kTRUE;
        if (IsAntiXiCandidate(cascade,positive_track,negative_track,bachelor_track,m_antixi,p_antixi)) isAntiXi = kTRUE;
        if ((!isXi)&&(!isAntiXi)) continue;
        
        //Pointing Angle
        Double_t cos_theta = cascade->GetCascadeCosineOfPointingAngle(vx,vy,vz);
        Double_t angle_rad = TMath::ACos(cos_theta);
        Double_t angle = (180.0/TMath::Pi())*angle_rad;

        //Xi Selection
        if (isXi)  {
            if (rC> 2.24) {containsRecXi_R0 = kTRUE; hMassXi_vs_P[0] -> Fill (p_xi,m_xi);}
            if (rC> 3.01) {containsRecXi_R1 = kTRUE; hMassXi_vs_P[1] -> Fill (p_xi,m_xi);}
            if (rC> 3.78) {containsRecXi_R2 = kTRUE; hMassXi_vs_P[2] -> Fill (p_xi,m_xi);}
            if (rC>19.44) {containsRecXi_R3 = kTRUE; hMassXi_vs_P[3] -> Fill (p_xi,m_xi);}
            if (rC>24.39) {containsRecXi_R4 = kTRUE; hMassXi_vs_P[4] -> Fill (p_xi,m_xi);}
            if (rC>34.23) {containsRecXi_R5 = kTRUE; hMassXi_vs_P[5] -> Fill (p_xi,m_xi);}
            if (rC>39.18) {containsRecXi_R6 = kTRUE; hMassXi_vs_P[6] -> Fill (p_xi,m_xi);}
        }
        
        //Pointing Angle Xi
        if (isXi && m_xi>1.320 && m_xi<1.327)  {
            if (rC> 2.24) {hXiPointingAngle_vs_P[0] -> Fill (p_xi,angle);}
            if (rC> 3.01) {hXiPointingAngle_vs_P[1] -> Fill (p_xi,angle);}
            if (rC> 3.78) {hXiPointingAngle_vs_P[2] -> Fill (p_xi,angle);}
            if (rC>19.44) {hXiPointingAngle_vs_P[3] -> Fill (p_xi,angle);}
            if (rC>24.39) {hXiPointingAngle_vs_P[4] -> Fill (p_xi,angle);}
            if (rC>34.23) {hXiPointingAngle_vs_P[5] -> Fill (p_xi,angle);}
            if (rC>39.18) {hXiPointingAngle_vs_P[6] -> Fill (p_xi,angle);}
        }
        
        //AntiXi Selection
        if (isAntiXi)  {
            if (rC> 2.24) {containsRecAntiXi_R0=kTRUE;hMassAntiXi_vs_P[0]->Fill(p_antixi,m_antixi);}
            if (rC> 3.01) {containsRecAntiXi_R1=kTRUE;hMassAntiXi_vs_P[1]->Fill(p_antixi,m_antixi);}
            if (rC> 3.78) {containsRecAntiXi_R2=kTRUE;hMassAntiXi_vs_P[2]->Fill(p_antixi,m_antixi);}
            if (rC>19.44) {containsRecAntiXi_R3=kTRUE;hMassAntiXi_vs_P[3]->Fill(p_antixi,m_antixi);}
            if (rC>24.39) {containsRecAntiXi_R4=kTRUE;hMassAntiXi_vs_P[4]->Fill(p_antixi,m_antixi);}
            if (rC>34.23) {containsRecAntiXi_R5=kTRUE;hMassAntiXi_vs_P[5]->Fill(p_antixi,m_antixi);}
            if (rC>39.18) {containsRecAntiXi_R6=kTRUE;hMassAntiXi_vs_P[6]->Fill(p_antixi,m_antixi);}
        }
        
        //Pointing Angle AntiXi
        if (isAntiXi && m_antixi>1.320 && m_antixi<1.327)  {
            if (rC> 2.24) {hAntiXiPointingAngle_vs_P[0] -> Fill (p_antixi,angle);}
            if (rC> 3.01) {hAntiXiPointingAngle_vs_P[1] -> Fill (p_antixi,angle);}
            if (rC> 3.78) {hAntiXiPointingAngle_vs_P[2] -> Fill (p_antixi,angle);}
            if (rC>19.44) {hAntiXiPointingAngle_vs_P[3] -> Fill (p_antixi,angle);}
            if (rC>24.39) {hAntiXiPointingAngle_vs_P[4] -> Fill (p_antixi,angle);}
            if (rC>34.23) {hAntiXiPointingAngle_vs_P[5] -> Fill (p_antixi,angle);}
            if (rC>39.18) {hAntiXiPointingAngle_vs_P[6] -> Fill (p_antixi,angle);}
        }
    }
    
    //Fill Reconstructed Distributions (Xi)
    if (containsRecXi_R0) hNumberOfXi->Fill(0.5);
    if (containsRecXi_R1) hNumberOfXi->Fill(1.5);
    if (containsRecXi_R2) hNumberOfXi->Fill(2.5);
    if (containsRecXi_R3) hNumberOfXi->Fill(3.5);
    if (containsRecXi_R4) hNumberOfXi->Fill(4.5);
    if (containsRecXi_R5) hNumberOfXi->Fill(5.5);
    if (containsRecXi_R6) hNumberOfXi->Fill(6.5);
    
    //Fill Reconstructed Distributions (AntiXi)
    if (containsRecAntiXi_R0) hNumberOfAntiXi->Fill(0.5);
    if (containsRecAntiXi_R1) hNumberOfAntiXi->Fill(1.5);
    if (containsRecAntiXi_R2) hNumberOfAntiXi->Fill(2.5);
    if (containsRecAntiXi_R3) hNumberOfAntiXi->Fill(3.5);
    if (containsRecAntiXi_R4) hNumberOfAntiXi->Fill(4.5);
    if (containsRecAntiXi_R5) hNumberOfAntiXi->Fill(5.5);
    if (containsRecAntiXi_R6) hNumberOfAntiXi->Fill(6.5);


    PostData(1, fOutputList);
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskXiNucleusInteraction::GetEvent ()  {
    
    //Get Input Event
    fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDevent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);
            
    //Standard Event Cuts
    if (!fESDeventCuts.AcceptEvent(fESDevent)) {
        PostData(2, fQAList);
        return kFALSE;
    }
    hNumberOfEvents -> Fill(1.5);

    //Reject Events with Incomplete DAQ
    if (fESDevent->IsIncompleteDAQ()) return kFALSE;
    hNumberOfEvents -> Fill(2.5);
        
    //V0 Timing Decision
    AliVVZERO *vzeroData = fESDevent->GetVZEROData();
    if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision())) return kFALSE;
    hNumberOfEvents -> Fill(3.5);
        
    //Pileup Rejection
    Int_t nClustersLayer0 = fESDevent->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = fESDevent->GetNumberOfITSClusters(1);
    Int_t nTracklets      = fESDevent->GetMultiplicity()->GetNumberOfTracklets();
    if ((nClustersLayer0 + nClustersLayer1) > 65.0 + (Double_t)nTracklets*4.0) return kFALSE;
    hNumberOfEvents -> Fill(4.5);

    //Primary Vertex Tracks
    AliESDVertex *vertex_tracks = (AliESDVertex*) fESDevent->GetPrimaryVertexTracks();
    if (!vertex_tracks) return kFALSE;
    hNumberOfEvents -> Fill(5.5);
        
    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(6.5);
        
    //Primary Vertex SPD
    AliESDVertex *vertex_SPD = (AliESDVertex*) fESDevent->GetPrimaryVertexSPD();
    if (!vertex_SPD) return kFALSE;
    hNumberOfEvents -> Fill(7.5);
        
    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(8.5);
        
    //SPD Pile-up in Mult Bins
    if (fESDevent->IsPileupFromSPDInMultBins()) return kFALSE;
    hNumberOfEvents -> Fill(9.5);
        
    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3) return kFALSE;
    hNumberOfEvents -> Fill(10.5);

    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -10.0 ) return kFALSE;
    if ( vertex_tracks->GetZ() > +10.0 ) return kFALSE;
    hNumberOfEvents -> Fill(11.5);
               
    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return kFALSE;
    hNumberOfEvents -> Fill(12.5);
                
    //Selection of Multiplicity Range
    Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    if (mult_percentile < fMultLow)  return kFALSE;
    if (mult_percentile > fMultHigh) return kFALSE;
    hNumberOfEvents -> Fill(13.5);
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    return kTRUE;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskXiNucleusInteraction::PassedTrackQualityCuts (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    //Track Selection Cuts
    if ( track->GetTPCsignalN() < 50 ) return passedTrkSelection;
    if ( !fESDtrackCuts->AcceptTrack (track) ) return passedTrkSelection;
    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskXiNucleusInteraction::PassedCascadeSelectionCuts (AliESDcascade *cascade)  {

    //Initialization
    Bool_t passedCascadeSelection=(kFALSE);

    //Primary Vertex
    AliESDVertex *primaryVertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    Double_t vx = primaryVertex->GetX();
    Double_t vy = primaryVertex->GetY();
    Double_t vz = primaryVertex->GetZ();
    
    //Topological Selections
    if (cascade->GetV0CosineOfPointingAngle(vx,vy,vz)<0.97)      return passedCascadeSelection;
    if (cascade->GetDcaV0Daughters()>1.6) return passedCascadeSelection;
    if (cascade->GetDcaXiDaughters()>1.6) return passedCascadeSelection;
    
    //Cascade Radius
    Double_t xCasc(0),yCasc(0),zCasc(0);
    cascade->GetXYZcascade (xCasc,yCasc,zCasc);
    Double_t rC = TMath::Sqrt(xCasc*xCasc + yCasc*yCasc);
    if (rC>34.0) return passedCascadeSelection;
    
    //V0 Radius
    Double_t xV0(0),yV0(0),zV0(0);
    cascade->GetXYZ(xV0,yV0,zV0);
    Double_t radiusV0 = TMath::Sqrt(xV0*xV0 + yV0*yV0);
    if (radiusV0>34.0) return passedCascadeSelection;
    
    passedCascadeSelection = kTRUE;
    return passedCascadeSelection;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskXiNucleusInteraction::IsXiCandidate (AliESDcascade *casc, AliESDtrack *pos, AliESDtrack *neg, AliESDtrack *bac, Double_t &m, Double_t &p)  {
    
    //PID Daughters
    Bool_t passedPID=(kFALSE);
    if (PassedPIDSelection(pos,AliPID::kProton) && PassedPIDSelection(neg,AliPID::kPion) && PassedPIDSelection(bac,AliPID::kPion)) passedPID=kTRUE;
    if (!passedPID) return kFALSE;
        
    //Reconstructed Momentum Components Cascade
    Double_t momentum_positive_Rec[3]={0,0,0};
    Double_t momentum_negative_Rec[3]={0,0,0};
    Double_t momentum_bachelor_Rec[3]={0,0,0};
    casc->GetPPxPyPz(momentum_positive_Rec[0], momentum_positive_Rec[1], momentum_positive_Rec[2]);
    casc->GetNPxPyPz(momentum_negative_Rec[0], momentum_negative_Rec[1], momentum_negative_Rec[2]);
    casc->GetBPxPyPz(momentum_bachelor_Rec[0], momentum_bachelor_Rec[1], momentum_bachelor_Rec[2]);
    
    //Proton and Pion Masses
    Double_t mass_pion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    Double_t mass_prot = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

    //4-Momentum Vectors
    TLorentzVector Ppos,Pneg,Pbac,Ptot;
    Ppos.SetXYZM(momentum_positive_Rec[0], momentum_positive_Rec[1], momentum_positive_Rec[2],mass_prot);
    Pneg.SetXYZM(momentum_negative_Rec[0], momentum_negative_Rec[1], momentum_negative_Rec[2],mass_pion);
    Pbac.SetXYZM(momentum_bachelor_Rec[0], momentum_bachelor_Rec[1], momentum_bachelor_Rec[2],mass_pion);
    Ptot=Ppos+Pneg+Pbac;
    Double_t mass = Ptot.M();
    
    //Mass Selection
    if (mass<1.30) return kFALSE;
    if (mass>1.35) return kFALSE;
    
    //Assignments
    m  = mass;
    TVector3 momentum = Ptot.Vect();
    p = momentum.Mag();

    return kTRUE;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskXiNucleusInteraction::IsAntiXiCandidate (AliESDcascade *casc, AliESDtrack *pos, AliESDtrack *neg, AliESDtrack *bac, Double_t &m, Double_t &p)  {
    
    //PID Daughters
    Bool_t passedPID=(kFALSE);
    if (PassedPIDSelection(pos,AliPID::kPion) && PassedPIDSelection(neg,AliPID::kProton) && PassedPIDSelection(bac,AliPID::kPion)) passedPID=kTRUE;
    if (!passedPID) return kFALSE;
        
    //Reconstructed Momentum Components Cascade
    Double_t momentum_positive_Rec[3]={0,0,0};
    Double_t momentum_negative_Rec[3]={0,0,0};
    Double_t momentum_bachelor_Rec[3]={0,0,0};
    casc->GetPPxPyPz(momentum_positive_Rec[0], momentum_positive_Rec[1], momentum_positive_Rec[2]);
    casc->GetNPxPyPz(momentum_negative_Rec[0], momentum_negative_Rec[1], momentum_negative_Rec[2]);
    casc->GetBPxPyPz(momentum_bachelor_Rec[0], momentum_bachelor_Rec[1], momentum_bachelor_Rec[2]);
    
    //Proton and Pion Masses
    Double_t mass_pion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    Double_t mass_prot = TDatabasePDG::Instance()->GetParticle(2212)->Mass();

    //4-Momentum Vectors
    TLorentzVector Ppos,Pneg,Pbac,Ptot;
    Ppos.SetXYZM(momentum_positive_Rec[0], momentum_positive_Rec[1], momentum_positive_Rec[2],mass_pion);
    Pneg.SetXYZM(momentum_negative_Rec[0], momentum_negative_Rec[1], momentum_negative_Rec[2],mass_prot);
    Pbac.SetXYZM(momentum_bachelor_Rec[0], momentum_bachelor_Rec[1], momentum_bachelor_Rec[2],mass_pion);
    Ptot=Ppos+Pneg+Pbac;
    Double_t mass = Ptot.M();
    
    //Mass Selection
    if (mass<1.30) return kFALSE;
    if (mass>1.35) return kFALSE;

    //Assignments
    m  = mass;
    TVector3 momentum = Ptot.Vect();
    p = momentum.Mag();

    return kTRUE;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskXiNucleusInteraction::PassedPIDSelection (AliESDtrack *track, AliPID::EParticleType type)  {
    
    //Initialization
    Bool_t passedPIDSelection=(kFALSE);
    
    //TPC Particle Identification
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,type);
    if (nsigmaTPC < -3.0) return passedPIDSelection;
    if (nsigmaTPC > +3.0) return passedPIDSelection;

    passedPIDSelection = kTRUE;
    return passedPIDSelection;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskXiNucleusInteraction::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_____________________________________________________________________________________________________________________________________
