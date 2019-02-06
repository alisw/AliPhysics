#include "AliAnalysisTaskHelium3EllipticFlow.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"

ClassImp(AliAnalysisTaskHelium3EllipticFlow)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHelium3EllipticFlow::AliAnalysisTaskHelium3EllipticFlow():
AliAnalysisTaskSE(),
fAODevent(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
fSignMagField(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(NULL),
fFilterBitMask(0),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fNumberClustersITSMin(0),
fNumberClustersTPCMin(0),
fCrossedRowsFindableClsMin(0),
fNumberClustersTPCdEdxMin(0),
fITSrequirement(NULL),
fDCAzMax(0),
fDCAxyMax(0),
fNumberSigmaTPCmax(0),
fHistoPtTransformationMatrix(NULL),
fCentralityMin(0),
fCentralityMax(0),
fHistoCentralityEvtPlaneWeight(NULL)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHelium3EllipticFlow::AliAnalysisTaskHelium3EllipticFlow(const char *name):
AliAnalysisTaskSE(name),
fAODevent(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
fSignMagField(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(NULL),
fFilterBitMask(0),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fNumberClustersITSMin(0),
fNumberClustersTPCMin(0),
fCrossedRowsFindableClsMin(0),
fNumberClustersTPCdEdxMin(0),
fITSrequirement(NULL),
fDCAzMax(0),
fDCAxyMax(0),
fNumberSigmaTPCmax(0),
fHistoPtTransformationMatrix(NULL),
fCentralityMin(0),
fCentralityMax(0),
fHistoCentralityEvtPlaneWeight(NULL)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHelium3EllipticFlow::~AliAnalysisTaskHelium3EllipticFlow()
{
    
    fOutputList->Clear();
    delete fAODevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    delete fQAList;
    delete fHistoPtTransformationMatrix;
    delete fHistoCentralityEvtPlaneWeight;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHelium3EllipticFlow::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
    fQAList = new TList();
    fQAList -> SetOwner();
    
    //QA Plots of Event Selection
    fAODeventCuts.AddQAplotsToList(fQAList);
    
    
    
    
    //Number of Events
    fHistoEvents = new TH1F ("fHistoEvents","",10,0,10);
    fOutputList -> Add(fHistoEvents);
    
    //Centrality Distribution 
    fHistoCentralityDistribution = new TH1F ("fHistoCentralityDistribution","",200,0,100);
    fOutputList -> Add(fHistoCentralityDistribution);
    
    //Centrality & Event Plane Angle
    fHistoCentrality_EvtPlaneAngle_NoFlattening = new TH2F ("fHistoCentrality_EvtPlaneAngle_NoFlattening","",200,-TMath::Pi()/2.0,TMath::Pi()/2.0,100,0,100);
    fHistoCentrality_EvtPlaneAngle_Flattening   = new TH2F ("fHistoCentrality_EvtPlaneAngle_Flattening","",200,-TMath::Pi()/2.0,TMath::Pi()/2.0,100,0,100);
    fHistoCentrality_EvtPlaneAngle_NoFlattening -> Sumw2();
    fHistoCentrality_EvtPlaneAngle_Flattening   -> Sumw2();
    fOutputList -> Add(fHistoCentrality_EvtPlaneAngle_NoFlattening);
    fOutputList -> Add(fHistoCentrality_EvtPlaneAngle_Flattening);
    
    //Event Plane Resolution Histograms
    fHistoCosn_PsiA_PsiB_vs_Centrality = new TH2F ("fHistoCosn_PsiA_PsiB_vs_Centrality","",20,0,100,2002,-1.001,1.001);
    fHistoCosn_PsiA_PsiC_vs_Centrality = new TH2F ("fHistoCosn_PsiA_PsiC_vs_Centrality","",20,0,100,2002,-1.001,1.001);
    fHistoCosn_PsiB_PsiC_vs_Centrality = new TH2F ("fHistoCosn_PsiB_PsiC_vs_Centrality","",20,0,100,2002,-1.001,1.001);
    fHistoCosn_PsiA_PsiB_vs_Centrality -> Sumw2();
    fHistoCosn_PsiA_PsiC_vs_Centrality -> Sumw2();
    fHistoCosn_PsiB_PsiC_vs_Centrality -> Sumw2();
    fOutputList -> Add(fHistoCosn_PsiA_PsiB_vs_Centrality);
    fOutputList -> Add(fHistoCosn_PsiA_PsiC_vs_Centrality);
    fOutputList -> Add(fHistoCosn_PsiB_PsiC_vs_Centrality);
    
    //He3 Centrality Distribution
    fHistoHe3_CentralityDistribution_Pt     = new TH2F ("fHistoHe3_CentralityDistribution_Pt","",100,0,10,20,0,100);
    fHistoAntiHe3_CentralityDistribution_Pt = new TH2F ("fHistoAntiHe3_CentralityDistribution_Pt","",100,0,10,20,0,100);
    fHistoHe3_CentralityDistribution_Pt     -> Sumw2();
    fHistoAntiHe3_CentralityDistribution_Pt -> Sumw2();
    fOutputList -> Add(fHistoHe3_CentralityDistribution_Pt);
    fOutputList -> Add(fHistoAntiHe3_CentralityDistribution_Pt);

    
    //DCAxy Distributions of Pure He3 Candidates
    fHistoHe3_DCAxy_vs_Pt_InPlane      = new TH2F ("fHistoHe3_DCAxy_vs_Pt_InPlane","",100,0,10,1000,-5,5);
    fHistoHe3_DCAxy_vs_Pt_OutPlane     = new TH2F ("fHistoHe3_DCAxy_vs_Pt_OutPlane","",100,0,10,1000,-5,5);
    fHistoAntiHe3_DCAxy_vs_Pt_InPlane  = new TH2F ("fHistoAntiHe3_DCAxy_vs_Pt_InPlane","",100,0,10,1000,-5,5);
    fHistoAntiHe3_DCAxy_vs_Pt_OutPlane = new TH2F ("fHistoAntiHe3_DCAxy_vs_Pt_OutPlane","",100,0,10,1000,-5,5);
    fHistoHe3_DCAxy_vs_Pt_InPlane      -> Sumw2();
    fHistoHe3_DCAxy_vs_Pt_OutPlane     -> Sumw2();
    fHistoAntiHe3_DCAxy_vs_Pt_InPlane  -> Sumw2();
    fHistoAntiHe3_DCAxy_vs_Pt_OutPlane -> Sumw2();
    fOutputList -> Add(fHistoHe3_DCAxy_vs_Pt_InPlane);
    fOutputList -> Add(fHistoHe3_DCAxy_vs_Pt_OutPlane);
    fOutputList -> Add(fHistoAntiHe3_DCAxy_vs_Pt_InPlane);
    fOutputList -> Add(fHistoAntiHe3_DCAxy_vs_Pt_OutPlane);
    
    
    //In-Plane & Out-Of-Plane He3 Candidates
    fHistoHe3_nsigmaTPC_vs_Pt_InPlane      = new TH2F ("fHistoHe3_nsigmaTPC_vs_Pt_InPlane","",100,0,10,200,-10,10);
    fHistoHe3_nsigmaTPC_vs_Pt_OutPlane     = new TH2F ("fHistoHe3_nsigmaTPC_vs_Pt_OutPlane","",100,0,10,200,-10,10);
    fHistoAntiHe3_nsigmaTPC_vs_Pt_InPlane  = new TH2F ("fHistoAntiHe3_nsigmaTPC_vs_Pt_InPlane","",100,0,10,200,-10,10);
    fHistoAntiHe3_nsigmaTPC_vs_Pt_OutPlane = new TH2F ("fHistoAntiHe3_nsigmaTPC_vs_Pt_OutPlane","",100,0,10,200,-10,10);
    fHistoHe3_nsigmaTPC_vs_Pt_InPlane      -> Sumw2();
    fHistoHe3_nsigmaTPC_vs_Pt_OutPlane     -> Sumw2();
    fHistoAntiHe3_nsigmaTPC_vs_Pt_InPlane  -> Sumw2();
    fHistoAntiHe3_nsigmaTPC_vs_Pt_OutPlane -> Sumw2();
    fOutputList -> Add(fHistoHe3_nsigmaTPC_vs_Pt_InPlane);
    fOutputList -> Add(fHistoHe3_nsigmaTPC_vs_Pt_OutPlane);
    fOutputList -> Add(fHistoAntiHe3_nsigmaTPC_vs_Pt_InPlane);
    fOutputList -> Add(fHistoAntiHe3_nsigmaTPC_vs_Pt_OutPlane);
    
    fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane      = new TH2F ("fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane","",300,0,3,200,-10,10);
    fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane     = new TH2F ("fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane","",300,0,3,200,-10,10);
    fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane  = new TH2F ("fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane","",300,0,3,200,-10,10);
    fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane = new TH2F ("fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane","",300,0,3,200,-10,10);
    fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane      -> Sumw2();
    fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane     -> Sumw2();
    fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane  -> Sumw2();
    fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane -> Sumw2();
    fOutputList -> Add(fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane);
    fOutputList -> Add(fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane);
    fOutputList -> Add(fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane);
    fOutputList -> Add(fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane);
    
    
    //Pt Spectra
    fHistoHe3_Pt     = new TH1F ("fHistoHe3_Pt","",100,0,10);
    fHistoAntiHe3_Pt = new TH1F ("fHistoAntiHe3_Pt","",100,0,10);
    fHistoHe3_Pt     -> Sumw2();
    fHistoAntiHe3_Pt -> Sumw2();
    fOutputList -> Add (fHistoHe3_Pt);
    fOutputList -> Add (fHistoAntiHe3_Pt);

    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHelium3EllipticFlow::UserExec(Option_t *)
{
    //Get Input Event
    Double_t centralityEvtPlaneWeight(1);
    if ( !GetEvent (centralityEvtPlaneWeight)) return;

    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    
    //Centrality & Event Plane
    fHistoCentrality_EvtPlaneAngle_NoFlattening -> Fill (GetEventPlaneAngle(),GetEventCentrality());
    fHistoCentrality_EvtPlaneAngle_Flattening   -> Fill (GetEventPlaneAngle(),GetEventCentrality(),centralityEvtPlaneWeight);
    
    
    //Calculate Event Plane Resolution
    GetEventPlaneResolution (centralityEvtPlaneWeight);
   
    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fAODevent->GetNumberOfTracks() ; i++)  {
        
        //Track Selection
        AliAODTrack *track = (AliAODTrack*) fAODevent -> GetTrack(i);
        if ( !track ) continue;
        if ( !PassedTrackQualityCuts (track)) continue;
       
        //Correct pt
        Double_t pt = GetCorrectedPt (track);
        
        //nsigmaTPC
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);

        //DCA
        Double_t DCAz  = GetDCAz (track);
        Double_t DCAxy = GetDCAxy (track);

        
        //DCA Distributions (Separation of Secondary He3 from Material)
        if ( PassedHelium3IDCuts (track) && TMath::Abs(DCAz) < fDCAzMax )  {
            if (track->Charge()>0 && IsTrackInPlane(track))  fHistoHe3_DCAxy_vs_Pt_InPlane      -> Fill(pt,DCAxy,centralityEvtPlaneWeight);
            if (track->Charge()>0 && !IsTrackInPlane(track)) fHistoHe3_DCAxy_vs_Pt_OutPlane     -> Fill(pt,DCAxy,centralityEvtPlaneWeight);
            if (track->Charge()<0 && IsTrackInPlane(track))  fHistoAntiHe3_DCAxy_vs_Pt_InPlane  -> Fill(pt,DCAxy,centralityEvtPlaneWeight);
            if (track->Charge()<0 && !IsTrackInPlane(track)) fHistoAntiHe3_DCAxy_vs_Pt_OutPlane -> Fill(pt,DCAxy,centralityEvtPlaneWeight);
        }

        
        //DCA Cuts
        if ( TMath::Abs(DCAxy) > fDCAxyMax ) continue;
        if ( TMath::Abs(DCAz)  > fDCAzMax )  continue;
        
        
        //Pt & Centrality Spectra
        if ( TMath::Abs(nsigmaTPC) < fNumberSigmaTPCmax )  {
            
            //Pt Spectra
            if ( track->Charge()>0) fHistoHe3_Pt     -> Fill (pt,centralityEvtPlaneWeight);
            if ( track->Charge()<0) fHistoAntiHe3_Pt -> Fill (pt,centralityEvtPlaneWeight);
            
            //Fill Centrality Distributions
            if ( track->Charge()>0) fHistoHe3_CentralityDistribution_Pt     -> Fill(pt,GetEventCentrality(),centralityEvtPlaneWeight);
            if ( track->Charge()<0) fHistoAntiHe3_CentralityDistribution_Pt -> Fill(pt,GetEventCentrality(),centralityEvtPlaneWeight);
        }
        
        
        //Fill In-Plane & Out-of-Plane Histograms (Primary He3 Candidates)
        if (track->Charge()>0)  {
            if (  IsTrackInPlane(track))  {
                fHistoHe3_nsigmaTPC_vs_Pt_InPlane             -> Fill(pt,nsigmaTPC,centralityEvtPlaneWeight);
                fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane  -> Fill(pt/3.0,nsigmaTPC,centralityEvtPlaneWeight);
            }
            
            if ( !IsTrackInPlane(track))  {
                fHistoHe3_nsigmaTPC_vs_Pt_OutPlane            -> Fill(pt,nsigmaTPC,centralityEvtPlaneWeight);
                fHistoHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane -> Fill(pt/3.0,nsigmaTPC,centralityEvtPlaneWeight);
            }
        }
        
        //Fill In-Plane & Out-of-Plane Histograms (Primary Anti-He3 Candidates)
        if (track->Charge()<0)  {
            if (  IsTrackInPlane(track))  {
                fHistoAntiHe3_nsigmaTPC_vs_Pt_InPlane             -> Fill(pt,nsigmaTPC,centralityEvtPlaneWeight);
                fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_InPlane  -> Fill(pt/3.0,nsigmaTPC,centralityEvtPlaneWeight);
            }
            
            if ( !IsTrackInPlane(track))  {
                fHistoAntiHe3_nsigmaTPC_vs_Pt_OutPlane            -> Fill(pt,nsigmaTPC,centralityEvtPlaneWeight);
                fHistoAntiHe3_nsigmaTPC_vs_Pt_PerNucleon_OutPlane -> Fill(pt/3.0,nsigmaTPC,centralityEvtPlaneWeight);
            }
        }
    }
    
    
    PostData(1, fOutputList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHelium3EllipticFlow::GetEvent (Double_t &centralityEvtPlaneWeight)  {
    
    //Get Input Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return false;
    fHistoEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fAODeventCuts.AcceptEvent(fAODevent)) {
        PostData(2, fQAList);
        return false;
    }
    fHistoEvents -> Fill(1.5);

    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    fHistoEvents -> Fill(2.5);
    Double_t centr = multiplicitySelection->GetMultiplicityPercentile(fCentralityEstimator);
    
    //Centrality Distributions
    fHistoCentralityDistribution -> Fill(centr);
    
    //Selection of Centrality Range
    if (centr<fCentralityMin || centr>=fCentralityMax ) return false;
    fHistoEvents -> Fill(3.5);

    //Primary Vertex
    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    fHistoEvents -> Fill(4.5);
    if ( vertex->GetZ() < fVertexZmin ) return false;
    if ( vertex->GetZ() > fVertexZmax ) return false;
    fHistoEvents -> Fill(5.5);
    if ( vertex->GetNContributors() < fNumberVertexContributorsMin ) return false;
    fHistoEvents -> Fill(6.5);

    //Event Plane
    AliEventplane *eventPlane =  fAODevent->GetEventplane();
    if (!eventPlane) return false;
    fHistoEvents -> Fill(7.5);

    //Magnetic Field Configuration
    Double_t signB = fAODevent->GetMagneticField()/TMath::Abs(fAODevent->GetMagneticField());
    if (fSignMagField != 0 && signB != fSignMagField) return false;
    fHistoEvents -> Fill(8.5);
    
    //Centrality Weight
    Int_t ix = fHistoCentralityEvtPlaneWeight -> GetXaxis() -> FindBin(GetEventPlaneAngle());
    Int_t iy = fHistoCentralityEvtPlaneWeight -> GetYaxis() -> FindBin(GetEventCentrality());
    Double_t weight = fHistoCentralityEvtPlaneWeight -> GetBinContent(ix,iy);
    
    if (weight!=0) centralityEvtPlaneWeight = weight;
    if (weight==0) centralityEvtPlaneWeight = 1;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHelium3EllipticFlow::GetEventCentrality()  {
    
    Double_t centralityPerc = 0;
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    centralityPerc = multiplicitySelection->GetMultiplicityPercentile(fCentralityEstimator);
    
    return centralityPerc;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHelium3EllipticFlow::GetEventPlaneAngle()  {
    
    AliEventplane *eventPlane =  fAODevent->GetEventplane();
    Double_t psiEP = eventPlane -> GetEventplane("V0",fAODevent, 2);
    return psiEP;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

void AliAnalysisTaskHelium3EllipticFlow::GetEventPlaneResolution (Double_t centralityEvtPlaneWeight)  {
    
    //Initialization
    Double_t qxEPa(0),qyEPa(0),qxEPc(0),qyEPc(0),qxEP(0),qyEP(0);
    Double_t Qx2(0),Qy2(0),Qx2p(0),Qy2p(0),Qx2n(0),Qy2n(0);
    
    //Loop on Tracks
    for (Int_t i=0 ; i<fAODevent->GetNumberOfTracks() ; i++)  {
        
        //Track Selection
        AliAODTrack *track = (AliAODTrack*) fAODevent -> GetTrack(i);
        if ( !track ) continue;
        if ( !PassedEvtPlaneResTrackQualityCuts(track)) continue;
        
        //Positive & Negative Eta Regions
        if (track->Eta()>0)  { Qx2p = Qx2p + TMath::Cos(2.0*track->Phi()); Qy2p = Qy2p + TMath::Sin(2.0*track->Phi()); }
        if (track->Eta()<0)  { Qx2n = Qx2n + TMath::Cos(2.0*track->Phi()); Qy2n = Qy2n + TMath::Sin(2.0*track->Phi()); }
        
        //Full TPC
        Qx2 = Qx2 + TMath::Cos(2.0*track->Phi());
        Qy2 = Qy2 + TMath::Sin(2.0*track->Phi());
    }
    
    //Event Plane Angles
    Double_t evPlAngV0A    = fAODevent->GetEventplane()->CalculateVZEROEventPlane(fAODevent, 8, 2, qxEPa, qyEPa);
    Double_t evPlAngV0C    = fAODevent->GetEventplane()->CalculateVZEROEventPlane(fAODevent, 9, 2, qxEPc, qyEPc);
    Double_t evPlAngV0     = fAODevent->GetEventplane()->CalculateVZEROEventPlane(fAODevent, 10, 2, qxEP, qyEP);
    Double_t evPlAngTPC    = TMath::ATan2(Qy2,Qx2)/2.0;
    Double_t evPlAngTPCneg = TMath::ATan2(Qy2n,Qx2n)/2.0;
    Double_t evPlAngTPCpos = TMath::ATan2(Qy2p,Qx2p)/2.0;
   
    //Fill Event Plane Resolution Histograms
    fHistoCosn_PsiA_PsiB_vs_Centrality -> Fill(GetEventCentrality(),TMath::Cos(2.0*(evPlAngV0-evPlAngTPCpos)),centralityEvtPlaneWeight);// VZERO - TPCpos
    fHistoCosn_PsiA_PsiC_vs_Centrality -> Fill(GetEventCentrality(),TMath::Cos(2.0*(evPlAngV0-evPlAngTPCneg)),centralityEvtPlaneWeight);// VZERO - TPCneg
    fHistoCosn_PsiB_PsiC_vs_Centrality -> Fill(GetEventCentrality(),TMath::Cos(2.0*(evPlAngTPCpos-evPlAngTPCneg)),centralityEvtPlaneWeight);// TPCpos - TPCneg
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

Bool_t AliAnalysisTaskHelium3EllipticFlow::PassedEvtPlaneResTrackQualityCuts (AliAODTrack *track)  {
    
    if ( !track->TestFilterBit(128))               return false;
    if ( TMath::Abs(track->Eta()) > 0.8)           return false;
    if ( track->Pt() < 0.2 || track->Pt() > 20.0)  return false;
    if ( track->GetTPCNcls() < 70)                 return false;
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHelium3EllipticFlow::PassedTrackQualityCuts (AliAODTrack* track)  {
    
    //Filterbit
    if(!track->TestFilterMask(fFilterBitMask)) return false;
    
    //Kinematic Cuts & Acceptance
    if ( track->Pt()<fPtMin || track->Pt()>fPtMax ) return false;
    if ( TMath::Abs(track->Eta()) > fEtaMax )       return false;
    
    //Track Quality Cuts
    Bool_t hitInITSLayer0 = track->HasPointOnITSLayer(0);
    Bool_t hitInITSLayer1 = track->HasPointOnITSLayer(1);
    
    
    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kFirst")==0  && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kSecond")==0 && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kAny")==0    && !hitInITSLayer0 && !hitInITSLayer1) return false;

    if ( track->GetITSNcls() < fNumberClustersITSMin ) return false;
    if ( track->GetTPCNcls() < fNumberClustersTPCMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCNCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberClustersTPCdEdxMin ) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHelium3EllipticFlow::GetDCAxy (AliAODTrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAxy = impactParameter[0];
    
    return DCAxy;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHelium3EllipticFlow::GetDCAz (AliAODTrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAz = impactParameter[1];
    
    return DCAz;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHelium3EllipticFlow::GetCorrectedPt (AliAODTrack *track)  {
    
    Double_t pt_rec  = 2.0*track->Pt();//rigidity -> pt
    Double_t pt_true = pt_rec;//Initialization
    Double_t eta = track->Eta();
    
    Int_t ix = fHistoPtTransformationMatrix->GetXaxis()->FindBin(pt_rec);
    Int_t iy = fHistoPtTransformationMatrix->GetYaxis()->FindBin(eta);
    
    Double_t Dpt = fHistoPtTransformationMatrix->GetBinContent(ix,iy);
    pt_true = Dpt + pt_rec;
    return pt_true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHelium3EllipticFlow::PassedHelium3IDCuts (AliAODTrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);
    
    if ( TMath::Abs(nsigmaTPC) > 3.0) return false;
    if ( TMath::Abs(nsigmaTOF) > 3.0) return false;
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHelium3EllipticFlow::IsTrackInPlane (AliAODTrack *track)  {
    
    Bool_t isInPlane = false;
    
    Double_t psiEP    = TVector2::Phi_0_2pi(GetEventPlaneAngle());
    Double_t deltaphi = TMath::Abs(track->Phi() - psiEP);
    
    if (deltaphi<=TMath::Pi()/4.0) isInPlane = true;
    if (deltaphi > 3.0*TMath::Pi()/4.0 && deltaphi <= 5.0*TMath::Pi()/4.0) isInPlane = true;
    if (deltaphi>7.0*TMath::Pi()/4.0) isInPlane = true;

    return isInPlane;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHelium3EllipticFlow::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

