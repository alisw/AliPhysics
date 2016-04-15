#include "AliAnalysisTaskDielectronsPbPb_Data.h"
#include "AliDielectronReducedTrack.h"
#include "AliInputEventHandler.h"
#include "AliEventPoolManager.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliAnalysisTaskDielectronsPbPb_Data)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDielectronsPbPb_Data::AliAnalysisTaskDielectronsPbPb_Data():
AliAnalysisTaskSE(),
fESDevent(0x0),
fPoolMgr(0x0),
fESDTrackCuts(0x0),
fPIDResponse(0x0),
fOutputList(0x0),
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
fMassLim(0),
fPhivLim(0),
fPtMin(0),
fPtMax(0),
fEtaLim(0),
fNcentralityBins(0),
fNvertexBins(0),
fNeventPlaneBins(0),
fMaxNumberEvts(0),
fMaxNumberTrks(0),
fNumberEvtsToMix(0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDielectronsPbPb_Data::AliAnalysisTaskDielectronsPbPb_Data(const char *name):
AliAnalysisTaskSE(name),
fESDevent(0x0),
fPoolMgr(0x0),
fESDTrackCuts(0x0),
fPIDResponse(0x0),
fOutputList(0x0),
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
fMassLim(0),
fPhivLim(0),
fPtMin(0),
fPtMax(0),
fEtaLim(0),
fNcentralityBins(0),
fNvertexBins(0),
fNeventPlaneBins(0),
fMaxNumberEvts(0),
fMaxNumberTrks(0),
fNumberEvtsToMix(0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDielectronsPbPb_Data::~AliAnalysisTaskDielectronsPbPb_Data()
{
    fOutputList->Clear();
    delete fESDevent;
    delete fPoolMgr;
    delete fESDTrackCuts;
    delete fPIDResponse;
    delete fOutputList;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDielectronsPbPb_Data::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
    //PID Response & Track Cuts
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    fESDTrackCuts = new AliESDtrackCuts ("AliESDtrackCuts");

    
    //Event Statistics
    fHistoEvents = new TH1F ("fHistoEvents","",2,0,2);
    fOutputList -> Add(fHistoEvents);
    
    //Eta Bins
    Double_t eta[] = { -0.81,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.81 };
    const Int_t nEtaBins = sizeof(eta)/sizeof(Double_t)-1;
    fHistoEtaBins = new TH1F ("fHistoEtaBins","",nEtaBins,eta);

    
    //Invariant Mass Spectra (Same Event)
    fHistoInvariantMass_ULS =      new TH2F ("fHistoInvariantMass_ULS","",500,0,5,20,0,10);
    fHistoInvariantMass_ULS_Pref = new TH2F ("fHistoInvariantMass_ULS_Pref","",500,0,5,20,0,10);
    fHistoInvariantMass_PLS =      new TH2F ("fHistoInvariantMass_PLS","",500,0,5,20,0,10);
    fHistoInvariantMass_PLS_Pref = new TH2F ("fHistoInvariantMass_PLS_Pref","",500,0,5,20,0,10);
    fHistoInvariantMass_NLS =      new TH2F ("fHistoInvariantMass_NLS","",500,0,5,20,0,10);
    fHistoInvariantMass_NLS_Pref = new TH2F ("fHistoInvariantMass_NLS_Pref","",500,0,5,20,0,10);
    
    fOutputList -> Add (fHistoInvariantMass_ULS);
    fOutputList -> Add (fHistoInvariantMass_ULS_Pref);
    fOutputList -> Add (fHistoInvariantMass_PLS);
    fOutputList -> Add (fHistoInvariantMass_PLS_Pref);
    fOutputList -> Add (fHistoInvariantMass_NLS);
    fOutputList -> Add (fHistoInvariantMass_NLS_Pref);
    
    
    //Invariant Mass Spectra (Mixed Event)
    fHistoInvariantMass_EvtMixing_ULS =      new TH2F ("fHistoInvariantMass_EvtMixing_ULS","",500,0,5,20,0,10);
    fHistoInvariantMass_EvtMixing_ULS_Pref = new TH2F ("fHistoInvariantMass_EvtMixing_ULS_Pref","",500,0,5,20,0,10);
    fHistoInvariantMass_EvtMixing_PLS =      new TH2F ("fHistoInvariantMass_EvtMixing_PLS","",500,0,5,20,0,10);
    fHistoInvariantMass_EvtMixing_PLS_Pref = new TH2F ("fHistoInvariantMass_EvtMixing_PLS_Pref","",500,0,5,20,0,10);
    fHistoInvariantMass_EvtMixing_NLS =      new TH2F ("fHistoInvariantMass_EvtMixing_NLS","",500,0,5,20,0,10);
    fHistoInvariantMass_EvtMixing_NLS_Pref = new TH2F ("fHistoInvariantMass_EvtMixing_NLS_Pref","",500,0,5,20,0,10);
    
    fOutputList->Add (fHistoInvariantMass_EvtMixing_ULS);
    fOutputList->Add (fHistoInvariantMass_EvtMixing_ULS_Pref);
    fOutputList->Add (fHistoInvariantMass_EvtMixing_PLS);
    fOutputList->Add (fHistoInvariantMass_EvtMixing_PLS_Pref);
    fOutputList->Add (fHistoInvariantMass_EvtMixing_NLS);
    fOutputList->Add (fHistoInvariantMass_EvtMixing_NLS_Pref);
    
    
    //Event Mixing Setting
    SetMixingPools();
    
    
    PostData(1, fOutputList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDielectronsPbPb_Data::SetMixingPools()  {
    
    //Centrality Bins
    Double_t *CentralityBins = new Double_t[fNcentralityBins+1];
    CentralityBins[0] = fCentralityMin;
    for ( Int_t ibin=0 ; ibin<fNcentralityBins ; ibin++ )
        CentralityBins[ibin+1] = CentralityBins[ibin] + (fCentralityMax-fCentralityMin)/(Double_t)fNcentralityBins;
    
    //Vertex Bins
    Double_t *VertexBins = new Double_t[fNvertexBins+1];
    VertexBins[0] = -10.0;
    for ( Int_t ibin=0 ; ibin<fNvertexBins ; ibin++ )
        VertexBins[ibin+1] = VertexBins[ibin] + 20.0/(Double_t)fNvertexBins;
    
    //Event Plane Bins
    Double_t *Psibins = new Double_t[fNeventPlaneBins+1];
    Psibins[0] = -TMath::Pi()/2.0;
    for ( Int_t ibin=0 ; ibin<fNeventPlaneBins ; ibin++ )
        Psibins[ibin+1] = Psibins[ibin] + TMath::Pi()/(Double_t)fNeventPlaneBins;
    
    //Event Pool Manager
    fPoolMgr = new AliEventPoolManager (fMaxNumberEvts,fMaxNumberTrks,fNcentralityBins,(Double_t*)CentralityBins,fNvertexBins,(Double_t*)VertexBins,fNeventPlaneBins,(Double_t*)Psibins);
    
    delete[] CentralityBins;
    delete[] VertexBins;
    delete[] Psibins;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDielectronsPbPb_Data::UserExec(Option_t *)
{
    
    //Get Input Event
    if ( !GetEvent()) return;
    
    
    //Event Multiplicity
    Double_t Ntrk = GetEventMultiplicity();
    
    
    //Track Arrays
    Int_t Nelec_before(0);
    Int_t Nelec_after(0);
    AliESDtrack *electron_before[200];
    AliESDtrack *electron_after[200];
    
    //Track Selection
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {
        
        AliESDtrack *esdtrack1 = (AliESDtrack*) fESDevent -> GetTrack(i);
        if ( !esdtrack1 ) continue;
        if ( !PassedTrackQualityCuts (esdtrack1)) continue;
        if ( !PassedPIDcuts (esdtrack1,Ntrk)) continue;

        //Save Track
        electron_before[Nelec_before] = esdtrack1;
        Nelec_before++;
        
        //Pre-filtering
        Bool_t IsConversionCandidate = false;
        for (Int_t j=0 ; j<fESDevent->GetNumberOfTracks() ; j++)  {
            
            if (i==j) continue;//Skip Self-pairing

            AliESDtrack *esdtrack2 = (AliESDtrack*) fESDevent -> GetTrack(j);
            if ( !esdtrack2 ) continue;
            if ( !PassedLooseTrackQualityCuts (esdtrack2)) continue;
            
            Double_t mass = GetMass (esdtrack1,esdtrack2);
            Double_t phiV = GetPhiV (esdtrack1,esdtrack2);
            
            //Pre-filtering Cuts
            if ( mass<=fMassLim && phiV<=fPhivLim)          { IsConversionCandidate = true; break; }
            if ( mass<=fMassLim && phiV>=(180.0-fPhivLim))  { IsConversionCandidate = true; break; }
        }
        
        if (IsConversionCandidate) continue;//Skip Conversion Candidates
        
        electron_after[Nelec_after] = esdtrack1;
        Nelec_after++;
    }
    
    
    //Dielectron Spectra Before Pre-filtering
    for ( Int_t i=0 ; i<Nelec_before ; i++ )  {
        for ( Int_t j=i+1 ; j<Nelec_before ; j++ )  {
        
            Short_t q1 = electron_before[i]->Charge();
            Short_t q2 = electron_before[j]->Charge();
            
            Double_t mass = GetMass (electron_before[i],electron_before[j]);
            Double_t ptee = GetPtee (electron_before[i],electron_before[j]);
            
            //Invariant Mass Spectra
            if ( q1*q2<0 )      fHistoInvariantMass_ULS -> Fill(mass,ptee);
            if ( q1>0 && q2>0 ) fHistoInvariantMass_PLS -> Fill(mass,ptee);
            if ( q1<0 && q2<0 ) fHistoInvariantMass_NLS -> Fill(mass,ptee);
        }
    }
    
    
    //Dielectron Spectra After Pre-filtering
    for ( Int_t i=0 ; i<Nelec_after ; i++ )  {
        for ( Int_t j=i+1 ; j<Nelec_after ; j++ )  {
            
            Short_t q1 = electron_after[i]->Charge();
            Short_t q2 = electron_after[j]->Charge();
            
            Double_t mass = GetMass (electron_after[i],electron_after[j]);
            Double_t ptee = GetPtee (electron_after[i],electron_after[j]);

            //Invariant Mass Spectra
            if ( q1*q2<0 )      fHistoInvariantMass_ULS_Pref -> Fill(mass,ptee);
            if ( q1>0 && q2>0 ) fHistoInvariantMass_PLS_Pref -> Fill(mass,ptee);
            if ( q1<0 && q2<0 ) fHistoInvariantMass_NLS_Pref -> Fill(mass,ptee);
        }
    }

    
    
    
    //Event Mixing Settings
    AliEventPool *pool = fPoolMgr->GetEventPool(Centrality(),ZVertex(),EventPlane());
    pool -> SetTargetEvents(fNumberEvtsToMix);
    TObjArray *TrackArray = new TObjArray();
    TrackArray->SetOwner(true);
    
    
    //Loop over i-th event (Current)
    for (Int_t i=0 ; i<Nelec_before ; i++)  {
        
        //Tag Conversion Candidates
        Bool_t IsConversionCandidate = TaggedByPrefilter (electron_before[i],electron_after,Nelec_after);
        
        //Save Array
        TrackArray->Add (new AliDielectronReducedTrack (electron_before[i]->Px(),electron_before[i]->Py(),electron_before[i]->Pz(),electron_before[i]->Charge(),IsConversionCandidate));
        if(!pool->IsReady()) continue;
        TObjArray *MixedEventsTrackArray = new TObjArray();
        
        //Loop over j-th event
        for (Int_t iEvent=0; iEvent<pool->GetCurrentNEvents(); iEvent++)  {
            
            MixedEventsTrackArray = pool->GetEvent(iEvent);//Track Array from j-th event
            
            for ( Int_t j=0 ; j<MixedEventsTrackArray->GetEntriesFast() ; j++ )  {
                
                AliDielectronReducedTrack *track_mix = (AliDielectronReducedTrack*) MixedEventsTrackArray->At(j);
                
                Double_t mass = GetMass (electron_before[i],track_mix);
                Double_t ptee = GetPtee (electron_before[i],track_mix);
               
                Short_t q1 = electron_before[i] -> Charge();
                Short_t q2 = track_mix -> Charge();
                
                //Dielectron Spectra Before Pre-filtering
                if ( q1*q2<0 )       fHistoInvariantMass_EvtMixing_ULS -> Fill(mass,ptee);
                if ( q1>0 && q2>0 )  fHistoInvariantMass_EvtMixing_PLS -> Fill(mass,ptee);
                if ( q1<0 && q2<0 )  fHistoInvariantMass_EvtMixing_NLS -> Fill(mass,ptee);
                
                //Skip Conversion Candidates
                if (IsConversionCandidate) continue;
                if (track_mix->IsConversionCandidate()) continue;
                
                //Dielectron Spectra After pre-filtering
                if ( q1*q2<0 )       fHistoInvariantMass_EvtMixing_ULS_Pref -> Fill(mass,ptee);
                if ( q1>0 && q2>0 )  fHistoInvariantMass_EvtMixing_PLS_Pref -> Fill(mass,ptee);
                if ( q1<0 && q2<0 )  fHistoInvariantMass_EvtMixing_NLS_Pref -> Fill(mass,ptee);
            }
        }
    }
    
    
    pool->UpdatePool(TrackArray);
    PostData(1, fOutputList);

}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Data::GetEvent ()  {
    
    //Get Input Event
    fESDevent = dynamic_cast <AliESDEvent*> (InputEvent());
    if (!fESDevent) return false;
    
    fHistoEvents -> Fill(0.5);
    
    //Trigger
    if (fCentralityMax <= 10 && !(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral)) return false;
    if (fCentralityMin >= 10 && !(fESDevent->GetTriggerMask() & (ULong64_t(1)<<7))) return false;
    
    //Centrality Range
    AliCentrality *centrality = fESDevent->GetCentrality();
    if (!centrality) return false;
    Double_t centralityPerc = centrality->GetCentralityPercentile("V0M");//V0A + V0C (Best Resolution)
    if (centralityPerc<fCentralityMin || centralityPerc>=fCentralityMax ) return false;
    
    //Primary Vertex
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    if ( TMath::Abs(vertex->GetZ() ) > 10.0 ) return false;
    if ( vertex->GetNContributors() < 1 ) return false;
    
    //Event Plane
    AliEventplane *eventPlane =  fESDevent->GetEventplane();
    if (!eventPlane) return false;
    
    fHistoEvents -> Fill(1.5);
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetEventMultiplicity()  {
    
    Double_t mult(0);
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {
        
        AliESDtrack *esdtrack = (AliESDtrack*) fESDevent->GetTrack(i);
        if ( !esdtrack ) continue;
        if ( !PassedTrackQualityCutsMultiplicity (esdtrack)) continue;
        mult++;
    }
    
    return mult;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::Centrality()  {
    
    AliCentrality *centrality = fESDevent->GetCentrality();
    Double_t centralityPerc = centrality->GetCentralityPercentile("V0M");
    return centralityPerc;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::EventPlane()  {
    
    AliEventplane *eventPlane =  fESDevent->GetEventplane();
    Double_t Psi_EP = eventPlane -> GetEventplane("V0",fESDevent, 2);
    return Psi_EP;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::ZVertex()  {
    
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    Double_t zVtx = vertex->GetZ();
    return zVtx;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Data::TaggedByPrefilter (AliESDtrack *track_before, AliESDtrack *track_after[], Int_t Nelec_After)  {
    
    Bool_t IsConversionCandidate = true;
    
    for ( Int_t i=0 ; i<Nelec_After ; i++ )
        if (track_before->GetID() == track_after[i]->GetID()) IsConversionCandidate = false;
    
    return IsConversionCandidate;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Data::PassedTrackQualityCutsMultiplicity (AliESDtrack* track)  {
    
    fESDTrackCuts -> SetMaxDCAToVertexZ(0.1);
    fESDTrackCuts -> SetMaxDCAToVertexXY(0.00515869 + 0.0101668/TMath::Power(track->Pt(),1.34489));
    fESDTrackCuts -> SetAcceptKinkDaughters(false);
    fESDTrackCuts -> SetMinNCrossedRowsTPC(100);
    fESDTrackCuts -> SetMinNClustersITS(4);
    fESDTrackCuts -> SetRequireTPCRefit(true);
    fESDTrackCuts -> SetRequireITSRefit(true);
    fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    fESDTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fESDTrackCuts -> SetMaxChi2TPCConstrainedGlobal(36);
    fESDTrackCuts -> SetMaxChi2PerClusterTPC(4.0);
    fESDTrackCuts -> SetMaxChi2PerClusterITS(36);
    fESDTrackCuts -> SetEtaRange(-0.8, 0.8);
    fESDTrackCuts -> SetPtRange(0.4,5.0);
    
    if ( track -> GetTPCsignalN()<50 ) return false;
    if ( fESDTrackCuts->AcceptTrack (track) ) return true;
    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Data::PassedLooseTrackQualityCuts (AliESDtrack* track)  {
    
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
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Data::PassedTrackQualityCuts (AliESDtrack* track)  {
    
    fESDTrackCuts -> SetMaxDCAToVertexZ(fDCAz_max);
    fESDTrackCuts -> SetMaxDCAToVertexXY(fDCAxy_param0 + fDCAxy_param1/TMath::Power(track->Pt(),fDCAxy_param2));
    fESDTrackCuts -> SetAcceptKinkDaughters(false);
    fESDTrackCuts -> SetMinNCrossedRowsTPC(fTPC_minCr);
    fESDTrackCuts -> SetMinNClustersITS(fITS_minNcls);
    fESDTrackCuts -> SetMinNClustersTPC(fTPC_minNcls);
    fESDTrackCuts -> SetRequireTPCRefit(true);
    fESDTrackCuts -> SetRequireITSRefit(true);
    
    //ITS Requirement
    if (strcmp(fITSreq,"kOff")==0)        fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
    if (strcmp(fITSreq,"kNone")==0)       fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
    if (strcmp(fITSreq,"kAny")==0)        fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    if (strcmp(fITSreq,"kFirst")==0)      fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kFirst);
    if (strcmp(fITSreq,"kOnlyFirst")==0)  fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOnlyFirst);
    if (strcmp(fITSreq,"kSecond")==0)     fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kSecond);
    if (strcmp(fITSreq,"kOnlySecond")==0) fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOnlySecond);
    if (strcmp(fITSreq,"kBoth")==0)       fESDTrackCuts -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kBoth);
    
    fESDTrackCuts -> SetMinRatioCrossedRowsOverFindableClustersTPC(fMinCrOverFindableCls);
    fESDTrackCuts -> SetMaxChi2TPCConstrainedGlobal(fMaxGoldenChi2);
    fESDTrackCuts -> SetMaxChi2PerClusterTPC(fMaxTPCchi2);
    fESDTrackCuts -> SetMaxChi2PerClusterITS(fMaxITSchi2);
    fESDTrackCuts -> SetEtaRange(-fEtaLim,fEtaLim);
    fESDTrackCuts -> SetPtRange(fPtMin,fPtMax);
    
    if ( track -> GetTPCsignalN() < fTPC_nClsdEdx ) return false;
    if ( FractionSharedClsITS(track) > fMaxFracSharedCls ) return false;
    if ( fESDTrackCuts->AcceptTrack (track) ) return true;
    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::FractionSharedClsITS (AliESDtrack *track)  {
    
    Double_t nSharedCls(0);
    Double_t fSharedCls(0);
    
    for ( Int_t i=0 ; i<6 ; i++ )
        if ( track->HasPointOnITSLayer(i) && track->HasSharedPointOnITSLayer(i) ) nSharedCls++;
    
    fSharedCls = nSharedCls/(Double_t)track->GetNcls(0);
    
    return fSharedCls;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDielectronsPbPb_Data::PassedPIDcuts (AliESDtrack* track, Double_t Ntrk )  {
        
    //TOF Response
    if ( !fPIDResponse -> CheckPIDStatus(AliPIDResponse::kTOF,track)) return false;
    
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kElectron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kElectron);
    Double_t nsigmaTPC = GetCorrectedTPCResponse (track,Ntrk);//Eta & Multiplicity Correction
    
    if ( TMath::Abs(nsigmaTOF ) > fnsigmaTOF_max )            return false;
    if ( nsigmaITS > fnsigmaITS_max )                         return false;
    if ( nsigmaTPC < fnsigmaTPC_min*TMath::Exp(-track->P()) ) return false;
    if ( nsigmaTPC > fnsigmaTPC_max )                         return false;
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetCorrectedTPCResponse (AliESDtrack *track, Double_t Ntrk)  {
    
    Double_t nsigmaTPC(0);
    nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kElectron);
    
    //Get Mean And Width
    Double_t mean  = GetMeanTPCnsigma  (track->Eta(),Ntrk);
    Double_t width = GetWidthTPCnsigma (track->Eta(),Ntrk);
    
    //Re-calibration
    nsigmaTPC = (nsigmaTPC-mean)/width;
    
    return nsigmaTPC;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetMeanTPCnsigma (Double_t eta,Double_t Ntrk)  {
    
    //Parametrizations
    Double_t mean[16];
    
    mean[0] =  -0.1855740 -0.000109412 * Ntrk;
    mean[1] =   0.0466150 -0.000753731 * Ntrk;
    mean[2] =  -0.1390870 -0.000620659 * Ntrk;
    mean[3] =  -0.3199010 -0.000427379 * Ntrk;
    mean[4] =  -0.0194713 -0.000772120 * Ntrk;
    mean[5] =   0.1154030 -0.000770983 * Ntrk;
    mean[6] =   0.2086960 -0.000767416 * Ntrk;
    mean[7] =   0.3132710 -0.000756067 * Ntrk;
    mean[8] =   0.3680350 -0.000745087 * Ntrk;
    mean[9] =   0.4594220 -0.000996809 * Ntrk;
    mean[10] =  0.3020870 -0.000938410 * Ntrk;
    mean[11] =  0.0699430 -0.000796447 * Ntrk;
    mean[12] = -0.1645030 -0.000537301 * Ntrk;
    mean[13] = -0.0862225 -0.000619860 * Ntrk;
    mean[14] = -0.1579340 -0.000294009 * Ntrk;
    mean[15] =  0.2577600 -0.000644332 * Ntrk;
    
    Int_t ieta = fHistoEtaBins->FindBin(eta)-1;
    
    return mean[ieta];
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetWidthTPCnsigma (Double_t eta,Double_t Ntrk)  {
    
    //Parametrizations
    Double_t width[16];
    
    width[0] =  1.014180 + 0.0003428740 * Ntrk;
    width[1] =  1.013650 + 0.0003123050 * Ntrk;
    width[2] =  1.075720 + 0.0001827030 * Ntrk;
    width[3] =  1.078130 + 0.0001200240 * Ntrk;
    width[4] =  0.977925 + 0.0003591680 * Ntrk;
    width[5] =  1.100640 + 0.0002447050 * Ntrk;
    width[6] =  0.813449 + 0.0007106460 * Ntrk;
    width[7] =  1.020480 + 0.0004094130 * Ntrk;
    width[8] =  1.209780 + 0.0000999875 * Ntrk;
    width[9] =  1.135020 + 0.0002038660 * Ntrk;
    width[10] = 1.066210 + 0.0002953440 * Ntrk;
    width[11] = 1.033010 + 0.0002973250 * Ntrk;
    width[12] = 0.936047 + 0.0004247790 * Ntrk;
    width[13] = 0.812811 + 0.0006544050 * Ntrk;
    width[14] = 1.117030 + 0.0001438510 * Ntrk;
    width[15] = 1.200970 + 0.0000283310 * Ntrk;
    
    Int_t ieta = fHistoEtaBins->FindBin(eta)-1;
    
    return width[ieta];
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetMass (AliESDtrack *track1, AliESDtrack *track2) {
    
    TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
    TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
    
    Double_t electron_mass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    
    Double_t E1 = TMath::Sqrt( TMath::Power(electron_mass,2) + P1.Mag2());
    Double_t E2 = TMath::Sqrt( TMath::Power(electron_mass,2) + P2.Mag2());
    Double_t mass = TMath::Sqrt( TMath::Power(E1+E2,2) - (P1+P2).Mag2());
    
    return mass;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetMass( AliESDtrack *track1, AliDielectronReducedTrack *track2 ) {
    
    TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
    TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
    
    Double_t electron_mass = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    
    Double_t E1 = TMath::Sqrt( TMath::Power(electron_mass,2) + P1.Mag2());
    Double_t E2 = TMath::Sqrt( TMath::Power(electron_mass,2) + P2.Mag2());
    Double_t mass = TMath::Sqrt( TMath::Power(E1+E2,2) - (P1+P2).Mag2());
    
    return mass;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetPtee (AliESDtrack *track1,AliESDtrack *track2)  {
    
    TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
    TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
    
    Double_t pt = (P1+P2).Pt();
    return pt;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetPtee (AliESDtrack *track1,AliDielectronReducedTrack *track2)  {
    
    TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
    TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
    
    Double_t pt = (P1+P2).Pt();
    return pt;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDielectronsPbPb_Data::GetPhiV (AliESDtrack *track1, AliESDtrack *track2 ) {
    
    TVector3 P1,P2;
    Double_t r(0);
    
    //Randomization of pair ordering (Symmetric Peaks)
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
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDielectronsPbPb_Data::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

