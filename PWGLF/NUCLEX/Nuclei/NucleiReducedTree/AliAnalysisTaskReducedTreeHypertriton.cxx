#include "AliAnalysisTaskReducedTreeHypertriton.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "TDatabasePDG.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliAODv0.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"

ClassImp(AliAnalysisTaskReducedTreeHypertriton)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeHypertriton::AliAnalysisTaskReducedTreeHypertriton():
AliAnalysisTaskSE(),
fAODevent(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
reducedTree_HyperTriton(NULL),
fcentralityMin(0),
fcentralityMax(0),
centrality(0),
px_Daughter1(0),
py_Daughter1(0),
pz_Daughter1(0),
q_Daughter1(0),
dcaxy_Daughter1(0),
nTPC_Clusters_Daughter1(0),
nTPC_FindableClusters_Daughter1(0),
nTPC_CrossedRows_Daughter1(0),
nTPC_Clusters_dEdx_Daughter1(0),
chi2_TPC_Daughter1(0),
nSigmaTPC_He3_Daughter1(0),
nSigmaTPC_Pion_Daughter1(0),
px_Daughter2(0),
py_Daughter2(0),
pz_Daughter2(0),
q_Daughter2(0),
dcaxy_Daughter2(0),
nTPC_Clusters_Daughter2(0),
nTPC_FindableClusters_Daughter2(0),
nTPC_CrossedRows_Daughter2(0),
nTPC_Clusters_dEdx_Daughter2(0),
chi2_TPC_Daughter2(0),
nSigmaTPC_He3_Daughter2(0),
nSigmaTPC_Pion_Daughter2(0),
isOnTheFlyV0(0),
cosPointingAngle(0),
dcaV0Daughters(0),
dcaV0ToVertex(0),
radius(0),
decayLength(0)
{}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeHypertriton::AliAnalysisTaskReducedTreeHypertriton(const char *name):
AliAnalysisTaskSE(name),
fAODevent(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
reducedTree_HyperTriton(NULL),
fcentralityMin(0),
fcentralityMax(0),
centrality(0),
px_Daughter1(0),
py_Daughter1(0),
pz_Daughter1(0),
q_Daughter1(0),
dcaxy_Daughter1(0),
nTPC_Clusters_Daughter1(0),
nTPC_FindableClusters_Daughter1(0),
nTPC_CrossedRows_Daughter1(0),
nTPC_Clusters_dEdx_Daughter1(0),
chi2_TPC_Daughter1(0),
nSigmaTPC_He3_Daughter1(0),
nSigmaTPC_Pion_Daughter1(0),
px_Daughter2(0),
py_Daughter2(0),
pz_Daughter2(0),
q_Daughter2(0),
dcaxy_Daughter2(0),
nTPC_Clusters_Daughter2(0),
nTPC_FindableClusters_Daughter2(0),
nTPC_CrossedRows_Daughter2(0),
nTPC_Clusters_dEdx_Daughter2(0),
chi2_TPC_Daughter2(0),
nSigmaTPC_He3_Daughter2(0),
nSigmaTPC_Pion_Daughter2(0),
isOnTheFlyV0(0),
cosPointingAngle(0),
dcaV0Daughters(0),
dcaV0ToVertex(0),
radius(0),
decayLength(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1,  TList::Class());
    DefineOutput(2,  TList::Class());
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeHypertriton::~AliAnalysisTaskReducedTreeHypertriton()
{
    fOutputList->Clear();
    delete fAODevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    delete fQAList;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeHypertriton::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
    fQAList = new TList();
    fQAList -> SetOwner();
    
    fAODeventCuts.AddQAplotsToList(fQAList); ///Add event selection QA plots
    
    
    //Histogram with number of events
    hEvents = new TH1F ("hEvents","",10,0,10);
    fOutputList->Add(hEvents);
    
    
    //Reduced Tree HyperTriton
    reducedTree_HyperTriton = new TTree("reducedTree_HyperTriton","reducedTree_HyperTriton");
    reducedTree_HyperTriton -> Branch("centrality",&centrality,"centrality/D");
    reducedTree_HyperTriton -> Branch("px_Daughter1",&px_Daughter1,"px_Daughter1/D");
    reducedTree_HyperTriton -> Branch("py_Daughter1",&py_Daughter1,"py_Daughter1/D");
    reducedTree_HyperTriton -> Branch("pz_Daughter1",&pz_Daughter1,"pz_Daughter1/D");
    reducedTree_HyperTriton -> Branch("q_Daughter1",&q_Daughter1,"q_Daughter1/I");
    reducedTree_HyperTriton -> Branch("dcaxy_Daughter1",&dcaxy_Daughter1,"dcaxy_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_Daughter1",&nTPC_Clusters_Daughter1,"nTPC_Clusters_Daughter1/I");
    reducedTree_HyperTriton -> Branch("nTPC_FindableClusters_Daughter1",&nTPC_FindableClusters_Daughter1,"nTPC_FindableClusters_Daughter1/I");
    reducedTree_HyperTriton -> Branch("nTPC_CrossedRows_Daughter1",&nTPC_CrossedRows_Daughter1,"nTPC_CrossedRows_Daughter1/I");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_dEdx_Daughter1",&nTPC_Clusters_dEdx_Daughter1,"nTPC_Clusters_dEdx_Daughter1/I");
    reducedTree_HyperTriton -> Branch("chi2_TPC_Daughter1",&chi2_TPC_Daughter1,"chi2_TPC_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_He3_Daughter1",&nSigmaTPC_He3_Daughter1,"nSigmaTPC_He3_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_Pion_Daughter1",&nSigmaTPC_Pion_Daughter1,"nSigmaTPC_Pion_Daughter1/D");
    reducedTree_HyperTriton -> Branch("px_Daughter2",&px_Daughter2,"px_Daughter2/D");
    reducedTree_HyperTriton -> Branch("py_Daughter2",&py_Daughter2,"py_Daughter2/D");
    reducedTree_HyperTriton -> Branch("pz_Daughter2",&pz_Daughter2,"pz_Daughter2/D");
    reducedTree_HyperTriton -> Branch("q_Daughter2",&q_Daughter2,"q_Daughter2/I");
    reducedTree_HyperTriton -> Branch("dcaxy_Daughter2",&dcaxy_Daughter2,"dcaxy_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_Daughter2",&nTPC_Clusters_Daughter2,"nTPC_Clusters_Daughter2/I");
    reducedTree_HyperTriton -> Branch("nTPC_FindableClusters_Daughter2",&nTPC_FindableClusters_Daughter2,"nTPC_FindableClusters_Daughter2/I");
    reducedTree_HyperTriton -> Branch("nTPC_CrossedRows_Daughter2",&nTPC_CrossedRows_Daughter2,"nTPC_CrossedRows_Daughter2/I");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_dEdx_Daughter2",&nTPC_Clusters_dEdx_Daughter2,"nTPC_Clusters_dEdx_Daughter2/I");
    reducedTree_HyperTriton -> Branch("chi2_TPC_Daughter2",&chi2_TPC_Daughter2,"chi2_TPC_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_He3_Daughter2",&nSigmaTPC_He3_Daughter2,"nSigmaTPC_He3_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_Pion_Daughter2",&nSigmaTPC_Pion_Daughter2,"nSigmaTPC_Pion_Daughter2/D");
    reducedTree_HyperTriton -> Branch("isOnTheFlyV0",&isOnTheFlyV0,"isOnTheFlyV0/I");
    reducedTree_HyperTriton -> Branch("cosPointingAngle",&cosPointingAngle,"cosPointingAngle/D");
    reducedTree_HyperTriton -> Branch("dcaV0Daughters",&dcaV0Daughters,"dcaV0Daughters/D");
    reducedTree_HyperTriton -> Branch("dcaV0ToVertex",&dcaV0ToVertex,"dcaV0ToVertex/D");
    reducedTree_HyperTriton -> Branch("radius",&radius,"radius/D");
    reducedTree_HyperTriton -> Branch("decayLength",&decayLength,"decayLength/D");
    fOutputList -> Add(reducedTree_HyperTriton);
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeHypertriton::UserExec(Option_t *)
{
    //Get Input Event
    if ( !GetInputEvent ()) return;
    
   
    //Primary Vertex
    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    Double_t vertexPosition[3] = { 0.0, 0.0, 0.0 };
    vertex -> GetXYZ (vertexPosition);
    
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    centrality = multiplicitySelection->GetMultiplicityPercentile("V0M");
    
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse) {
        AliError("No PID Response found");
        return;
    }
    
    
    //Magnetic Field
    //if(fAODevent->GetMagneticField() < 0) magFieldSign = -1;
    //if(fAODevent->GetMagneticField() > 0) magFieldSign = 1;
    
    
    for ( Int_t iV0=0 ; iV0<fAODevent->GetNumberOfV0s() ; iV0++ ) {
        
        //Get V0 candidate
        AliAODv0 *V0 = (AliAODv0*)fAODevent->GetV0(iV0);
        if (!V0) continue;
        if ( V0->GetOnFlyStatus()) isOnTheFlyV0=1; //V0-Online
        if (!V0->GetOnFlyStatus()) isOnTheFlyV0=0; //V0-Offline

        
        //Get V0 Daughters
        AliAODTrack *posTrack = (AliAODTrack*)(V0->GetDaughter(0));
        AliAODTrack *negTrack = (AliAODTrack*)(V0->GetDaughter(1));
        if (!posTrack) continue;
        if (!negTrack) continue;
        if (posTrack->Charge() == negTrack->Charge()) continue;

        //Quality Requirements
        if (!PassedBasicTrackQualityCuts (posTrack)) continue;
        if (!PassedBasicTrackQualityCuts (negTrack)) continue;
        if (!PassedV0QualityCuts(V0)) continue;
        if (!IsHyperTritonCandidate(V0)) continue;

        //Daughter1 (Positive Charge)
        px_Daughter1                    = V0->MomPosX();
        py_Daughter1                    = V0->MomPosY();
        pz_Daughter1                    = V0->MomPosZ();
        q_Daughter1                     = (Int_t) posTrack -> Charge();
        dcaxy_Daughter1                 = V0->DcaPosToPrimVertex();
        nTPC_Clusters_Daughter1         = posTrack->GetTPCNcls();
        nTPC_FindableClusters_Daughter1 = posTrack->GetTPCNclsF();
        nTPC_CrossedRows_Daughter1      = posTrack->GetTPCNCrossedRows();
        nTPC_Clusters_dEdx_Daughter1    = posTrack -> GetTPCsignalN();
        chi2_TPC_Daughter1              = posTrack -> GetTPCchi2();
        nSigmaTPC_He3_Daughter1         = fPIDResponse -> NumberOfSigmasTPC (posTrack,AliPID::kHe3);
        nSigmaTPC_Pion_Daughter1        = fPIDResponse -> NumberOfSigmasTPC (posTrack,AliPID::kPion);

        //Daughter2  (Negative Charge)
        px_Daughter2                    = V0->MomNegX();
        py_Daughter2                    = V0->MomNegY();
        pz_Daughter2                    = V0->MomNegZ();
        q_Daughter2                     = (Int_t) negTrack -> Charge();
        dcaxy_Daughter2                 = V0->DcaNegToPrimVertex();
        nTPC_Clusters_Daughter2         = negTrack->GetTPCNcls();
        nTPC_FindableClusters_Daughter2 = negTrack->GetTPCNclsF();
        nTPC_CrossedRows_Daughter2      = negTrack->GetTPCNCrossedRows();
        nTPC_Clusters_dEdx_Daughter2    = negTrack -> GetTPCsignalN();
        chi2_TPC_Daughter2              = negTrack -> GetTPCchi2();
        nSigmaTPC_He3_Daughter2         = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kHe3);
        nSigmaTPC_Pion_Daughter2        = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kPion);

        //Pair Variables
        cosPointingAngle = V0->CosPointingAngle(vertexPosition);
        dcaV0Daughters   = V0->DcaV0Daughters();
        dcaV0ToVertex    = V0->DcaV0ToPrimVertex();
        radius           = V0->RadiusV0();
        decayLength      = V0->DecayLengthV0(vertexPosition);
        
        //Fill Reduced Tree
        reducedTree_HyperTriton -> Fill();
    }
    
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::GetInputEvent ()  {
    
    //Get Input Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return false;
    hEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fAODeventCuts.AcceptEvent(fAODevent)) {
        PostData(2, fQAList);
        return false;
    }
    hEvents -> Fill(1.5);
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    hEvents -> Fill(2.5);
    Double_t centralityPercentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    
    //Selection of Centrality Range
    if (centralityPercentile<fcentralityMin || centralityPercentile>=fcentralityMax ) return false;
    hEvents -> Fill(3.5);
    
    //Primary Vertex
    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    hEvents -> Fill(4.5);
    
    //Primary Vertex Selection
    if ( vertex->GetZ() < -10.0 ) return false;
    if ( vertex->GetZ() > +10.0 ) return false;
    hEvents -> Fill(5.5);
    
    if ( vertex->GetNContributors() < 1 ) return false;
    hEvents -> Fill(6.5);
    
    //Event Plane
    AliEventplane *eventPlane =  fAODevent->GetEventplane();
    if (!eventPlane) return false;
    hEvents -> Fill(7.5);
    
    return true;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::PassedBasicTrackQualityCuts (AliAODTrack *track)  {
  
    //Filter Bit
    if(!track->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return false;

    //Kinematic Cuts & Acceptance
    if ( track->Pt()<0.2 || track->Pt()>100.0 ) return false;
    if ( TMath::Abs(track->Eta()) > 1.2 )       return false;
    
    //Track Selection Cuts
    if ( track->GetTPCNcls() < 60 ) return false;
    if ( track->GetTPCNclsF() == 0) return false;
    if ( static_cast<Double_t>(track->GetTPCNCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < 0.5) return false;
    if ( track->GetTPCsignalN() < 40 ) return false;
    
    return true;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::PassedV0QualityCuts (AliAODv0 *V0)  {
 
    
    //Position of Primary Vertex
    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    Double_t vertexPosition[3] = { 0.0, 0.0, 0.0 };
    vertex->GetXYZ(vertexPosition);
    
    //if (V0->Chi2V0()>10.0) return false;
    if (V0->DcaV0Daughters()>2.0) return false;
    if (V0->RadiusV0()<3.0) return false;
    if (V0->CosPointingAngle(vertexPosition)<0.5) return false;
    
    return true;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::IsHyperTritonCandidate (AliAODv0 *V0)  {
    
    //Get V0 Daughters
    AliAODTrack *track0 = (AliAODTrack*)(V0->GetDaughter(0));
    AliAODTrack *track1 = (AliAODTrack*)(V0->GetDaughter(1));
    
    //Pair Requirements
    if ( IsPionCandidate (track0) && (!Is3HeCandidate (track1)))  return false;
    if ( Is3HeCandidate  (track0) && (!IsPionCandidate (track1))) return false;

    return true;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::IsPionCandidate (AliAODTrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
    if (TMath::Abs(nsigmaTPC) > 4.0)        return false;
    if (track->Pt()<0.1 || track->Pt()>1.2) return false;
    
    return true;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::Is3HeCandidate (AliAODTrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTPC) > 4.0) return false;
    if (2.0*track->Pt()<1.5)         return false;
    
    return true;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeHypertriton::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_____________________________________________________________________________________________________________________________________

