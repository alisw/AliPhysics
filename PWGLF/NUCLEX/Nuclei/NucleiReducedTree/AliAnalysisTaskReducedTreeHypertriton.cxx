#include "AliAnalysisTaskReducedTreeHypertriton.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliEventCuts.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliESDv0.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"

ClassImp(AliAnalysisTaskReducedTreeHypertriton)

//_____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeHypertriton::AliAnalysisTaskReducedTreeHypertriton():
AliAnalysisTaskSE(),
fESDevent(NULL),
fPIDResponse(NULL),
fESDtrackCuts_Pos(NULL),
fESDtrackCuts_Neg(NULL),
fESDeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
hEvents(NULL),
reducedTree_HyperTriton(NULL),
fCentralityMin(0),
fCentralityMax(90),
centrality(0),
px_Daughter1(0),
py_Daughter1(0),
pz_Daughter1(0),
q_Daughter1(0),
dcaxy_Daughter1(0),
nTPC_Clusters_Daughter1(0),
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
nTPC_Clusters_dEdx_Daughter2(0),
chi2_TPC_Daughter2(0),
nSigmaTPC_He3_Daughter2(0),
nSigmaTPC_Pion_Daughter2(0),
isOnTheFlyV0(0),
cosPointingAngle(0),
dcaV0Daughters(0),
radius(0),
chi2V0(0),
decayLength(0)
{}
//_____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeHypertriton::AliAnalysisTaskReducedTreeHypertriton(const char *name):
AliAnalysisTaskSE(name),
fESDevent(NULL),
fPIDResponse(NULL),
fESDtrackCuts_Pos(NULL),
fESDtrackCuts_Neg(NULL),
fESDeventCuts(),
fUtils(NULL),
fOutputList(NULL),
fQAList(NULL),
hEvents(NULL),
reducedTree_HyperTriton(NULL),
fCentralityMin(0),
fCentralityMax(90),
centrality(0),
px_Daughter1(0),
py_Daughter1(0),
pz_Daughter1(0),
q_Daughter1(0),
dcaxy_Daughter1(0),
nTPC_Clusters_Daughter1(0),
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
nTPC_Clusters_dEdx_Daughter2(0),
chi2_TPC_Daughter2(0),
nSigmaTPC_He3_Daughter2(0),
nSigmaTPC_Pion_Daughter2(0),
isOnTheFlyV0(0),
cosPointingAngle(0),
dcaV0Daughters(0),
radius(0),
chi2V0(0),
decayLength(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1,  TList::Class());
    DefineOutput(2,  TList::Class());
}
//_____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeHypertriton::~AliAnalysisTaskReducedTreeHypertriton()
{
    fOutputList->Clear();
    delete fESDevent;
    delete fPIDResponse;
    delete fESDtrackCuts_Pos;
    delete fESDtrackCuts_Neg;
    delete fUtils;
    delete fOutputList;
    delete fQAList;
}
//_____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeHypertriton::UserCreateOutputObjects()
{
    //Output List
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
    //QA List
    fQAList = new TList();
    fQAList -> SetOwner();
    
    //Add Event Selection QA Plots
    fESDeventCuts.AddQAplotsToList(fQAList);
    
    
    //Histogram with number of events
    hEvents = new TH1F ("hEvents","",10,0,10);
    fOutputList -> Add(hEvents);
    
    
    //Reduced Tree HyperTriton
    reducedTree_HyperTriton = new TTree("reducedTree_HyperTriton","reducedTree_HyperTriton");
    reducedTree_HyperTriton -> Branch("centrality",&centrality,"centrality/D");
    reducedTree_HyperTriton -> Branch("px_Daughter1",&px_Daughter1,"px_Daughter1/D");
    reducedTree_HyperTriton -> Branch("py_Daughter1",&py_Daughter1,"py_Daughter1/D");
    reducedTree_HyperTriton -> Branch("pz_Daughter1",&pz_Daughter1,"pz_Daughter1/D");
    reducedTree_HyperTriton -> Branch("q_Daughter1",&q_Daughter1,"q_Daughter1/I");
    reducedTree_HyperTriton -> Branch("dcaxy_Daughter1",&dcaxy_Daughter1,"dcaxy_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_Daughter1",&nTPC_Clusters_Daughter1,"nTPC_Clusters_Daughter1/I");
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
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_dEdx_Daughter2",&nTPC_Clusters_dEdx_Daughter2,"nTPC_Clusters_dEdx_Daughter2/I");
    reducedTree_HyperTriton -> Branch("chi2_TPC_Daughter2",&chi2_TPC_Daughter2,"chi2_TPC_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_He3_Daughter2",&nSigmaTPC_He3_Daughter2,"nSigmaTPC_He3_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_Pion_Daughter2",&nSigmaTPC_Pion_Daughter2,"nSigmaTPC_Pion_Daughter2/D");
    reducedTree_HyperTriton -> Branch("isOnTheFlyV0",&isOnTheFlyV0,"isOnTheFlyV0/I");
    reducedTree_HyperTriton -> Branch("cosPointingAngle",&cosPointingAngle,"cosPointingAngle/D");
    reducedTree_HyperTriton -> Branch("dcaV0Daughters",&dcaV0Daughters,"dcaV0Daughters/D");
    reducedTree_HyperTriton -> Branch("radius",&radius,"radius/D");
    reducedTree_HyperTriton -> Branch("chi2V0",&chi2V0,"chi2V0/D");
    reducedTree_HyperTriton -> Branch("decayLength",&decayLength,"decayLength/D");
    fOutputList -> Add(reducedTree_HyperTriton);
    
    
    //Track Cuts Object
    fESDtrackCuts_Pos = new AliESDtrackCuts ("fESDtrackCuts_Pos");
    fESDtrackCuts_Neg = new AliESDtrackCuts ("fESDtrackCuts_Neg");

    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeHypertriton::UserExec(Option_t *)  {
    
    //Get Input Event
    if ( !GetInputEvent ()) return;
    
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    centrality = multiplicitySelection->GetMultiplicityPercentile("V0M");

    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse) {
        AliError("No PID Response found");
        return;
    }
    
    
    //Loop Over Reconstructed V0s
    for ( Int_t iV0=0 ; iV0<fESDevent->GetNumberOfV0s() ; iV0++ ) {
        
        //Get V0 Candidate
        AliESDv0 *V0 = (AliESDv0*)fESDevent->GetV0(iV0);
        if (!V0) continue;
        if ( V0->GetOnFlyStatus()) isOnTheFlyV0=1;//V0-Online
        if (!V0->GetOnFlyStatus()) isOnTheFlyV0=0;//V0-Offline
        if (!PassedMinimalQualityCutsV0(V0)) continue;

        
        //Get V0 Daughters
        AliESDtrack *posTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
        AliESDtrack *negTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
        if (!posTrack) continue;
        if (!negTrack) continue;
        if (posTrack->Charge() == negTrack->Charge()) continue;
        if (posTrack->GetID()  == negTrack->GetID() ) continue;

        //Quality Requirements
        if (!PassedBasicTrackQualityCuts_Pos (posTrack)) continue;
        if (!PassedBasicTrackQualityCuts_Neg (negTrack)) continue;
        
        //Hypertriton Candidate Selection
        if (!IsHyperTritonCandidate(V0)) continue;

        
        //Momentum Components of V0 Daughters
        Double_t posMomentum[3] = { 0.0, 0.0, 0.0 };
        Double_t negMomentum[3] = { 0.0, 0.0, 0.0 };
        V0->GetPPxPyPz(posMomentum[0],posMomentum[1],posMomentum[2]);
        V0->GetNPxPyPz(negMomentum[0],negMomentum[1],negMomentum[2]);
        
        
        //Daughter1 (Positive Charge)
        px_Daughter1                    = posMomentum[0];
        py_Daughter1                    = posMomentum[1];
        pz_Daughter1                    = posMomentum[2];
        q_Daughter1                     = (Int_t) posTrack -> Charge();
        dcaxy_Daughter1                 = GetTransverseDCA (posTrack);
        nTPC_Clusters_Daughter1         = posTrack -> GetTPCNcls();
        nTPC_Clusters_dEdx_Daughter1    = posTrack -> GetTPCsignalN();
        chi2_TPC_Daughter1              = posTrack -> GetTPCchi2();
        nSigmaTPC_He3_Daughter1         = fPIDResponse -> NumberOfSigmasTPC (posTrack,AliPID::kHe3);
        nSigmaTPC_Pion_Daughter1        = fPIDResponse -> NumberOfSigmasTPC (posTrack,AliPID::kPion);

        //Daughter2  (Negative Charge)
        px_Daughter2                    = negMomentum[0];
        py_Daughter2                    = negMomentum[1];
        pz_Daughter2                    = negMomentum[2];
        q_Daughter2                     = (Int_t) negTrack -> Charge();
        dcaxy_Daughter2                 = GetTransverseDCA (negTrack);
        nTPC_Clusters_Daughter2         = negTrack -> GetTPCNcls();
        nTPC_Clusters_dEdx_Daughter2    = negTrack -> GetTPCsignalN();
        chi2_TPC_Daughter2              = negTrack -> GetTPCchi2();
        nSigmaTPC_He3_Daughter2         = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kHe3);
        nSigmaTPC_Pion_Daughter2        = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kPion);

        
        //Pair Variables
        cosPointingAngle = V0->GetV0CosineOfPointingAngle();
        dcaV0Daughters   = V0->GetDcaV0Daughters();
        radius           = V0->GetRr();
        chi2V0           = V0->GetChi2V0();
        decayLength      = GetDecayLengthV0 (V0);
        
        //Fill Reduced Tree
        reducedTree_HyperTriton -> Fill();
    }
    
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::GetInputEvent ()  {
    
    //Get Input Event
    fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDevent) return false;
    hEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fESDeventCuts.AcceptEvent(fESDevent)) {
        PostData(2, fQAList);
        return false;
    }
    hEvents -> Fill(1.5);
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    hEvents -> Fill(2.5);
    Double_t centralityPerc = multiplicitySelection->GetMultiplicityPercentile("V0M");
    
    //Selection of Centrality Range
    if (centralityPerc<fCentralityMin || centralityPerc>=fCentralityMax ) return false;
    hEvents -> Fill(3.5);
    
    //Primary Vertex
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    hEvents -> Fill(4.5);
    
    //Primary Vertex Selection
    if ( vertex->GetZ() < -10.0 ) return false;
    if ( vertex->GetZ() > +10.0 ) return false;
    hEvents -> Fill(5.5);
    
    //Vertex Contributors
    if ( vertex->GetNContributors() < 1 ) return false;
    hEvents -> Fill(6.5);
    
    //Event Plane
    AliEventplane *eventPlane =  fESDevent->GetEventplane();
    if (!eventPlane) return false;
    hEvents -> Fill(7.5);
    
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::PassedBasicTrackQualityCuts_Pos (AliESDtrack *track)  {
  
    fESDtrackCuts_Pos -> SetAcceptKinkDaughters(false);
    fESDtrackCuts_Pos -> SetMinNClustersTPC(50);
    fESDtrackCuts_Pos -> SetRequireTPCRefit(true);
    fESDtrackCuts_Pos -> SetMaxChi2PerClusterTPC(10);
    fESDtrackCuts_Pos -> SetEtaRange (-1.0,1.0);
    if ( !fESDtrackCuts_Pos->AcceptTrack (track) ) return false;
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::PassedBasicTrackQualityCuts_Neg (AliESDtrack *track)  {
    
    fESDtrackCuts_Neg -> SetAcceptKinkDaughters(false);
    fESDtrackCuts_Neg -> SetMinNClustersTPC(50);
    fESDtrackCuts_Neg -> SetRequireTPCRefit(true);
    fESDtrackCuts_Neg -> SetMaxChi2PerClusterTPC(10);
    fESDtrackCuts_Neg -> SetEtaRange (-1.0,1.0);
    if ( !fESDtrackCuts_Neg->AcceptTrack (track) ) return false;
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::PassedMinimalQualityCutsV0 (AliESDv0 *V0)  {
 
    //Basic Cuts
    if (V0->GetDcaV0Daughters()>2.0) return false;
    if (V0->GetRr()<3.0) return false;
    if (V0->GetV0CosineOfPointingAngle()<0.8) return false;
    
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::IsHyperTritonCandidate (AliESDv0 *V0)  {
    
    //Get V0 Daughters
    AliESDtrack *track0 = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
    AliESDtrack *track1 = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
    
    //Pair Requirements
    if ( IsPionCandidate (track0) && (!Is3HeCandidate  (track1))) return false;
    if ( Is3HeCandidate  (track0) && (!IsPionCandidate (track1))) return false;

    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::IsPionCandidate (AliESDtrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
    if (TMath::Abs(nsigmaTPC) > 4.0) return false;
    if (track->Pt()>1.5) return false;
    
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeHypertriton::Is3HeCandidate (AliESDtrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTPC) > 4.0) return false;
    
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeHypertriton::GetDecayLengthV0 (AliESDv0 *V0)  {
    
    
    //Initialization
    Double_t decayLengthV0 = 0;
    
    //Secondary Vertex Position
    Double_t secVertex[3] = { 0.0, 0.0, 0.0 };
    V0->GetXYZ(secVertex[0],secVertex[1],secVertex[2]);
    
    //Primary Vertex Position
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    Double_t primVertex[3] = { 0.0, 0.0, 0.0 };
    vertex->GetXYZ(primVertex);

    //Decay Length
    Double_t Dx = primVertex[0]-secVertex[0];
    Double_t Dy = primVertex[1]-secVertex[1];
    Double_t Dz = primVertex[2]-secVertex[2];
    decayLengthV0 = TMath::Sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
    
    return decayLengthV0;
}
//_____________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeHypertriton::GetTransverseDCA (AliESDtrack *track)  {
    
    /*
    Double_t impactParameter[2];
    track -> GetImpactParameters(impactParameter[0],impactParameter[1]);
    */
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAxy = impactParameter[0];
    
    return DCAxy;
}
//_____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeHypertriton::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_____________________________________________________________________________________________________________________________________________________

