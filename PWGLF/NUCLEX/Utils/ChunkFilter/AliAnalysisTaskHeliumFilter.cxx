#include "AliAnalysisTaskHeliumFilter.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliTimeRangeCut.h"
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
#include "TString.h"
#include "TRandom.h"
#include "AliPID.h"
#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

ClassImp(AliAnalysisTaskHeliumFilter)

//_____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHeliumFilter::AliAnalysisTaskHeliumFilter():
AliAnalysisTaskSE(),
fESDevent(NULL),
fPIDResponse(NULL),
fESDtrackCuts(NULL),
fESDeventCuts(),
fTimeRangeCut(),
fUtils(NULL),
fOutputList(NULL),
hEvents(NULL),
hdEdx_vs_p(NULL),
hdEdx_vs_p_Helium(NULL),
tree_ListOfFiles(NULL),
fEventIdFile(0),
fFileName(0)
{}
//_____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHeliumFilter::AliAnalysisTaskHeliumFilter(const char *name):
AliAnalysisTaskSE(name),
fESDevent(NULL),
fPIDResponse(NULL),
fESDtrackCuts(NULL),
fESDeventCuts(),
fTimeRangeCut(),
fUtils(NULL),
fOutputList(NULL),
hEvents(NULL),
hdEdx_vs_p(NULL),
hdEdx_vs_p_Helium(NULL),
tree_ListOfFiles(NULL),
fEventIdFile(0),
fFileName(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1,  TList::Class());
    DefineOutput(2,  TTree::Class());
    
    //Parametrization of Bethe-Block
    Double_t paramDandTdata[5] = { 6.70549, 6.11866, 8.86205e-15, 2.34059, 1.07029};
    Double_t paramHe3data[5]   = { 1.74962, 27.4992, 4.00313e-15, 2.48485, 8.31768};
    Double_t paramDandTmc[5]   = { 20.1533, 2.58127, 0.00114169,  2.03730, 0.502123};
    Double_t paramHe3mc[5]     = { 20.1533, 2.58127, 0.00114169,  2.03730, 0.502123};
    
    for(Int_t iParam=0; iParam < 5; iParam++)  {
        
       fParamDeuteron[iParam]   = paramDandTdata[iParam];
       fParamTriton[iParam]     = paramDandTdata[iParam];
       fParamHe3[iParam]        = paramHe3data[iParam];
       fParamDeuteronMC[iParam] = paramDandTmc[iParam];
       fParamTritonMC[iParam]   = paramDandTmc[iParam];
       fParamHe3MC[iParam]      = paramHe3mc[iParam];
     }
}
//_____________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHeliumFilter::~AliAnalysisTaskHeliumFilter()
{
    fOutputList->Clear();
    delete fESDevent;
    delete fPIDResponse;
    delete fESDtrackCuts;
    delete fUtils;
    delete fOutputList;
    delete hdEdx_vs_p;
    delete hdEdx_vs_p_Helium;
}
//_____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHeliumFilter::UserCreateOutputObjects()
{
    //Output List
    fOutputList = new TList();
    fOutputList -> SetOwner();
    
   
    //Histogram with number of events
    hEvents = new TH1F ("hEvents","",10,0,10);
    fOutputList -> Add(hEvents);
    
    //Tree
    tree_ListOfFiles = new TTree("tree_ListOfFiles","TreeHelium3Candidates");
    tree_ListOfFiles->Branch("fEventIdFile",&fEventIdFile,"fEventIdFile/I");
    tree_ListOfFiles->Branch("fFileName",&fFileName,16000,0);
    fOutputList -> Add (tree_ListOfFiles);
    
    //dE/dx in the TPC
    hdEdx_vs_p        = new TH2F ("hdEdx_vs_p", "", 2000, -5.0, 5.0, 1200, 0.0, 1200.0);
    hdEdx_vs_p_Helium = new TH2F ("hdEdx_vs_p_Helium", "", 2000, -5.0, 5.0, 1200, 0.0, 1200.0);
    fOutputList -> Add (hdEdx_vs_p);
    fOutputList -> Add (hdEdx_vs_p_Helium);

    //Track Cuts Object
    fESDtrackCuts = new AliESDtrackCuts ("fESDtrackCuts");
    
    PostData(1, fOutputList);
}
//_____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHeliumFilter::UserExec(Option_t *)  {
    
    //Get Input Event
    if ( !GetInputEvent ()) return;
    
    fTimeRangeCut.InitFromEvent(InputEvent());
    Bool_t cutThisEvent = fTimeRangeCut.CutEvent(InputEvent());
    if (cutThisEvent) return;
    
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse) {
        AliError("No PID Response found");
        return;
    }
    
    //Flag
    Bool_t isEventSelected = false;
    
    
    //Loop Over Reconstructed Tracks
    for ( Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++ ) {
          
        //Get Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
        if (!track) continue;
        if (!PassedBasicTrackQualityCuts(track)) continue;
        if (!track->GetInnerParam()) continue;
        
        //Variables
        Double_t dEdx_au = track->GetTPCsignal();
        Double_t p = track->GetInnerParam()->GetP();
        Int_t q = track->Charge();
        Double_t mass = AliPID::ParticleMass (AliPID::kHe3);

        //Fill dE/dx
        hdEdx_vs_p -> Fill (q*p,dEdx_au);
        
        //Selection of 3He Candidates
        if (!IsHeliumCandidate (track)) continue;
        hdEdx_vs_p_Helium -> Fill (q*p,dEdx_au);

        isEventSelected = true;
    }
    
    
    //Fill List of Files
    if (isEventSelected) {
      
        fFileName = inputHandler->GetTree()->GetCurrentFile()->GetName();
        fEventIdFile = fESDevent->GetHeader()->GetEventNumberInFile();
        tree_ListOfFiles->Fill();
    }
    
    PostData(1, fOutputList);
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHeliumFilter::GetInputEvent ()  {
    
    //Get Input Event
    fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDevent) return false;
    hEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fESDeventCuts.AcceptEvent(fESDevent)) {
        return false;
    }
    hEvents -> Fill(1.5);
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    hEvents -> Fill(2.5);
    
    //Primary Vertex
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    hEvents -> Fill(3.5);
    
    //Primary Vertex Selection
    if ( vertex->GetZ() < -10.0 ) return false;
    if ( vertex->GetZ() > +10.0 ) return false;
    hEvents -> Fill(4.5);
    
    //Vertex Contributors
    if ( vertex->GetNContributors() < 1 ) return false;
    hEvents -> Fill(5.5);
    
    //Event Plane
    AliEventplane *eventPlane =  fESDevent->GetEventplane();
    if (!eventPlane) return false;
    hEvents -> Fill(6.5);
    
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHeliumFilter::PassedBasicTrackQualityCuts (AliESDtrack *track)  {
   
    fESDtrackCuts -> SetAcceptKinkDaughters(false);
    fESDtrackCuts -> SetMinNClustersTPC(50);
    fESDtrackCuts -> SetRequireTPCRefit(true);
    fESDtrackCuts -> SetMaxChi2PerClusterTPC(10.0);
    fESDtrackCuts -> SetEtaRange (-0.9,0.9);
    if ( track -> GetTPCsignalN() < 40 ) return false;
    if ( !fESDtrackCuts->AcceptTrack (track) ) return false;
    return true;
}
//_____________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHeliumFilter::IsHeliumCandidate (AliESDtrack *track)  {
    
    //Initialization
    Bool_t IsHeliumCandidate = false;
    
    //Variables
    Double_t p = track->GetInnerParam()->GetP();
    Double_t mass = AliPID::ParticleMass (AliPID::kHe3);
    Double_t dEdx_au = track->GetTPCsignal();

    //Expected dE/dx for 3He
    Float_t hel3Exp = 4.0*AliExternalTrackParam::BetheBlochAleph(2.0*p/mass,fParamHe3[0],fParamHe3[1],fParamHe3[2],fParamHe3[3],fParamHe3[4]);
    Double_t sigma = 0.07;//dE/dx Resolution for 3He (7%)
    Double_t nSigmaHe3  = (dEdx_au - hel3Exp)/(sigma*hel3Exp);
    
    if (nSigmaHe3>-4.0) IsHeliumCandidate = true;
    
    return IsHeliumCandidate;
}
//_____________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHeliumFilter::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_____________________________________________________________________________________________________________________________________________________

