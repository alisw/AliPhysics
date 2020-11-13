#include "AliAnalysisTaskNdimFitsTPCsplines.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliTPCPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliTimeRangeCut.h"
#include "TLorentzVector.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliEventCuts.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "TObjArray.h"
#include "THnSparse.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliESDv0.h"
#include "TString.h"
#include "TRandom.h"
#include "AliPID.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskNdimFitsTPCsplines)

//______________________________________________________________________________________________________________________________________
AliAnalysisTaskNdimFitsTPCsplines::AliAnalysisTaskNdimFitsTPCsplines():
AliAnalysisTaskSE(),
fESDevent(nullptr),
fPIDResponse(nullptr),
fESDtrackCuts_V0daugh(nullptr),
fESDtrackCuts_Primary(nullptr),
fESDtrackCuts_Nuclei(nullptr),
fESDeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fTimeRangeCut(),
hNumberOfEvents(nullptr),
hTPCdEdx_Electrons(nullptr),
hTPCdEdx_Pions(nullptr),
hTPCdEdx_Kaons(nullptr),
hTPCdEdx_Protons(nullptr),
hTPCdEdx_Tritons(nullptr),
hTPCdEdx_Helium3(nullptr)
{}
//______________________________________________________________________________________________________________________________________
AliAnalysisTaskNdimFitsTPCsplines::AliAnalysisTaskNdimFitsTPCsplines(const char *name):
AliAnalysisTaskSE(name),
fESDevent(nullptr),
fPIDResponse(nullptr),
fESDtrackCuts_V0daugh(nullptr),
fESDtrackCuts_Primary(nullptr),
fESDtrackCuts_Nuclei(nullptr),
fESDeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fTimeRangeCut(),
hNumberOfEvents(nullptr),
hTPCdEdx_Electrons(nullptr),
hTPCdEdx_Pions(nullptr),
hTPCdEdx_Kaons(nullptr),
hTPCdEdx_Protons(nullptr),
hTPCdEdx_Tritons(nullptr),
hTPCdEdx_Helium3(nullptr)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//______________________________________________________________________________________________________________________________________
AliAnalysisTaskNdimFitsTPCsplines::~AliAnalysisTaskNdimFitsTPCsplines()  {
    
    fOutputList->Clear();
    delete fESDevent;
    delete fPIDResponse;
    delete fESDtrackCuts_V0daugh;
    delete fESDtrackCuts_Primary;
    delete fESDtrackCuts_Nuclei;
    delete fUtils;
    delete fOutputList;
    delete fQAList;
    delete hNumberOfEvents;
    delete hTPCdEdx_Electrons;
    delete hTPCdEdx_Pions;
    delete hTPCdEdx_Kaons;
    delete hTPCdEdx_Protons;
    delete hTPCdEdx_Tritons;
    delete hTPCdEdx_Helium3;

}
//______________________________________________________________________________________________________________________________________
void AliAnalysisTaskNdimFitsTPCsplines::UserCreateOutputObjects()  {
    
    //Output Lists
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();

    //QA Plots
    fESDeventCuts.AddQAplotsToList(fQAList);
    
    //Histogram with Number of Events
    hNumberOfEvents = new TH1F ("hNumberOfEvents","",10,0,10);
    fOutputList -> Add(hNumberOfEvents);
    
    
    //THnSparse Binning: Pions, Kaons, Protons, Electrons
    const Int_t nDim = 4;
    //*********************  dE/dx,   p,  eta, centr
    Int_t    bins[nDim] = {    500,  100,   18,   18 };
    Double_t xmin[nDim] = {    0.0,  0.1, -0.9,  0.0 };
    Double_t xmax[nDim] = { 2000.0, 10.0,  0.9, 90.0 };
    
    //THnSparse Binning: Tritons & Helium3
    //***************************  dE/dx,    p,  eta, centr
    Int_t    binsNuclei[nDim] = {    500,   19,    4,    9 };
    Double_t xminNuclei[nDim] = {    0.0,  0.5, -0.9,  0.0 };
    Double_t xmaxNuclei[nDim] = { 2000.0, 10.0,  0.9, 90.0 };

    
    //nDimensional Histograms: Raw TPC dE/dx
    hTPCdEdx_Electrons = new THnSparseF ("hTPCdEdx_Electrons","",nDim, bins, xmin, xmax);
    hTPCdEdx_Pions     = new THnSparseF ("hTPCdEdx_Pions","    ",nDim, bins, xmin, xmax);
    hTPCdEdx_Kaons     = new THnSparseF ("hTPCdEdx_Kaons","    ",nDim, bins, xmin, xmax);
    hTPCdEdx_Protons   = new THnSparseF ("hTPCdEdx_Protons","  ",nDim, bins, xmin, xmax);
    hTPCdEdx_Tritons   = new THnSparseF ("hTPCdEdx_Tritons","  ",nDim, binsNuclei, xminNuclei, xmaxNuclei);
    hTPCdEdx_Helium3   = new THnSparseF ("hTPCdEdx_Helium3","  ",nDim, binsNuclei, xminNuclei, xmaxNuclei);

    BinLogAxis (hTPCdEdx_Electrons,1);
    BinLogAxis (hTPCdEdx_Pions,1);
    BinLogAxis (hTPCdEdx_Kaons,1);
    BinLogAxis (hTPCdEdx_Protons,1);
    //BinLogAxis (hTPCdEdx_Tritons,1);
    //BinLogAxis (hTPCdEdx_Helium3,1);

    hTPCdEdx_Electrons -> Sumw2();
    hTPCdEdx_Pions     -> Sumw2();
    hTPCdEdx_Kaons     -> Sumw2();
    hTPCdEdx_Protons   -> Sumw2();
    hTPCdEdx_Tritons   -> Sumw2();
    hTPCdEdx_Helium3   -> Sumw2();

    fOutputList -> Add (hTPCdEdx_Electrons);
    fOutputList -> Add (hTPCdEdx_Pions);
    fOutputList -> Add (hTPCdEdx_Kaons);
    fOutputList -> Add (hTPCdEdx_Protons);
    fOutputList -> Add (hTPCdEdx_Tritons);
    fOutputList -> Add (hTPCdEdx_Helium3);

    
    //Track Cuts Objects
    fESDtrackCuts_V0daugh = new AliESDtrackCuts ("fESDtrackCuts_V0daugh");
    fESDtrackCuts_Primary = new AliESDtrackCuts ("fESDtrackCuts_Primary");
    fESDtrackCuts_Nuclei  = new AliESDtrackCuts ("fESDtrackCuts_Nuclei");

    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//______________________________________________________________________________________________________________________________________
void AliAnalysisTaskNdimFitsTPCsplines::UserExec(Option_t *)  {
    
    
    //Get Input Event
    if ( !GetInputEvent ()) return;
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    //Get Centrality
    AliMultSelection *multiplicitySel = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    Double_t centrality = multiplicitySel->GetMultiplicityPercentile("V0M");
    
    //ID Candidates
    vector<Int_t> pion_ID;
    vector<Int_t> prot_ID;
    vector<Int_t> elec_ID;
    vector<Int_t> primaryElec_ID;
    vector<Int_t> trit_ID;
    vector<Int_t> he3_ID;


    //Loop Over Reconstructed V0s
    for (Int_t i=0 ; i<fESDevent->GetNumberOfV0s() ; i++)  {
        
        //Get V0 Candidate
        AliESDv0 *V0 = (AliESDv0*)fESDevent->GetV0(i);
        if (!V0) continue;
        if (V0->GetOnFlyStatus()) continue;//Offline V0s Selection

        //Get V0 Daughters
        AliESDtrack *posTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
        AliESDtrack *negTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
        if (!posTrack) continue;
        if (!negTrack) continue;
        if (posTrack->Charge() == negTrack->Charge()) continue;
        if (posTrack->GetID()  == negTrack->GetID() ) continue;

        //Track Quality Cuts
        if (!PassedTrackSelectionV0daugh (posTrack)) continue;
        if (!PassedTrackSelectionV0daugh (negTrack)) continue;
        if (!PassedV0Selection(V0)) continue;
        
        //Store Candidate IDs
        if (PassedLambdaSelection(V0))          { prot_ID.push_back(V0->GetPindex()); pion_ID.push_back(V0->GetNindex()); }
        if (PassedAntiLambdaSelection(V0))      { prot_ID.push_back(V0->GetNindex()); pion_ID.push_back(V0->GetPindex()); }
        if (PassedK0shortSelection(V0))         { pion_ID.push_back(V0->GetPindex()); pion_ID.push_back(V0->GetNindex()); }
        if (PassedGammaConversionSelection(V0)) { elec_ID.push_back(V0->GetPindex()); elec_ID.push_back(V0->GetNindex()); }
    }
    

    //Loop Over Reconstructed Tracks
    for (Int_t i=0; i<fESDevent->GetNumberOfTracks(); i++)  {
        
        //Get ESD Track
        AliESDtrack *track = static_cast<AliESDtrack*>(fESDevent->GetTrack(i));
        if (!track) continue;
        
        //Store Nuclei Candidates
        if (PassedNucleiTrackSelection(track))  {
            
            Double_t nsigmaTPC_triton  = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
            Double_t nsigmaTPC_helium3 = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);

            if  (TMath::Abs(nsigmaTPC_triton)<5.0) trit_ID.push_back(i);
            if  (nsigmaTPC_helium3 > -4.0)         he3_ID.push_back(i);
        }
        
        //Primary Track Selection
        if (!PassedPrimaryTrackSelection (track)) continue;
        
        //TOF Requirements
        if(!track->IsOn(AliESDtrack::kTOFout)) continue;
        if(!track->IsOn(AliESDtrack::kTIME))   continue;
               
        //Kaon Selection
        Double_t nsigmaTOF_kaon = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kKaon);
        Double_t nsigmaITS_kaon = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kKaon);
        if (TMath::Abs(nsigmaTOF_kaon) < 2.0 && TMath::Abs(nsigmaITS_kaon) < 3.0)  {
            
            //Variables
            Double_t dEdx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track);
            Double_t p    = track->GetInnerParam()->GetP();
            Double_t eta  = track->Eta();
                              
            //Fill Histogram
            Double_t x[4] = {dEdx,p,eta,centrality};
            hTPCdEdx_Kaons -> Fill (x);
        }

        //Electron Selection
        Double_t nsigmaTOF_elec = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kElectron);
        if (TMath::Abs(nsigmaTOF_elec) < 3.0) primaryElec_ID.push_back(i);
    }
    
    
    
    
    //Electron Selection
    Int_t nPrimElec = (Int_t)primaryElec_ID.size();
    for (Int_t i=0 ; i<nPrimElec ; i++)  {
        for (Int_t j=i+1 ; j<nPrimElec ; j++)  {
            
            AliESDtrack *track1 = static_cast<AliESDtrack*>(fESDevent->GetTrack(primaryElec_ID[i]));
            AliESDtrack *track2 = static_cast<AliESDtrack*>(fESDevent->GetTrack(primaryElec_ID[j]));
            if (track1->Charge() == track2->Charge()) continue;
            if (track1->GetID()  == track2->GetID())  continue;

            TVector3 P1 (track1->Px(),track1->Py(),track1->Pz());
            TVector3 P2 (track2->Px(),track2->Py(),track2->Pz());
            Double_t mass = MassDielectron (P1,P2);
            if  (mass>0.1) continue;
            
            elec_ID.push_back (primaryElec_ID[i]);
            elec_ID.push_back (primaryElec_ID[j]);
        }
    }
    
    
    
    
    //Sort Vectors
    sort( pion_ID.begin(), pion_ID.end() );
    sort( prot_ID.begin(), prot_ID.end() );
    sort( elec_ID.begin(), elec_ID.end() );

    //Remove Replicas
    pion_ID.erase( unique( pion_ID.begin(), pion_ID.end() ), pion_ID.end() );
    prot_ID.erase( unique( prot_ID.begin(), prot_ID.end() ), prot_ID.end() );
    elec_ID.erase( unique( elec_ID.begin(), elec_ID.end() ), elec_ID.end() );

    
    //Fill Histogram for Pions
    Int_t nPions = (Int_t)pion_ID.size();
    for (Int_t i=0 ; i<nPions ; i++)  {
        
        //Get Track
        AliESDtrack *track = static_cast<AliESDtrack*>(fESDevent->GetTrack(pion_ID[i]));

        //Variables
        Double_t dEdx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track);
        Double_t p    = track->GetInnerParam()->GetP();
        Double_t eta  = track->Eta();
        
        //Fill Histogram
        Double_t x[4] = {dEdx,p,eta,centrality};
        hTPCdEdx_Pions -> Fill (x);
    }
    
    
    //Fill Histogram for Protons
    Int_t nProtons = (Int_t)prot_ID.size();
    for (Int_t i=0 ; i<nProtons ; i++)  {
        
        //Get Track
        AliESDtrack *track = static_cast<AliESDtrack*>(fESDevent->GetTrack(prot_ID[i]));

        //Variables
        Double_t dEdx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track);
        Double_t p    = track->GetInnerParam()->GetP();
        Double_t eta  = track->Eta();
        
        //Fill Histogram
        Double_t x[4] = {dEdx,p,eta,centrality};
        hTPCdEdx_Protons -> Fill (x);
    }
    
    
    //Fill Histogram for Electrons
    Int_t nElectrons = (Int_t)elec_ID.size();
    for (Int_t i=0 ; i<nElectrons ; i++)  {
        
        //Get Track
        AliESDtrack *track = static_cast<AliESDtrack*>(fESDevent->GetTrack(elec_ID[i]));

        //Variables
        Double_t dEdx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track);
        Double_t p    = track->GetInnerParam()->GetP();
        Double_t eta  = track->Eta();
        
        //Fill Histogram
        Double_t x[4] = {dEdx,p,eta,centrality};
        hTPCdEdx_Electrons -> Fill (x);
    }
    
    
    
    //Fill Histogram for Tritons
    Int_t nTritons = (Int_t)trit_ID.size();
    for (Int_t i=0 ; i<nTritons ; i++)  {
           
        //Get Track
        AliESDtrack *track = static_cast<AliESDtrack*>(fESDevent->GetTrack(trit_ID[i]));

        //Variables
        Double_t dEdx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track);
        Double_t p    = track->GetInnerParam()->GetP();
        Double_t eta  = track->Eta();
           
        //Fill Histogram
        Double_t x[4] = {dEdx,p,eta,centrality};
        hTPCdEdx_Tritons -> Fill (x);
    }
    
    
    //Fill Histogram for Helium3
    Int_t nHelium3 = (Int_t)he3_ID.size();
    for (Int_t i=0 ; i<nHelium3 ; i++)  {
        
        //Get Track
        AliESDtrack *track = static_cast<AliESDtrack*>(fESDevent->GetTrack(he3_ID[i]));
           
        //Variables
        Double_t dEdx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track);
        Double_t p    = 2.0*track->GetInnerParam()->GetP();
        Double_t eta  = track->Eta();
        
        //Fill Histogram
        Double_t x[4] = {dEdx,p,eta,centrality};
        hTPCdEdx_Helium3 -> Fill (x);
    }
       
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::GetInputEvent ()  {
    
    //Get Input Event
    fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDevent) return false;
    hNumberOfEvents -> Fill(0.5);
     
    //Primary Vertex
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    if (!vertex) return false;
    hNumberOfEvents -> Fill(1.5);
     
    //Standard Event Cuts
    if (!fESDeventCuts.AcceptEvent(fESDevent)) { PostData(2, fQAList); return false; }
    hNumberOfEvents -> Fill(2.5);
     
    //Primary Vertex Selection
    if ( vertex->GetZ() < -10.0 ) return false;
    if ( vertex->GetZ() > +10.0 ) return false;
    hNumberOfEvents -> Fill(3.5);
        
    //Vertex Contributors
    if ( vertex->GetNContributors() < 2 ) return false;
    hNumberOfEvents -> Fill(4.5);

    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    hNumberOfEvents -> Fill(5.5);
     
    //Selection of Centrality Range
    Double_t centrality = multiplicitySelection->GetMultiplicityPercentile("V0M");
    if (centrality <   0.0) return false;
    if (centrality >= 90.0) return false;
    hNumberOfEvents -> Fill(6.5);

    //Time-Range Selection (for LHC18r)
    fTimeRangeCut.InitFromEvent(InputEvent());
    const Bool_t cutThisEvent = fTimeRangeCut.CutEvent(InputEvent());
    if (cutThisEvent) return false;
    hNumberOfEvents -> Fill(7.5);

    return true;
}
//______________________________________________________________________________________________________________________________________
void AliAnalysisTaskNdimFitsTPCsplines::BinLogAxis (THnSparseF *h, Int_t axisNumber)  {
    
    TAxis *axis = h->GetAxis(axisNumber);
    Int_t bins = axis->GetNbins();

    Double_t min = axis->GetXmin();
    Double_t max = axis->GetXmax();
    Double_t *newBins = new Double_t[bins + 1];
   
    newBins[0] = min;
    Double_t factor = TMath::Power (max/min, 1.0/(Double_t)bins);
  
    for (Int_t i=1; i <= bins; i++) {
        newBins[i] = factor * newBins[i-1];
    }
    axis->Set(bins, newBins);
    delete [] newBins;
    
    return;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedTrackSelectionV0daugh (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
       
    //Track Selection
    fESDtrackCuts_V0daugh -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts_V0daugh -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts_V0daugh -> SetMinNCrossedRowsTPC(80);
    fESDtrackCuts_V0daugh -> SetMinNClustersTPC(70);
    fESDtrackCuts_V0daugh -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fESDtrackCuts_V0daugh -> SetMaxChi2PerClusterTPC(4.0);
    fESDtrackCuts_V0daugh -> SetEtaRange (-0.9,0.9);
    if ( track->GetTPCsignalN() < 70 )                 return passedTrkSelection;
    if ( !fESDtrackCuts_V0daugh->AcceptTrack (track) ) return passedTrkSelection;

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedPrimaryTrackSelection (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
       
    //Track Selection Cuts
    fESDtrackCuts_Primary -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts_Primary -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts_Primary -> SetRequireITSRefit(kTRUE);
    fESDtrackCuts_Primary -> SetMinNClustersTPC(70);
    fESDtrackCuts_Primary -> SetMinNClustersITS(4);
    fESDtrackCuts_Primary -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    fESDtrackCuts_Primary -> SetMinNCrossedRowsTPC(80);
    fESDtrackCuts_Primary -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    fESDtrackCuts_Primary -> SetMaxChi2PerClusterTPC(4.0);
    fESDtrackCuts_Primary -> SetEtaRange (-0.9,0.9);
    fESDtrackCuts_Primary -> SetMaxDCAToVertexXY(0.1);
    fESDtrackCuts_Primary -> SetMaxDCAToVertexZ(0.1);
    fESDtrackCuts_Primary -> SetDCAToVertex2D(kTRUE);
    if ( track->GetTPCsignalN() < 70 )                 return passedTrkSelection;
    if ( !fESDtrackCuts_Primary->AcceptTrack (track) ) return passedTrkSelection;

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedNucleiTrackSelection (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
       
    //Track Selection Cuts
    fESDtrackCuts_Nuclei -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts_Nuclei -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts_Nuclei -> SetMinNClustersTPC(70);
    fESDtrackCuts_Nuclei -> SetMinNCrossedRowsTPC(80);
    fESDtrackCuts_Nuclei -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
    fESDtrackCuts_Nuclei -> SetMaxChi2PerClusterTPC(5.0);
    fESDtrackCuts_Nuclei -> SetEtaRange (-0.9,0.9);
    if ( track->GetTPCsignalN() < 50 )                return passedTrkSelection;
    if ( !fESDtrackCuts_Nuclei->AcceptTrack (track) ) return passedTrkSelection;

    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedV0Selection (AliESDv0 *V0)  {
    
    //Initialization
    Bool_t passedV0Selection=(kFALSE);
    
    //Primary Vertex
    AliESDVertex *primaryVertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    Double_t vx = primaryVertex->GetX();
    Double_t vy = primaryVertex->GetY();
    Double_t vz = primaryVertex->GetZ();
    
    //Pair Cuts
    if (V0->GetD(vx,vy,vz) > 1.0 )               return passedV0Selection;
    if (V0->GetDcaV0Daughters() > 0.5)           return passedV0Selection;
    if (V0->GetV0CosineOfPointingAngle() < 0.99) return passedV0Selection;
    if (V0->GetChi2V0() > 50.0)                  return passedV0Selection;
    
    passedV0Selection=kTRUE;
    return passedV0Selection;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedLambdaSelection (AliESDv0 *V0)  {
    
    //Initialization
    Bool_t passedLambdaSelection=(kFALSE);
    
    //Get V0 Daughters
    AliESDtrack *posTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
    AliESDtrack *negTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
    
    //Momenta of the Daughters
    Double_t posMomentum[3] = { 0.0, 0.0, 0.0 };
    Double_t negMomentum[3] = { 0.0, 0.0, 0.0 };
    V0->GetPPxPyPz(posMomentum[0],posMomentum[1],posMomentum[2]);
    V0->GetNPxPyPz(negMomentum[0],negMomentum[1],negMomentum[2]);
    TVector3 Ppos (posMomentum[0],posMomentum[1],posMomentum[2]);
    TVector3 Pneg (negMomentum[0],negMomentum[1],negMomentum[2]);
    
    //Selection on Armenteros Plot
    if ( V0->AlphaV0()<-0.9 || V0->AlphaV0()>0.9 ) return passedLambdaSelection;
    if ( V0->AlphaV0()>-0.4 && V0->AlphaV0()<0.4 ) return passedLambdaSelection;
    if ( V0->PtArmV0() > 0.12) return passedLambdaSelection;
    
    //Invariant-Mass Selection: Lambda^{0}
    Double_t mass = MassLambda(Pneg,Ppos);//Pion,Proton
    if (mass<1.105) return passedLambdaSelection;
    if (mass>1.125) return passedLambdaSelection;
    
    
    passedLambdaSelection=kTRUE;
    return passedLambdaSelection;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedAntiLambdaSelection (AliESDv0 *V0)  {
    
    //Initialization
    Bool_t passedLambdaSelection=(kFALSE);
    
    //Get V0 Daughters
    AliESDtrack *posTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
    AliESDtrack *negTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
    
    //Momenta of the Daughters
    Double_t posMomentum[3] = { 0.0, 0.0, 0.0 };
    Double_t negMomentum[3] = { 0.0, 0.0, 0.0 };
    V0->GetPPxPyPz(posMomentum[0],posMomentum[1],posMomentum[2]);
    V0->GetNPxPyPz(negMomentum[0],negMomentum[1],negMomentum[2]);
    TVector3 Ppos (posMomentum[0],posMomentum[1],posMomentum[2]);
    TVector3 Pneg (negMomentum[0],negMomentum[1],negMomentum[2]);
    
    //Selection on Armenteros Plot
    if ( V0->AlphaV0()<-0.9 || V0->AlphaV0()>0.9 ) return passedLambdaSelection;
    if ( V0->AlphaV0()>-0.4 && V0->AlphaV0()<0.4 ) return passedLambdaSelection;
    if ( V0->PtArmV0() > 0.12) return passedLambdaSelection;
    
    //Invariant-Mass Selection: Anti-Lambda^{0}
    Double_t mass = MassLambda(Ppos,Pneg);//Pion,Proton
    if (mass<1.105) return passedLambdaSelection;
    if (mass>1.125) return passedLambdaSelection;
    
    
    passedLambdaSelection=kTRUE;
    return passedLambdaSelection;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedK0shortSelection (AliESDv0 *V0)  {
    
    //Initialization
    Bool_t passedK0shortSelection=(kFALSE);
    
    //Get V0 Daughters
    AliESDtrack *posTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
    AliESDtrack *negTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
    
    //Momenta of the Daughters
    Double_t posMomentum[3] = { 0.0, 0.0, 0.0 };
    Double_t negMomentum[3] = { 0.0, 0.0, 0.0 };
    V0->GetPPxPyPz(posMomentum[0],posMomentum[1],posMomentum[2]);
    V0->GetNPxPyPz(negMomentum[0],negMomentum[1],negMomentum[2]);
    TVector3 Ppos (posMomentum[0],posMomentum[1],posMomentum[2]);
    TVector3 Pneg (negMomentum[0],negMomentum[1],negMomentum[2]);
    
    //Selection on Armenteros Plot
    if ( V0->PtArmV0() < 0.10) return passedK0shortSelection;
    if ( V0->PtArmV0() > 0.25) return passedK0shortSelection;
    
    //Invariant-Mass Selection: K^{0}
    Double_t mass = MassK0short(Pneg,Ppos);//Pion,Proton
    if (mass<0.45) return passedK0shortSelection;
    if (mass>0.55) return passedK0shortSelection;
    
    
    passedK0shortSelection=kTRUE;
    return passedK0shortSelection;
}
//______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskNdimFitsTPCsplines::PassedGammaConversionSelection (AliESDv0 *V0)  {
    
    //Initialization
    Bool_t passedGammaSelection=(kFALSE);
    
    //Get V0 Daughters
    AliESDtrack *posTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetPindex());
    AliESDtrack *negTrack = (AliESDtrack*) fESDevent->GetTrack(V0->GetNindex());
    
    //Momenta of the Daughters
    Double_t posMomentum[3] = { 0.0, 0.0, 0.0 };
    Double_t negMomentum[3] = { 0.0, 0.0, 0.0 };
    V0->GetPPxPyPz(posMomentum[0],posMomentum[1],posMomentum[2]);
    V0->GetNPxPyPz(negMomentum[0],negMomentum[1],negMomentum[2]);
    TVector3 Ppos (posMomentum[0],posMomentum[1],posMomentum[2]);
    TVector3 Pneg (negMomentum[0],negMomentum[1],negMomentum[2]);
    
    //Selection on Armenteros Plot
    if ( V0->PtArmV0() > 0.02) return passedGammaSelection;
    
    //Invariant-Mass & PhiV Selection
    Double_t mass = MassDielectron(Ppos,Pneg);
    Double_t phiV = GetPhiV (posTrack,negTrack);
    
    if (mass>0.1)                return passedGammaSelection;
    if (phiV>60.0 && phiV<120.0) return passedGammaSelection;
    
    passedGammaSelection=kTRUE;
    return passedGammaSelection;
}
//______________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskNdimFitsTPCsplines::GetPhiV (AliESDtrack *track1, AliESDtrack *track2 ) {
    
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
//______________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskNdimFitsTPCsplines::MassLambda (TVector3 Ppion, TVector3 Pprot)  {
    
    //Initialization
    Double_t mass(0);
    
    //Particle Masses
    Double_t mPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    Double_t mProt = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    
    //4-Momentum Vectors
    TLorentzVector P1; P1.SetXYZM(Ppion.Px(),Ppion.Py(),Ppion.Pz(),mPion);
    TLorentzVector P2; P2.SetXYZM(Pprot.Px(),Pprot.Py(),Pprot.Pz(),mProt);
    
    //Invariant Mass
    mass = (P1 + P2).M();
    
    return mass;
}
//______________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskNdimFitsTPCsplines::MassDielectron (TVector3 Pelec1, TVector3 Pelec2)  {
    
    //Initialization
    Double_t mass(0);
    
    //Particle Mass
    Double_t mElec = TDatabasePDG::Instance()->GetParticle(11)->Mass();
    
    //4-Momentum Vectors
    TLorentzVector P1; P1.SetXYZM(Pelec1.Px(),Pelec1.Py(),Pelec1.Pz(),mElec);
    TLorentzVector P2; P2.SetXYZM(Pelec2.Px(),Pelec2.Py(),Pelec2.Pz(),mElec);

    //Invariant Mass
    mass = (P1 + P2).M();
    
    return mass;
}
//______________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskNdimFitsTPCsplines::MassK0short (TVector3 Ppion1, TVector3 Ppion2)  {
    
    //Initialization
    Double_t mass(0);
    
    //Particle Mass
    Double_t mPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
    
    //4-Momentum Vectors
    TLorentzVector P1; P1.SetXYZM(Ppion1.Px(),Ppion1.Py(),Ppion1.Pz(),mPion);
    TLorentzVector P2; P2.SetXYZM(Ppion2.Px(),Ppion2.Py(),Ppion2.Pz(),mPion);

    //Invariant Mass
    mass = (P1 + P2).M();
    
    return mass;
}
//______________________________________________________________________________________________________________________________________
void AliAnalysisTaskNdimFitsTPCsplines::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//______________________________________________________________________________________________________________________________________

