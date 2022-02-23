#include "AliAnalysisTaskTritonESD_PbPb.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "TDatabasePDG.h"
#include "AliEventCuts.h"
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliAnalysisTaskTritonESD_PbPb)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonESD_PbPb::AliAnalysisTaskTritonESD_PbPb():
AliAnalysisTaskSE(),
fESDevent(NULL),
fPIDResponse(NULL),
fESDeventCuts(),
fUtils(NULL),
fOutputList(NULL),
//fQAList(NULL),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(NULL),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fYMax(0),
fNumberClustersITSMin(0),
fNumberClustersTPCMin(0),
fNumberCrossedRowsTPCMin(0),
fCrossedRowsFindableClsMin(0),
fNumberClustersTPCdEdxMin(0),
fChiSquarePerNDFMax(0),
fITSrequirement(NULL),
fDCAzMax(0),
fDCAxyMax(0),
fnSigmaTOFmax(0),
fnSigmaTPCmax(0),
fTRDntracklets(0),
fpar0_mean_TPC(0),
fpar1_mean_TPC(0),
fpar0_sigma_TPC(0),
fpar0_mean_TOF(0),
fpar1_mean_TOF(0),
fpar0_sigma_TOF(0),
fpar1_sigma_TOF(0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonESD_PbPb::AliAnalysisTaskTritonESD_PbPb(const char *name):
AliAnalysisTaskSE(name),
fESDevent(NULL),
fPIDResponse(NULL),
fESDeventCuts(),
fUtils(NULL),
fOutputList(NULL),
//fQAList(NULL),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(NULL),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fYMax(0),
fNumberClustersITSMin(0),
fNumberClustersTPCMin(0),
fNumberCrossedRowsTPCMin(0),
fCrossedRowsFindableClsMin(0),
fNumberClustersTPCdEdxMin(0),
fChiSquarePerNDFMax(0),
fITSrequirement(NULL),
fDCAzMax(0),
fDCAxyMax(0),
fnSigmaTOFmax(0),
fnSigmaTPCmax(0),
fTRDntracklets(0),
fpar0_mean_TPC(0),
fpar1_mean_TPC(0),
fpar0_sigma_TPC(0),
fpar0_mean_TOF(0),
fpar1_mean_TOF(0),
fpar0_sigma_TOF(0),
fpar1_sigma_TOF(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(2, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonESD_PbPb::~AliAnalysisTaskTritonESD_PbPb()
{
    fOutputList->Clear();
    delete fESDevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    //delete fQAList;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonESD_PbPb::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList -> SetOwner();
   
    //fQAList = new TList();
    //fQAList -> SetOwner();
    
    //fESDeventCuts.AddQAplotsToList(fQAList);//Add event selection QA plots

    
    //Number of Events
    histoNumberOfEvents = new TH1F("histoNumberOfEvents","Events after selection steps",10,0,10);
    fOutputList -> Add (histoNumberOfEvents);
    
    
    //Signal Extraction
    histoNsigmaTPCtriton_vs_pt     = new TH2F ("histoNsigmaTPCtriton_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFtriton_vs_pt     = new TH2F ("histoNsigmaTOFtriton_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_pt = new TH2F ("histoNsigmaTPCantitriton_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_pt = new TH2F ("histoNsigmaTOFantitriton_vs_pt","",500,0,5,1000,-20,20);
 
    histoNsigmaTPCtriton_vs_pt_centered     = new TH2F ("histoNsigmaTPCtriton_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_pt_centered = new TH2F ("histoNsigmaTPCantitriton_vs_pt_centered","",500,0,5,1000,-20,20);
    
    histoNsigmaTOFtriton_vs_pt_centered     = new TH2F ("histoNsigmaTOFtriton_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_pt_centered = new TH2F ("histoNsigmaTOFantitriton_vs_pt_centered","",500,0,5,1000,-20,20);
    
    histoNsigmaTOFtriton_vs_pt_trd     = new TH2F ("histoNsigmaTOFtriton_vs_pt_trd","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_pt_trd = new TH2F ("histoNsigmaTOFantitriton_vs_pt_trd","",500,0,5,1000,-20,20);
        
    histoNsigmaTPCtriton_vs_p         = new TH2F ("histoNsigmaTPCtriton_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_p     = new TH2F ("histoNsigmaTPCantitriton_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFtriton_vs_p         = new TH2F ("histoNsigmaTOFtriton_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_p     = new TH2F ("histoNsigmaTOFantitriton_vs_p","",500,0,5,1000,-20,20);
    
    histoNsigmaTPCtriton_vs_p_notof     = new TH2F ("histoNsigmaTPCtriton_vs_p_notof","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_p_notof = new TH2F ("histoNsigmaTPCantitriton_vs_p_notof","",500,0,5,1000,-20,20);

    histoNsigmaTPCtriton_vs_pt     -> Sumw2();
    histoNsigmaTOFtriton_vs_pt     -> Sumw2();
    histoNsigmaTPCantitriton_vs_pt -> Sumw2();
    histoNsigmaTOFantitriton_vs_pt -> Sumw2();
  
    histoNsigmaTPCtriton_vs_pt_centered           -> Sumw2();
    histoNsigmaTPCantitriton_vs_pt_centered       -> Sumw2();
    
    histoNsigmaTOFtriton_vs_pt_centered           -> Sumw2();
    histoNsigmaTOFantitriton_vs_pt_centered       -> Sumw2();
    
    histoNsigmaTOFtriton_vs_pt_trd                -> Sumw2();
    histoNsigmaTOFantitriton_vs_pt_trd            -> Sumw2();
    
    histoNsigmaTPCtriton_vs_p                     -> Sumw2();
    histoNsigmaTPCantitriton_vs_p                 -> Sumw2();
    histoNsigmaTOFtriton_vs_p                     -> Sumw2();
    histoNsigmaTOFantitriton_vs_p                 -> Sumw2();
    
    histoNsigmaTPCtriton_vs_p_notof              -> Sumw2();
	histoNsigmaTPCantitriton_vs_p_notof          -> Sumw2();
    
    fOutputList -> Add(histoNsigmaTPCtriton_vs_pt);
    fOutputList -> Add(histoNsigmaTOFtriton_vs_pt);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_pt);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_pt);
  
    fOutputList -> Add(histoNsigmaTPCtriton_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTOFtriton_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_pt_centered);
    
    fOutputList -> Add(histoNsigmaTOFtriton_vs_pt_trd);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_pt_trd);
        
    fOutputList -> Add(histoNsigmaTPCtriton_vs_p);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_p);
    fOutputList -> Add(histoNsigmaTOFtriton_vs_p);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_p);
    
    fOutputList -> Add(histoNsigmaTPCtriton_vs_p_notof);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_p_notof);

    
    
    
    //DCA Distributions
    histoDCAxyTriton_vs_pt     = new TH2F ("histoDCAxyTriton_vs_pt","",500,0,5,500,-5,5);
    histoDCAxyAntiTriton_vs_pt = new TH2F ("histoDCAxyAntiTriton_vs_pt","",500,0,5,500,-5,5);

    histoDCAxyTriton_vs_pt     -> Sumw2();
    histoDCAxyAntiTriton_vs_pt -> Sumw2();
    
    fOutputList -> Add (histoDCAxyTriton_vs_pt);
    fOutputList -> Add (histoDCAxyAntiTriton_vs_pt);

    
    
    
    PostData(1, fOutputList);
    //PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonESD_PbPb::UserExec(Option_t *)
{
    
    //Get Input Event
    if ( !GetInputEvent ()) return;
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {
        
        //Track Selection
        AliESDtrack *track = (AliESDtrack*) fESDevent -> GetTrack(i);
        if ( !track ) continue;
        if ( PassedTrackQualityCutsNoDCA (track)) {

        if (IsCleanTritonCandidate(track))  {
            if (track->Charge()>0) histoDCAxyTriton_vs_pt     -> Fill (track->Pt(),GetDCAxy(track));
            if (track->Charge()<0) histoDCAxyAntiTriton_vs_pt -> Fill (track->Pt(),GetDCAxy(track));
        }
        
        }

        if ( !PassedTrackQualityCuts (track)) continue;

        
        //Variables
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);

        
        //TPC Signal vs. pT
        
        
      	if (track->Charge()>0) histoNsigmaTPCtriton_vs_p_notof	    -> Fill (track->P(),nsigmaTPC);
        if (track->Charge()<0) histoNsigmaTPCantitriton_vs_p_notof -> Fill (track->P(),nsigmaTPC);
        
        
        if (PassedTOFSelection(track))  {
            
            if (track->Charge()>0) histoNsigmaTPCtriton_vs_pt     -> Fill (track->Pt(),nsigmaTPC);
            if (track->Charge()<0) histoNsigmaTPCantitriton_vs_pt -> Fill (track->Pt(),nsigmaTPC);
	        if (track->Charge()>0) histoNsigmaTPCtriton_vs_p      -> Fill (track->P(),nsigmaTPC);
            if (track->Charge()<0) histoNsigmaTPCantitriton_vs_p  -> Fill (track->P(),nsigmaTPC);
        
            if (track->Charge()>0) histoNsigmaTPCtriton_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTPC(track));
            if (track->Charge()<0) histoNsigmaTPCantitriton_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTPC(track));
        
        
        }
        
        //TOF Signal vs. pT
        if (PassedTPCSelection(track))  {
            
            if (track->Charge()>0) histoNsigmaTOFtriton_vs_pt     -> Fill (track->Pt(),nsigmaTOF);
            if (track->Charge()<0) histoNsigmaTOFantitriton_vs_pt -> Fill (track->Pt(),nsigmaTOF);
            if (track->Charge()>0) histoNsigmaTOFtriton_vs_p      -> Fill (track->P(),nsigmaTOF);
            if (track->Charge()<0) histoNsigmaTOFantitriton_vs_p  -> Fill (track->P(),nsigmaTOF);
        
            
            if (track->Charge()>0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFtriton_vs_pt_trd     -> Fill (track->Pt(),nsigmaTOF);
            if (track->Charge()<0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFantitriton_vs_pt_trd -> Fill (track->Pt(),nsigmaTOF);
            
            if (track->Charge()>0) histoNsigmaTOFtriton_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTOF(track));
            if (track->Charge()<0) histoNsigmaTOFantitriton_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTOF(track));
        
            
            
        }
    }
    
    
    PostData(1, fOutputList);
    //PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::GetInputEvent ()  {
    
    //Get Input Event
    fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDevent) return false;
    histoNumberOfEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fESDeventCuts.AcceptEvent(fESDevent)) {
        //PostData(2, fQAList);
        return false;
    }
    histoNumberOfEvents -> Fill(1.5);
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    histoNumberOfEvents -> Fill(2.5);
    Double_t centrality = multiplicitySelection->GetMultiplicityPercentile(fCentralityEstimator);
    
   
    //Selection of Centrality Range
    if (centrality<fCentralityMin || centrality>=fCentralityMax ) return false;
    histoNumberOfEvents -> Fill(3.5);
    
    //Primary Vertex
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    histoNumberOfEvents -> Fill(4.5);
    
    //Primary Vertex Selection
    if ( vertex->GetZ() < fVertexZmin ) return false;
    if ( vertex->GetZ() > fVertexZmax ) return false;
    histoNumberOfEvents -> Fill(5.5);
    
    if ( vertex->GetNContributors() < fNumberVertexContributorsMin ) return false;
    histoNumberOfEvents -> Fill(6.5);
    
   
    return true;
    
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTrackQualityCuts (AliESDtrack* track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    if ( track->GetTPCsignalN() < 50 )         return passedTrkSelection;
    
    //Rapidity Calculation
    Double_t m  = AliPID::ParticleMass(AliPID::kTriton);
    Double_t px = track -> Px();
    Double_t py = track -> Py();
    Double_t pz = track -> Pz();
    Double_t E = TMath::Sqrt(m*m + px*px + py*py + pz*pz);
    TLorentzVector P (px,py,pz,E);
    Double_t y = P.Rapidity();
    
    
    //Kinematic Cuts & Acceptance
    if ( track->Pt()<fPtMin || track->Pt()>fPtMax ) return false;
    if ( TMath::Abs(track->Eta()) > fEtaMax )       return false;
    if ( TMath::Abs(y) > fYMax )       return false;
    
    //Track Selection Cuts
    if ( track->GetITSNcls() < fNumberClustersITSMin ) return false;
    if ( track->GetTPCNcls() < fNumberClustersTPCMin ) return false;
    if ( track->GetTPCCrossedRows() < fNumberCrossedRowsTPCMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberClustersTPCdEdxMin ) return false;
    if ( (((Double_t)track->GetTPCchi2())/((Double_t)track->GetTPCNcls())) > fChiSquarePerNDFMax) return false;
    
    //ITS Requirement
    Bool_t hitInITSLayer0 = track->HasPointOnITSLayer(0);
    Bool_t hitInITSLayer1 = track->HasPointOnITSLayer(1);
    

    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kFirst")==0  && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kSecond")==0 && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kAny")==0    && (!hitInITSLayer0) && (!hitInITSLayer1)) return false;

    
    //DCA Cuts
    Double_t dcaxy = GetDCAxy (track);
    Double_t dcaz  = GetDCAz (track);
    if (TMath::Abs(dcaxy) > fDCAxyMax) return false;
    if (TMath::Abs(dcaz)  > fDCAzMax)  return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTrackQualityCutsNoDCA (AliESDtrack* track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    if ( track->GetTPCsignalN() < 50 )         return passedTrkSelection;
    
    
    //Rapidity Calculation
    Double_t m  = AliPID::ParticleMass(AliPID::kTriton);
    Double_t px = track -> Px();
    Double_t py = track -> Py();
    Double_t pz = track -> Pz();
    Double_t E = TMath::Sqrt(m*m + px*px + py*py + pz*pz);
    TLorentzVector P (px,py,pz,E);
    Double_t y = P.Rapidity();
    
    
    //Kinematic Cuts & Acceptance
    if ( track->Pt()<fPtMin || track->Pt()>fPtMax ) return false;
    if ( TMath::Abs(track->Eta()) > fEtaMax )       return false;
    if ( TMath::Abs(y) > fYMax )       return false;
    
    //Track Selection Cuts
    if ( track->GetITSNcls() < fNumberClustersITSMin ) return false;
    if ( track->GetTPCNcls() < fNumberClustersTPCMin ) return false;
    if ( track->GetTPCCrossedRows() < fNumberCrossedRowsTPCMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberClustersTPCdEdxMin ) return false;
    if ( (((Double_t)track->GetTPCchi2())/((Double_t)track->GetTPCNcls())) > fChiSquarePerNDFMax) return false;
    
    //ITS Requirement
    Bool_t hitInITSLayer0 = track->HasPointOnITSLayer(0);
    Bool_t hitInITSLayer1 = track->HasPointOnITSLayer(1);
    

    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kFirst")==0  && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kSecond")==0 && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kAny")==0    && (!hitInITSLayer0) && (!hitInITSLayer1)) return false;

    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::IsCleanTritonCandidate (AliESDtrack *track)  {
    
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::Centered_nsigmaTPC (AliESDtrack *track)  {
   
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
   
    Double_t mean_fitted  = fpar0_mean_TPC*exp(fpar1_mean_TPC*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TPC*(track->P());
    
    nsigmaTPC = (nsigmaTPC - mean_fitted)/sigma_fitted;
   
    return nsigmaTPC;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::Centered_nsigmaTOF (AliESDtrack *track)  {
   
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
   
    Double_t mean_fitted  = fpar0_mean_TOF*exp(fpar1_mean_TOF*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TOF*exp(fpar1_sigma_TOF*(track->P()));
    
    nsigmaTOF = (nsigmaTOF - mean_fitted)/sigma_fitted;
   
    return nsigmaTOF;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTOFSelection (AliESDtrack *track)  {
    
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTPCSelection (AliESDtrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;
    if ((track->GetStatus()&AliESDtrack::kTOFout)==0) return false;

    //TPC-TOF Matching
//   if ((track->GetStatus()&AliESDtrack::kTOFout)==0) hasTOFhit=0;//Track with no TOF hit
//   if ((track->GetStatus()&AliESDtrack::kTOFout)!=0) hasTOFhit=1;//Track with TOF hit
    
    return true;
}

//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::GetDCAxy (AliESDtrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAxy = impactParameter[0];
    
    return DCAxy;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::GetDCAz (AliESDtrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAz = impactParameter[1];
    
    return DCAz;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonESD_PbPb::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

