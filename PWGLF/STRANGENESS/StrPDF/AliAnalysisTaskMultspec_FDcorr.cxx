class TTree;
class TParticle;
class TVector3;

class AliESDVertex;
class AliESDv0;
//class AliESDcascade;
class AliAODVertex;
class AliAODv0;
class AliAODcascade;

#include <Riostream.h>
#include "TFile.h"
#include "TCanvas.h"
#include "THistManager.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliPID.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliESDcascade.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisTaskMultspec_FDcorr.h"


ClassImp(AliAnalysisTaskMultspec_FDcorr)

AliAnalysisTaskMultspec_FDcorr::AliAnalysisTaskMultspec_FDcorr() : AliAnalysisTaskSE(),
fHistos_misc(nullptr), 
fTree(nullptr), 
fPIDResponse(0), 
fEventCuts(0),
ffillV0(nullptr)
{
  //default constructor
}

AliAnalysisTaskMultspec_FDcorr::AliAnalysisTaskMultspec_FDcorr(const char *name, TString lExtraOptions) : AliAnalysisTaskSE(name),
fHistos_misc(nullptr),
fTree(nullptr),
fPIDResponse(0),
fEventCuts(0),
ffillV0(nullptr)
{

  //Standard output
  DefineOutput(1, TList::Class()); // Miscellaneous Histograms
  DefineOutput(2, TTree::Class()); // Output Tree
}

AliAnalysisTaskMultspec_FDcorr::~AliAnalysisTaskMultspec_FDcorr()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    if (fHistos_misc) {
        delete fHistos_misc;
        fHistos_misc = 0x0;
    }
    if (fTree) {
        delete fTree;
        fTree = 0x0;
    }
 
 
}

//________________________________________________________________________
void AliAnalysisTaskMultspec_FDcorr::UserCreateOutputObjects()
{

  //miscellaneous histograms
  fHistos_misc = new THistManager("fHistos_misc");
  fHistos_misc->CreateTH1("hNevts", "", 1, 0, 1);  // #events histo
  fHistos_misc->CreateTH1("hEvSelCode", "", 1000, 0, 1000, "s");  //event selection code, for debug
  fHistos_misc->CreateTH1("hPtGen_Lam_prim", "", 400, 0, 20, "s");
  fHistos_misc->CreateTH1("hPtGen_ALam_prim", "", 400, 0, 20, "s");
  fHistos_misc->CreateTH1("hPtGen_Lam_seco", "", 400, 0, 20, "s");
  fHistos_misc->CreateTH1("hPtGen_ALam_seco", "", 400, 0, 20, "s");
  fHistos_misc->CreateTH1("hPtGen_Xi", "", 400, 0, 20, "s");
  fHistos_misc->CreateTH1("hPtGen_AXi", "", 400, 0, 20, "s");
  
  //tree
  ffillV0 = new LamFiller_FD;
  fTree = new TTree("fTree","Lambdas");
  fTree->Branch("LamFiller_FD", ffillV0);
  
  // PID Setup
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  inputHandler->SetNeedField();
  
  //fEventCuts Setup
  fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);

  //Output posting
  DataPosting();

}// end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskMultspec_FDcorr::UserExec(Option_t *)
{
  //get event from the input handler and cast it into the desired type of event
  AliAODEvent *lAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!lAODevent) { 
    AliWarning("ERROR: event not available \n"); 
    DataPosting(); 
    return; 
  }

  //get MC event and fill histo
  AliMCEvent  *lMCev  = MCEvent();
  if (!lMCev) {
    Printf("ERROR: Could not retrieve MC event in file %s\n",fInputHandler->GetTree()->GetCurrentFile()->GetName());
    DataPosting();
    return;
  }

 //get MC header and MC array
  AliAODMCHeader* header = static_cast<AliAODMCHeader*>(lAODevent->FindListObject(AliAODMCHeader::StdBranchName()));
  if (!header) {
    AliWarning("No header found.");
    DataPosting();
    return;
  }
  TClonesArray* MCTrackArray = dynamic_cast<TClonesArray*>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (MCTrackArray == NULL){
    AliWarning("No MC track array found.");
    DataPosting();
    return;
  }
  
  // take track of number of analyzed events
  fHistos_misc->FillTH1("hNevts", 0.5);
  
  // accepting events that pass internal selections in "AliEventCuts" 
  bool isEvtAccepted = fEventCuts.AcceptEvent(lAODevent);

  if (!isEvtAccepted) {
    DataPosting(); 
    return;
  }

  // multiplicity Information 
  int lEvSelCode = -666;
  AliMultSelection *MultSelection = (AliMultSelection*) lAODevent->FindListObject("MultSelection");
  if (!MultSelection) { 
    AliWarning("AliMultSelection object not found!"); 
    DataPosting(); 
    return; 
  } else {
    lEvSelCode = MultSelection->GetEvSelCode(); //==0 means event is good. Set by AliMultSelectionTask
    fHistos_misc->FillTH1("hEvSelCode", lEvSelCode);
  }
  if(lEvSelCode != 0) {   //skip everything if event selection code !=0
    DataPosting(); 
    return; 
  }

  //acquire best PV
  const AliVVertex *lBestPrimVtx = lAODevent->GetPrimaryVertex();
  double lBestPV[3] = {-666., -666., -666.};
  lBestPrimVtx->GetXYZ(lBestPV);
  float lfBestPV[3];
  for (int i=0; i<3; i++) lfBestPV[i] = (float)lBestPV[i];

  //acquire magnetic field
  double lMagField = -666;
  lMagField = lAODevent->GetMagneticField();
  
  int MCtracks = 0;
  
  // start v0 part
  int nv0s = 0;
  nv0s = (int) lAODevent->GetNumberOfV0s();

  //loop over generated MC
  MCtracks = lMCev->GetNumberOfTracks();

  int genv0 = 0;
  for (int i_MCtrk = 0;  i_MCtrk < MCtracks; i_MCtrk++){
    AliVParticle *lPart = (AliAODMCParticle*) lMCev->GetTrack(i_MCtrk);
    if(!lPart || lPart->Y()<-0.5 || lPart->Y()>0.5) continue;
    if( lPart->PdgCode()==3122 && lPart->IsPhysicalPrimary() )  {fHistos_misc->FillTH1("hPtGen_Lam_prim", lPart->Pt()); genv0++;}
    if( lPart->PdgCode()==-3122 && lPart->IsPhysicalPrimary() ) {fHistos_misc->FillTH1("hPtGen_ALam_prim", lPart->Pt()); genv0++;}
    if( lPart->PdgCode()==3312 && lPart->IsPhysicalPrimary() )  fHistos_misc->FillTH1("hPtGen_Xi", lPart->Pt());
    if( lPart->PdgCode()==-3312 && lPart->IsPhysicalPrimary() ) fHistos_misc->FillTH1("hPtGen_AXi", lPart->Pt());
    if( lPart->PdgCode()==3122 && !lPart->IsPhysicalPrimary() && ((AliAODMCParticle*) lMCev->GetTrack(lPart->GetMother()))->PdgCode()==3312 )  fHistos_misc->FillTH1("hPtGen_Lam_seco", lPart->Pt());
    if( lPart->PdgCode()==-3122 && !lPart->IsPhysicalPrimary() && ((AliAODMCParticle*) lMCev->GetTrack(lPart->GetMother()))->PdgCode()==-3312 )  fHistos_misc->FillTH1("hPtGen_ALam_seco", lPart->Pt());
   }

   //start loop over V0 candidates
   for (int iV0 = 0; iV0 < nv0s; iV0++) {

     //get the iV0_th candidate in the event
     AliAODv0 *v0 = lAODevent->GetV0(iV0);
     if (!v0 || v0->GetOnFlyStatus()) continue;
     
     //retrieve daughter AODTracks
     AliAODTrack *pTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(0);
     AliAODTrack *nTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(1);
     if (!pTrack || !nTrack) { AliWarning("ERROR: Could not retrieve one of the daughter tracks"); continue; }
     
     //labels
     int pTrack_lbl = (int) TMath::Abs(pTrack->GetLabel());
     int nTrack_lbl = (int) TMath::Abs(nTrack->GetLabel());
     AliMCParticle *mc_pTrack = (AliMCParticle*) lMCev->GetTrack(pTrack_lbl);
     AliMCParticle *mc_nTrack = (AliMCParticle*) lMCev->GetTrack(nTrack_lbl);
     if (!mc_pTrack || !mc_nTrack) { AliWarning("ERROR: Could not retrieve one of the daughter MC tracks"); continue; }
     int pTrack_moth_lbl = mc_pTrack->GetMother();
     int nTrack_moth_lbl = mc_nTrack->GetMother();
     int v0num = -666;
     if(pTrack_moth_lbl == nTrack_moth_lbl) v0num = pTrack_moth_lbl;
     else continue;

     bool islam  = kFALSE;
     AliMCParticle *mc_v0 = (AliMCParticle*) lMCev->GetTrack(v0num);
     if (!mc_v0) { AliWarning("ERROR: Could not retrieve v0 MC track"); continue; }
     if(TMath::Abs(mc_v0->PdgCode())==3122) islam = kTRUE;
     if(!islam) continue;

     //primary status
     ffillV0->IsPrimary = mc_v0->IsPhysicalPrimary() ? kTRUE : kFALSE ;

     //GetMother
     int pdgmother = 0;
     AliMCParticle *mc_v0_mother = 0x0;
     if(!ffillV0->IsPrimary){
        mc_v0_mother =  (AliMCParticle*) lMCev->GetTrack(mc_v0->GetMother());
        if (!mc_v0_mother) { AliWarning("ERROR: Could not retrieve v0 MC mother track"); continue; }
        pdgmother = mc_v0_mother->PdgCode();
        if(!ffillV0->IsPrimary && TMath::Abs(pdgmother)!= 3312) continue;
     }

     //get the 2D radius of the decay vertex
     ffillV0->V0Rad = v0->RadiusV0();
  
     //pt and rapidity (cut the unwanted rapidities)
     ffillV0->Pt_lam = v0->Pt();
     if( !ffillV0->IsPrimary ) ffillV0->Pt_xi = mc_v0_mother->Pt();
     if( TMath::Abs(v0->RapLambda())>0.5 ) continue;
  
     //Daughters' Eta must be inside 0.8 to be true
     bool lV0_etaPos = TMath::Abs(pTrack->Eta())<0.8 ? kTRUE : kFALSE ;
     bool lV0_etaNeg = TMath::Abs(nTrack->Eta())<0.8 ? kTRUE : kFALSE ;
     if(!lV0_etaPos || !lV0_etaNeg) continue;
  
     //like-sign V0s (if any)
     if(((int) pTrack->GetSign()) == ((int) nTrack->GetSign()))  continue;
  
     //crossed raws
     double_t lCrosRawsPos = pTrack->GetTPCClusterInfo(2, 1);
     double_t lCrosRawsNeg = nTrack->GetTPCClusterInfo(2, 1);
     ffillV0->LeastCRaws = (int) TMath::Min(lCrosRawsPos, lCrosRawsNeg);

     //crossed raws / Findable clusters
     if(pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0) ffillV0->LeastCRawsOvF = 666;
     else{
       double_t lCrosRawsOvFPos = lCrosRawsPos/((double) (pTrack->GetTPCNclsF()));
       double_t lCrosRawsOvFNeg = lCrosRawsNeg/((double) (nTrack->GetTPCNclsF()));
       ffillV0->LeastCRawsOvF = TMath::Min(lCrosRawsOvFPos, lCrosRawsOvFNeg);
     }

     //dca info
     ffillV0->DcaPosToPV  = v0->DcaPosToPrimVertex();
     ffillV0->DcaNegToPV  = v0->DcaNegToPrimVertex();
     ffillV0->DcaV0Daught = v0->DcaV0Daughters();
  
     //cosPA
     ffillV0->V0CosPA = v0->CosPointingAngle(lBestPV);
  
     //particle antiparticle
     printf("pdgcode: %d\n",mc_v0->PdgCode());
     ffillV0->IsParticle = mc_v0->PdgCode() > 0 ? kTRUE : kFALSE ;

     //PID
     if( ffillV0->IsParticle ) {
         ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
         ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
     }
     else {
         ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
         ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton);
     }
  
     //distance over total momentum
     ffillV0->DistOverTotP = v0->DecayLengthV0(lBestPV)/(v0->P()+1e-10);//avoid division by zero
  
     //TOF matching
     bool lNegTOFmatching = ( nTrack->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE );
     bool lPosTOFmatching = ( pTrack->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE );
     ffillV0->TOFmatch = (!lNegTOFmatching && !lPosTOFmatching) ? kFALSE : kTRUE;

     //ITS matching
     bool lNegITSmatching = ( !(nTrack->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE );
     bool lPosITSmatching = ( !(pTrack->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE );
     ffillV0->ITSmatch = (!lNegITSmatching && !lPosITSmatching) ? kFALSE : kTRUE;

     //tree filling
     fTree->Fill();

   } // end of V0 loop

  DataPosting();

}

//________________________________________________________________________
void AliAnalysisTaskMultspec_FDcorr::Terminate(Option_t *)
{

}

//________________________________________________________________________
void AliAnalysisTaskMultspec_FDcorr::DataPosting() {

  PostData(1, fHistos_misc->GetListOfHistograms());
  PostData(2, fTree);
  
}
