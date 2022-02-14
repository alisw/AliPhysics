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

#include "AliAnalysisTaskV0multspec.h"


ClassImp(AliAnalysisTaskV0multspec)

AliAnalysisTaskV0multspec::AliAnalysisTaskV0multspec() : AliAnalysisTaskSE(),
fHistos_misc(nullptr), 
fTree(nullptr), 
fpart(666),
fisMC(kTRUE),
fPIDResponse(0), 
ffillV0(nullptr),
ffillCasc(nullptr)
{
  //default constructor
}

AliAnalysisTaskV0multspec::AliAnalysisTaskV0multspec(const char *name, int particle, TString lExtraOptions, bool ismc) : AliAnalysisTaskSE(name),
fHistos_misc(nullptr), 
fTree(nullptr), 
fpart(particle),
fisMC(ismc),
fPIDResponse(0), 
ffillV0(nullptr),
ffillCasc(nullptr)
{

  //SetAddCasc(AddCasc);
  SetIsMC(ismc);
  SetParticle(particle);
    
  //Standard output
  DefineOutput(1, TList::Class()); // Miscellaneous Histograms
  DefineOutput(2, TTree::Class()); // Output Tree
  
}

AliAnalysisTaskV0multspec::~AliAnalysisTaskV0multspec()
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
void AliAnalysisTaskV0multspec::UserCreateOutputObjects()
{

  //miscellaneous histograms
  fHistos_misc = new THistManager("fHistos_misc");
  fHistos_misc->CreateTH1("hNevts", "", 1, 0, 1);  // #events histo
  fHistos_misc->CreateTH1("hmultspec_nocut", "", 1000, 0, 1000, "s");  // multiplicity spectrum for checking
  fHistos_misc->CreateTH2("hmultVSEvmult_nocut", "", 1000, 0, 1000, 100,0,100, "s");  //particle multiplicity VS Event multiplicity spectrum for checking
  fHistos_misc->CreateTH1("hEvSelCode", "", 1000, 0, 1000, "s");  //event selection code, for debug
  
  //MC-related TH
  if(fisMC) {
      fHistos_misc->CreateTH2("hR0gen", "", 200, -100, 100, 100, 0, 100, "s");  //histogram Gen_distrib for zero-Rec events for given part VS mult
      fHistos_misc->CreateTH2("hR0genBoth", "", 200, -100, 100, 100, 0, 100, "s");  //same as before, but for particle and anti-particle (for cascades)
      fHistos_misc->CreateTH2("hmultspec_Gen", "", 200, -100, 100, 100, 0, 100, "s");  //histogram Gen_distrib given part VS mult
  }
  
  //histos for "zero" issue
  fHistos_misc->CreateTH1("hZeros", "", 200, -100, 100, "s");  //events with zero positive or negative rec cand (no cand at all or not passing rough selections)
  fHistos_misc->CreateTH1("hZerosBoth", "", 100, 0, 100, "s");  //events with both zero pos and neg rec cascade (no casc at all or not passing rough selections)
  fHistos_misc->CreateTH1("hPlusAndMinFilled", "", 100, 0, 100, "s");  //histo for debug
  
  //Events passing different cuts
  //TString labtext[5] = {"PassEvSelCode",">0Cand","FillTree","FillTreeP","FillTreeM"};
  fHistos_misc->CreateTH1("TH1EventMonitor", "", 5, 0, 5);
  //((TH1D*)fHistos_misc->FindObject("TH1EventMonitor"))->GetXaxis()->SetLabelOffset(0.05);
  //((TH1D*)fHistos_misc->FindObject("TH1EventMonitor"))->GetXaxis()->SetNdivisions(5);
  fHistos_misc->CreateTH2("EventMonitor", "", 5, 0, 5, 101, 0, 101);
  //((TH2D*)fHistos_misc->FindObject("EventMonitor"))->GetXaxis()->SetLabelOffset(0.05);
  //((TH2D*)fHistos_misc->FindObject("EventMonitor"))->GetXaxis()->SetNdivisions(5);
  //for (int i=0; i<5; i++) {
  //    ((TH1D*)fHistos_misc->FindObject("TH1EventMonitor"))->GetXaxis()->ChangeLabel(i+1,30,-1,13,-1,-1,labtext[i].Data());
  //    ((TH2D*)fHistos_misc->FindObject("EventMonitor"))->GetXaxis()->ChangeLabel(i+1,30,-1,13,-1,-1,labtext[i].Data());
  //}
  
  //tree
  if(fpart<kxi){
    ffillV0 = new V0filler;
    fTree = new TTree("fTree","V0s");
    fTree->Branch("V0filler", ffillV0);
  }
  else {
    ffillCasc = new Cascfiller;
    fTree = new TTree("fTree","Cascades");
    fTree->Branch("Cascfiller",ffillCasc);
  }
  
  // PID Setup
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  inputHandler->SetNeedField();

  //Output posting
  DataPosting();

}// end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisTaskV0multspec::UserExec(Option_t *)
{

  //get event from the input handler and cast it into the desired type of event
  AliAODEvent *lAODevent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!lAODevent) { 
    AliWarning("ERROR: event not available \n"); 
    DataPosting(); 
    return; 
  }

  //get MC event and fill histo
  AliMCEvent  *lMCev  = 0x0;
  if(fisMC){
    lMCev = MCEvent();
    if (!lMCev) {
      Printf("ERROR: Could not retrieve MC event in file %s\n",fInputHandler->GetTree()->GetCurrentFile()->GetName());
      DataPosting(); 
      return;
    }
  }
  //get MC header and MC array
  AliAODMCHeader* header = 0x0;
  TClonesArray* MCTrackArray = 0x0;
  if(fisMC) {
      header = static_cast<AliAODMCHeader*>(lAODevent->FindListObject(AliAODMCHeader::StdBranchName()));
      if (!header) {
        AliWarning("No header found.");
        DataPosting(); 
        return;
      } 
      MCTrackArray = dynamic_cast<TClonesArray*>(lAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (MCTrackArray == NULL){
        AliWarning("No MC track array found.");
        DataPosting(); 
        return;
      }  
  }
  
  // take track of number of analyzed events
  fHistos_misc->FillTH1("hNevts", 0.5);

  //get trigger information
  //fTriggerMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  // multiplicity Information 
  int lEvSelCode = -666;
  double lMultiplicity = -666.;
  AliMultSelection *MultSelection = (AliMultSelection*) lAODevent->FindListObject("MultSelection");
  if (!MultSelection) { 
    AliWarning("AliMultSelection object not found!"); 
    DataPosting(); 
    return; 
  } else {
    lMultiplicity = MultSelection->GetMultiplicityPercentile("V0M");
    if(fpart<kxi) ffillV0->Mult = lMultiplicity;
    else ffillCasc->Mult = lMultiplicity;
    lEvSelCode = MultSelection->GetEvSelCode(); //==0 means event is good. Set by AliMultSelectionTask
    fHistos_misc->FillTH1("hEvSelCode", lEvSelCode);
  }
  if(lEvSelCode != 0) {   //skip everything if event selection code !=0
    DataPosting(); 
    return; 
  }

  //fill event monitor with all good events
  fHistos_misc->FillTH2("EventMonitor", 0.5, 100.5);
  fHistos_misc->FillTH2("EventMonitor", 0.5, lMultiplicity);
  fHistos_misc->FillTH1("TH1EventMonitor", 0.5);
  
  //acquire best PV
  const AliVVertex *lBestPrimVtx = lAODevent->GetPrimaryVertex();
  double lBestPV[3] = {-666., -666., -666.};
  lBestPrimVtx->GetXYZ(lBestPV);
  float lfBestPV[3];
  for (int i=0; i<3; i++) lfBestPV[i] = (float)lBestPV[i];

  //acquire magnetic field
  double lMagField = -666;
  lMagField = lAODevent->GetMagneticField();

  // start v0 part
  if(fpart<kxi){
    
    int nv0s = 0;
    nv0s = (int) lAODevent->GetNumberOfV0s();
    fHistos_misc->FillTH1("hmultspec_nocut",nv0s);
    fHistos_misc->FillTH2("hmultVSEvmult_nocut",nv0s,lMultiplicity);
    if(nv0s!=0)   {
      fHistos_misc->FillTH2("EventMonitor", 1.5, 100.5);
      fHistos_misc->FillTH2("EventMonitor", 1.5, lMultiplicity);
      fHistos_misc->FillTH1("TH1EventMonitor", 1.5);
    }
    bool newV0 = kTRUE;
    bool isZeroV0 = kTRUE; //if this flag remains true, either all V0s in the loop have been discarded by "continue" or there were no candidates at all
    
    //loop over generated MC 
    int initgens = 253;
    ffillV0->nGen  = initgens;
    int ngen = 0;
    if(fisMC){
      for (int i_MCtrk = 0;  i_MCtrk < lMCev->GetNumberOfTracks(); i_MCtrk++){
        AliVParticle *lPart = (AliAODMCParticle*) lMCev->GetTrack(i_MCtrk);
        if(!lPart || lPart->Y()<-0.5 || lPart->Y()>0.5 || !lPart->IsPhysicalPrimary()) continue;
        if( fpart==kk0s && lPart->PdgCode()==310   ) ngen++;
        if( fpart==klam && lPart->PdgCode()==3122  ) ngen++;
        if( fpart==kalam && lPart->PdgCode()==-3122 ) ngen++;
      }
      ffillV0->nGen  = ngen;
      fHistos_misc->FillTH2("hmultspec_Gen", lMultiplicity, ngen);      //fill generated histo
    }
    
    //start loop over V0 candidates
    for (int iV0 = 0; iV0 < nv0s; iV0++) {
        
      //set a fV0_IsNewEvt flag if this is the first V0 of the loop
      ffillV0->IsNewEvt = newV0;
  
      //get the iV0_th candidate in the event
      AliAODv0 *v0 = lAODevent->GetV0(iV0);
      if (!v0 || v0->GetOnFlyStatus()) continue;
  
      //get the 2D radius of the decay vertex
      ffillV0->V0Rad = v0->RadiusV0();
  
      //pt and rapidity (cut the unwanted rapidities)
      ffillV0->Pt = v0->Pt();
      if( fpart==kk0s && TMath::Abs(v0->RapK0Short())>0.5 ) continue;
      if( (fpart==klam || fpart==kalam) && TMath::Abs(v0->RapLambda())>0.5 ) continue;
  
      //retrieve daughter AODTracks
      AliAODTrack *pTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(0);
      AliAODTrack *nTrack = (AliAODTrack*) v0->GetSecondaryVtx()->GetDaughter(1);
      if (!pTrack || !nTrack) { AliWarning("ERROR: Could not retrieve one of the daughter tracks"); continue; }
  
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
  
      // invariant mass with different hypotheses
      if( fpart==kk0s ) {
          if ( v0->MassK0Short()<0.4||v0->MassK0Short()>0.6 ) continue;
          else ffillV0->InvMass = v0->MassK0Short();
      }
      else if( fpart==klam ){
          if ( v0->MassLambda()<1.08||v0->MassLambda()>1.16 ) continue;
          else ffillV0->InvMass = v0->MassLambda();
      }
      else if( fpart==kalam ){
          if ( v0->MassAntiLambda()<1.08||v0->MassAntiLambda()>1.16 ) continue;
          else ffillV0->InvMass = v0->MassAntiLambda();
      }
  
      //PID
      if( fpart==kk0s ) {
          ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion);
          ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
      }
      else if( fpart==klam ) {
          ffillV0->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton);
          ffillV0->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion);
      }
      else if( fpart==kalam ) {
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
      
      //MC-related PDG matching flag
      if(fisMC) {
        int labMothPosDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrack->GetLabel())))->GetMother();
        int labMothNegDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrack->GetLabel())))->GetMother();
        int pdgV0 = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(labMothPosDaught)))-> PdgCode();
        if((labMothPosDaught==labMothNegDaught) && (pdgV0==PDGcodes[fpart])) ffillV0->IsPDGmatched = kTRUE;
        else ffillV0->IsPDGmatched = kFALSE;
      }
      
      //tree filling 
      isZeroV0 = kFALSE; //if I fill at least once in this event, this becomes kFALSE
      fTree->Fill();
      newV0 = kFALSE;
  
    } // end of V0 loop

    if(!isZeroV0) {
        fHistos_misc->FillTH2("EventMonitor", 2.5, 100.5);
        fHistos_misc->FillTH2("EventMonitor", 2.5, lMultiplicity);
        fHistos_misc->FillTH1("TH1EventMonitor", 2.5);
    }
    else{
        //fill histos with zero candidates
        fHistos_misc->FillTH1("hZeros",lMultiplicity);
        if(fisMC) fHistos_misc->FillTH2("hR0gen", lMultiplicity, ngen);
    }    
    
  } //end V0 part
  else{  // start cascade part

    int ncasc = 0;
    ncasc = lAODevent->GetNumberOfCascades();
    fHistos_misc->FillTH1("hmultspec_nocut",ncasc);
    fHistos_misc->FillTH2("hmultVSEvmult_nocut",ncasc,lMultiplicity);
    if(ncasc!=0)   {
      fHistos_misc->FillTH2("EventMonitor", 1.5, 100.5);
      fHistos_misc->FillTH2("EventMonitor", 1.5, lMultiplicity);
      fHistos_misc->FillTH1("TH1EventMonitor", 1.5);
    }
    bool newCasc = kTRUE;
    bool isZeroCascM = kTRUE; //if this flag remains true, either all negative Casc in the loop have been discarded by "continue" or there were no candidates at all
    bool isZeroCascP = kTRUE;
    
    //loop over generated MC 
    int initgens = 253;
    ffillCasc->nGenP  = initgens;
    ffillCasc->nGenM  = initgens;
    int ngenP = 0;
    int ngenM = 0;
    if(fisMC){
      for (int i_MCtrk = 0;  i_MCtrk < lMCev->GetNumberOfTracks(); i_MCtrk++){
        AliVParticle *lPart = (AliAODMCParticle*) lMCev->GetTrack(i_MCtrk);
        if(!lPart || lPart->Y()<-0.5 || lPart->Y()>0.5 || !lPart->IsPhysicalPrimary()) continue;
        if( fpart==kxi && lPart->PdgCode()==-3312  ) ngenP++;
        if( fpart==kxi && lPart->PdgCode()==3312   ) ngenM++;
        if( fpart==kom && lPart->PdgCode()==-3334  ) ngenP++;
        if( fpart==kom && lPart->PdgCode()==3334   ) ngenM++;
      }
      ffillCasc->nGenP  = ngenP;
      ffillCasc->nGenM  = ngenM;
      //fill generated histo
      fHistos_misc->FillTH2("hmultspec_Gen", lMultiplicity, ngenP);
      fHistos_misc->FillTH2("hmultspec_Gen", -lMultiplicity, ngenM);
    }  
  
    //start casc loop
    for (int i_casc = 0; i_casc < ncasc; i_casc++) {

      //set a fCasc_IsNewEvt flag if this is the first V0 of the loop
      ffillCasc->IsNewEvt = newCasc;
    
      //get the iCasc_th candidate in the event
      AliAODcascade *casc = lAODevent->GetCascade(i_casc);
      if (!casc || casc->GetOnFlyStatus()) continue;
    
      //charge
      ffillCasc->charge = casc->ChargeXi()<0 ? kTRUE : kFALSE ;
    
      //cascade and V0 2D radii
      double lVtxCasc[3];
      lVtxCasc[0] = casc->DecayVertexXiX();
      lVtxCasc[1] = casc->DecayVertexXiY();
      lVtxCasc[2] = casc->DecayVertexXiZ();
      ffillCasc->CascRad = TMath::Sqrt(lVtxCasc[0]*lVtxCasc[0]+lVtxCasc[1]*lVtxCasc[1]);
      ffillCasc->V0Rad = casc->RadiusSecVtx();
    
      //transverse momentum
      ffillCasc->Pt = TMath::Sqrt(casc->Pt2Xi());
    
      //candidate's rapidity (mass hypothesis dependent)
      if(fpart==kxi && (TMath::Abs(casc->RapXi())>0.5)) continue;
      if(fpart==kom && (TMath::Abs(casc->RapOmega())>0.5)) continue;
    
      //get daughter tracks (positive, negative and bachelor)
      AliAODTrack *pTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(0));
      AliAODTrack *nTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDaughter(1));
      AliAODTrack *bTrackCasc = dynamic_cast<AliAODTrack*> (casc->GetDecayVertexXi()->GetDaughter(0));
      if (!pTrackCasc || !nTrackCasc || !bTrackCasc) { 
        AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...\n"); 
        continue; 
      }
    
      //daughters' etas must be inside 0.8 to be true
      bool letaPos = TMath::Abs(pTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      bool letaNeg = TMath::Abs(nTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      bool letaBac = TMath::Abs(bTrackCasc->Eta())<0.8 ? kTRUE : kFALSE ;
      if(!letaPos || !letaNeg || !letaBac) continue;
      
      //crossed raws
      double lCrosRawsPos = pTrackCasc->GetTPCClusterInfo(2, 1);
      double lCrosRawsNeg = nTrackCasc->GetTPCClusterInfo(2, 1);
      double lCrosRawsBac = bTrackCasc->GetTPCClusterInfo(2, 1);
      ffillCasc->LeastCRaws = (int) (lCrosRawsPos<lCrosRawsNeg ? std::min(lCrosRawsPos, lCrosRawsBac) : std::min(lCrosRawsNeg, lCrosRawsBac));
      
      //crossed raws / Findable clusters
      if(pTrackCasc->GetTPCNclsF()<=0 || nTrackCasc->GetTPCNclsF()<=0 || bTrackCasc->GetTPCNclsF()<=0) ffillCasc->LeastCRawsOvF = 0; //avoid division by zero
      else{    
        double lCrosRawsOvFPos = lCrosRawsPos / ((double)(pTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFNeg = lCrosRawsNeg / ((double)(nTrackCasc->GetTPCNclsF()));
        double lCrosRawsOvFBac = lCrosRawsBac / ((double)(bTrackCasc->GetTPCNclsF()));
        ffillCasc->LeastCRawsOvF = lCrosRawsOvFPos<lCrosRawsOvFNeg ? std::min(lCrosRawsOvFPos, lCrosRawsOvFBac) : std::min(lCrosRawsOvFNeg, lCrosRawsOvFBac);
      }
      
      //DCA info
      ffillCasc->DcaCascDaught = casc->DcaXiDaughters();
      ffillCasc->DcaBachToPV = casc->DcaBachToPrimVertex();
      ffillCasc->DcaPosToPV = casc->DcaPosToPrimVertex();
      ffillCasc->DcaNegToPV = casc->DcaNegToPrimVertex();
      ffillCasc->DcaV0Daught = casc->DcaV0Daughters();
      ffillCasc->DcaV0ToPV = casc->DcaV0ToPrimVertex();
    
      //cascade and V0 cosine of pointing angle
      ffillCasc->CascCosPA = casc->CosPointingAngleXi((const Double_t&) lBestPV[0], (const Double_t&) lBestPV[1], (const Double_t&) lBestPV[2]);
      ffillCasc->V0CosPA = casc->CosPointingAngle(lBestPV);
    
      //candidate's invariant mass
      if( fpart==kxi ) {
        if ( casc->MassXi()<1.26||casc->MassXi()>1.38 ) continue;
        else ffillCasc->InvMass = casc->MassXi();
      }
      else if( fpart==kom ){
        if ( casc->MassOmega()<1.62||casc->MassOmega()>1.73 ) continue;
        else ffillCasc->InvMass = casc->MassOmega();
      }
    
      //PID
      if(casc->ChargeXi()<0){
        ffillCasc->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kProton);
        ffillCasc->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kPion);
      }
      else {
        ffillCasc->NSigPos = fPIDResponse->NumberOfSigmasTPC(pTrackCasc, AliPID::kPion);
        ffillCasc->NSigNeg = fPIDResponse->NumberOfSigmasTPC(nTrackCasc, AliPID::kProton);
      }
      //check if bachelor's sign is always checked in compliance with Lambda/antilambda nature
      if(fpart==kxi) ffillCasc->NSigBac = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kPion);
      else ffillCasc->NSigBac = fPIDResponse->NumberOfSigmasTPC(bTrackCasc, AliPID::kKaon);
      if( casc->ChargeXi()<0 && bTrackCasc->Charge()>0 ) printf("Wrong charge association\n");
    
      //distance over total momentum
      ffillCasc->DistOverTotP = (TMath::Sqrt(TMath::Power(lVtxCasc[0]-lBestPV[0], 2)+TMath::Power(lVtxCasc[1]-lBestPV[1], 2)+TMath::Power(lVtxCasc[2]-lBestPV[2], 2)))/(TMath::Sqrt(casc->Ptot2Xi())+1e-10);
    
      //TOF matching
      bool lNegTOFmatching = nTrackCasc->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE ;
      bool lPosTOFmatching = pTrackCasc->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE ;
      bool lBacTOFmatching = bTrackCasc->GetTOFBunchCrossing(lMagField)<-95. ? kFALSE : kTRUE ;
      ffillCasc->TOFmatch = (!lNegTOFmatching && !lPosTOFmatching && !lBacTOFmatching) ? kFALSE : kTRUE;
      
      //ITS matching
      bool lNegITSmatching = !(nTrackCasc->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE ;
      bool lPosITSmatching = !(pTrackCasc->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE ;
      bool lBacITSmatching = !(bTrackCasc->GetStatus() & AliESDtrack::kITSrefit) ? kFALSE : kTRUE ;
      ffillCasc->ITSmatch = (!lNegITSmatching && !lPosITSmatching && !lBacITSmatching) ? kFALSE : kTRUE;
      
      //V0 daughter mass
      double_t imassla = 0.;
      (casc->ChargeXi()<0 ? imassla=casc->MassLambda() : imassla=casc->MassAntiLambda());
      if (TMath::Abs(imassla-1.115683)<0.008) ffillCasc->GoodInvMassLam = kTRUE;
      else ffillCasc->GoodInvMassLam = kFALSE;
    
      //calculate DCA Bachelor-Baryon to remove "bump" structure in InvMass
      ffillCasc->BacBarCosPA = casc->BachBaryonCosPA();
    
      //MC-related PDG matching flag
      if(fisMC) {
        int labMothPosDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(pTrackCasc->GetLabel())))->GetMother();
        int labMothNegDaught = ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(nTrackCasc->GetLabel())))->GetMother();
        int labelV0 = (int)TMath::Abs(labMothPosDaught);
        int pdgV0 = ((AliMCParticle*) lMCev->GetTrack(labelV0))-> PdgCode();
        int mothV0 = ((AliMCParticle*) lMCev->GetTrack(labelV0))->GetMother();
        int mothBac= ((AliMCParticle*) lMCev->GetTrack((int)TMath::Abs(bTrackCasc->GetLabel())))->GetMother();
        int labelCasc = (int)TMath::Abs(mothV0);
        int pdgCasc = ((AliMCParticle*) lMCev->GetTrack(labelCasc))-> PdgCode();
        if((labMothPosDaught!=labMothNegDaught) || (pdgV0!=PDGcodes[(casc->ChargeXi()<0 ? klam : kalam)])) ffillCasc->IsPDGmatched = kFALSE;
        else if((mothV0!=mothBac) || (pdgCasc!=(casc->ChargeXi()<0 ? PDGcodes[fpart] : -PDGcodes[fpart]))) ffillCasc->IsPDGmatched = kFALSE;
        else ffillCasc->IsPDGmatched = kTRUE;
      }
    
      //tree filling 
      if(casc->ChargeXi()<0) isZeroCascM = kFALSE; //if I fill at least once, this becomes kFALSE
      else isZeroCascP = kFALSE; //if I fill at least once, this becomes kFALSE
      fTree->Fill();
      newCasc = kFALSE;
      
    } //end casc loop
    
    //fill EventMonitor with events where the tree is filled (since isZeroCascP or isZeroCascM switched to kFALSE)
    if(!isZeroCascM || !isZeroCascP) {
      fHistos_misc->FillTH2("EventMonitor", 2.5, 100.5);
      fHistos_misc->FillTH2("EventMonitor", 2.5, lMultiplicity);
      fHistos_misc->FillTH1("TH1EventMonitor", 2.5);
      if(!isZeroCascP) {fHistos_misc->FillTH2("EventMonitor", 3.5, 100.5);
        fHistos_misc->FillTH2("EventMonitor", 3.5, lMultiplicity);
        fHistos_misc->FillTH1("TH1EventMonitor", 3.5);
      }
      if(!isZeroCascM) {fHistos_misc->FillTH2("EventMonitor", 4.5, lMultiplicity);
        fHistos_misc->FillTH2("EventMonitor", 4.5, 100.5);
        fHistos_misc->FillTH1("TH1EventMonitor", 4.5);
      }
    }
    
    //fill debug histo
    if(!isZeroCascM && !isZeroCascP) fHistos_misc->FillTH1("hPlusAndMinFilled",lMultiplicity);
    
    //remove from zero-plus and zero-minus those events which are zero for both
    bool isZeroBoth=kFALSE;
    if(isZeroCascM && isZeroCascP) {
        isZeroBoth=kTRUE;
        isZeroCascM=kFALSE;
        isZeroCascP=kFALSE;
    }
      
    //fill histos with zero candidates
    if(isZeroCascM) {
      fHistos_misc->FillTH1("hZeros",-lMultiplicity-1e-6);
      if(fisMC) fHistos_misc->FillTH2("hR0gen", -lMultiplicity-1e-6, ngenM);
    }
    if(isZeroCascP) {
      fHistos_misc->FillTH1("hZeros",lMultiplicity+1e-6);
      if(fisMC) fHistos_misc->FillTH2("hR0gen", lMultiplicity+1e-6, ngenP);
    }
    if(isZeroBoth) {
      fHistos_misc->FillTH1("hZerosBoth",lMultiplicity);
      if(fisMC){
        fHistos_misc->FillTH2("hR0genBoth", -lMultiplicity-1e-6, ngenM);
        fHistos_misc->FillTH2("hR0genBoth", lMultiplicity+1e-6, ngenP);
      }
    }
    
  }
  
  DataPosting();

}

//________________________________________________________________________
void AliAnalysisTaskV0multspec::Terminate(Option_t *)
{

}

//________________________________________________________________________
void AliAnalysisTaskV0multspec::DataPosting() {

  PostData(1, fHistos_misc->GetListOfHistograms());
  PostData(2, fTree);
  
}
